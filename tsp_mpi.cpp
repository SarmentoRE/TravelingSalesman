#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <iterator>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <mpi.h>
#include <stddef.h>
using namespace std;

#define NUMCITYTAG 1
#define NUMBLOCKTAG 2
#define GRIDSIZETAG 3
#define BLOCKSTOSENDTAG 4
#define INITIALBLOCKTAG 5

struct City
{
	int id;
	double x;
	double y;
};
struct Salesman
{
	double minCost;
	vector<int> minPath;
};
struct Block
{
	int id;
	double minCost;
	vector<City> minPath;
};
vector<Block> allBlocks;
int procNum;
int numProc;
double blockSpace;
int numBlocks;
int numCities;
int coords[2];

inline double distance(City city1, City city2) { return sqrt((pow((city2.x - city1.x),2)) + (pow((city2.y - city1.y),2))); }
inline bool isSquare(int x) { return (sqrt(x) - floor(sqrt(x)) == 0); }

double fRand(double min, double max) 
{
	double f = (double)rand() / RAND_MAX;
	return min + f * (max - min);
}

double** distanceMatrix(vector<City> loc)
{
	int i, j, size = loc.size();
	double **dMatrix = (double **)malloc(size * sizeof(double)); 
	for(i = 0; i < size; i++)
	{
		dMatrix[i] = (double *)malloc(size * sizeof(double));
	}

	for(i = 0; i < size; i++)
	{
		for(j = 0; j < size; j++)
		{
			dMatrix[i][j] = distance(loc[i], loc[j]);
		}
	}
	return dMatrix;
}

vector<vector<int>> generateSubsets(int size, int n)
{
	int count = 0;
	vector<vector<int> > allSets;
	vector<int> set;
	vector<bool> v((unsigned long)n);
	fill(v.begin(), v.begin() + size, true);

	do {
		for (int i = 0; i < n; ++i) {
			if (v[i]) {
				count++;
				set.push_back(i + 1);
				if (count == size) {
					allSets.push_back(set);
					set.clear();
					count = 0;
				}
			}
		}
	} while (prev_permutation(v.begin(), v.end()));
	return allSets;
}
	
vector<City> sortCities(vector<int> minPath, vector<City> cities){
	vector<City> sortedCities;
	for(int i = 0; i < minPath.size(); i++)
	{
		sortedCities.push_back(cities[minPath[i]]);
	}
	return sortedCities;
}
/*
void convertPath(Salesman salesman, vector<City> cities, int blockId)
{
	vector<City> sortedCities;
	Block block;
	vector<int> path = salesman.minPath;

	for(int i = 0; i < path.size(); i++)
	{
		sortedCities.push_back(cities[path[i]]);
	}

	block.minCost = salesman.minCost;
	block.minPath = sortedCities;
	block.id = blockId;
	block.start = true;
	block.end = true;
	allBlocks.push_back(block);
}
*/
void genKey(vector<int> &set, int z, long long &key)
{

	//we set our key up in a way that the lowest byte is the z in C(z,S) and the remaining 56 bits are a mask
	//where if set S = {1,2,3} bits 1, 2, and 3 are 1
	key = 0;
	key |= z;
	for(int j = 0; j < set.size(); j++)
	{
		key |= (1 << set[j]+8);
	}
}

Block travel(Block myBlock)
{
	vector<City> cityList = myBlock.minPath;
	map <long long, Salesman> shortPath;
	double **dMatrix = distanceMatrix(cityList);
	int size = (int)cityList.size();
	long long key;
	Salesman finalSalesman;
	double currentCost = 0;
	vector<int> minPath;
	double minCost = 0;

	if(size < 3)
	{
		finalSalesman.minCost = 0;
		finalSalesman.minPath.push_back(0);
		if(size == 2)
		{
			finalSalesman.minPath.push_back(1);
			finalSalesman.minCost = dMatrix[0][1];
		}
		Block returnBlock;
		vector<City> sorted = sortCities(finalSalesman.minPath, cityList);
		returnBlock.minPath = sorted;
		returnBlock.minCost = finalSalesman.minCost;
		return returnBlock;
	}
	// find C(1,{}), C(2,{}) = dMatrix[0][1] dMatrix[0][2] etc.
	// and the ones with set cardnality = 1. these are special cases for our dynamic approach
	for(int i = 1; i < size; i++)
	{
		for(int j = 1; j < size; j++) 
		{
			if(i == j)
			{
				continue;
			}
			vector<int> iSet {i};
			genKey(iSet, j, key);
			Salesman s;
			s.minCost = dMatrix[i][j] + dMatrix[0][i];
			vector<int> minPath{0,i};
			s.minPath = minPath;
			shortPath.insert(pair <long long, Salesman> (key, s));
		}	
	}

	for(int i = 2; i < size; i++)//from subset size 2 to n-1
	{
		vector<vector<int>> subsets = generateSubsets(i, size-1); //generate subsets of size i
		for(int s = 0; s < subsets.size(); s++) //iterate through every subset of size i
		{	
			vector<int> currentSet = subsets[s];
			for(int k = 0; k < currentSet.size(); k++) //in order to find C(z,S) we need to find the min C(k,{S}-{k}) for all k
			{
				vector<int> kSet;
				kSet.push_back(currentSet[k]);
				vector<int> setWithoutK;
				set_difference(currentSet.begin(), currentSet.end(), kSet.begin(), kSet.end(), inserter(setWithoutK, setWithoutK.begin()));
				minCost = DBL_MAX;
				int minM = 0;
				for(int m = 0; m < setWithoutK.size(); m++) //find min C(m, S-{k}) + dMatrix[m][k]
				{	
					vector<int> mSet;
					mSet.push_back(setWithoutK[m]);
					vector<int> setWithoutM;	
					set_difference(setWithoutK.begin(), setWithoutK.end(), mSet.begin(), mSet.end(), inserter(setWithoutM, setWithoutM.begin()));

					genKey(setWithoutM, setWithoutK[m], key);
					currentCost = shortPath[key].minCost + dMatrix[setWithoutK[m]][currentSet[k]];
					if(currentCost < minCost)
					{
						minCost = currentCost;
						minPath = shortPath[key].minPath;
						minM = m;
					}
				}
				//add the path to the min path for this route
				minPath.push_back(setWithoutK[minM]);
				// now store the answer for the next iteration 
				genKey(setWithoutK, currentSet[k], key);
				Salesman s;
				s.minCost = minCost;
				s.minPath = minPath;
				finalSalesman = s;
				shortPath.insert(pair <long long, Salesman> (key, s));
			}
		}
	}
	
	vector<int> allCities;
	for(int m = 1; m < size; m++)
	{
		allCities.push_back(m);
	}	
	minCost = DBL_MAX;
	int minM = 0;
	vector<int> mSet;
	vector<int> setWithoutM;
	
	for(int m = 0; m < allCities.size(); m++)
	{
		setWithoutM.clear();
		mSet.clear();
		mSet.push_back(allCities[m]);
		
		set_difference(allCities.begin(), allCities.end(), mSet.begin(), mSet.end(), inserter(setWithoutM, setWithoutM.begin()));
		
		genKey(setWithoutM, allCities[m], key);
		currentCost = shortPath[key].minCost + dMatrix[allCities[m]][0];
		if(currentCost < minCost)
		{
			finalSalesman.minCost = currentCost;
			finalSalesman.minPath = shortPath[key].minPath;
			minM = m;
		}
	}
	
	finalSalesman.minPath.push_back(allCities[minM]);	
	finalSalesman.minPath.push_back(0);

	//Print data
		for(int q = 0; q < finalSalesman.minPath.size(); q++)
		{
			printf("%i,",finalSalesman.minPath[q]);
		}	
		printf("\n");
		printf("MinCost: %f\n",finalSalesman.minCost);
	//	convertPath(finalSalesman, cityList, blockId);
	Block returnBlock;
	vector<City> sorted = sortCities(finalSalesman.minPath, cityList);
	returnBlock.id = myBlock.id;
	returnBlock.minPath = sorted;
	//printf("min cost is %f\n",finalSalesman.minCost);
	returnBlock.minCost = finalSalesman.minCost;
	return returnBlock;
}

bool compareCitiesy(City c1, City c2)
{
	return (c1.y > c2.y);
}

bool compareCitiesx(City c1, City c2)
{
	return (c1.x < c2.x);
}
/*
vector<vector<City>> blockBreaker(vector<City> cities){	
	vector<City> cityList = cities;	
	int numberOfBlocks = ceil(cityList.size() / pow(BLOCK_SIZE,2));
	sort(cityList.begin(), cityList.end(), compareCitiesy);
	vector<City> row;	
	vector<vector<City>> blocks;

	int counter = 0;
	unsigned int leftOver;
	vector<City> block1;
	vector<City> block2;
	for(int i = 0; i < cityList.size(); i++)
	{
		row.clear();
		for(int j = 0; j < BLOCK_SIZE*2; j++)
		{
			row.push_back(cityList[counter]);
			counter++;
			if(counter >= cityList.size())
				break;
		}
		sort(row.begin(), row.end(), compareCitiesx);

		if(row.size() > BLOCK_SIZE)
		{
			block1.insert(block1.end(),row.begin(), row.begin()+(BLOCK_SIZE));
			block2.insert(block2.end(),row.begin()+BLOCK_SIZE, row.end());
		}
		else
		{
			leftOver = (row.size() / 2);
			for(int j = 0; j < row.size(); j++)
			{
				if(j < leftOver)
				{
					block1.push_back(row[j]);
				}
				else
				{
					block2.push_back(row[j]);
				}
			}
		}
		if(counter >= cityList.size())
		{
			blocks.push_back(block1);
			blocks.push_back(block2);
			break;
		}
		if((i+1)%BLOCK_SIZE == 0 && i != 0)
		{
			if(cityList.size() - counter < 5)
			{
				for(int j = 0; j < BLOCK_SIZE; j++)
				{
					block1.pop_back();
					block2.pop_back();
					counter -= 2;
				}
			}
			blocks.push_back(block1);
			blocks.push_back(block2);
			block1.clear();
			block2.clear();
		}
	}	

	//	for(int i = 0; i < blocks.size(); i++)
	//	{
	//		printf("Block %i = \n[",i);
	//		for(int j = 0; j < block1.size(); j ++)
	//		{
	//			printf("(%.2f,%.2f)",block1[j].x,block1[j].y);
	//		}
	//		printf("\n");
	//	}	
	return blocks;
}
*/
bool compareId(Block b1, Block b2){
	return (b1.id < b2.id);
}
/*
void stitch(vector<Block> block)
{
	vector<Block> blocks = block;
	sort(blocks.begin(), blocks.end(), compareId);

	vector<City> totalPath;
	double totalCost = 0;
	double minDistance = INT_MAX;
	double currDistance;
	int minTry;
	int currBlock = 1;
	bool reverse = false;
	bool done = false;
	int startCase;
	//lets create a circle :)
	for(int i = 1; i <= blocks.size(); i++)
	{
		minDistance = INT_MAX;
		Block *block1;
		Block *block2;

		if(i == ((int)(blocks.size()/2)))
		{
			reverse = true;
			block1 = &blocks[currBlock];
			block2 = &blocks[currBlock-1];
			currBlock -= 1;
		}
		else if(reverse && currBlock == 0)
		{
			block1 = &blocks[currBlock];
			block2 = &blocks[currBlock + 1];
			done = true;
		}
		else if(reverse)
		{
			block1 = &blocks[currBlock];
			block2 = &blocks[currBlock-2];
			currBlock -= 2;
		}
		else
		{
			block1 = &blocks[currBlock];
			block2 = &blocks[currBlock+2];
			currBlock += 2;
		}

		if(block1->start)
		{
			if(block2->start)
			{
				//block1 start to block2 start
				currDistance = distance(block1->minPath[0], block2->minPath[0]);
				if(currDistance < minDistance)
				{
					minDistance = currDistance;
					minTry = 0;
				}
			}

			if(block2->end == true)
			{
				//block1 start to block2 end
				currDistance = distance(block1->minPath[0], block2->minPath.back());
				if(currDistance < minDistance)
				{
					minDistance = currDistance;
					minTry = 1;
				}
			}
		}
		if(block1->end == true)
		{
			if(block2->start == true)
			{
				//block1 end to block2 start
				currDistance = distance(block1->minPath.back(), block2->minPath[0]);
				if(currDistance < minDistance)
				{
					minDistance = currDistance;
					minTry = 2;
				}
			}
			if(block2->end == true)
			{
				
				//block1 end to block2 end
				currDistance = distance(block1->minPath.back(), block2->minPath.back());
				if(currDistance < minDistance)
				{
					minDistance = currDistance;
					minTry = 3;
				}
			}
		}

		switch(minTry)
		{
			case 0:
				if(i == 1)
				{
					startCase = 0;
				}
				if(done)
				{
					totalCost += block1->minCost + minDistance;
                                        totalPath.insert(totalPath.end(), block1->minPath.rbegin(), block1->minPath.rend());
					break;
				}
				block1->start = false;
				block2->start = false;
				totalCost += block1->minCost + minDistance;
				totalPath.insert(totalPath.end(), block1->minPath.rbegin(), block1->minPath.rend());
				break;	
			case 1:	
				if(i == 1)
				{
					startCase = 1;
				}
				if(done)
				{
					totalCost += block1->minCost + minDistance;
                                        totalPath.insert(totalPath.end(), block1->minPath.rbegin(), block1->minPath.rend());
					break;
				}
				block1->start = false;
				block2->end = false;
				totalCost += block1->minCost + minDistance;
				totalPath.insert(totalPath.end(), block1->minPath.rbegin(), block1->minPath.rend());
				break;	
			case 2:
				if(i == 1)
				{
					startCase = 2;
				}
				if(done)
				{
					totalCost += block1->minCost + minDistance;
					totalPath.insert(totalPath.end(), block1->minPath.begin(), block1->minPath.end());
					break;
				}
				block1->end = false;
				block2->start = false;
				totalCost += block1->minCost + minDistance;
				totalPath.insert(totalPath.end(), block1->minPath.begin(), block1->minPath.end());
				break;
			case 3:
				if(i == 1)
				{
					startCase = 3;
				}
				if(done)
				{
					totalCost += block1->minCost + minDistance;
					totalPath.insert(totalPath.end(), block1->minPath.begin(), block1->minPath.end());
					break;
				}
				block1->end = false;
				block2->end = false;
				totalCost += block1->minCost + minDistance;
				totalPath.insert(totalPath.end(), block1->minPath.begin(), block1->minPath.end());
				break;
		}
		if(done)
		{

			switch(startCase)
			{
				case 0:
					totalPath.push_back(block2->minPath.back());
					break;
				case 1:
					totalPath.push_back(block2->minPath.back());
					break;
				case 2:
					totalPath.push_back(block2->minPath[0]);
					break;
				case 3:
					totalPath.push_back(block2->minPath[0]);
					break;
			}	
		}
	}
	printf("total cost was: %.2f\n",totalCost);
	printf("total path:\n");
	for(int i = 0; i < totalPath.size(); i ++)
	{
		printf("%i,",totalPath[i].id);
	}
	printf("\n");
}
*/
vector<City> generateCities(int blockId, int rowSize, int colSize)
{
	time_t t;
	srand(t);
	vector<City> cities(numCities);
	City city;
	
	//Find what row and column the block should be in
	//I am creating a matrix that snakes so all blocks are guarenteed to be somewhat close
	int row = ((blockId - (blockId % rowSize)) / rowSize);
	int col;
	if((row + 1) % 2)
	{
		col = blockId % colSize;
	}
	else
	{
		col = (colSize - (blockId % colSize)) - 1;
	}

	for(int i = 0; i < numCities; i++)
	{
		city.id = i;
		city.x = fRand(col * blockSpace, (col + 1) * blockSpace);	
		city.y = fRand(row * blockSpace, (row + 1) * blockSpace);
		cities[i] = city;
	}
	return cities;
}
	
vector<Block> generateBlocks(int blocksPerProc, int blockId, vector<int> dimentions)
{
	vector<Block> blocks(blocksPerProc);
	Block block;
	for(int i = 0; i < blocksPerProc; i++)
	{
		block.id = blockId+i;
		block.minCost = 0;
		block.minPath = generateCities(block.id, dimentions[0], dimentions[1]);
		blocks[i] = block;
	}
	return blocks;
}

vector<int> cartesianProcDims()
{
	int row;
	int col;
	
	if(isSquare(numProc))
	{
		row = sqrt(numProc);
		col = sqrt(numProc);
	}
	else
	{
		int divisor = 2;
		while(numProc % divisor)
		{
			divisor++;
		}
		row = divisor;
		col = numProc / divisor;
	}
	return vector<int>({row,col});
}

int main(int argc, char* argv[])
{
	
	int blocksPerProcess;
	double gridSize;
	int initialBlock;
	vector<Block> myBlocks;
	clock_t before = clock();
	
	//initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procNum);
	MPI_Comm_size(MPI_COMM_WORLD, &numProc);
	MPI_Request send_request;
	
	const int nitems=3;
	int blocklengths[3] = {1,1,1};
	MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Datatype mpi_city_type;
	MPI_Aint offsets[3];
	offsets[0] = offsetof(City,id);
	offsets[1] = offsetof(City,x);
	offsets[2] = offsetof(City,y);

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_city_type);
	MPI_Type_commit(&mpi_city_type);
	MPI_Status status;

	if(argc == 4)
	{
		numCities = atoi(argv[1]);
		numBlocks = atoi(argv[2]);
		gridSize = atoi(argv[3]);
	}
	else
	{
		if(!procNum)
		{
			printf("usage: ./tsp_mpi citiesPerBlock numberOfBlocks, GridSize\n");
		}
		MPI_Finalize();
		exit(1);
	}

	if(procNum > numBlocks)
	{
		MPI_Finalize();
	}
	
	//how big a single block is
	blockSpace = (double)gridSize / (double)sqrt(numBlocks);
	//the dimentions of the martrix	
	vector<int> dimentions = cartesianProcDims();
	
	if(!procNum)
	{	
		vector<int> blocksPerProc(numProc, 0);
		int blocksLeft = numBlocks;
		int currentProc = 0;
		int blocksToSend;
		initialBlock = 0;


		//find how many blocks each process is going to handle
		while(blocksLeft)
		{
			blocksPerProc[currentProc] += 1;
			currentProc++;
			blocksLeft--;
			if(currentProc == numProc) 
			{
				currentProc = 0;
			}
		}
		
		//send each process how many blocks it is working on and the first id it is responsible for
		for(int i = 1; i < numProc; i++)
		{
			blocksToSend = blocksPerProc[i];
			MPI_Isend(&blocksToSend, 1, MPI_INT, i, BLOCKSTOSENDTAG, MPI_COMM_WORLD, &send_request);
			
			//the initial block index is the addition of how many blocks each previous process has already taken care of
			initialBlock += blocksPerProc[i-1];
			MPI_Isend(&initialBlock, 1, MPI_INT, i, INITIALBLOCKTAG, MPI_COMM_WORLD, &send_request);
		}
		//set proc 0 to its correct stats
		blocksPerProcess = blocksPerProc[0];
		initialBlock = 0;
	}
	else
	{
		MPI_Recv(&blocksPerProcess, 1, MPI_INT, 0, BLOCKSTOSENDTAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&initialBlock, 1, MPI_INT, 0, INITIALBLOCKTAG, MPI_COMM_WORLD, &status);
	}
	myBlocks = generateBlocks(blocksPerProcess, initialBlock, dimentions);
	
	//I could spawn threads here and run each tsp this process is responsible in parrallel but it seems out of the scope?
	for(int i = 0; i < blocksPerProcess; i++)
	{
		allBlocks.push_back(travel(myBlocks[i]));
	}
		
	//set up the communication grid
	MPI_Comm COMM_MATRIX;
	int wrap[2] = {0,0};
	MPI_Cart_create(MPI_COMM_WORLD, 2, &dimentions[0], wrap, 0, &COMM_MATRIX);
	MPI_Comm_rank(COMM_MATRIX, &numProc);
	MPI_Cart_coords(COMM_MATRIX, procNum, 2, coords);

//	clock_t after = clock() - before;
//	int msec = after * 1000 / CLOCKS_PER_SEC;
//	printf("Total time elapsed for tsp: %d seconds %d milliseconds\n",msec/1000,msec%1000);
	MPI_Finalize();
	return 0;
}
