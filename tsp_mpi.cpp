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
#include <cmath>
#include <stddef.h>
#include <utility>
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
MPI_Datatype mpi_city_type;

inline double distance(City city1, City city2) { return sqrt((pow((city2.x - city1.x),2)) + (pow((city2.y - city1.y),2))); }
inline bool isSquare(int x) { return (sqrt(x) - floor(sqrt(x)) == 0); }
inline double swapFTS(pair<City, City> left, pair<City, City> right)
{
	return (distance(left.first, right.second) + distance(left.second, right.first) - distance(left.first, left.second) - distance(right.first, right.second));
}
inline double swapFTF(pair<City, City> left, pair<City, City> right)
{
	return (distance(left.first, right.first) + distance(left.second, right.second) - distance(left.first, left.second) - distance(right.first, right.second));
}

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

Block convertPath(Salesman salesman, vector<City> cities, int blockId)
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
	return block;
}

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

	return convertPath(finalSalesman, cityList, myBlock.id);
}

bool compareCitiesy(City c1, City c2)
{
	return (c1.y > c2.y);
}

bool compareCitiesx(City c1, City c2)
{
	return (c1.x < c2.x);
}

Block merge(Block block1, Block block2)
{
	double swapCostFTS;
	double swapCostFTF;
	double swapCost;
	double minCost = DBL_MAX;
	bool FTFswap;
	bool minSwap;
	vector<City> left = block1.minPath;
	vector<City> right = block2.minPath;
	pair<City, City> p1;
	pair<City, City> p2;
	pair<pair<City, City>, pair<City, City>> minPair;
	Block returnBlock;

	//this function alone makes me love c++	
	for(int i = 0; i < left.size(); i++)
	{
		p1 = make_pair(left[0], left[1]);
		for(int j = 0; j < right.size(); j++)
		{
			p2 = make_pair(right[0], right[1]);
			swapCostFTS = swapFTS(p1,p2);
			swapCostFTF = swapFTF(p1,p2);
			
			if(swapCostFTF < swapCostFTS)
			{
				swapCost = swapCostFTF;
				FTFswap = true;	
			}
			else
			{
				swapCost = swapCostFTS;
				FTFswap = false;
			}
				
			if(swapCost < minCost)
			{
				minCost = swapCost;
				minPair = make_pair(p1,p2);
				minSwap = FTFswap;
			}
			rotate(right.begin(), right.begin() + 1, right.end());
		}
		rotate(left.begin(), left.begin() + 1, left.end());
	}

	//now that you have found the min edge swap connect them together
	returnBlock.id = block1.id;
	returnBlock.minCost = block1.minCost + block2.minCost + minCost; 
	
	
//	printf("the best swap is from city1 %i to city2 %i and city1 %i to city2 %i\n",minPair.first.first.id, minPair.second.first.id, minPair.first.second.id, minPair.second.second.id);
	//disconnect block2 circle
	block2.minPath.pop_back();
	
	//rotate until block2's new edge is at the end
	while(minPair.second.first.id != block2.minPath[0].id)
	{
		rotate(block2.minPath.begin(), block2.minPath.begin() + 1, block2.minPath.end());
	}
	rotate(block2.minPath.begin(), block2.minPath.begin() + 1, block2.minPath.end());
	
	//insert until you hit the new edge in block1
	int counter = -1;
	do
	{
		counter++;
		returnBlock.minPath.push_back(block1.minPath[counter]);
	}while(minPair.first.first.id != block1.minPath[counter].id);
	
	//insert block2
	if(minSwap)
	{
		for(int i = block2.minPath.size()-1; i >= 0; i--)
		{
			returnBlock.minPath.push_back(block2.minPath[i]);
		}
	}
	else
	{
		for(int i = 0; i < block2.minPath.size(); i++)
		{
			returnBlock.minPath.push_back(block2.minPath[i]);
		}
	}

	//insert the rest of block1
	for(int i = counter+1; i < block1.minPath.size(); i++)
	{
		returnBlock.minPath.push_back(block1.minPath[i]);
	}

	return returnBlock;
}

//Gotta give credit where credit is due
// https://gist.github.com/rmcgibbo/7178576
// rmcgibbo made this beautiful custom reduce that uses only sends and recvs
// I customized it to work for me of course
Block MPI_ManualReduce(MPI_Comm comm)
{
	int recbuffer;
	Block mergeBlock;
	int tag = 0;
	const int lastpower = 1 << ((int)log2(numProc));
	double minCost;
	int size;
	vector<City> cityList;
	MPI_Status status;

	//for all the ones greator than a power of 2
	//do procs/2 communications
	for(int i = lastpower; i < numProc; i++)
	{
		if(procNum == i)
		{
			size = allBlocks[0].minPath.size();
			minCost = allBlocks[0].minCost;
			City *cityArr = (City *)malloc(size * sizeof(City));
			copy(allBlocks[0].minPath.begin(), allBlocks[0].minPath.end(), cityArr);
			MPI_Send(&size, 1, MPI_INT, i - lastpower, tag, comm);
			MPI_Send(cityArr, size, mpi_city_type, i - lastpower, tag, comm);
			MPI_Send(&minCost, 1, MPI_DOUBLE, i - lastpower, tag, comm);
		}
	}
	for(int i = 0; i < numProc - lastpower; i++)
	{
		if(procNum == i)
		{
			MPI_Recv(&size, 1, MPI_INT, i + lastpower, tag, comm, &status);
	
			City *cityArr2 = (City *)malloc(size * sizeof(City));
			MPI_Recv(cityArr2, size, mpi_city_type, i + lastpower, tag, comm, &status);
			MPI_Recv(&minCost, 1, MPI_DOUBLE, i + lastpower, tag, comm, &status);
			
			cityList.clear();
			for(int j = 0; j < size; j++)
			{
				cityList.push_back(cityArr2[j]);
			}
			mergeBlock.id = 0;
			mergeBlock.minPath = cityList;
			mergeBlock.minCost = minCost;
			allBlocks[0] = merge(allBlocks[0], mergeBlock);	
		}
	}
	///////////
	
	//for the rest = power of 2
	//do perfect logp reduction
	for(int i = 0; i < ((int)log2(lastpower)); i++)
	{
		for(int j = 0; j < lastpower; j += 1 << (i+1))
		{
			const int receiver = j;
			const int sender = j + (1 << i);
			if(procNum == receiver)
			{
				MPI_Recv(&size, 1, MPI_INT, sender, tag, comm, &status);
		
				City *cityArr3 = (City *)malloc(size * sizeof(City));
				MPI_Recv(cityArr3, size, mpi_city_type, sender, tag, comm, &status);
				MPI_Recv(&minCost, 1, MPI_DOUBLE, sender, tag, comm, &status);
				
				cityList.clear();
				for(int j = 0; j < size; j++)
				{
					cityList.push_back(cityArr3[j]);
				}
				mergeBlock.id = 0;
				mergeBlock.minPath = cityList;
				mergeBlock.minCost = minCost;
				allBlocks[0] = merge(allBlocks[0], mergeBlock);	
			}
			else if(procNum == sender)
			{
				size = allBlocks[0].minPath.size();
				minCost = allBlocks[0].minCost;
				City *cityArr4 = (City *)malloc(size * sizeof(City));
				copy(allBlocks[0].minPath.begin(), allBlocks[0].minPath.end(), cityArr4);
				MPI_Send(&size, 1, MPI_INT, receiver, tag, comm);
				MPI_Send(cityArr4, size, mpi_city_type,receiver, tag, comm);
				MPI_Send(&minCost, 1, MPI_DOUBLE, receiver, tag, comm);
			}
		}
	}
	return allBlocks[0];
}	

vector<City> generateCities(int blockId, int rowSize, int colSize)
{
	time_t t;
	srand(t);
	vector<City> cities(numCities);
	City city;
	string str;
	
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
		str = to_string(i+1);
		str += to_string(blockId);
		str += to_string(procNum);
		city.id = atoi(str.c_str());
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

void printData()
{
	printf("\n");
	for(int i = 0; i < allBlocks.size(); i++)
	{
		printf("MinCost: %.2f \n",allBlocks[i].minCost);
	//	printf("MinPath:\n");
	//	for(int j = 0; j < allBlocks[i].minPath.size(); j++)
	//	{
	//		printf("%i,",allBlocks[i].minPath[j].id);
	//	}
	//	printf("\n");
	}
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
	
	MPI_Comm_size(MPI_COMM_WORLD, &numProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &procNum);	
	//the dimentions of the martrix	
	vector<int> dimentions = cartesianProcDims();
	//set up the communication grid
	MPI_Comm COMM_MATRIX;
	int wrap[2] = {1,1};
	MPI_Cart_create(MPI_COMM_WORLD, 2, &dimentions[0], wrap, 0, &COMM_MATRIX);
	MPI_Comm_rank(COMM_MATRIX, &procNum);
	MPI_Comm_size(COMM_MATRIX, &numProc);
	MPI_Cart_coords(COMM_MATRIX, procNum, 2, coords);
	MPI_Request send_request;
	
	const int nitems=3;
	int blocklengths[3] = {1,1,1};
	MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
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
		
	//how big a single block is
	blockSpace = (double)gridSize / (double)sqrt(numBlocks);
	
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
			MPI_Isend(&blocksToSend, 1, MPI_INT, i, BLOCKSTOSENDTAG, COMM_MATRIX, &send_request);
			
			//the initial block index is the addition of how many blocks each previous process has already taken care of
			initialBlock += blocksPerProc[i-1];
			MPI_Isend(&initialBlock, 1, MPI_INT, i, INITIALBLOCKTAG, COMM_MATRIX, &send_request);
		}
		//set proc 0 to its correct stats
		blocksPerProcess = blocksPerProc[0];
		initialBlock = 0;
	}
	else
	{
		MPI_Recv(&blocksPerProcess, 1, MPI_INT, 0, BLOCKSTOSENDTAG, COMM_MATRIX, &status);
		MPI_Recv(&initialBlock, 1, MPI_INT, 0, INITIALBLOCKTAG, COMM_MATRIX, &status);
	}
	myBlocks = generateBlocks(blocksPerProcess, initialBlock, dimentions);

	//I could spawn threads here and run each tsp this process is responsible in parrallel but it seems out of the scope?
	for(int i = 0; i < blocksPerProcess; i++)
	{
		allBlocks.push_back(travel(myBlocks[i]));
	}

	while(allBlocks.size() > 1)
	{
		allBlocks[0] = merge(allBlocks[0],allBlocks[1]);
		allBlocks.erase(allBlocks.begin()+1);
	}
			
	allBlocks[0] = MPI_ManualReduce(COMM_MATRIX);

	if(!procNum)
	{
		printData();
	}
	
/*
	if(!procNum)
	{
		vector<City> finalPath = allBlocks[0].minPath;
		int pathSize = finalPath.size();
		double pathCost = 0;
		City previous = finalPath[0];
		for(int i = 1; i < pathSize; i++)
		{
			pathCost += distance(previous, finalPath[i]);
			previous = finalPath[i];
		}
		printf("Actual cost = %.2f\n",pathCost);
	}
*/
	
	if(!procNum)
	{
		clock_t after = clock() - before;
		int msec = after * 1000 / CLOCKS_PER_SEC;
		printf("Total time elapsed for tsp: %d seconds %d milliseconds\n",msec/1000,msec%1000);
	}

	MPI_Finalize();
	return(0);
}
