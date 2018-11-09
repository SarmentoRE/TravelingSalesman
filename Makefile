

tsp: tsp_seq.cpp 
	g++ -Wall -g -std=c++11 $^ -o $@
tsp_t: tsp_thread.cpp
	g++ -Wall -pthread -g -std=c++11 $^ -o $@
tsp_mpi: tsp_mpi.cpp
	mpic++ -Wall -g -std=c++11 $^ -o $@
