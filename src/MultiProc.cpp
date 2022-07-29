#include<iostream>
#include "MultiProc.hpp"

void MultiProc::Init(){
    MPI_Init(NULL, NULL);
}

double MultiProc::getTime(){
    return MPI_Wtime();
}

int MultiProc::getMyRank(){
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    return myrank;
}

int MultiProc::getTotalRanks(){
    int numranks;
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    return numranks;
}

void MultiProc::End(){
    MPI_Finalize();
}
