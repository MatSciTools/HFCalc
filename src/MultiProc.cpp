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

void MultiProc::synchronize(){
    MPI_Barrier(MPI_COMM_WORLD);
}

void MultiProc::collectMatrix(Eigen::MatrixXd &m, Eigen::MatrixXd &mbase, int sendsize){
    MPI_Gather(m.data(), sendsize, MPI_DOUBLE, mbase.data(), sendsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MultiProc::distributeMatrix(Eigen::MatrixXd &m, Eigen::MatrixXd &mbase, int sendsize){
    MPI_Scatter(mbase.data(), sendsize, MPI_DOUBLE, m.data(), sendsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MultiProc::sendVectorEverywhere(Eigen::VectorXd &vec, int sendsize){
    MPI_Bcast(vec.data(), sendsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
void MultiProc::sendMatrixEverywhere(Eigen::MatrixXd &m, int sendsize){
    MPI_Bcast(m.data(), sendsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MultiProc::sumOverProcesses(long double &val, long double &temp){
    MPI_Allreduce(&val, &temp, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
void MultiProc::abortMultiProc(){
    MPI_Abort(MPI_COMM_WORLD,0);
}
void MultiProc::End(){
    MPI_Finalize();
}
