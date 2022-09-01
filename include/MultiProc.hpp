#include<iostream>
#include<Eigen/Dense>
#include<mpi.h>

class MultiProc {
public:
     MultiProc(){};
     static void Init();
     static int getMyRank();
     static int getTotalRanks();
     static double getTime();
     static void synchronize();
     static void abortMultiProc();
     static void collectMatrix(Eigen::MatrixXd &m, Eigen::MatrixXd &mbase, int sendsize);
     static void distributeMatrix(Eigen::MatrixXd &m, Eigen::MatrixXd &mbase, int sendsize);
     static void sendMatrixEverywhere(Eigen::MatrixXd &m, int sendsize);
     static void sendVectorEverywhere(Eigen::VectorXd &vec, int sendsize);
     static void sumOverProcesses(long double &val);
     static void End();
};
