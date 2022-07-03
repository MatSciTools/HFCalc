#include <iostream>
#include <vector>
#include "VectorMath.hpp"


class MatrixElement {
public:
    static long double E(int i, int j, int t, long double dist, 
    long double alpha, long double beta);
    static long double boys(int n, long double T);
    static long double R(long double p, std::vector <long double> posC, 
    std::vector <long double> posP, int t, int u, int v, int n);
    static long double calOverlapMatrixElement(int xi, int yi, int zi, int xj, int yj, int zj, 
    std::vector <long double> Qpos, long double alpha, long double beta);
    static long double calKineticMatrixElement(int xi, int yi, int zi, int xj, int yj, int zj,
    std::vector <long double> Qpos, long double alpha, long double beta);
    static long double calPotentialMatrixElement(int xi, int yi, int zi, int xj, int yj, int zj,
    std::vector <long double> posA, std::vector <long double> posB, std::vector <long double> posC,
    long double alpha, long double beta);
    static long double calElectronRepulsionElement(int xi, int yi, int zi, int xj, int yj, int zj, 
    int xk, int yk, int zk, int xl, int yl, int zl, std::vector <long double> posA, std::vector <long double> posB, 
    std::vector <long double> posC, std::vector <long double> posD, long double alpha, long double beta, 
    long double gamma, long double delta);
};