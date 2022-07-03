#include <math.h>
#include <iostream>
#include <vector>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include "MatrixElement.hpp"

long double MatrixElement::E(int i, int j, int t, long double dir, long double alpha, long double beta) {
    long double p = alpha + beta;
    long double q = alpha*beta/p;
    long double KAB = exp(-q*pow(dir, 2));
    if (t < 0 || t > i + j){
        return 0.0;
    }
    else if (j > 0){
        return (0.5/p)*E(i,j-1,t-1,dir,alpha,beta) + 
        (q*dir/beta)*E(i,j-1,t,dir,alpha,beta) + 
        (t + 1)*E(i,j-1,t+1,dir,alpha,beta);
    }
    else if (i > 0){
        return (0.5/p)*E(i-1,j,t-1,dir,alpha,beta) -
        (q*dir/alpha)*E(i-1,j,t,dir,alpha,beta) +
        (t + 1)*E(i-1,j,t+1,dir,alpha,beta);
    }
    else if (i == 0 && j == 0 && t == 0) {
        return KAB;
    }
    else {
        return 0.0;
    }
}

long double MatrixElement::boys(int n, long double T){
    return boost::math::hypergeometric_1F1(n + 0.5, n+1.5, -1.0*T)/(2.0*n + 1.0);
    /*
    if (n >= 1000){
        return 0.0;
    }
    else {
        return (2.0*T*MatrixElement::boys(n+1,T) + exp(-1.0*T))/(2.0*n + 1.0);
    }
    */
}

long double MatrixElement::R(long double p, std::vector <long double> posC, 
    std::vector <long double> posP, int t, int u, int v, int n){
        long double XPC = posP[0] - posC[0];
        long double YPC = posP[1] - posC[1];
        long double ZPC = posP[2] - posC[2];
        long double RPC = VectorMath::calDistance(posP, posC);
        if (t == 0 && u == 0 && v == 0) {
            return pow(-2.0*p, n)*MatrixElement::boys(n, p*pow(RPC, 2));
        }
        else if (t > 1) {
            return XPC*MatrixElement::R(p, posC, posP, t-1, u, v, n+1) + 
            (t-1.0)*MatrixElement::R(p, posC, posP, t-2, u, v, n+1);
        }
        else if (t == 1) {
            return XPC*MatrixElement::R(p, posC, posP, t-1, u, v, n+1);
        }
        else if (u > 1) {
            return YPC*MatrixElement::R(p, posC, posP, t, u-1, v, n+1) +
            (u-1.0)*MatrixElement::R(p, posC, posP, t, u-2, v, n+1);
        }
        else if (u == 1) {
            return YPC*MatrixElement::R(p, posC, posP, t, u-1, v, n+1);
        }
        else if (v > 1){
            return ZPC*MatrixElement::R(p, posC, posP, t, u, v-1, n+1) + 
            (v-1.0)*MatrixElement::R(p, posC, posP, t, u, v-2, n+1);
        }
        else if (v == 1){
            return ZPC*MatrixElement::R(p, posC, posP, t, u, v-1, n+1);
        }
}

long double MatrixElement::calOverlapMatrixElement(int xi, int yi, int zi, int xj, int yj, int zj, 
std::vector <long double> Qpos, long double alpha, long double beta) {
    long double prefactor = pow(M_PI/(alpha + beta), 1.5);
    long double Ex = E(xi, xj, 0, Qpos[0], alpha, beta);
    long double Ey = E(yi, yj, 0, Qpos[1], alpha, beta);
    long double Ez = E(zi, zj, 0, Qpos[2], alpha, beta);
    return prefactor*Ex*Ey*Ez;
}

long double MatrixElement::calKineticMatrixElement(int xi, int yi, int zi, int xj, int yj, int zj, 
std::vector <long double> Qpos, long double alpha, long double beta){
    long double term1 = beta*(2*(xj + yj + zj) + 3)*
    MatrixElement::calOverlapMatrixElement(xi,yi,zi,xj,yj,zj,Qpos,alpha,beta);
    long double term2 = -2.0*pow(beta, 2)*
    (MatrixElement::calOverlapMatrixElement(xi,yi,zi,xj+2,yj,zj,Qpos,alpha,beta) +
    MatrixElement::calOverlapMatrixElement(xi,yi,zi,xj,yj+2,zj,Qpos,alpha,beta) +
    MatrixElement::calOverlapMatrixElement(xi,yi,zi,xj,yj,zj+2,Qpos,alpha,beta));
    long double term3 = -0.5*(xj*(xj-1)*MatrixElement::calOverlapMatrixElement(xi,yi,zi,xj-2,yj,zj,Qpos,alpha,beta)
    + yj*(yj-1)*MatrixElement::calOverlapMatrixElement(xi,yi,zi,xj,yj-2,zj,Qpos,alpha,beta) +
    + zj*(zj-1)*MatrixElement::calOverlapMatrixElement(xi,yi,zi,xj,yj,zj-2,Qpos,alpha,beta));
    return term1+term2+term3;
}

long double MatrixElement::calPotentialMatrixElement(int xi, int yi, int zi, int xj, int yj, int zj,
    std::vector <long double> posA, std::vector <long double> posB, std::vector <long double> posC,
    long double alpha, long double beta) {
        long double val = 0.0; int t,u,v;
        std::vector <long double> posQ = VectorMath::subtractVectors(posA, posB);
        std::vector <long double> posP = VectorMath::calMeanVector(posA, alpha, posB, beta);
        long double p = alpha + beta;
        for (t=0;t<=xi+xj;t++){
            for(u=0;u<=yi+yj;u++){
                for(v=0;v<=zi+zj;v++){
                    val = val + MatrixElement::E(xi, xj, t, posQ[0], alpha, beta)*
                    MatrixElement::E(yi, yj, u, posQ[1],alpha,beta)*
                    MatrixElement::E(zi, zj, v, posQ[2],alpha,beta)*MatrixElement::R(p, posC, posP, t, u, v, 0);
                }
            }
        }
        return (2.0*M_PI/(1.0*p))*val;
    }

long double MatrixElement::calElectronRepulsionElement(int xi, int yi, int zi, int xj, int yj, int zj, 
    int xk, int yk, int zk, int xl, int yl, int zl, std::vector <long double> posA, std::vector <long double> posB, 
    std::vector <long double> posC, std::vector <long double> posD, long double alpha, long double beta, 
    long double gamma, long double delta){
        long double val = 0.0; long double term1,term2,term3;
        long double p = alpha + beta;
        long double q = gamma + delta;
        std::vector <long double> posP = VectorMath::calMeanVector(posA, alpha, posB, beta);
        std::vector <long double> posQ = VectorMath::calMeanVector(posC, gamma, posD, delta);
        long double pq = p*q/(p + q);
        for(int t=0;t<=xi+xj;t++){
            for (int u=0;u<=yi+yj;u++){
                for (int v=0;v<=zi+zj;v++){
                    term1 = MatrixElement::E(xi, xj, t, posA[0]-posB[0],alpha,beta)*
                    MatrixElement::E(yi, yj, u, posA[1]-posB[1], alpha, beta)*
                    MatrixElement::E(zi, zj, v, posA[2]-posB[2],alpha,beta);
                    for (int tau=0;tau<=xk+xl;tau++){
                        for (int nu=0;nu<=yk+yl;nu++){
                            for (int phi=0;phi<=zk+zl;phi++){
                                term2 = MatrixElement::E(xk, xl, tau, posC[0]-posD[0], gamma, delta)*
                                MatrixElement::E(yk, yl, nu, posC[1]-posD[1], gamma, delta)*
                                MatrixElement::E(zk, zl, phi, posC[2]-posD[2], gamma, delta);
                                term3 = pow(-1.0, tau+nu+phi)*MatrixElement::R(pq, posQ, posP, tau+t,nu+u,phi+v, 0);
                                val = val + term1*term2*term3;
                            }
                        }
                    }
                }
            }
        }
        val = val*(2*pow(M_PI, 2.5)/(p*q*sqrt(p+q)));
        return val;
    }