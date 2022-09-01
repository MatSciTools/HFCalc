#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include "ElectronIntegrals.hpp"

void ElectronIntegrals::Init(Molecule m){
    mol = m;
    matsize = mol.totalContractedGaussians;
    divideMatrixAmongRanks();
    S.resize(matsize, localmatsize[MultiProc::getMyRank()]);
    T.resize(matsize, localmatsize[MultiProc::getMyRank()]);
    V.resize(matsize, localmatsize[MultiProc::getMyRank()]);
    g.resize(matsize*matsize, matsize*matsize);
    calSingleElectronMatrices();
    calElectronRepulsionMatrix();
    MultiProc::synchronize();

}

void ElectronIntegrals::divideMatrixAmongRanks(){
    int lmat = matsize/MultiProc::getTotalRanks();
    for(int j=0;j<MultiProc::getTotalRanks();j++){
        localmatsize.push_back(lmat);
        if (matsize % MultiProc::getTotalRanks() != 0){
            if (j < matsize % MultiProc::getTotalRanks()){
                localmatsize[j] = localmatsize[j] + 1;
            }
        }
    }
}

int ElectronIntegrals::returnGlobalIndex(int j){
    if (MultiProc::getMyRank() == 0){
        return j;
    }
    else {
        int matcumulative = 0;
        for (int k=0;k<MultiProc::getMyRank();k++){
            matcumulative = matcumulative + localmatsize[k];
        }
        return matcumulative + j;
    }
}

void ElectronIntegrals::calSingleElectronMatrices(){
    int i; int j,jl; std::vector <int> A, B;
    std::vector <long double> posA, posB, posC, alphaA, alphaB, dmatA, dmatB;
    int xi, yi, zi, xj, yj, zj; int m, n, l; int start,end;
    std::vector <long double> normA, normB, Q;
    long double element, element2, element3, Zl;
    for (i = 0; i < matsize; i++){
        for (jl = 0; jl < localmatsize[MultiProc::getMyRank()]; jl++){
            j = returnGlobalIndex(jl);
            A = mol.getAtomOfContractedGaussian(i);
            posA = mol.atoms[A[0]].position;
            alphaA = mol.atoms[A[0]].orbitals[A[1]].alpha;
            dmatA = mol.atoms[A[0]].orbitals[A[1]].dmat;
            normA = mol.atoms[A[0]].orbitals[A[1]].norm;
            xi = mol.atoms[A[0]].orbitals[A[1]].xi;
            yi = mol.atoms[A[0]].orbitals[A[1]].yi;
            zi = mol.atoms[A[0]].orbitals[A[1]].zi;

            B = mol.getAtomOfContractedGaussian(j);
            posB = mol.atoms[B[0]].position;
            alphaB = mol.atoms[B[0]].orbitals[B[1]].alpha;
            dmatB = mol.atoms[B[0]].orbitals[B[1]].dmat;
            normB = mol.atoms[B[0]].orbitals[B[1]].norm;
            xj = mol.atoms[B[0]].orbitals[B[1]].xi;
            yj = mol.atoms[B[0]].orbitals[B[1]].yi;
            zj = mol.atoms[B[0]].orbitals[B[1]].zi;

            Q = VectorMath::subtractVectors(posA, posB); 
            element = 0.0; element2 = 0.0; element3 = 0.0;
            for (m=0;m<mol.ngauss;m++){
                for(n=0;n<mol.ngauss;n++){
                    element = element + dmatA[m]*dmatB[n]*normA[m]*normB[n]*
                    MatrixElement::calOverlapMatrixElement(xi,yi,zi,xj,yj,zj,Q,alphaA[m], alphaB[n]);
                    element2 = element2 + dmatA[m]*dmatB[n]*normA[m]*normB[n]*
                    MatrixElement::calKineticMatrixElement(xi,yi,zi,xj,yj,zj,Q,alphaA[m], alphaB[n]);
                    for (l=0;l<mol.atoms.size();l++){
                            Zl = mol.atoms[l].Z;
                        posC = mol.atoms[l].position;
                        element3 = element3 - Zl*dmatA[m]*dmatB[n]*normA[m]*normB[n]*
                          MatrixElement::calPotentialMatrixElement(xi, yi, zi, xj, yj, zj, posA, posB, posC, alphaA[m], alphaB[n]);
                        }
                }
            }

            S(i, jl) = element;
            T(i, jl) = element2;
            V(i, jl) = element3;
        }
    }
}

void ElectronIntegrals::calElectronRepulsionMatrix(){
    int i, j, k ,l; std::vector <int> A, B, C, D;
    std::vector <long double> posA, posB, posC, posD, alphaA, alphaB, alphaC, alphaD;
    std::vector <long double> dmatA, dmatB, dmatC, dmatD, normA, normB, normC, normD;
    int xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl;
    int m, n, o, p, start, end;
    long double element = 0.0;
    for (i=0;i<matsize;i++){
        for(j=0;j<matsize;j++){
            for(k=0;k<matsize;k++){
                for (l=0;l<matsize;l++){
                A = mol.getAtomOfContractedGaussian(i);
                posA = mol.atoms[A[0]].position;
                alphaA = mol.atoms[A[0]].orbitals[A[1]].alpha;
                dmatA = mol.atoms[A[0]].orbitals[A[1]].dmat;
                normA = mol.atoms[A[0]].orbitals[A[1]].norm;
                xi = mol.atoms[A[0]].orbitals[A[1]].xi;
                yi = mol.atoms[A[0]].orbitals[A[1]].yi;
                zi = mol.atoms[A[0]].orbitals[A[1]].zi;

                B = mol.getAtomOfContractedGaussian(j);
                posB = mol.atoms[B[0]].position;
                alphaB = mol.atoms[B[0]].orbitals[B[1]].alpha;
                dmatB = mol.atoms[B[0]].orbitals[B[1]].dmat;
                normB = mol.atoms[B[0]].orbitals[B[1]].norm;
                xj = mol.atoms[B[0]].orbitals[B[1]].xi;
                yj = mol.atoms[B[0]].orbitals[B[1]].yi;
                zj = mol.atoms[B[0]].orbitals[B[1]].zi;

                C = mol.getAtomOfContractedGaussian(k);
                posC = mol.atoms[C[0]].position;
                alphaC = mol.atoms[C[0]].orbitals[C[1]].alpha;
                dmatC = mol.atoms[C[0]].orbitals[C[1]].dmat;
                normC = mol.atoms[C[0]].orbitals[C[1]].norm;
                xk = mol.atoms[C[0]].orbitals[C[1]].xi;
                yk = mol.atoms[C[0]].orbitals[C[1]].yi;
                zk = mol.atoms[C[0]].orbitals[C[1]].zi;

                D = mol.getAtomOfContractedGaussian(l);
                posD = mol.atoms[D[0]].position;
                alphaD = mol.atoms[D[0]].orbitals[D[1]].alpha;
                dmatD = mol.atoms[D[0]].orbitals[D[1]].dmat;
                normD = mol.atoms[D[0]].orbitals[D[1]].norm;
                xl = mol.atoms[D[0]].orbitals[D[1]].xi;
                yl = mol.atoms[D[0]].orbitals[D[1]].yi;
                zl = mol.atoms[D[0]].orbitals[D[1]].zi;

                element = 0.0;
                for (m=0;m<mol.ngauss;m++){
                    for (n=0;n<mol.ngauss;n++){
                        for (o=0;o<mol.ngauss;o++){
                            for(p=0;p<mol.ngauss;p++){
                                element = element + dmatA[m]*dmatB[n]*dmatC[o]*dmatD[p]*
                                normA[m]*normB[n]*normC[o]*normD[p]*MatrixElement::calElectronRepulsionElement(xi, yi, zi, 
                                xj, yj, zj, xk, yk, zk, xl, yl, zl, posA, posB, posC, posD, alphaA[m],alphaB[n], alphaC[o],alphaD[p]);
                            }
                        }
                    }
                }
                g(i*matsize + k, j*matsize+l) = element;
                }
            }
        }
    }
}