#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "RHF.hpp"

RHF::RHF(int ngauss, int iterations, long double tolerance, int mixing_type, 
long double mixing_parameter, std::string posdata, 
std::string basisfile, int wm, OutputHandler& outhandle){
    maxiter = iterations;
    tol = tolerance;
    writemode = wm;
    mixtype = mixing_type;
    mixparam = mixing_parameter;
    out = &outhandle;
    mol.Init(ngauss, posdata, basisfile);
    out->writeNewline();
    for (int i=0;i<mol.atoms.size();i++){
        out->writeString("Atom "+std::to_string(i)+": "+mol.atoms[i].symbol);
        out->writeStringVector("Position:", mol.atoms[i].position);
        out->writeStringInt("Number of Atomic Orbitals: ", mol.atoms[i].numAtomicOrbitals);
        out->writeNewline();
    }
    out->writeNewline();
    out->writeStringInt("Total Atomic Orbitals/Contracted Gaussians: ", mol.totalContractedGaussians);
    out->writeStringInt("Number of Molecular Orbitals: ", mol.numMolecularOrbitals);
    out->writeNewline();
    ei.Init(mol);
    C.resize(mol.numMolecularOrbitals, ei.matsize);
    P.resize(ei.matsize, ei.matsize);
    Pen.resize(ei.matsize, ei.localmatsize[MultiProc::getMyRank()]);
    newP.resize(ei.matsize, ei.matsize);
    F.resize(ei.matsize, ei.localmatsize[MultiProc::getMyRank()]);
    if (MultiProc::getMyRank() == 0){
        FFull.resize(ei.matsize, ei.matsize);
        SFull.resize(ei.matsize, ei.matsize);
    }
    fockvals.resize(mol.numMolecularOrbitals);
}

std::vector <long double> RHF::getOrbitalVector(int k){
    int i;
    std::vector <long double> Cvec;
    for (i=0;i<ei.matsize;i++){
        Cvec.push_back(C(k,i));
    }
    return Cvec;

}
void RHF::initializeDensity(){
    int i,j,start,end;
    for (i=0;i<mol.numMolecularOrbitals;i++){
        for (j=0;j<ei.matsize;j++){
            C(i, j) = 0.0;
        }
    }
    for (i=0;i<ei.matsize;i++){
        for (j=0;j<ei.matsize;j++){
            P(i, j) = 0.0;
            newP(i, j) = 0.0;
        }
    }
}

void RHF::updateP(){
    int i, j, k, start, end;
    for (i=0;i<ei.matsize;i++){
        for (j=0;j<ei.matsize;j++){
            newP(i,j) = 0.0;
            for (k=0;k<mol.numMolecularOrbitals;k++){
                newP(i, j) = newP(i, j) + 2*C(k, i)*C(k, j);
            }
        }
    }
}

void RHF::mixP(){
    int i, j, start, end;
    if (mixtype == 0){
        for (i=0;i<ei.matsize;i++){
            for (j=0;j<ei.matsize;j++){
                P(i,j) = mixparam*newP(i,j) + (1.0-mixparam)*P(i,j);
            }
        }
    }
}

void RHF::normalizeC(){
    int k,i; long double norm = 0.0;
    long double normreduce = 0.0; 
    std::vector <long double> vec;
    if (MultiProc::getMyRank() == 0) {
    for(k=0;k<mol.numMolecularOrbitals;k++){
        vec = RHF::getOrbitalVector(k);
        norm = VectorMath::calExpectationValue(vec, SFull, ei.matsize);
        for (i=0;i<ei.matsize;i++){
            C(k, i) = (1.0/sqrt(norm))*C(k, i);
        }
    }
    }
    MultiProc::sendMatrixEverywhere(C, mol.numMolecularOrbitals*ei.matsize);
}

void RHF::constructF(){
    int i, j, k, l, jl;
    for(i=0;i<ei.matsize;i++){
        for (jl=0;jl<ei.localmatsize[MultiProc::getMyRank()];jl++){
            j = ei.returnGlobalIndex(jl);
            F(i, jl) = ei.T(i, jl) + ei.V(i, jl);
            for (k=0;k<ei.matsize;k++){
                for(l=0;l<ei.matsize;l++){
                    F(i, jl) = F(i, jl) + (ei.g(i*ei.matsize+k, j*ei.matsize+l) 
                    - 0.5*ei.g(i*ei.matsize+k, l*ei.matsize+j))*P(k, l);
                }
            }
        }
    }
}

long double RHF::solveSystem(){
    long double focksum = 0.0;
    Eigen::VectorXd w(ei.matsize);
    Eigen::MatrixXd v(mol.numMolecularOrbitals, ei.matsize);
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges;
    if (MultiProc::getMyRank() == 0){
       ges.compute(FFull, SFull);
       w = ges.eigenvalues().transpose().real();
       v = ges.eigenvectors().transpose().real();
    }
    MultiProc::sendMatrixEverywhere(v, mol.numMolecularOrbitals*ei.matsize);
    MultiProc::sendVectorEverywhere(w, ei.matsize);
    for (int k=0;k<mol.numMolecularOrbitals;k++){
        for (int i=0;i<ei.matsize;i++){
            C(k, i) = v(k, i);
        }
        fockvals(k) = w(k);
        focksum = focksum + w(k);
    }
    return focksum;
}

void RHF::calEnergy(long double focksum){
    long double nuclen = mol.nuclearEnergy;
    long double total_energy = 0.0;
    long double total_energy_sum = 0.0;
    
    for(int i=0;i<ei.matsize;i++){
        for (int j=0;j<ei.localmatsize[MultiProc::getMyRank()];j++){
            total_energy = total_energy + 0.5*(ei.T(i,j) + ei.V(i,j))*Pen(i,j);
        }
    }
    MultiProc::sumOverProcesses(total_energy, total_energy_sum);
    single_en = total_energy_sum;
    fock_en = focksum;
    total_energy =  single_en + focksum + nuclen;
    E.push_back(total_energy);
}

long double RHF::getTotalEnergy(){
    return E[E.size()-1];
}

void RHF::writeC(){
    out->writeString("Contracted Gaussian coefficients");
    for(int k=0;k<mol.numMolecularOrbitals;k++){
        out->writeNewline();
        std::vector <long double> Cvec = getOrbitalVector(k);
        out->writeStringVector("MO"+std::to_string(k),Cvec);
        out->writeNewline();
    }
}

int RHF::doSCF(){
    long double error;
    out->writeNewline();
    if (mixtype == 0){
        out->writeString("Mixing Type : Simple");
    }
    else {
        out->writeString("Error: Mixing type chosen not yet implemented!");
        return 1;
    }
    out->writeStringFloat("Mixing Parameter :", mixparam);
    out->writeString("Starting SCF Loop.....");
    out->writeNewline();
    initializeDensity();
    for (int i=0;i<maxiter;i++){
        out->writeNewline();
        out->writeString("#############################");
        out->writeStringInt("Iteration ", i);
        out->writeString("#############################");
        out->writeNewline();
        constructF();
        MultiProc::collectMatrix(F, FFull, ei.matsize*ei.localmatsize[MultiProc::getMyRank()]);
        MultiProc::collectMatrix(ei.S, SFull, ei.matsize*ei.localmatsize[MultiProc::getMyRank()]);
        long double focksum = solveSystem();
        normalizeC();
        updateP();
        mixP();
        MultiProc::distributeMatrix(Pen, P, ei.matsize*ei.localmatsize[MultiProc::getMyRank()]);
        calEnergy(focksum);
        if (writemode == 1){
            writeC();
        }
        out->writeStringFloat("Energy : ", E[E.size()-1]);
        if (E.size() >= 2){
            error = abs(E[E.size()-1] - E[E.size()-2]);
            out->writeStringFloat("Error : ", error);
            if (error < tol){
                out->writeNewline();
                out->writeString("#######################");
                out->writeString("SCF CONVERGED!");
                out->writeString("#######################");
                out->writeNewline();
                break;
            }
        }
    }
    return 0;
}