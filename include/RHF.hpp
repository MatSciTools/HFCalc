#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include "ElectronIntegrals.hpp"
#include "OutputHandler.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

class RHF {
public:
   long double energy = 0.0;
   Molecule mol; ElectronIntegrals ei;
   OutputHandler *out;
   int maxiter, mixtype, writemode;
   long double tol, mixparam, single_en, fock_en;
   Eigen::MatrixXd C, P, newP, F;
   Eigen::VectorXf fockvals;
   std::vector <long double> E;
   RHF(int ngauss, int iterations, long double tolerance, int mixing_type, 
   long double mixing_parameter, std::string posdata, std::string basisfile, int wm, OutputHandler &outhandle);
   int doSCF();
   std::vector <long double> getOrbitalVector(int k);
   void initializeDensity();
   void constructF();
   void normalizeC();
   long double solveSystem();
   void updateP();
   void mixP();
   void writeC();
   void calEnergy(long double focksum);
   long double getTotalEnergy();
};