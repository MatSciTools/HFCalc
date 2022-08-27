#include <iostream>
#include <vector>
#include <string>
#include "MatrixElement.hpp"
#include "Molecule.hpp"
#include "OutputHandler.hpp"
#include <Eigen/Dense>

class ElectronIntegrals{
public:
    int matsize; 
    std::vector <int> localmatsize;
    Eigen::MatrixXd S, T, V, g;
    Molecule mol;
    ElectronIntegrals(){};
    void Init(Molecule m);
    void divideMatrixAmongRanks();
    void calSingleElectronMatrices();
    void calElectronRepulsionMatrix();
};
