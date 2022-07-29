#include <iostream>
#include <vector>
#include <string>
#include "MatrixElement.hpp"
#include "Molecule.hpp"
#include "MultiProc.hpp"
#include <Eigen/Dense>

class ElectronIntegrals{
public:
    int matsize;
    Eigen::MatrixXd S, T, V, g;
    Molecule mol;
    ElectronIntegrals(){};
    void Init(Molecule m);
    void calSingleElectronMatrices();
    void calElectronRepulsionMatrix();
};
