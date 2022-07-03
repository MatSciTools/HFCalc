#include <vector>
#include <string>
#include "Orbital.hpp"
#include "Basis.hpp"

class Atom {
public:
    std::string symbol;
    int Z,ne,ngauss,numAtomicOrbitals;
    std::vector <long double> position;
    std::vector <std::string> atom_basis;
    std::vector <Orbital> orbitals;
    Basis mol_basis;
    Atom(std::string asym, int AtomicNumber, int numElectrons, 
    int numGauss, std::vector<long double> pos, std::string bfilename);
};