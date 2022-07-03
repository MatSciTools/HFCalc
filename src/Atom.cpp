#include <vector>
#include <string>
#include "Atom.hpp"

Atom::Atom(std::string asym, int AtomicNumber, int numElectrons, 
    int numGauss, std::vector <long double> pos, std::string bfilename){
        symbol = asym;
        Z = AtomicNumber;
        ne = numElectrons;
        ngauss = numGauss;
        position = pos;
        mol_basis.Init(bfilename);
        atom_basis = mol_basis.getAtomSpecificBasis(symbol);
        orbitals.push_back(Orbital("1s", atom_basis, ngauss));
        numAtomicOrbitals = 1;
        if (ne > 2) {
            numAtomicOrbitals = 5;
            orbitals.push_back(Orbital("2s", atom_basis, ngauss));
            orbitals.push_back(Orbital("2px", atom_basis, ngauss));
            orbitals.push_back(Orbital("2py", atom_basis, ngauss));
            orbitals.push_back(Orbital("2pz", atom_basis, ngauss));
        }
        if (ne > 10) {
            numAtomicOrbitals = 9;
            orbitals.push_back(Orbital("3s", atom_basis, ngauss));
            orbitals.push_back(Orbital("3px", atom_basis, ngauss));
            orbitals.push_back(Orbital("3py", atom_basis, ngauss));
            orbitals.push_back(Orbital("3pz", atom_basis, ngauss));
        }
        if (ne > 18) {
            numAtomicOrbitals = 13;
        }
}