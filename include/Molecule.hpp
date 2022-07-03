#include <string>
#include <vector>
#include "Atom.hpp"

class Molecule {
public:
    int ngauss, natoms; 
    int totalElectrons, numMolecularOrbitals, totalContractedGaussians;
    long double nuclearEnergy;
    std::string bfilename;
    std::vector <int> MolecularOrbitalOccupancy; 
    std::vector <std::string> geodata;
    std::vector <Atom> atoms;
    Molecule(){};
    void Init(int ng, std::string posfile, std::string basisfile);
    void readGeodata();
    void calNuclearRepulsionEnergy();
    void calTotalElectrons();
    void calTotalAtomicOrbitals();
    int calCumulativeAtomicOrbitals(int index);
    void calTotalMolecularOrbitals();
    std::vector <int> getAtomOfContractedGaussian(int cgf);
};