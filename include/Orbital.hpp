#include <vector>
#include <iostream>
#include <string>

class Orbital {
public:
    int xi;
    int yi;
    int zi;
    std::string otype;
    std::vector <long double> alpha;
    std::vector <long double> dmat;
    std::vector <long double> norm;
    Orbital(std::string orbital_type, std::vector <std::string> atom_data, int ngauss);
    void populateAlphaDmat(std::vector <std::string> atom_data, int ngauss);
};