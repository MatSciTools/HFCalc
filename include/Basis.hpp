#include <iostream>
#include <string>
#include <vector>

class Basis {
    std::vector <std::string> filedata;
    int numlines;
public:
    Basis(){};
    void Init(std::string filename);
    std::vector <std::string> getAtomSpecificBasis(std::string atom);
};