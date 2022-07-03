#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "Basis.hpp"

void Basis::Init (std::string filename) {
    std::ifstream file(filename);
    std::string temp;
    if (file.is_open()) {
        while(std::getline(file,temp)) {
            std::istringstream ss(temp);
            filedata.push_back(ss.str());
        }
    }
    numlines = filedata.size();
}

std::vector <std::string> Basis::getAtomSpecificBasis(std::string atom) {
    int index=0;
    int startline=0;
    int endline=0;
    std::vector <std::string> atom_spec_basis;

    // Find startline
    size_t found;
    for(index=13;index<numlines;index++){
        found = filedata[index].find(atom);
        if (found != std::string::npos){
            startline=index;
            break;
        }
    }

    // Find endline
    for(index=startline;index<numlines;index++){
        found = filedata[index].find("#BASIS SET:");
        if (found != std::string::npos){
            endline=index;
            break;
        }
    }

    // Populate new vector
    for(index=startline;index<endline;index++){
        atom_spec_basis.push_back(filedata[index]);
    }

    return atom_spec_basis;
}