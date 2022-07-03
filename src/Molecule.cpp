#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include "VectorMath.hpp"
#include "Molecule.hpp"

void Molecule::Init(int ng, std::string posfile, std::string basisfile){
    ngauss = ng;
    std::ifstream file(posfile);
    std::string temp;
    if (file.is_open()) {
        while(std::getline(file,temp)) {
            std::istringstream ss(temp);
            geodata.push_back(ss.str());
        }
    }
    natoms = geodata.size();
    bfilename = basisfile;
    readGeodata();
    calNuclearRepulsionEnergy();
    calTotalElectrons();
    calTotalAtomicOrbitals();
    calTotalMolecularOrbitals();
}

void Molecule::readGeodata(){
    size_t found1, found2;
    std::string temp,temp2; int i,Z,ne; std::string sym;
    std::vector <long double> pos;
    for (i=0;i<natoms;i++){
        found1 = geodata[i].find_first_not_of(" ");
        temp2 = geodata[i].substr(found1,geodata[i].size()-found1+1);
        found1 = temp2.find_first_of(" ");
        sym = temp2.substr(0,found1);
        temp = temp2.substr(found1+1,temp2.size()-found1 + 1);
        found2 = temp.find_first_not_of(" ");
        Z = (int) temp[found2] - (int) '0';
        temp2 = temp.substr(found2+1,temp.size()-found2+1);
        found1 = temp2.find_first_not_of(" ");
        ne = (int) temp2[found1] - (int) '0';

        temp = temp2.substr(found1+1,temp2.size()-found1+1);
        found2 = temp.find_first_not_of(" ");
        temp2 = temp.substr(found2,temp.size()-found2+1);
        found1 = temp2.find_first_of(" ");
        pos.push_back(std::stold(temp2.substr(0,found1)));
        
        temp = temp2.substr(found1,temp2.size()-found1+1);
        found2 = temp.find_first_not_of(" ");
        temp2 = temp.substr(found2,temp.size()-found2+1);
        found1 = temp2.find_first_of(" ");
        pos.push_back(std::stold(temp2.substr(0,found1)));

        temp = temp2.substr(found1,temp2.size()-found1+1);
        found2 = temp.find_first_not_of(" ");
        temp2 = temp.substr(found2,temp.size()-found2+1);
        found1 = temp2.find_first_of(" ");
        pos.push_back(std::stold(temp2.substr(0,found1)));
        
        atoms.push_back(Atom(sym, Z, ne, ngauss, pos, bfilename));
        
        pos.clear();
    }
}

void Molecule::calNuclearRepulsionEnergy(){
    int i,j; long double dist;
    nuclearEnergy=0.0;
    for(i=0;i<atoms.size();i++){
        for(j=0;j<atoms.size();j++){
            if (i != j) {
                dist = VectorMath::calDistance(atoms[i].position, atoms[j].position);
                nuclearEnergy = nuclearEnergy + 0.5*atoms[i].Z*atoms[j].Z/dist;
            }
        }
    }
}

void Molecule::calTotalElectrons(){
    int i;
    totalElectrons = 0;
    for(i=0;i<atoms.size();i++){
        totalElectrons = totalElectrons + atoms[i].ne;
    }
}

void Molecule::calTotalAtomicOrbitals(){
    int i;
    totalContractedGaussians = 0;
    for(i=0;i<atoms.size();i++){
        totalContractedGaussians = totalContractedGaussians + atoms[i].numAtomicOrbitals;
    }
}

int Molecule::calCumulativeAtomicOrbitals(int index){
    int torbitals = 0; int i;
    for (i=0;i<=index;i++){
        torbitals = torbitals + atoms[i].numAtomicOrbitals;
    }
    return torbitals;
}

void Molecule::calTotalMolecularOrbitals(){
    int i;
    numMolecularOrbitals = totalElectrons/2 + totalElectrons % 2;
    for (i=0;i<numMolecularOrbitals-1;i++){
        MolecularOrbitalOccupancy.push_back(2);
    }
    if (numMolecularOrbitals % 2 == 0){
        MolecularOrbitalOccupancy.push_back(2);
    }
    else {
        MolecularOrbitalOccupancy.push_back(1);
    }
}

std::vector <int> Molecule::getAtomOfContractedGaussian(int cgf){
    std::vector <int> atom_info; int i;
    for(i=0;i<atoms.size();i++){
        if(cgf < Molecule::calCumulativeAtomicOrbitals(i)) {
            atom_info.push_back(i);
            break;
        }
    }
    if (atom_info[0] == 0){
        atom_info.push_back(cgf);
    }
    else {
        atom_info.push_back(cgf - Molecule::calCumulativeAtomicOrbitals(atom_info[0]-1));
    }
    return atom_info;
}