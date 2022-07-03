#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include "Orbital.hpp"

Orbital::Orbital(std::string orbital_type, 
std::vector <std::string> atom_data, int ngauss){
        int i; long double n;
        otype = orbital_type;
        populateAlphaDmat(atom_data, ngauss);
        if (otype[1] == 's'){
            xi = 0; yi = 0; zi = 0;
            for(i=0;i<=alpha.size();i++){
                n = pow(2.0*alpha[i]/M_PI, 0.75);
                norm.push_back(n);
            }
        }
        else if (otype[1] == 'p'){
            if (otype[2] == 'x'){
                xi = 1; yi = 0; zi = 0;
            }
            else if (otype[2] == 'y'){
                xi = 0; yi = 1; zi = 0;
            }
            else if (otype[2] == 'z'){
                xi = 0; yi = 0; zi = 1;
            }
            else {
                std::cerr << "Incorrect orbital direction input!" << otype[2] << std::endl;
            }
            for(i=0;i<=alpha.size();i++){
                n = pow(128.0*pow(alpha[i], 5)/pow(M_PI, 3), 0.25);
                norm.push_back(n);
            }
        }
}

void Orbital::populateAlphaDmat(std::vector <std::string> atom_data, int ngauss){
    std::string temp,temp2,temp3,temp4; int i;
    size_t start,end,found,len1,second,len2;
    long double alpha_single,dmat_single;
    int onumber = (int) otype[0] - (int) '0';
    int line = (onumber-1)*(ngauss+1) + 1;
    for(i=line;i<line+ngauss;i++){
        start = atom_data[i].find_first_not_of(" ");
        end = atom_data[i].find_last_not_of(" ");
        found = atom_data[i].find(" ", start+1);
        len1 = found - start;
        alpha_single = std::stold(atom_data[i].substr(start,len1));
        alpha.push_back(alpha_single);
        temp = atom_data[i].substr(found,end-found+1);
        second = temp.find_first_not_of(" ");
        temp2 = temp.substr(second,temp.size());
        if (otype[1] == 's') {
            dmat_single = std::stold(temp2);
            dmat.push_back(dmat_single);
        }
        if (otype[1] == 'p') {
            found = temp2.find(" ",0);
            temp3 = temp2.substr(found,temp2.size());
            second = temp3.find_first_not_of(" ");
            dmat_single = std::stold(temp3.substr(second,temp3.size()));
            dmat.push_back(dmat_single);
        }
    }
}