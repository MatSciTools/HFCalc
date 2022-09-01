#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include "OutputHandler.hpp"

void OutputHandler::Init(std::string ofile){
    outfile.open(ofile);
    if (MultiProc::getMyRank() == 0){
     writeString("###################################################################################");
     writeString("                                                                                   ");
     writeString("***       ***    ***********   **********          *          ***        **********");
     writeString("***       ***    ***********   **********         ***         ***        **********");
     writeString("***       ***    ***           **                *****        ***        **        ");
     writeString("*************    *******       **               *** ***       ***        **        ");
     writeString("*************    *******       **              ***   ***      ***        **        "); 
     writeString("***       ***    ***           **             ***********     ***        **        ");
     writeString("***       ***    ***           **********    ***       ***    *********  **********");
     writeString("***       ***    ***           **********   ***         ***   *********  **********");
     writeString("                                                                                   ");
     writeString("###################################################################################");
     writeString("Version 1.1 (September 2022)");
     writeString("###################################################################################");
     writeStringInt("Number of MPI Ranks: ", MultiProc::getTotalRanks());
     std::cout << "###################################################################################" << std::endl;
     std::cout << "                                                                                   " << std::endl;
     std::cout << "***       ***    ***********   **********          *          ***        **********" << std::endl;
     std::cout << "***       ***    ***********   **********         ***         ***        **********" << std::endl;
     std::cout << "***       ***    ***           **                *****        ***        **        " << std::endl;
     std::cout << "*************    *******       **               *** ***       ***        **        " << std::endl;
     std::cout << "*************    *******       **              ***   ***      ***        **        " << std::endl;
     std::cout << "***       ***    ***           **             ***********     ***        **        " << std::endl;
     std::cout << "***       ***    ***           **********    ***       ***    *********  **********" << std::endl;
     std::cout << "***       ***    ***           **********   ***         ***   *********  **********" << std::endl;
     std::cout << "                                                                                   " << std::endl;
     std::cout << "###################################################################################" << std::endl;
     std::cout << "Version 1.1 (September 2022)" << std::endl;
     std::cout << "###################################################################################" << std::endl;
     std::cout << "Using " << MultiProc::getTotalRanks() << " MPI Ranks" << std::endl;
    }
}

void OutputHandler::writeNewline(){
    if (MultiProc::getMyRank() == 0) {
       outfile << std::endl;
    }
}

void OutputHandler::writeString(std::string str){
    if (MultiProc::getMyRank() == 0){
       outfile << str << "\n";
    }
}

void OutputHandler::writeStringFloat(std::string str, long double num){
    if (MultiProc::getMyRank() == 0){
       outfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << str << num << "\n";
    }
}

void OutputHandler::writeStringInt(std::string str, int num){
    if (MultiProc::getMyRank() == 0){
       outfile << str << num << std::endl;
    }
}

void OutputHandler::writeStringVector(std::string str, std::vector <long double> vec){
    if (MultiProc::getMyRank() == 0){
      outfile << str << "  ";
      for (int i=0;i<vec.size();i++){
         outfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << vec[i] << "  ";
      }
      outfile << std::endl;
    }
}

void OutputHandler::closeOutputFile(){
    outfile.close();
}