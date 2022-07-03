#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include "OutputHandler.hpp"

void OutputHandler::Init(std::string ofile){
    outfile.open(ofile);
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
    writeString("Version 1.0 (July 2022)");
    writeString("###################################################################################");
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
    std::cout << "Version 1.0 (July 2022)" << std::endl;
    std::cout << "###################################################################################" << std::endl;
}

void OutputHandler::writeNewline(){
    outfile << std::endl;
}

void OutputHandler::writeString(std::string str){
    outfile << str << "\n";
}

void OutputHandler::writeStringFloat(std::string str, long double num){
    outfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << str << num << "\n";
}

void OutputHandler::writeStringInt(std::string str, int num){
    outfile << str << num << std::endl;
}

void OutputHandler::writeStringVector(std::string str, std::vector <long double> vec){
    outfile << str << "  ";
    for (int i=0;i<vec.size();i++){
        outfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << vec[i] << "  ";
    }
    outfile << std::endl;
}

void OutputHandler::closeOutputFile(){
    outfile.close();
}