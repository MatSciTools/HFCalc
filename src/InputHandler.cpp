#include <iostream>
#include <fstream>
#include <sstream>
#include "InputHandler.hpp"

InputHandler::InputHandler(){
    std::ifstream file("config.in");
    std::string temp;
    if (file.is_open()) {
        while(std::getline(file,temp)) {
            std::istringstream ss(temp);
            indata.push_back(ss.str());
        }
    }
    getPositionFile();
    getBasisFolder();
    getBasisFile();
    determineBasisFile();
    getNumIterations();
    getEnergyTolerance();
    getMode();
    getMixParameter();
    getMixType();
    getRelaxationMode();
    getWriteMode();
    getOutputFileName();
}

void InputHandler::getPositionFile(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("POS");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+3,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            posfile = temp2.substr(0,found2);
            break;
        }
    }
}

void InputHandler::getBasisFile(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("BASIS");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+5,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            basischoice = temp2.substr(0,found2);
            break;
        }
    }
}

void InputHandler::getBasisFolder(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("BFOLDER");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+7,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            bfolder = temp2.substr(0,found2);
            break;
        }
    }
}

void InputHandler::getNumIterations(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("NITER");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+5,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            niter = std::stoi(temp2.substr(0,found2));
            break;
        }
    }
}

void InputHandler::getEnergyTolerance(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("ETOL");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+4,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            etol = std::stold(temp2.substr(0,found2));
            break;
        }
    }
}

void InputHandler::getMixType(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("MIXTYPE");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+7,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            mixtype = std::stoi(temp2.substr(0,found2));
            break;
        }
    }
}

void InputHandler::getMixParameter(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("MIXPARAM");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+8,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            mixparam = std::stold(temp2.substr(0,found2));
            break;
        }
    }
}

void InputHandler::getMode(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("MODE");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+4,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            caltype = std::stoi(temp2.substr(0,found2));
            break;
        }
    }
}

void InputHandler::getWriteMode(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("WRITE");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+5,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            write = std::stoi(temp2.substr(0,found2));
            break;
        }
    }    
}

void InputHandler::getRelaxationMode(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("RELAX");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+5,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            relax = std::stoi(temp2.substr(0,found2));
            break;
        }
    }
}

void InputHandler::getOutputFileName(){
    std::string temp, temp2; int found1, found2;
    for (int i=0;i<indata.size();i++){
        int startpos = indata[i].find("OUT");
        if (startpos != std::string::npos){
            temp = indata[i].substr(startpos+3,indata[i].size()-startpos+1);
            found1 = temp.find_first_not_of(" ");
            std::string temp2 = temp.substr(found1,temp.size()-found1+1);
            found2 = temp2.find_first_of("  ");
            outfile = temp2.substr(0,found2);
            break;
        }
    }   
}
void InputHandler::determineBasisFile(){
    if (basischoice == "STO2G"){
        basisfile = bfolder+"sto2g.txt";
        ngauss = 2;
    }
    else if (basischoice == "STO3G"){
        basisfile = bfolder+"sto3g.txt";
        ngauss = 3;
    }
    else if (basischoice == "STO4G"){
        basisfile = bfolder+"sto4g.txt";
        ngauss = 4;
    }
    else if (basischoice == "STO5G"){
        basisfile = bfolder+"sto5g.txt";
        ngauss = 5;
    }
    else if (basischoice == "STO6G"){
        basisfile = bfolder+"sto6g.txt";
        ngauss = 6;
    }
}