#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "MultiProc.hpp"

class OutputHandler {
public:
     std::ofstream outfile;
     OutputHandler(){};
     void Init(std::string ofile);
     void writeNewline();
     void writeString(std::string str);
     void writeStringFloat(std::string str, long double num);
     void writeStringInt(std::string str, int num);
     void writeStringVector(std::string str, std::vector <long double> vec);
     void closeOutputFile();
};