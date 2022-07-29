#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include "RHF.hpp"
#include "InputHandler.hpp"

int main(){
  // Initializing timekeeping
    std::chrono::time_point<std::chrono::system_clock> tstart;
    tstart = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(tstart);
    std::cout << "Execution started at " << std::ctime(&start_time);
  
  // Initializing MPI
    MultiProc::Init();
    double wtstart = MultiProc::getTime();

  // Handling I/O
    InputHandler in;
    OutputHandler out;
    out.Init(in.outfile);
  
  // Determining mode of calculation
    out.writeNewline();
    out.writeString("-------------------------------------------");
    if (in.caltype == 0){
      out.writeString("Calculation Mode: Restricted Hartree-Fock  ");
      if (in.relax == 1){
        out.writeString("Structural Relaxation not yet implemented");
        out.writeString("Execution failed!");
        return 0;
      }
      out.writeString("-------------------------------------------");
      out.writeNewline();
      RHF rf(in.ngauss, in.niter, in.etol, in.mixtype, in.mixparam, in.posfile, in.basisfile, in.write, out);
      int status = rf.doSCF();
      if (status == 1){
        std::cout << "Execution error, see output file: " << in.outfile << std::endl;
        out.writeString("Execution error!");
        return 0;
      }
      // Print and write final results
      out.writeStringFloat("Total Single-Electron Energy (in Hartrees): ", rf.single_en);
      out.writeStringFloat("Total Interaction Matrix Energy (in Hartrees): ", rf.fock_en);
      out.writeStringFloat("Nuclear Energy (in  Hartrees): ", rf.mol.nuclearEnergy);
      out.writeStringFloat("Total Converged Energy (in Hartrees): ", rf.getTotalEnergy());
      out.writeString("#################################################");
      out.writeNewline();
      std::cout << "Total Converged Energy : " << rf.getTotalEnergy() << " Hartrees" << std::endl;
      std::cout << "Data written to : " << in.outfile << std::endl;
    }
    else {
      out.writeString("Error: Only RHF implemented in Version 1.0! Stay tuned for more.");
      return 0;
    }
  
  // Calculating walltime elapsed and ending MPI
    double wtend = MultiProc::getTime();
    MultiProc::End();

  // Printing elapsed time and end I/O
    out.writeString("Execution Time : "+std::to_string(wtend-wtstart)+"s");
    out.closeOutputFile();
}
