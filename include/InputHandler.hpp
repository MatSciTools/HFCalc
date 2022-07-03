#include <vector>
#include <string>

class InputHandler {
public:
    InputHandler();
    std::vector <std::string> indata;
    std::string posfile,basischoice,basisfile,bfolder;
    int caltype; int relax;
    int niter; long double etol;
    int write;
    int ngauss; int mixtype; long double mixparam;
    std::string outfile;
    void getPositionFile();
    void getBasisFile();
    void getBasisFolder();
    void getNumIterations();
    void getEnergyTolerance();
    void getMixType();
    void getMixParameter();
    void getMode();
    void getRelaxationMode();
    void getOutputFileName();
    void getWriteMode();
    void determineBasisFile();
};