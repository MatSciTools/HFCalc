#include<iostream>
#include<mpi.h>

class MultiProc {
public:
     MultiProc(){};
     static void Init();
     static int getMyRank();
     static int getTotalRanks();
     static double getTime();
     static void End();
};
