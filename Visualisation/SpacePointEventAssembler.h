#ifndef _SpacePointEventAssembler_h_
#define _SpacePointEventAssembler_h_

#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

using std::vector;
using std::string;

namespace hyperonreco {

   struct Event {
       int run;
       int subrun;
       int event;
       std::vector<std::vector<double>> PFPSpacePoint_X;
       std::vector<std::vector<double>> PFPSpacePoint_Y;
       std::vector<std::vector<double>> PFPSpacePoint_Z;
       std::vector<std::vector<std::vector<double>>> PFPHit_Channel;
       std::vector<std::vector<std::vector<double>>> PFPHit_Tick;
       std::vector<std::vector<std::vector<double>>> PFPHit_Width;
   };

   class SpacePointEventAssembler {

      public:

         // Setters and getters
         void SetFile(std::string infilename,std::string treename);
         void Close();

         Event GetEvent(int i);
         Long64_t GetNEvents(){ return nEvents; }

      private:

         // Input file and event tree
         TFile * f_in=nullptr;
         TTree * t_in=nullptr;
         int nEvents;

         int run,subrun,event;
         std::vector<std::vector<double>> *PFPSpacePoint_X=0;
         std::vector<std::vector<double>> *PFPSpacePoint_Y=0;
         std::vector<std::vector<double>> *PFPSpacePoint_Z=0;
         std::vector<std::vector<std::vector<double>>> *PFPHit_Channel=0;
         std::vector<std::vector<std::vector<double>>> *PFPHit_Tick=0;
         std::vector<std::vector<std::vector<double>>> *PFPHit_Width=0;

   };

}

#endif
