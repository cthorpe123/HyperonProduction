#ifndef _SpacePointEventAssembler_cxx_
#define _SpacePointEventAssembler_cxx_

#include "SpacePointEventAssembler.h"

namespace hyperonreco { 

//////////////////////////////////////////////////////////////////////////////////////////////

void SpacePointEventAssembler::SetFile(std::string infilename,std::string treename){

  f_in = TFile::Open(infilename.c_str());
  f_in->GetObject(("spacepoints/"+treename).c_str(),t_in);
  nEvents = t_in->GetEntries();

  PFPSpacePoint_X=0;
  PFPSpacePoint_Y=0;
  PFPSpacePoint_Z=0;
  PFPHit_Channel=0;
  PFPHit_Tick=0;
  PFPHit_Width=0;

  t_in->SetBranchStatus("run",1);
  t_in->SetBranchStatus("subrun",1);
  t_in->SetBranchStatus("event",1);
  t_in->SetBranchStatus("PFPSpacePoint_X",1);
  t_in->SetBranchStatus("PFPSpacePoint_Y",1);
  t_in->SetBranchStatus("PFPSpacePoint_Z",1);
  t_in->SetBranchStatus("PFPHit_Channel",1);
  t_in->SetBranchStatus("PFPHit_Tick",1);
  t_in->SetBranchStatus("PFPHit_Width",1);

  t_in->SetBranchAddress("run",&run);
  t_in->SetBranchAddress("subrun",&subrun);
  t_in->SetBranchAddress("event",&event);
  t_in->SetBranchAddress("PFPSpacePoint_X",&PFPSpacePoint_X);
  t_in->SetBranchAddress("PFPSpacePoint_Y",&PFPSpacePoint_Y);
  t_in->SetBranchAddress("PFPSpacePoint_Z",&PFPSpacePoint_Z);
  t_in->SetBranchAddress("PFPHit_Channel",&PFPHit_Channel);
  t_in->SetBranchAddress("PFPHit_Tick",&PFPHit_Tick);
  t_in->SetBranchAddress("PFPHit_Width",&PFPHit_Width);

}

//////////////////////////////////////////////////////////////////////////////////////////////

void SpacePointEventAssembler::Close(){

  if(f_in != nullptr) f_in->Close();

}

//////////////////////////////////////////////////////////////////////////////////////////////

Event SpacePointEventAssembler::GetEvent(int i){

  t_in->GetEntry(i);

  Event e;

  e.run = run;
  e.subrun = subrun;
  e.event = event;

  e.PFPSpacePoint_X = *PFPSpacePoint_X;   
  e.PFPSpacePoint_Y = *PFPSpacePoint_Y;   
  e.PFPSpacePoint_Z = *PFPSpacePoint_Z;   
  e.PFPHit_Channel = *PFPHit_Channel;
  e.PFPHit_Tick = *PFPHit_Tick;
  e.PFPHit_Width = *PFPHit_Width;

  return e;
}

//////////////////////////////////////////////////////////////////////////////////////////////

}

#endif
