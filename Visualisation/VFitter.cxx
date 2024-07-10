#ifndef VFitter_cxx_
#define VFitter_cxx_

#include "VFitter.h"

namespace hyperonreco {

VFitter::VFitter(){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::SetROI(double roi){

ROI = roi;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::LoadData(const std::vector<std::vector<double>>& channel_v,const std::vector<std::vector<double>>& tick_v,const std::vector<std::vector<double>>& width_v){
  Channel_v = channel_v;
  Tick_v = tick_v;
  Width_v = width_v;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _InsideLArSoft_
FittedV VFitter::DoFit(TVector3 start,std::map<art::Ptr<recob::SpacePoint>,art::Ptr<recob::Hit>> hitspacepointhap) const {

  FittedV result;

  std::cout << "Fitting V to shower with starting position " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;

  return result;

}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance between a hit and a 2D line divided by width of hit 

double VFitter::HitLineSeparation2(const double& channel,const double& tick,const double& width,const LineWireTick& line){

  // Tick line passes through at same channel as hit
  double line_tick = line.StartTick + (channel - line.StartChannel)*line.Gradient;

  return pow((tick-line_tick)/width,2);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance between a hit and a 2D line divided by width of hit 

#ifdef _InsideLArSoft_
double VFitter::HitLineSeparation(art::Ptr<recob::Hit> hit,LineWireTick line){

  int hit_channel = hit->Channel();
  float hit_peak = hit->PeakTime();
  float hit_width = hit->EndTick() - hit->StartTick();

  return HitLineSeparation(hit_channel,hit_peak,hit_width,line);

}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fit metric between 3D V and collection of hits 

std::pair<double,int> VFitter::FitScore(const std::vector<std::vector<double>>& channel_v,const std::vector<std::vector<double>>& tick_v,const std::vector<std::vector<double>>& width_v,FittedV fittedv){

  double score = 0.0;

  int ndof = 0;
  for(int i_pl=0;i_pl<kPlaneInvalid;i_pl++){
  
    LineWireTick line_1 = fittedv.GetArm1_2D(i_pl); 
    LineWireTick line_2 = fittedv.GetArm2_2D(i_pl); 

    for(size_t i_h=0;i_h<channel_v.at(i_pl).size();i_h++){
      double channel = channel_v.at(i_pl).at(i_h);
      double tick = tick_v.at(i_pl).at(i_h);
      double width = width_v.at(i_pl).at(i_h);
      double arm1_sep = HitLineSeparation2(channel,tick,width,line_1);
      double arm2_sep = HitLineSeparation2(channel,tick,width,line_2);
      score += std::min(arm1_sep,arm2_sep);
      ndof++;
    }
  }

  return std::make_pair(score,ndof);

} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fit metric between 3D V and collection of hits 

#ifdef _InsideLArSoft_
std::pair<double,int> VFitter::FitScore(std::vector<art::Ptr<recob::Hit>> hit_v,FittedV fittedv){

  double score = 0.0;

  for(int i_pl=0;i_pl<kInvalid;i_pl++){
    fittedv.Arm1_2D.at(i_pl) = ProjectXYZWireTick(fittedv.Vertex,fittedv.Arm1Dir,i_pl);
    fittedv.Arm2_2D.at(i_pl) = ProjectXYZWireTick(fittedv.Vertex,fittedv.Arm2Dir,i_pl);
  }

  int ndof = 0;
  for(art::Ptr<recob::Hit> hit : hit_v){
    double arm1_sep = HitLineSeparation(hit,fittedv.Arm1_2D.at(hit->View()));
    double arm2_sep = HitLineSeparation(hit,fittedv.Arm2_2D.at(hit->View()));
    score += std::min(arm1_sep,arm2_sep);
    ndof++;
  }

  return std::make_pair(score,ndof);

} 
#endif


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Do Fit, returns bool indicating if fit was successful 

bool VFitter::DoFit(FittedV& fittedv){

   std::cout << "Setting up fit" << std::endl;

   ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff ){

       //std::cout << "Evaluating fit function" << std::endl;

       FittedV v;
       v.Vertex = TVector3(coeff[0],coeff[1],coeff[2]);
       v.Arm1Theta = coeff[3];
       v.Arm1Phi = coeff[4];
       v.Arm2Theta = coeff[5];
       v.Arm2Phi = coeff[6];
/*
       std::cout << "Vertex: " << v.Vertex.X() << "  " << v.Vertex.Y() << "  " << v.Vertex.Z() << std::endl;
       std::cout << "Arm 1: " << v.GetArm1Dir().X() << "  " << v.GetArm1Dir().Y() << "  " << v.GetArm1Dir().Z() << std::endl;
       std::cout << "Arm 1: " << v.GetArm2Dir().X() << "  " << v.GetArm2Dir().Y() << "  " << v.GetArm2Dir().Z() << std::endl;
*/
       std::pair<double,int> FitVal = FitScore(Channel_v,Tick_v,Width_v,v);

       //std::cout << FitVal.first/FitVal.second << std::endl;

       return FitVal.first/FitVal.second;

       } , 7);

   std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
     ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

   fMinimizer->SetMaxFunctionCalls(10000);
   fMinimizer->SetTolerance( 0.1 );

   fMinimizer->SetVariable(0,"Vertex X",fittedv.Vertex.X(),0.1);
   fMinimizer->SetVariable(1,"Vertex Y",fittedv.Vertex.Y(),0.1);
   fMinimizer->SetVariable(2,"Vertex Z",fittedv.Vertex.Z(),0.1);
   fMinimizer->SetVariable(3,"Arm 1 Theta",fittedv.Arm1Theta,0.1);
   fMinimizer->SetVariable(4,"Arm 1 Phi",fittedv.Arm1Phi,0.1);
   fMinimizer->SetVariable(5,"Arm 2 Theta",fittedv.Arm2Theta,0.1);
   fMinimizer->SetVariable(6,"Arm 2 Phi",fittedv.Arm2Phi,0.1);

   fMinimizer->SetFunction(min);
   fMinimizer->Minimize();

   fittedv.Vertex = TVector3(fMinimizer->X()[0],fMinimizer->X()[1],fMinimizer->X()[2]);
   fittedv.Arm1Theta = fMinimizer->X()[3];
   fittedv.Arm1Phi = fMinimizer->X()[4];
   fittedv.Arm2Theta = fMinimizer->X()[5];
   fittedv.Arm2Phi = fMinimizer->X()[6];

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool VFitter::InROI(TVector3 vertex,double channel,double tick,int plane) const {

if(plane == 0) return PointHitDistanceU(vertex,channel,tick) < ROI; 
else return true;

}


}

#endif
