#ifndef VFitter_cxx_
#define VFitter_cxx_

#include "VFitter.h"

namespace hyperon {

VFitter::VFitter(){

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
// Project a 3D line onto one of the 2D views in wire/tick space

LineWireTick VFitter::ProjectXYZWireTick(TVector3 start,TVector3 dir,int plane){

  if(plane == kPlane0) return ProjectXYZWireTick_Plane0(start,dir);
  else if(plane == kPlane1) return ProjectXYZWireTick_Plane1(start,dir);
  else if(plane == kPlane2) return ProjectXYZWireTick_Plane2(start,dir);
  else throw std::invalid_argument("VFitter: Invalid plane number");

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Project a 3D line onto the 2D view of Plane0

LineWireTick VFitter::ProjectXYZWireTick_Plane0(TVector3 start,TVector3 dir){

  LineWireTick line;
  line.Plane = kPlane0;
  line.StartChannel = U_wire(start); 
  line.StartTick = tick(start); 
  line.Gradient = (tick(start+dir)-tick(start))/(U_wire(start+dir)-U_wire(start)); 

  return line;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Project a 3D line onto the 2D view of Plane1

LineWireTick VFitter::ProjectXYZWireTick_Plane1(TVector3 start,TVector3 dir){

  LineWireTick line;
  line.Plane = kPlane1;
  line.StartChannel = V_wire(start); 
  line.StartTick = tick(start); 
  line.Gradient = (tick(start+dir)-tick(start))/(V_wire(start+dir)-V_wire(start)); 

  return line;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Project a 3D line onto the 2D view of Plane2

LineWireTick VFitter::ProjectXYZWireTick_Plane2(TVector3 start,TVector3 dir){

  LineWireTick line;
  line.Plane = kPlane2;
  line.StartChannel = Y_wire(start); 
  line.StartTick = tick(start); 
  line.Gradient = (tick(start+dir)-tick(start))/(Y_wire(start+dir)-Y_wire(start)); 

  return line;

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
// Distance between a hit and a 2D line divided by width of hit 
double VFitter::HitLineSeparation(int hit_channel,float hit_peak,float hit_width,LineWireTick line){

  // Tick line passes through at same channel as hit
  double line_tick = line.StartTick + (hit_channel - line.StartChannel)*line.Gradient;

  return (hit_peak-line_tick)/hit_width;

}

}

#endif
