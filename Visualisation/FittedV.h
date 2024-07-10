#ifndef _FittedV_h_
#define _FittedV_h_

// C++ STL includes
#include <vector>
#include <stdexcept>

#ifdef _InsideLArSoft_
// larsoft includes
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#endif

// root includes
#include "TVector3.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

// Local includes
#ifdef _InsideLArSoft_ 
#include "ubana/HyperonProduction/Alg/Position_To_Wire.h"
#else 
#include "Position_To_Wire.h" 
#endif

namespace hyperonreco {

  enum e_Planes {kPlane0,kPlane1,kPlane2,kPlaneInvalid};

  // Make a unit TVector3 from theta/phi values
  inline TVector3 MakeTVector3(double theta,double phi){
    return TVector3(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  }

  // Line in Wire tick space
  struct LineWireTick {
    int StartChannel;
    int StartTick;
    double Gradient; // Gradient in ticks/channel
    int Plane;
  };

  inline LineWireTick ProjectXYZWireTick_Plane0(TVector3 start,TVector3 dir){
    LineWireTick line;
    line.Plane = kPlane0;
    line.StartChannel = U_wire(start); 
    line.StartTick = tick(start); 
    line.Gradient = (tick(start+dir)-tick(start))/(U_wire(start+dir)-U_wire(start)); 
    return line;
  }

  inline LineWireTick ProjectXYZWireTick_Plane1(TVector3 start,TVector3 dir){
    LineWireTick line;
    line.Plane = kPlane1;
    line.StartChannel = V_wire(start); 
    line.StartTick = tick(start); 
    line.Gradient = (tick(start+dir)-tick(start))/(V_wire(start+dir)-V_wire(start)); 
    return line;
  }

  inline LineWireTick ProjectXYZWireTick_Plane2(TVector3 start,TVector3 dir){
    LineWireTick line;
    line.Plane = kPlane2;
    line.StartChannel = Y_wire(start); 
    line.StartTick = tick(start); 
    line.Gradient = (tick(start+dir*1000)-tick(start))/(Y_wire(start+dir*1000)-Y_wire(start)); 
    return line;
  }

  inline LineWireTick ProjectXYZWireTick(TVector3 start,TVector3 dir,int plane){
    if(plane == kPlane0) return ProjectXYZWireTick_Plane0(start,dir);
    else if(plane == kPlane1) return ProjectXYZWireTick_Plane1(start,dir);
    else if(plane == kPlane2) return ProjectXYZWireTick_Plane2(start,dir);
    else throw std::invalid_argument("VFitter: Invalid plane number");
  }


  // 3D fitted V object
  class FittedV {

    public: 

      TVector3 Vertex;
      //double Arm1Len,Arm2Len;     
      double OpeningAngle;
      double Arm1Theta,Arm2Theta;
      double Arm1Phi,Arm2Phi;

      LineWireTick GetArm1_2D(int i_pl) const;
      LineWireTick GetArm2_2D(int i_pl) const;

      TVector3 GetArm1Dir() const;
      TVector3 GetArm2Dir() const;

      double Chi2;
      int NDof;

      #ifdef _InsideLArSoft_
      std::vector<art::Ptr<recob::Hit>> Hits;
      std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>> HitSpacePointMap; 
      #endif

  };


}
#endif
