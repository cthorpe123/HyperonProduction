#ifndef _VFitter_h_
#define _VFitter_h_

// This object is compilable outside of LArSoft, comment out this line to do so
#define _InsideLArSoft_

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

// Local includes
#include "ubana/HyperonProduction/Alg/Position_To_Wire.h"
//#include "Position_To_Wire.h"

namespace hyperon {

  // 3D fitted V object
  struct FittedV {
    TVector3 Vertex;
    TVector3 Arm1Dir,Arm2Dir;
    double Arm1Len,Arm2Len;     
    double OpeningAngle;
    double Chi2;
    int NDof;
    #ifdef _InsideLArSoft_
    std::vector<art::Ptr<recob::Hit>> Hits;
    std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>> HitSpacePointMap; 
    #endif
  };

  // Line in Wire tick space
  struct LineWireTick {
    int StartChannel;
    int StartTick;
    double Gradient; // Gradient in ticks/channel
    int Plane;
  };

  class VFitter {

    public:

      VFitter();

      #ifdef _InsideLArSoft_
      FittedV DoFit(TVector3 start,std::map<art::Ptr<recob::SpacePoint>,art::Ptr<recob::Hit>> hitspacepointhap) const;
      #endif

    private:

      enum e_Planes {kPlane0,kPlane1,kPlane2,kPlaneInvalid};

      LineWireTick ProjectXYZWireTick(TVector3 start,TVector3 dir,int plane);
      LineWireTick ProjectXYZWireTick_Plane0(TVector3 start,TVector3 dir);
      LineWireTick ProjectXYZWireTick_Plane1(TVector3 start,TVector3 dir);
      LineWireTick ProjectXYZWireTick_Plane2(TVector3 start,TVector3 dir);

      #ifdef _InsideLArSoft_
      double HitLineSeparation(art::Ptr<recob::Hit> hit,LineWireTick line);
      #endif

      double HitLineSeparation(int hit_channel,float hit_peak,float hit_width,LineWireTick line);

  };

}

#endif
