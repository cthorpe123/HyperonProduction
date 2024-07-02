#ifndef _VFitter_h_
#define _VFitter_h_

// C++ STL includes
#include <vector>

// larsoft includes
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// root includes
#include "TVector3.h"

// Local includes

namespace hyperon {

  struct FittedV {

   TVector3 Vertex;
   TVector3 Arm1Dir,Arm2Dir;
   double Arm1Len,Arm2Len;     
   double OpeningAngle;
   double Chi2;
   int NDof;
   std::vector<art::Ptr<recob::Hit>> Hits;
 
  };

  class VFitter {

    public:

       VFitter();
     
       void LoadData(std::map<recob::SpacePoint,recob::Hit> hitspacepointhap);
 
    private:
     
        std::map<recob::SpacePoint,recob::Hit> HitSpacePointMap;
  };

}

#endif
