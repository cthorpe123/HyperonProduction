#ifndef _VFitter_h_
#define _VFitter_h_

// This object is compilable outside of LArSoft, comment out this line to do so
//#define _InsideLArSoft_

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
#include "FittedV.h"

namespace hyperonreco {

  class VFitter {

    public:

      VFitter();

      void LoadData(const std::vector<std::vector<double>>& channel_v,const std::vector<std::vector<double>>& tick_v,const std::vector<std::vector<double>>& width_v);

      #ifdef _InsideLArSoft_
      FittedV DoFit(TVector3 start,std::map<art::Ptr<recob::SpacePoint>,art::Ptr<recob::Hit>> hitspacepointhap) const;
      #endif

      bool DoFit(FittedV& fittedv);

      void SetROI(double roi);

    private:

      double ROI = 1e10;

      std::vector<std::vector<double>> Channel_v,Tick_v,Width_v;

      #ifdef _InsideLArSoft_
      double HitLineSeparation(art::Ptr<recob::Hit> hit,LineWireTick line);
      std::pair<double,int> FitScore(std::vector<art::Ptr<recob::Hit>> hit_v,FittedV fittedv);
      #endif

      double HitLineSeparation2(const double& channel,const double& tick,const double& width,const LineWireTick& line);
      std::pair<double,int> FitScore(const std::vector<std::vector<double>>& channel_v,const std::vector<std::vector<double>>& tick_v,const std::vector<std::vector<double>>& width_v,FittedV fittedv);

      bool InROI(TVector3 vertex,double channel,double tick,int plane) const;

  };


}

#endif
