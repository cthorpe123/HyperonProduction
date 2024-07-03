#ifndef _ShowerReanalyzer_h_
#define _ShowerReanalyzer_h_

#include "lardataobj/RecoBase/Shower.h"
#include "HitCollectionToolBase.h"
#include "VFitter.h"

namespace hyperon {

  class ShowerReanalyzer : public HitCollectionToolBase {

    public:

      ShowerReanalyzer(const fhicl::ParameterSet& p) : 
        HitCollectionToolBase(p),
        Fitter()
    {}

      void MakeHitCollections(std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
          std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
          std::vector<pandora::CartesianVector>& r_vertex) const;

      FittedV FitVToShower(const art::Ptr<recob::Shower> shower) const; 

    private:

      VFitter Fitter;    

  };

}


#endif
