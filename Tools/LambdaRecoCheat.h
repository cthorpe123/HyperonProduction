#ifndef _LambdaRecoCheat_h_
#define _LambdaRecoCheat_h_

#include "HitCollectionToolBase.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// Class to get hits from Lambda decay products via truth matching

namespace hyperon {

  class LambdaRecoCheat : public HitCollectionToolBase {

    public: 

      LambdaRecoCheat(const fhicl::ParameterSet& p) : 
        HitCollectionToolBase(p),
        Alg(p.get<int>("Alg"))
    {}

      void MakeHitCollections(std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
          std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
          std::vector<pandora::CartesianVector>& r_vertex) const;

    private: 

      const int Alg;

  };

}

#endif


