#ifndef _ChargedSigmaRecoCheat_h_
#define _ChargedSigmaRecoCheat_h_

#include "HitCollectionToolBase.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

namespace hyperon {

  class ChargedSigmaRecoCheat : public HitCollectionToolBase {

    public: 

      ChargedSigmaRecoCheat(const fhicl::ParameterSet& p) : 
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


