#ifndef _LambdaRecoCheat2_h_
#define _LambdaRecoCheat2_h_

#include "HitCollectionToolBase.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// Class to get hits from Lambda decay products via truth matching, gets
// each particle via a separate pass (unlike LambdaRecoCheat)

namespace hyperon {

  class LambdaRecoCheat2 : public HitCollectionToolBase {

  public: 

    LambdaRecoCheat2(const fhicl::ParameterSet& p) : 
      HitCollectionToolBase(p), 
      Particle(p.get<int>("ParticlePDG")),
      Alg(p.get<int>("Alg"))
    {
      if(Particle != 2212 && Particle != -211)
        throw cet::exception("LambdaRecoCheat2") << "Request either pdg 2212 or -211" << std::endl;
    }

    void MakeHitCollections(std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
        std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
        std::vector<pandora::CartesianVector>& r_vertex) const;

  private: 

    const int Particle;
    const int Alg;

  };

}

#endif


