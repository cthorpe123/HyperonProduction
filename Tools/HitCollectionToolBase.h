#ifndef _HitCollectionToolBase_h_
#define _HitCollectionToolBase_h_

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "ubana/KReco/TrackRebuilder/TrackRebuilder.h"
#include "ubana/HyperonProduction/Headers/FV.h"

// Base class for all hit manipulationtools

namespace hyperon {

  class HitCollectionToolBase {

    public: 

      HitCollectionToolBase(const fhicl::ParameterSet& p) : 
        params(p),
        IsData(p.get<bool>("IsData")),
        HitFraction(p.get<double>("HitFraction",0.1)),
        TrueRecoTolerance(p.get<double>("TrueRecoTolerance",1)), 
        SCE(lar::providerFrom<spacecharge::SpaceChargeService>())
    {}

      //virtual ~HitCollectionToolBase(); 

      // makes a vector of the collections required to generate a new track - each entry is for a separate track
      virtual void MakeHitCollections(std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
          std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
          std::vector<pandora::CartesianVector>& r_vertex) const = 0;

      // Load a new event
      void LoadEvent(art::Event const& e);

    protected:

      // Parameters 
      const fhicl::ParameterSet params;  
      const bool IsData;
      const double HitFraction;
      const double TrueRecoTolerance;

      // Services
      const spacecharge::SpaceCharge* SCE;

      // Handles and vectors 
      art::Handle<std::vector<simb::MCParticle>> Handle_G4;
      std::vector<art::Ptr<simb::MCParticle>> Vect_G4;
      art::Handle<std::vector<recob::Hit>> Handle_Hit;
      std::vector<art::Ptr<recob::Hit>> Vect_Hit;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* ParticlesPerHit;
      art::FindManyP<recob::SpacePoint>* Assoc_HitSpacePoint;

      // Maps between hits and spacepoints
      std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>> m_HitsSpacePoints;
      std::map<art::Ptr<recob::SpacePoint>,art::Ptr<recob::Hit>> m_SpacePointsHits;

      // MC Particle by ID lookup
      std::map<int,art::Ptr<simb::MCParticle>> partByID;

      bool WriteHit(std::vector<art::Ptr<recob::Hit>>& r_hits,
          std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap,
          art::Ptr<recob::Hit> hit) const;

      // Different algorithms for getting hits that belong to a single MC particle

      void GetTruthMatchedHits(const unsigned int& trackid,
          std::vector<art::Ptr<recob::Hit>>& r_hits,
          std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap) const;

      void GetTruthMatchedHits2(const unsigned int& trackid,
          std::vector<art::Ptr<recob::Hit>>& r_hits,
          std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap) const;

      void GetTruthMatchedHits3(const unsigned int& trackid,
          std::vector<art::Ptr<recob::Hit>>& r_hits,
          std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap) const;

  };

}

#endif
