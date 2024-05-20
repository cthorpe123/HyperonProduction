#ifndef _SubModuleRecoRepass_h_
#define _SubModuleRecoRepass_h_

#include <string>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "ubana/HyperonProduction/Headers/ParticleTypes.h"
#//include "ubana/HyperonProduction/Headers/LLR_PID.h"
#//include "ubana/HyperonProduction/Headers/LLRPID_proton_muon_lookup.h"
//#include "ubana/HyperonProduction/Headers/LLR_PID_K.h"
//#include "ubana/HyperonProduction/Headers/LLRPID_kaon_proton_lookup.h"
#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"
#include "ubana/HyperonProduction/Alg/PIDManager.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/HyperonProduction/Alg/BDTHandle.h"

#include "TVector3.h"

// Temporary submodule to collect info about regenerated tracks
// produced using new reco tools

using std::string;

namespace hyperon {

  struct RecoRepassData {

    TVector3 RecoPrimaryVertex = TVector3(-1000,-1000,-1000);

    int NPrimaryDaughters; 
    int NPrimaryTrackDaughters;
    int NPrimaryShowerDaughters;

    std::vector<RecoParticle> TrackPrimaryDaughters;
    std::vector<TVector3> TrackStarts;

    size_t TrueDecayProtonIndex = -1;
    size_t TrueDecayPionIndex = -1;

    bool GoodReco = false;
  };

  class SubModuleRecoRepass {

    public:

      SubModuleRecoRepass(art::Event const& e,bool isdata,string tracklabel,
          string pidlabel,string calolabel,string hitlabel,
          string hittruthassnlabel,string trackhitassnlabel,string genlabel,
          string g4label,fhicl::ParameterSet pidsettings,bool dogetpids);

      SubModuleRecoRepass(art::Event const& e,bool isdata,fhicl::ParameterSet pset,bool particlegunmode=false);

      RecoRepassData GetInfo();

    private:

      art::Handle<std::vector<recob::Track>> Handle_Track;
      std::vector<art::Ptr<recob::Track>> Vect_Track;

      art::Handle<std::vector<recob::Hit>> Handle_Hit;
      std::vector<art::Ptr<recob::Hit>> Vect_Hit;

      RecoParticle MakeRecoParticle(const art::Ptr<recob::Track> &trk);

      art::FindManyP<recob::Hit>* Assoc_TrackHit;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* Assoc_MCParticleBacktracker;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* ParticlesPerHit;
      art::FindManyP<anab::Calorimetry>* Assoc_TrackCalo;
      art::FindManyP<anab::ParticleID>* Assoc_TrackPID;

      SubModuleG4Truth* G4T = nullptr;
      PIDManager PIDCalc;      

      RecoRepassData theData;
      size_t neutrinoID = 99999;

      void GetTrackData(const art::Ptr<recob::Track> &trk,RecoParticle &P);
      void TruthMatch(const art::Ptr<recob::Track> &trk,RecoParticle &P);
      void GetPIDs(const art::Ptr<recob::Track> &trk,RecoParticle &P);

      bool IsData;
      bool DoGetPIDs=true;

  };

}

#endif


