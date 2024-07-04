////////////////////////////////////////////////////////////////////////
// Class:       SpacePointTreeMaker
// Plugin Type: analyzer (art v3_03_01)
// File:        SpacePointTreeMaker_module.cc
//
// Generated at Mon Jan 20 06:07:14 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

// cpp stl includes
#include <vector>
#include <string>

// art/larsoft includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

//uboonecode includes
#include "ubevt/Utilities/SignalShapingServiceMicroBooNE.h"

// local includes
#include "ubana/HyperonProduction/Headers/ParticleTypes.h"
#include "ubana/HyperonProduction/Alg/Position_To_Wire.h"

// Root includes
#include "TTree.h"
#include "TVector3.h"

namespace hyperon {
  class SpacePointTreeMaker;
}


class hyperon::SpacePointTreeMaker : public art::EDAnalyzer {

  public:

    explicit SpacePointTreeMaker(fhicl::ParameterSet const& p);

    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    SpacePointTreeMaker(SpacePointTreeMaker const&) = delete;
    SpacePointTreeMaker(SpacePointTreeMaker&&) = delete;
    SpacePointTreeMaker& operator=(SpacePointTreeMaker const&) = delete;
    SpacePointTreeMaker& operator=(SpacePointTreeMaker&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;


    void beginSubRun(const art::SubRun& sr);
    void endSubRun(const art::SubRun& sr);

  private:

    TTree *t_SpacePointTree;
    TTree *t_SpacePointTree_Tracks;
    TTree *t_SpacePointTree_Showers;

    // General event info
    unsigned int t_EventID;
    int t_run,t_subrun,t_event;

    // Spacepoint info
    std::vector<std::vector<double>> t_PFPSpacePoint_X;
    std::vector<std::vector<double>> t_PFPSpacePoint_Y;
    std::vector<std::vector<double>> t_PFPSpacePoint_Z;
    std::vector<std::vector<int>> t_PFPSpacePoint_PDG;
    std::vector<std::vector<double>> t_PFPSpacePoint_PDGPur;
    std::vector<std::vector<std::vector<double>>> t_PFPHit_Channel;
    std::vector<std::vector<std::vector<double>>> t_PFPHit_Tick;
    std::vector<std::vector<std::vector<double>>> t_PFPHit_Width;

    std::vector<std::vector<double>> t_PFPSpacePoint_X_Tracks;
    std::vector<std::vector<double>> t_PFPSpacePoint_Y_Tracks;
    std::vector<std::vector<double>> t_PFPSpacePoint_Z_Tracks;
    std::vector<std::vector<int>> t_PFPSpacePoint_PDG_Tracks;
    std::vector<std::vector<double>> t_PFPSpacePoint_PDGPur_Tracks;
    std::vector<std::vector<std::vector<double>>> t_PFPHit_Channel_Tracks;
    std::vector<std::vector<std::vector<double>>> t_PFPHit_Tick_Tracks;
    std::vector<std::vector<std::vector<double>>> t_PFPHit_Width_Tracks;
    
    std::vector<std::vector<double>> t_PFPSpacePoint_X_Showers;
    std::vector<std::vector<double>> t_PFPSpacePoint_Y_Showers;
    std::vector<std::vector<double>> t_PFPSpacePoint_Z_Showers;
    std::vector<std::vector<int>> t_PFPSpacePoint_PDG_Showers;
    std::vector<std::vector<double>> t_PFPSpacePoint_PDGPur_Showers;
    std::vector<std::vector<std::vector<double>>> t_PFPHit_Channel_Showers;
    std::vector<std::vector<std::vector<double>>> t_PFPHit_Tick_Showers;
    std::vector<std::vector<std::vector<double>>> t_PFPHit_Width_Showers;

    //////////////////////////
    //   FHICL PARAMETERS   //
    //////////////////////////

    bool fDebug;

    // Producer module labels
    const std::string fPFParticleModuleLabel;
    const std::string fSpacePointModuleLabel;
    const std::string fHitModuleLabel;
    const std::string fPFParticleTrackAssnLabel;
    const std::string fPFParticleShowerAssnLabel;
    const std::string fPFParticleSpacePointAssnLabel;
    const std::string fSpacePointHitAssnLabel;
    const std::string fHitTruthAssnLabel;
    //const std::string fTrackModuleLabel;
};

////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////

hyperon::SpacePointTreeMaker::SpacePointTreeMaker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fPFParticleModuleLabel(p.get<std::string>("PFParticleModuleLabel")),
  fSpacePointModuleLabel(p.get<std::string>("SpacePointModuleLabel")),
  fHitModuleLabel(p.get<std::string>("HitModuleLabel")),
  fPFParticleTrackAssnLabel(p.get<std::string>("PFParticleTrackAssnLabel")),
  fPFParticleShowerAssnLabel(p.get<std::string>("PFParticleShowerAssnLabel")),
  fPFParticleSpacePointAssnLabel(p.get<std::string>("PFParticleSpacePointAssnLabel")),
  fSpacePointHitAssnLabel(p.get<std::string>("SpacePointHitAssnLabel")),
  fHitTruthAssnLabel(p.get<std::string>("HitTruthAssnLabel"))
  // More initializers here.
{
  fDebug = p.get<bool>("Debug","false");
}

void hyperon::SpacePointTreeMaker::analyze(art::Event const& e)
{

  if(fDebug) std::cout << "New Event" << std::endl;

  // Get event ID info
  t_EventID = e.id().event();
  t_run = e.run();
  t_subrun = e.subRun();
  t_event = e.event();

  // Clear vectors
  t_PFPSpacePoint_X.clear();
  t_PFPSpacePoint_Y.clear();
  t_PFPSpacePoint_Z.clear();
  t_PFPSpacePoint_PDG.clear();
  t_PFPSpacePoint_PDGPur.clear();
  t_PFPHit_Channel.clear();
  t_PFPHit_Tick.clear();
  t_PFPHit_Width.clear();

  t_PFPSpacePoint_X_Tracks.clear();
  t_PFPSpacePoint_Y_Tracks.clear();
  t_PFPSpacePoint_Z_Tracks.clear();
  t_PFPSpacePoint_PDG_Tracks.clear();
  t_PFPSpacePoint_PDGPur_Tracks.clear();
  t_PFPHit_Channel_Tracks.clear();
  t_PFPHit_Tick_Tracks.clear();
  t_PFPHit_Width_Tracks.clear();

  t_PFPSpacePoint_X_Showers.clear();
  t_PFPSpacePoint_Y_Showers.clear();
  t_PFPSpacePoint_Z_Showers.clear();
  t_PFPSpacePoint_PDG_Showers.clear();
  t_PFPSpacePoint_PDGPur_Showers.clear();
  t_PFPHit_Channel_Showers.clear();
  t_PFPHit_Tick_Showers.clear();
  t_PFPHit_Width_Showers.clear();

  // Setup handles/vectors 
  art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
  std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;
  if(!e.getByLabel(fPFParticleModuleLabel,Handle_PFParticle)) 
    throw cet::exception("SpacePointTreeMaker") << "No PFParticle Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_PFParticle,Handle_PFParticle);

  art::Handle<std::vector<recob::SpacePoint>> Handle_SpacePoint;
  std::vector<art::Ptr<recob::SpacePoint>> Vect_SpacePoint;
  if(!e.getByLabel(fSpacePointModuleLabel,Handle_SpacePoint)) 
    throw cet::exception("SpacePointTreeMaker") << "No SpacePoint Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_SpacePoint,Handle_SpacePoint);

  art::Handle<std::vector<recob::Hit>> Handle_Hit;
  std::vector<art::Ptr<recob::Hit>> Vect_Hit;
  if(!e.getByLabel(fHitModuleLabel,Handle_Hit)) 
    throw cet::exception("SpacePointTreeMaker") << "No Hit Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_Hit,Handle_Hit);

  art::FindManyP<recob::SpacePoint> Assoc_PFParticleSpacePoint(Vect_PFParticle,e,fPFParticleSpacePointAssnLabel);
  art::FindManyP<recob::Hit> Assoc_SpacePointHit(Vect_SpacePoint,e,fSpacePointHitAssnLabel);
  art::FindManyP<recob::Track> Assoc_PFParticleTrack(Vect_PFParticle,e,fPFParticleTrackAssnLabel);
  art::FindManyP<recob::Shower> Assoc_PFParticleShower(Vect_PFParticle,e,fPFParticleShowerAssnLabel);
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> ParticlesPerHit(Handle_Hit,e,fHitTruthAssnLabel);

  // Find the neutrino PFP
  size_t neutrinoID = 99999;
  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle)
    if(pfp->IsPrimary() && isNeutrino(pfp->PdgCode()))
      neutrinoID = pfp->Self();

  if(neutrinoID == 99999){
    if(fDebug) std::cout << "No neutrino candidate in event" << std::endl;
    t_SpacePointTree->Fill();
    return;
  }

  // Get all PFPs that are children of the neutrino, and their respective spacepoints
  int ctr = 0; // TODO: This is a crude way to find the muon candidate - should improve
  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
    if(pfp->Parent() != neutrinoID) continue;
    ctr++;
    if(ctr == 1) continue; // cude way to ignore muon candidate 
    bool is_track=false,is_shower=false;
      
    std::vector<art::Ptr<recob::Track>> pfpTracks = Assoc_PFParticleTrack.at(pfp.key());
    std::vector<art::Ptr<recob::Shower>> pfpShowers = Assoc_PFParticleShower.at(pfp.key());

    //if(pfpTracks.size() == 1) is_track = true;
    //if(pfpShowers.size() == 1) is_shower = true;
     
    if(pfp->PdgCode() == 13) is_track = true; 
    else if(pfp->PdgCode() == 11) is_shower = true;

    if(!is_track && !is_shower) 
      throw cet::exception("SpacePointTreeMaker") << "PFParticle has reconstructed PDG code " << pfp->PdgCode() << std::endl;

    std::vector<art::Ptr<recob::SpacePoint>> pfpSpacePoints = Assoc_PFParticleSpacePoint.at(pfp.key());

    t_PFPSpacePoint_X.push_back(std::vector<double>());
    t_PFPSpacePoint_Y.push_back(std::vector<double>());
    t_PFPSpacePoint_Z.push_back(std::vector<double>());
    t_PFPSpacePoint_PDG.push_back(std::vector<int>());
    t_PFPSpacePoint_PDGPur.push_back(std::vector<double>());
    t_PFPHit_Channel.push_back(std::vector<std::vector<double>>(3,std::vector<double>()));
    t_PFPHit_Tick.push_back(std::vector<std::vector<double>>(3,std::vector<double>()));
    t_PFPHit_Width.push_back(std::vector<std::vector<double>>(3,std::vector<double>()));

    // If the spacepoint has an associated hit, truth match it and save
    for(const art::Ptr<recob::SpacePoint> &sp : pfpSpacePoints){
      std::vector<art::Ptr<recob::Hit>> hits = Assoc_SpacePointHit.at(sp.key());
      if(hits.size() != 1) continue;
      //std::cout << "Found a hit assoc with this spacepoint" << std::endl;

      std::vector<simb::MCParticle const*> particleVec;
      std::vector<anab::BackTrackerHitMatchingData const*> matchVec;
      ParticlesPerHit.get(hits.at(0).key(),particleVec,matchVec);

      int pdg = 0;
      double maxe = 0.0;
      double e = 0.0;
      for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){
        e += matchVec.at(i_particle)->energy;
        if(matchVec.at(i_particle)->energy > maxe){
          pdg = particleVec.at(i_particle)->PdgCode();
          maxe = matchVec.at(i_particle)->energy;
        }
      }
      //std::cout << pdg << std::endl;

      if(hits.at(0)->View() != 0 && hits.at(0)->View() != 1 && hits.at(0)->View() != 2)
        throw cet::exception("SpacePointTreeMaker") << "Hit has unrecognised View number " << hits.at(0)->View() << std::endl;
        
      t_PFPSpacePoint_X.back().push_back(sp->XYZ()[0]); 
      t_PFPSpacePoint_Y.back().push_back(sp->XYZ()[1]); 
      t_PFPSpacePoint_Z.back().push_back(sp->XYZ()[2]); 
      t_PFPHit_Channel.back().at(hits.at(0)->View()).push_back(hits.at(0)->Channel());
      t_PFPHit_Tick.back().at(hits.at(0)->View()).push_back(hits.at(0)->PeakTime());
      t_PFPHit_Width.back().at(hits.at(0)->View()).push_back(hits.at(0)->EndTick() - hits.at(0)->StartTick());
      t_PFPSpacePoint_PDG.back().push_back(pdg); 
      t_PFPSpacePoint_PDGPur.back().push_back(maxe/e); 

    }

    if(is_track){
      t_PFPSpacePoint_X_Tracks.push_back(t_PFPSpacePoint_X.back());
      t_PFPSpacePoint_Y_Tracks.push_back(t_PFPSpacePoint_Y.back());
      t_PFPSpacePoint_Z_Tracks.push_back(t_PFPSpacePoint_Z.back());
      t_PFPSpacePoint_PDG_Tracks.push_back(t_PFPSpacePoint_PDG.back());
      t_PFPSpacePoint_PDGPur_Tracks.push_back(t_PFPSpacePoint_PDGPur.back());
      t_PFPHit_Channel_Tracks.push_back(t_PFPHit_Channel.back());
      t_PFPHit_Tick_Tracks.push_back(t_PFPHit_Tick.back());
      t_PFPHit_Width_Tracks.push_back(t_PFPHit_Width.back());
    }

    if(is_shower){
      t_PFPSpacePoint_X_Showers.push_back(t_PFPSpacePoint_X.back());
      t_PFPSpacePoint_Y_Showers.push_back(t_PFPSpacePoint_Y.back());
      t_PFPSpacePoint_Z_Showers.push_back(t_PFPSpacePoint_Z.back());
      t_PFPSpacePoint_PDG_Showers.push_back(t_PFPSpacePoint_PDG.back());
      t_PFPSpacePoint_PDGPur_Showers.push_back(t_PFPSpacePoint_PDGPur.back());
      t_PFPHit_Channel_Showers.push_back(t_PFPHit_Channel.back());
      t_PFPHit_Tick_Showers.push_back(t_PFPHit_Tick.back());
      t_PFPHit_Width_Showers.push_back(t_PFPHit_Width.back());
    }
   
  } // pfp

  // Fill tree! 
  t_SpacePointTree->Fill();

}

//////////////////////////////////////////////////////////////////

void hyperon::SpacePointTreeMaker::beginJob(){

  if(fDebug) std::cout << "Begin job" << std::endl;

  art::ServiceHandle<art::TFileService> tfs;

  t_SpacePointTree=tfs->make<TTree>("SpacePointTree","Wire Tree");
  t_SpacePointTree_Tracks=tfs->make<TTree>("SpacePointTree_Tracks","Wire Tree");
  t_SpacePointTree_Showers=tfs->make<TTree>("SpacePointTree_Showers","Wire Tree");

  t_SpacePointTree->Branch("EventID",&t_EventID);
  t_SpacePointTree->Branch("run",&t_run);
  t_SpacePointTree->Branch("subrun",&t_subrun);
  t_SpacePointTree->Branch("event",&t_event);
  t_SpacePointTree->Branch("PFPSpacePoint_X",&t_PFPSpacePoint_X);
  t_SpacePointTree->Branch("PFPSpacePoint_Y",&t_PFPSpacePoint_Y);
  t_SpacePointTree->Branch("PFPSpacePoint_Z",&t_PFPSpacePoint_Z);
  t_SpacePointTree->Branch("PFPHit_Channel","vector<vector<vector<double>>>",&t_PFPHit_Channel);
  t_SpacePointTree->Branch("PFPHit_Tick","vector<vector<vector<double>>>",&t_PFPHit_Tick);
  t_SpacePointTree->Branch("PFPHit_Width","vector<vector<vector<double>>>",&t_PFPHit_Width);
  t_SpacePointTree->Branch("PFPSpacePoint_PDG",&t_PFPSpacePoint_PDG);
  t_SpacePointTree->Branch("PFPSpacePoint_PDGPur",&t_PFPSpacePoint_PDGPur);

  t_SpacePointTree_Tracks->Branch("EventID",&t_EventID);
  t_SpacePointTree_Tracks->Branch("run",&t_run);
  t_SpacePointTree_Tracks->Branch("subrun",&t_subrun);
  t_SpacePointTree_Tracks->Branch("event",&t_event);
  t_SpacePointTree_Tracks->Branch("PFPSpacePoint_X",&t_PFPSpacePoint_X_Tracks);
  t_SpacePointTree_Tracks->Branch("PFPSpacePoint_Y",&t_PFPSpacePoint_Y_Tracks);
  t_SpacePointTree_Tracks->Branch("PFPSpacePoint_Z",&t_PFPSpacePoint_Z_Tracks);
  t_SpacePointTree_Tracks->Branch("PFPHit_Channel","vector<vector<vector<double>>>",&t_PFPHit_Channel_Tracks);
  t_SpacePointTree_Tracks->Branch("PFPHit_Tick","vector<vector<vector<double>>>",&t_PFPHit_Tick_Tracks);
  t_SpacePointTree_Tracks->Branch("PFPHit_Width","vector<vector<vector<double>>>",&t_PFPHit_Width_Tracks);
  t_SpacePointTree_Tracks->Branch("PFPSpacePoint_PDG",&t_PFPSpacePoint_PDG_Tracks);
  t_SpacePointTree_Tracks->Branch("PFPSpacePoint_PDGPur",&t_PFPSpacePoint_PDGPur_Tracks);

  t_SpacePointTree_Showers->Branch("EventID",&t_EventID);
  t_SpacePointTree_Showers->Branch("run",&t_run);
  t_SpacePointTree_Showers->Branch("subrun",&t_subrun);
  t_SpacePointTree_Showers->Branch("event",&t_event);
  t_SpacePointTree_Showers->Branch("PFPSpacePoint_X",&t_PFPSpacePoint_X_Showers);
  t_SpacePointTree_Showers->Branch("PFPSpacePoint_Y",&t_PFPSpacePoint_Y_Showers);
  t_SpacePointTree_Showers->Branch("PFPSpacePoint_Z",&t_PFPSpacePoint_Z_Showers);
  t_SpacePointTree_Showers->Branch("PFPHit_Channel","vector<vector<vector<double>>>",&t_PFPHit_Channel_Showers);
  t_SpacePointTree_Showers->Branch("PFPHit_Tick","vector<vector<vector<double>>>",&t_PFPHit_Tick_Showers);
  t_SpacePointTree_Showers->Branch("PFPHit_Width","vector<vector<vector<double>>>",&t_PFPHit_Width_Showers);
  t_SpacePointTree_Showers->Branch("PFPSpacePoint_PDG",&t_PFPSpacePoint_PDG_Showers);
  t_SpacePointTree_Showers->Branch("PFPSpacePoint_PDGPur",&t_PFPSpacePoint_PDGPur_Showers);

  if(fDebug) std::cout << "Finished begin job" << std::endl;

}

void hyperon::SpacePointTreeMaker::endJob()
{

}

void hyperon::SpacePointTreeMaker::beginSubRun(const art::SubRun& sr)
{

}

void hyperon::SpacePointTreeMaker::endSubRun(const art::SubRun& sr)
{

}

DEFINE_ART_MODULE(hyperon::SpacePointTreeMaker)
