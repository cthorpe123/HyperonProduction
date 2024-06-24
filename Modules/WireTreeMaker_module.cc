////////////////////////////////////////////////////////////////////////
// Class:       WireTreeMaker
// Plugin Type: analyzer (art v3_03_01)
// File:        WireTreeMaker_module.cc
//
// Generated at Mon Jan 20 06:07:14 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <string>

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "ubevt/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// Root includes
#include "TTree.h"
#include "TVector3.h"

// Local includes
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleReco.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleRecoRepass.h"
#include "ubana/HyperonProduction/Alg/Position_To_Wire.h"


namespace hyperon {
  class WireTreeMaker;
}


class hyperon::WireTreeMaker : public art::EDAnalyzer {

  public:

    explicit WireTreeMaker(fhicl::ParameterSet const& p);

    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    WireTreeMaker(WireTreeMaker const&) = delete;
    WireTreeMaker(WireTreeMaker&&) = delete;
    WireTreeMaker& operator=(WireTreeMaker const&) = delete;
    WireTreeMaker& operator=(WireTreeMaker&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;


    void beginSubRun(const art::SubRun& sr);
    void endSubRun(const art::SubRun& sr);

  private:

    // General event info
    unsigned int t_EventID;
    int t_run,t_subrun,t_event;

    ///////////////////////////////
    //     Wire Signal Info      //
    ///////////////////////////////

    TTree *t_WireTree;

    std::vector<int> t_TrackIndex; // Index values used to identify tracks by HyperonProduction/Analysis code

    std::vector<int> t_TrackStart_Channel_Plane0;
    std::vector<int> t_TrackStart_Time_Plane0;
    std::vector<int> t_TrackEnd_Channel_Plane0;
    std::vector<int> t_TrackEnd_Time_Plane0;

    std::vector<int> t_TrackStart_Channel_Plane1;
    std::vector<int> t_TrackStart_Time_Plane1;
    std::vector<int> t_TrackEnd_Channel_Plane1;
    std::vector<int> t_TrackEnd_Time_Plane1;

    std::vector<int> t_TrackStart_Channel_Plane2;
    std::vector<int> t_TrackStart_Time_Plane2;
    std::vector<int> t_TrackEnd_Channel_Plane2;
    std::vector<int> t_TrackEnd_Time_Plane2;

    std::vector<double> t_TrackStart_X;
    std::vector<double> t_TrackStart_Y;
    std::vector<double> t_TrackStart_Z;
    std::vector<double> t_TrackEnd_X;
    std::vector<double> t_TrackEnd_Y;
    std::vector<double> t_TrackEnd_Z;

    std::vector<double> t_TrackDir_X;
    std::vector<double> t_TrackDir_Y;
    std::vector<double> t_TrackDir_Z;
    std::vector<double> t_TrackEndDir_X;
    std::vector<double> t_TrackEndDir_Y;
    std::vector<double> t_TrackEndDir_Z;

    // Wire Signals
    std::vector<int> t_Wire_Channel_Plane0;
    std::vector<int> t_Wire_Tick_Plane0;
    std::vector<double> t_Wire_Signal_Plane0;

    std::vector<int> t_Wire_Channel_Plane1;
    std::vector<int> t_Wire_Tick_Plane1;
    std::vector<double> t_Wire_Signal_Plane1;

    std::vector<int> t_Wire_Channel_Plane2;
    std::vector<int> t_Wire_Tick_Plane2;
    std::vector<double> t_Wire_Signal_Plane2;

    //////////////////////////
    //   FHICL PARAMETERS   //
    //////////////////////////

    bool fDebug;

    // Producer module labels
    //std::string fTrackLabel;
    //std::string fPFParticleLabel;


    fhicl::ParameterSet f_Reco;
    bool f_GetRecoRepassInfo;
    fhicl::ParameterSet f_RecoRepass;
    std::string f_WireLabel;
    bool f_IsData;
};

////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////

hyperon::WireTreeMaker::WireTreeMaker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  f_Reco(p.get<fhicl::ParameterSet>("Reco")),
  f_GetRecoRepassInfo(p.get<bool>("GetRecoRepassInfo",false)),   
  f_RecoRepass(p.get<fhicl::ParameterSet>("RecoRepass"))
  // More initializers here.
{
  fDebug = p.get<bool>("Debug","false");

  // Module labels
  f_WireLabel = p.get<std::string>("WireLabel");
}

void hyperon::WireTreeMaker::analyze(art::Event const& e)
{

  if(fDebug) std::cout << "New Event" << std::endl;

  ////////////////////////////////
  //  Reset All of the Vectors  //
  ////////////////////////////////

  t_TrackIndex.clear();

  t_TrackStart_Channel_Plane0.clear();
  t_TrackStart_Time_Plane0.clear();
  t_TrackEnd_Channel_Plane0.clear();
  t_TrackEnd_Time_Plane0.clear();

  t_TrackStart_Channel_Plane1.clear();
  t_TrackStart_Time_Plane1.clear();
  t_TrackEnd_Channel_Plane1.clear();
  t_TrackEnd_Time_Plane1.clear();

  t_TrackStart_Channel_Plane2.clear();
  t_TrackStart_Time_Plane2.clear();
  t_TrackEnd_Channel_Plane2.clear();
  t_TrackEnd_Time_Plane2.clear();

  t_TrackStart_X.clear();
  t_TrackStart_Y.clear();
  t_TrackStart_Z.clear();
  t_TrackEnd_X.clear();
  t_TrackEnd_Y.clear();
  t_TrackEnd_Z.clear();

  t_TrackDir_X.clear();
  t_TrackDir_Y.clear();
  t_TrackDir_Z.clear();
  t_TrackEndDir_X.clear();
  t_TrackEndDir_Y.clear();
  t_TrackEndDir_Z.clear();

  t_Wire_Channel_Plane0.clear();
  t_Wire_Tick_Plane0.clear();
  t_Wire_Signal_Plane0.clear();

  t_Wire_Channel_Plane1.clear();
  t_Wire_Tick_Plane1.clear();
  t_Wire_Signal_Plane1.clear();

  t_Wire_Channel_Plane2.clear();
  t_Wire_Tick_Plane2.clear();
  t_Wire_Signal_Plane2.clear();

  ////////////////////////////
  //  Event ID Information  //
  ////////////////////////////

  t_EventID = e.id().event();
  t_run = e.run();
  t_subrun = e.subRun();
  t_event = e.event();

  /////////////////////////////////////////////////////////
  //  Obtain Start Positions of Tracks in the Hierarchy  //
  /////////////////////////////////////////////////////////

  if(fDebug) std::cout << "Getting Reco'd Particles" << std::endl;

  SubModuleReco* Reco_SM = new SubModuleReco(e,f_IsData,f_Reco,false);
  Reco_SM->PrepareInfo();
  RecoData RecoD =  Reco_SM->GetInfo();   

  for(RecoParticle trk : RecoD.TrackPrimaryDaughters){

    t_TrackIndex.push_back(trk.Index);

    t_TrackStart_Channel_Plane0.push_back(trk.TrackStart_Channel_Plane0); 
    t_TrackStart_Time_Plane0.push_back(trk.TrackStart_Time); 
    t_TrackEnd_Channel_Plane0.push_back(trk.TrackEnd_Channel_Plane0); 
    t_TrackEnd_Time_Plane0.push_back(trk.TrackEnd_Time); 

    t_TrackStart_Channel_Plane1.push_back(trk.TrackStart_Channel_Plane1); 
    t_TrackStart_Time_Plane1.push_back(trk.TrackStart_Time); 
    t_TrackEnd_Channel_Plane1.push_back(trk.TrackEnd_Channel_Plane1); 
    t_TrackEnd_Time_Plane1.push_back(trk.TrackEnd_Time); 

    t_TrackStart_Channel_Plane2.push_back(trk.TrackStart_Channel_Plane2); 
    t_TrackStart_Time_Plane2.push_back(trk.TrackStart_Time); 
    t_TrackEnd_Channel_Plane2.push_back(trk.TrackEnd_Channel_Plane2); 
    t_TrackEnd_Time_Plane2.push_back(trk.TrackEnd_Time); 

    t_TrackStart_X.push_back(trk.TrackStartX);
    t_TrackStart_Y.push_back(trk.TrackStartY);
    t_TrackStart_Z.push_back(trk.TrackStartZ);
    t_TrackEnd_X.push_back(trk.TrackEndX);
    t_TrackEnd_Y.push_back(trk.TrackEndY);
    t_TrackEnd_Z.push_back(trk.TrackEndZ);

    t_TrackDir_X.push_back(trk.TrackDirectionX);
    t_TrackDir_Y.push_back(trk.TrackDirectionY);
    t_TrackDir_Z.push_back(trk.TrackDirectionZ);

  }

  delete Reco_SM;

  if(f_GetRecoRepassInfo){

    SubModuleRecoRepass* Reco_SM_Repass = new SubModuleRecoRepass(e,f_IsData,f_RecoRepass);
    RecoRepassData RecoD_Repass =  Reco_SM_Repass->GetInfo();   

    // Add the starts of these new tracks to the end of the track start vector 
    // so we can use them in the CT test

    for(RecoParticle trk : RecoD_Repass.TrackPrimaryDaughters){

      t_TrackIndex.push_back(trk.Index);

      t_TrackStart_Channel_Plane0.push_back(trk.TrackStart_Channel_Plane0); 
      t_TrackStart_Time_Plane0.push_back(trk.TrackStart_Time); 
      t_TrackEnd_Channel_Plane0.push_back(trk.TrackEnd_Channel_Plane0); 
      t_TrackEnd_Time_Plane0.push_back(trk.TrackEnd_Time); 

      t_TrackStart_Channel_Plane1.push_back(trk.TrackStart_Channel_Plane1); 
      t_TrackStart_Time_Plane1.push_back(trk.TrackStart_Time); 
      t_TrackEnd_Channel_Plane1.push_back(trk.TrackEnd_Channel_Plane1); 
      t_TrackEnd_Time_Plane1.push_back(trk.TrackEnd_Time); 

      t_TrackStart_Channel_Plane2.push_back(trk.TrackStart_Channel_Plane2); 
      t_TrackStart_Time_Plane2.push_back(trk.TrackStart_Time); 
      t_TrackEnd_Channel_Plane2.push_back(trk.TrackEnd_Channel_Plane2); 
      t_TrackEnd_Time_Plane2.push_back(trk.TrackEnd_Time); 

      t_TrackStart_X.push_back(trk.TrackStartX);
      t_TrackStart_Y.push_back(trk.TrackStartY);
      t_TrackStart_Z.push_back(trk.TrackStartZ);
      t_TrackEnd_X.push_back(trk.TrackEndX);
      t_TrackEnd_Y.push_back(trk.TrackEndY);
      t_TrackEnd_Z.push_back(trk.TrackEndZ);

      t_TrackDir_X.push_back(trk.TrackDirectionX);
      t_TrackDir_Y.push_back(trk.TrackDirectionY);
      t_TrackDir_Z.push_back(trk.TrackDirectionZ);

    }

    delete Reco_SM_Repass;
  }

  /////////////////////////////
  //   Obtain Wire Signals   //
  /////////////////////////////

  // Setup handles
  art::Handle<std::vector<recob::Wire>> wireHandle;
  std::vector<art::Ptr<recob::Wire>> wireVect;

  // Fill Wire vector
  if(e.getByLabel(f_WireLabel,wireHandle)) art::fill_ptr_vector(wireVect,wireHandle);
  else
    std::cout << "Wire handle not setup" << std::endl;

  // Iterate through all of the wires, record signal at every tick with nonzero signal
  for(const art::Ptr<recob::Wire> &wire : wireVect){

    // Get regions of interest        
    unsigned int NROI = wire->SignalROI().n_ranges();
    for(size_t i_roi=0; i_roi<NROI; ++i_roi){

      // Region of tick space with nonzero activity
      recob::Wire::RegionsOfInterest_t::datarange_t const& range = wire->SignalROI().range(i_roi);

      // Iterate through the ticks in this ROI, record signal
      unsigned int thisTick = range.begin_index();

      while(thisTick < range.end_index()){

        if(wire->View() == 0){
          t_Wire_Channel_Plane0.push_back(wire->Channel());
          t_Wire_Tick_Plane0.push_back(thisTick);
          t_Wire_Signal_Plane0.push_back(wire->Signal().at(thisTick));
        }

        if(wire->View() == 1){
          t_Wire_Channel_Plane1.push_back(wire->Channel());
          t_Wire_Tick_Plane1.push_back(thisTick);
          t_Wire_Signal_Plane1.push_back(wire->Signal().at(thisTick));
        }

        if(wire->View() == 2){
          t_Wire_Channel_Plane2.push_back(wire->Channel());
          t_Wire_Tick_Plane2.push_back(thisTick);
          t_Wire_Signal_Plane2.push_back(wire->Signal().at(thisTick));
        }

        thisTick++;

      } // while(thisTick < range.end_index

    } // loop over ROI

  } // loop over wires

  // Fill tree! 
  t_WireTree->Fill();
}

//////////////////////////////////////////////////////////////////

void hyperon::WireTreeMaker::beginJob(){

  if(fDebug) std::cout << "Begin job" << std::endl;

  art::ServiceHandle<art::TFileService> tfs;

  //////////////////////////////////////////
  //              Wire Tree		   //
  //////////////////////////////////////////

  t_WireTree=tfs->make<TTree>("WireTree","Wire Tree");

  t_WireTree->Branch("EventID",&t_EventID);
  t_WireTree->Branch("run",&t_run);
  t_WireTree->Branch("subrun",&t_subrun);
  t_WireTree->Branch("event",&t_event);

  t_WireTree->Branch("Wire_Channel_Plane0",&t_Wire_Channel_Plane0);
  t_WireTree->Branch("Wire_Tick_Plane0",&t_Wire_Tick_Plane0);
  t_WireTree->Branch("Wire_Signal_Plane0",&t_Wire_Signal_Plane0);

  t_WireTree->Branch("Wire_Channel_Plane1",&t_Wire_Channel_Plane1);
  t_WireTree->Branch("Wire_Tick_Plane1",&t_Wire_Tick_Plane1);
  t_WireTree->Branch("Wire_Signal_Plane1",&t_Wire_Signal_Plane1);

  t_WireTree->Branch("Wire_Channel_Plane2",&t_Wire_Channel_Plane2);
  t_WireTree->Branch("Wire_Tick_Plane2",&t_Wire_Tick_Plane2);
  t_WireTree->Branch("Wire_Signal_Plane2",&t_Wire_Signal_Plane2);

  t_WireTree->Branch("TrackIndex",&t_TrackIndex);

  t_WireTree->Branch("TrackStart_Channel_Plane0",&t_TrackStart_Channel_Plane0);
  t_WireTree->Branch("TrackStart_Time_Plane0",&t_TrackStart_Time_Plane0);
  t_WireTree->Branch("TrackEnd_Channel_Plane0",&t_TrackEnd_Channel_Plane0);
  t_WireTree->Branch("TrackEnd_Time_Plane0",&t_TrackEnd_Time_Plane0);

  t_WireTree->Branch("TrackStart_Channel_Plane1",&t_TrackStart_Channel_Plane1);
  t_WireTree->Branch("TrackStart_Time_Plane1",&t_TrackStart_Time_Plane1);
  t_WireTree->Branch("TrackEnd_Channel_Plane1",&t_TrackEnd_Channel_Plane1);
  t_WireTree->Branch("TrackEnd_Time_Plane1",&t_TrackEnd_Time_Plane1);

  t_WireTree->Branch("TrackStart_Channel_Plane2",&t_TrackStart_Channel_Plane2);
  t_WireTree->Branch("TrackStart_Time_Plane2",&t_TrackStart_Time_Plane2);
  t_WireTree->Branch("TrackEnd_Channel_Plane2",&t_TrackEnd_Channel_Plane2);
  t_WireTree->Branch("TrackEnd_Time_Plane2",&t_TrackEnd_Time_Plane2);

  t_WireTree->Branch("TrackStart_X",&t_TrackStart_X);
  t_WireTree->Branch("TrackStart_Y",&t_TrackStart_Y);
  t_WireTree->Branch("TrackStart_Z",&t_TrackStart_Z);
  t_WireTree->Branch("TrackEnd_X",&t_TrackEnd_X);
  t_WireTree->Branch("TrackEnd_Y",&t_TrackEnd_Y);
  t_WireTree->Branch("TrackEnd_Z",&t_TrackEnd_Z);

  t_WireTree->Branch("TrackDir_X",&t_TrackDir_X);
  t_WireTree->Branch("TrackDir_Y",&t_TrackDir_Y);
  t_WireTree->Branch("TrackDir_Z",&t_TrackDir_Z);
  t_WireTree->Branch("TrackEndDir_X",&t_TrackEndDir_X);
  t_WireTree->Branch("TrackEndDir_Y",&t_TrackEndDir_Y);
  t_WireTree->Branch("TrackEndDir_Z",&t_TrackEndDir_Z);

  if(fDebug) std::cout << "Finished begin job" << std::endl;
}

void hyperon::WireTreeMaker::endJob()
{

}

void hyperon::WireTreeMaker::beginSubRun(const art::SubRun& sr)
{

}

void hyperon::WireTreeMaker::endSubRun(const art::SubRun& sr)
{

}

DEFINE_ART_MODULE(hyperon::WireTreeMaker)
