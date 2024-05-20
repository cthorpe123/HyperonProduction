#ifndef _SubModuleRecoRepass_cxx_
#define _SubModuleRecoRepass_cxx_

#include "SubModuleRecoRepass.h"

using namespace hyperon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

//SubModuleRecoRepass::SubModuleRecoRepass(){}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleRecoRepass::SubModuleRecoRepass(art::Event const& e,bool isdata,fhicl::ParameterSet pset,bool particlegunmode) :
SubModuleRecoRepass(e,isdata,
                  pset.get<std::string>("TrackModuleLabel"),
                  pset.get<std::string>("PIDModuleLabel"),
                  pset.get<std::string>("CaloModuleLabel"),
                  pset.get<std::string>("HitModuleLabel"),
                  pset.get<std::string>("HitTruthAssnLabel"),
                  pset.get<std::string>("TrackHitAssnLabel"),
                  pset.get<std::string>("GeneratorModuleLabel"),
                  pset.get<std::string>("G4ModuleLabel"),
                  pset.get<fhicl::ParameterSet>("PIDSettings"),
                  pset.get<bool>("DoGetPIDs",true))
{

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleRecoRepass::SubModuleRecoRepass(art::Event const& e,bool isdata,string tracklabel,
                                     string pidlabel,string calolabel,string hitlabel,
                                     string hittruthassnlabel,string trackhitassnlabel,string genlabel,
                                     string g4label,fhicl::ParameterSet pidsettings,bool dogetpids) :
PIDCalc(pidsettings),
DoGetPIDs(dogetpids)
{

   IsData = isdata;

   if(!e.getByLabel(tracklabel,Handle_Track)) 
      throw cet::exception("SubModuleRecoRepass") << "No Track Data Products Found!" << std::endl;

   if(!e.getByLabel(hitlabel,Handle_Hit)) 
      throw cet::exception("SubModuleRecoRepass") << "No Hit Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_Track,Handle_Track);
   art::fill_ptr_vector(Vect_Hit,Handle_Hit);

   Assoc_TrackHit = new  art::FindManyP<recob::Hit>(Vect_Track,e,trackhitassnlabel);
   ParticlesPerHit = new art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>(Handle_Hit,e,hittruthassnlabel);

   if(DoGetPIDs){
      Assoc_TrackCalo = new art::FindManyP<anab::Calorimetry>(Vect_Track,e,calolabel);
      Assoc_TrackPID = new art::FindManyP<anab::ParticleID>(Vect_Track,e,pidlabel);
   }

   if(!IsData){
      G4T = new SubModuleG4Truth(e,genlabel,g4label);
      G4T->GetParticleLists();
   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

RecoRepassData SubModuleRecoRepass::GetInfo(){


   for(art::Ptr<recob::Track> trk : Vect_Track)
        theData.TrackPrimaryDaughters.push_back(MakeRecoParticle(trk));

   return theData;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO: This function is a bit pointless
RecoParticle SubModuleRecoRepass::MakeRecoParticle(const art::Ptr<recob::Track> &trk){

   RecoParticle P;
   GetTrackData(trk,P);

   return P;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleRecoRepass::GetTrackData(const art::Ptr<recob::Track> &trk,RecoParticle &P){

   // Sets track length/position related variables
   SetTrackVariables(P,trk);
   if(!IsData) TruthMatch(trk,P);
   if(DoGetPIDs) GetPIDs(trk,P);
   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleRecoRepass::TruthMatch(const art::Ptr<recob::Track> &trk,RecoParticle &P){

   std::vector<art::Ptr<recob::Hit>> hits = Assoc_TrackHit->at(trk.key());

   std::unordered_map<int,double>  trkide;
   int maxhits=-1;

   simb::MCParticle const* matchedParticle = NULL;

   std::vector<simb::MCParticle const*> particleVec;
   std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

   for(size_t i_hit=0;i_hit<hits.size();++i_hit){

      particleVec.clear();
      matchVec.clear();
      ParticlesPerHit->get(hits[i_hit].key(),particleVec,matchVec);

      for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){

         trkide[particleVec[i_particle]->TrackId()]++; 

         //new method - choose particle depositing energy in the most hits
         if(trkide[particleVec[i_particle]->TrackId()] > maxhits){
            maxhits = trkide[particleVec[i_particle]->TrackId()];
            matchedParticle = particleVec[i_particle];
         }
      }
   }

   if(matchedParticle != NULL){ 

      SimParticle SP = MakeSimParticle(*matchedParticle);
            
      SP.Origin = G4T->GetOrigin(matchedParticle->TrackId());
      G4T->MCTruthMatch(SP,matchedParticle->TrackId());
 
      P.HasTruth = true;
      P.MCTruthIndex = SP.MCTruthIndex;
      P.TrackTruePDG = SP.PDG;
      P.TrackTrueE = SP.E;
      P.TrackTruePx = SP.Px;
      P.TrackTruePy = SP.Py;
      P.TrackTruePz = SP.Pz;
      P.TrackTrueEndE = SP.E;
      P.TrackTrueEndPx = SP.EndPx;
      P.TrackTrueEndPy = SP.EndPy;
      P.TrackTrueEndPz = SP.EndPz;
      P.TrackTrueModMomentum = SP.ModMomentum;
      P.TrackTrueEndModMomentum = SP.EndModMomentum;
      P.TrackTrueKE = SP.KE;
      P.TrackTrueEndKE = SP.EndKE;
      P.TrackTrueLength = SP.Travel;
      P.TrackTrueOrigin = SP.Origin;
      P.TrackTruthPurity = (double)maxhits/hits.size();
   }
   else P.HasTruth = false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleRecoRepass::GetPIDs(const art::Ptr<recob::Track> &trk,RecoParticle &P){

   std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = Assoc_TrackCalo->at(trk.key());
   std::vector<art::Ptr<anab::ParticleID>> trackPID = Assoc_TrackPID->at(trk.key());
   std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

   PIDStore store = PIDCalc.GetPIDs(trk,caloFromTrack,AlgScoresVec);
   P.Track_LLR_PID = store.LLR;
   P.Track_LLR_PID_Kaon = store.LLR_Kaon;
   P.Track_LLR_PID_Kaon_Partial = store.LLR_Kaon_Partial;
   P.MeandEdX_Plane0 = store.MeandEdX_Plane0;
   P.MeandEdX_Plane1 = store.MeandEdX_Plane1;
   P.MeandEdX_Plane2 = store.MeandEdX_Plane2;
   P.MeandEdX_ThreePlane = store.MeandEdX_3Plane;
   P.Track_Bragg_PID_Kaon = store.Bragg_Kaon_3Plane;
   P.Track_LLR_PID_SigmaKaon = store.LLR_SigmaKaon;
   P.Track_LLR_PID_SigmaProton = store.LLR_SigmaProton;
   P.Track_LLR_PID_SigmaMuon = store.LLR_SigmaMuon;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
