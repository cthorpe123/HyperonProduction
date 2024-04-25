#ifndef _HitCollectionToolBase_cxx_
#define _HitCollectionToolBase_cxx_

#include "HitCollectionToolBase.h"

using namespace hyperon;

void HitCollectionToolBase::LoadEvent(art::Event const& e){

  Vect_G4.clear();
  Vect_Hit.clear();  

  if(!e.getByLabel(params.get<std::string>("HitModuleLabel"),Handle_Hit)) 
    throw cet::exception("HitCollectionToolBase") << "No Hit Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_Hit,Handle_Hit);

  Assoc_HitSpacePoint = new art::FindManyP<recob::SpacePoint>(Vect_Hit,e,params.get<std::string>("HitSpacePointAssnLabel"));

  if(!IsData){

    ParticlesPerHit = new art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>(Handle_Hit,e,params.get<std::string>("HitTruthAssnLabel"));

    if(!e.getByLabel(params.get<std::string>("G4ModuleLabel"),Handle_G4)) 
      throw cet::exception("HitCollectionToolBase") << "No Geant4 Data Products Found!" << std::endl;

    art::fill_ptr_vector(Vect_G4,Handle_G4);

    partByID.clear();

    for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4)
      partByID.insert(std::make_pair(g4p->TrackId(),g4p));

  }

}

bool HitCollectionToolBase::WriteHit(std::vector<art::Ptr<recob::Hit>>& r_hits,std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap,art::Ptr<recob::Hit> hit) const {

  std::vector<art::Ptr<recob::SpacePoint>> spacepoints = Assoc_HitSpacePoint->at(hit.key());
  if(spacepoints.size() != 1) return false;
  r_hits.push_back(hit);
  r_hitspacepointmap[hit] = spacepoints.at(0);

  return true;

}

void HitCollectionToolBase::GetTruthMatchedHits(const unsigned int& trackid,std::vector<art::Ptr<recob::Hit>>& r_hits,std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap) const {

  // Find all of the hits these particles deposit energy in
  std::vector<simb::MCParticle const*> particleVec;
  std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

  for(art::Ptr<recob::Hit> hit : Vect_Hit){

    particleVec.clear();
    matchVec.clear();
    ParticlesPerHit->get(hit.key(),particleVec,matchVec);
    std::unordered_map<int,double>  trkide;
 
    // Find the MCParticle that deposited the most energy in this hit
    unsigned int id = 0;
    double maxe = 0.0;
    for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){
      if(matchVec[i_particle]->energy > maxe){
        id = particleVec[i_particle]->TrackId(); 
        maxe = matchVec[i_particle]->energy;
      }
    }

    if(trackid == id)
      WriteHit(r_hits,r_hitspacepointmap,hit);

  } 

}

#endif
