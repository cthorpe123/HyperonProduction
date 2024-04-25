#ifndef _ChargedSigmaRecoCheat_cxx_
#define _ChargedSigmaRecoCheat_cxx_

#include "ubana/HyperonProduction/Tools/ChargedSigmaRecoCheat.h" 

using namespace hyperon;

void ChargedSigmaRecoCheat::MakeHitCollections(std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
    std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
    std::vector<pandora::CartesianVector>& r_vertex) const {

  // Iterate through the list of particles - find any ChargedSigmas and their decay products
  std::vector<unsigned int> sigma_children;
  unsigned int pion_trkid = 0, sigma_trkid = 0;
  for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4){
    if(g4p->Mother() != 0 || !(g4p->PdgCode() == 3222 || g4p->PdgCode() == 3112) || g4p->EndProcess() != "Decay") continue;
    sigma_trkid = g4p->TrackId();
    for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){
      if(partByID.find(g4p->Daughter(i_d)) == partByID.end()) continue;
      art::Ptr<simb::MCParticle> daughter = partByID.at(g4p->Daughter(i_d));
      if(daughter->PdgCode() > 10000) continue;
      if(abs(daughter->PdgCode()) == 211) pion_trkid = daughter->TrackId();
    }
  } 

  if(sigma_trkid == 0 || pion_trkid == 0) return;

  r_hits.resize(2);
  r_hitspacepointmap.resize(2);
  r_vertex.resize(2,pandora::CartesianVector(-1000,-1000,-1000)); 

  GetTruthMatchedHits(sigma_trkid,r_hits.at(0),r_hitspacepointmap.at(0));
  if(r_hits.at(0).size()){ 
    art::Ptr<recob::SpacePoint> point = r_hitspacepointmap.at(0).begin()->second;
    r_vertex.at(0) = pandora::CartesianVector(point->XYZ()[0],point->XYZ()[1],point->XYZ()[2]);
  }

  GetTruthMatchedHits(pion_trkid,r_hits.at(1),r_hitspacepointmap.at(1));
  if(r_hits.at(1).size()){ 
    art::Ptr<recob::SpacePoint> point = r_hitspacepointmap.at(1).begin()->second;
    r_vertex.at(1) = pandora::CartesianVector(point->XYZ()[0],point->XYZ()[1],point->XYZ()[2]);
  }

  std::cout << "Hits found: " << r_hits.at(0).size() << " " << r_hits.at(1).size() << std::endl;

}


#endif
