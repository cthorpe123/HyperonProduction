#ifndef _ChargedSigmaRecoCheat_cxx_
#define _ChargedSigmaRecoCheat_cxx_

#include "ubana/HyperonProduction/Tools/ChargedSigmaRecoCheat.h" 

using namespace hyperon;

void ChargedSigmaRecoCheat::MakeHitCollections(std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
    std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
    std::vector<pandora::CartesianVector>& r_vertex) const {

  // Iterate through the list of particles - find any ChargedSigmas and their decay products
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

  std::vector<unsigned int> ids = {sigma_trkid,pion_trkid}; 

  r_hits.resize(2);
  r_hitspacepointmap.resize(2);
  r_vertex.resize(2,pandora::CartesianVector(-1000,-1000,-1000)); 

  for(size_t i_p=0;i_p<ids.size();i_p++){

    art::Ptr<simb::MCParticle> part = partByID.at(ids.at(i_p));

    switch (Alg){
      case 1: GetTruthMatchedHits(ids.at(i_p),r_hits.at(i_p),r_hitspacepointmap.at(i_p)); break;
      case 2: GetTruthMatchedHits2(ids.at(i_p),r_hits.at(i_p),r_hitspacepointmap.at(i_p)); break;
      case 3: GetTruthMatchedHits3(ids.at(i_p),r_hits.at(i_p),r_hitspacepointmap.at(i_p)); break;
    }

    geo::Point_t point = {part->Vx(),part->Vy(),part->Vz()};
    geo::Vector_t sce_corr = SCE->GetPosOffsets(point);
    r_vertex.at(i_p) = pandora::CartesianVector(point.X()-sce_corr.X()+1.76,point.Y()+sce_corr.Y(),point.Z()+sce_corr.Z()+0.1);

  }

  std::cout << "Hits found: " << r_hits.at(0).size() << " " << r_hits.at(1).size() << std::endl;

}


#endif
