#ifndef _LambdaRecoCheat_cxx_
#define _LambdaRecoCheat_cxx_

#include "ubana/HyperonProduction/Tools/LambdaRecoCheat.h" 

using namespace hyperon;

void LambdaRecoCheat::MakeHitCollections(std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
    std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
    std::vector<pandora::CartesianVector>& r_vertex) const {

  // Iterate through the list of particles - find any Lambdas and their decay products
  std::vector<unsigned int> lambda_children;
  bool found_proton=false,found_pion=false;
  for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4){
    if(g4p->Mother() != 0 || g4p->PdgCode() != 3122 || g4p->EndProcess() != "Decay") continue;
    for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){
      if(partByID.find(g4p->Daughter(i_d)) == partByID.end()) continue;
      art::Ptr<simb::MCParticle> daughter = partByID.at(g4p->Daughter(i_d));
      if(daughter->PdgCode() > 10000) continue;
      lambda_children.push_back(daughter->TrackId());
      if(daughter->PdgCode() == 2212) found_proton = true;
      if(daughter->PdgCode() == -211) found_pion = true;
    }
  } 

  if(lambda_children.size() != 2 || !found_proton || !found_pion) return;

  r_hits.resize(2);
  r_hitspacepointmap.resize(2);
  r_vertex.resize(2,pandora::CartesianVector(-1000,-1000,-1000)); 
  std::vector<bool> got_start(2,false);

  for(size_t i_d=0;i_d<lambda_children.size();i_d++){

    unsigned int particle_id = lambda_children.at(i_d);

    switch (Alg){
      case 1: GetTruthMatchedHits(particle_id,r_hits.at(i_d),r_hitspacepointmap.at(i_d)); break;
      case 2: GetTruthMatchedHits2(particle_id,r_hits.at(i_d),r_hitspacepointmap.at(i_d)); break;
      case 3: GetTruthMatchedHits3(particle_id,r_hits.at(i_d),r_hitspacepointmap.at(i_d)); break;
    }

    art::Ptr<simb::MCParticle> daughter = partByID.at(particle_id);
    geo::Point_t point = {daughter->Vx(),daughter->Vy(),daughter->Vz()};
    geo::Vector_t sce_corr = SCE->GetPosOffsets(point);
    r_vertex.at(i_d) = pandora::CartesianVector(point.X()-sce_corr.X()+1.76,point.Y()+sce_corr.Y(),point.Z()+sce_corr.Z()+0.1);

  } 

  std::cout << "Hits found: " << r_hits.at(0).size() << " " << r_hits.at(1).size() << std::endl;

}


#endif
