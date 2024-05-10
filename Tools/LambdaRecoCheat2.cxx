#ifndef _LambdaRecoCheat2_cxx_
#define _LambdaRecoCheat2_cxx_

#include "ubana/HyperonProduction/Tools/LambdaRecoCheat2.h" 

using namespace hyperon;

void LambdaRecoCheat2::MakeHitCollections(std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
    std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
    std::vector<pandora::CartesianVector>& r_vertex) const {

  // Iterate through the list of particles - find any Lambdas and their decay products
  int lambda_children = 0;
  unsigned int particle_id = 0;   
  for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4){
    if(g4p->Mother() != 0 || g4p->PdgCode() != 3122 || g4p->EndProcess() != "Decay") continue;
    for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){
      if(partByID.find(g4p->Daughter(i_d)) == partByID.end()) continue;
      art::Ptr<simb::MCParticle> daughter = partByID.at(g4p->Daughter(i_d));
      if(daughter->PdgCode() > 10000) continue;
      lambda_children++;
      if(daughter->PdgCode() == Particle) particle_id = daughter->TrackId(); 
    }
  } 

  if(lambda_children != 2 || particle_id == 0) return;

  r_hits.resize(1);
  r_hitspacepointmap.resize(1);
  r_vertex.resize(1,pandora::CartesianVector(-1000,-1000,-1000)); 

  switch (Alg){
    case 1: GetTruthMatchedHits(particle_id,r_hits.at(0),r_hitspacepointmap.at(0)); break;
    case 2: GetTruthMatchedHits2(particle_id,r_hits.at(0),r_hitspacepointmap.at(0)); break;
    case 3: GetTruthMatchedHits3(particle_id,r_hits.at(0),r_hitspacepointmap.at(0)); break;
  }

  art::Ptr<simb::MCParticle> daughter = partByID.at(particle_id);
  geo::Point_t point = {daughter->Vx(),daughter->Vy(),daughter->Vz()};
  geo::Vector_t sce_corr = SCE->GetPosOffsets(point);
  r_vertex.at(0) = pandora::CartesianVector(point.X()-sce_corr.X()+1.76,point.Y()+sce_corr.Y(),point.Z()+sce_corr.Z()+0.1);

  std::cout << "Hits found for pdg " << Particle << ": " << r_hits.at(0).size() << std::endl;

}


#endif
