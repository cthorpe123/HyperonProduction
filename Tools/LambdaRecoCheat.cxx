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
    GetTruthMatchedHits(lambda_children.at(i_d),r_hits.at(i_d),r_hitspacepointmap.at(i_d));
    if(r_hits.at(i_d).size()){ 
      art::Ptr<recob::SpacePoint> point = r_hitspacepointmap.at(i_d).begin()->second;
      r_vertex.at(i_d) = pandora::CartesianVector(point->XYZ()[0],point->XYZ()[1],point->XYZ()[2]);
    }
  } 

  std::cout << "Hits found: " << r_hits.at(0).size() << " " << r_hits.at(1).size() << std::endl;

}


#endif
