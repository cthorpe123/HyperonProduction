#ifndef _ShowerReanalyzer_cxx_
#define _ShowerReanalyzer_cxx_

#include "ShowerReanalyzer.h"

using namespace hyperon;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ShowerReanalyzer::MakeHitCollections(std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
    std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
    std::vector<pandora::CartesianVector>& r_vertex) const {

  // Only consider showers that are children of neutrino for now 
  size_t neutrinoID = 99999;
  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
    if(pfp->IsPrimary() && (pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){
      neutrinoID = pfp->Self();
    } 
  }       
  if(neutrinoID == 99999) return;

  std::vector<unsigned int> neutrino_child_ids;
  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){

    // Skip particles not children/grandchildren of neutrino
    if(pfp->Parent() != neutrinoID && std::find(neutrino_child_ids.begin(),neutrino_child_ids.end(),pfp->Self()) == neutrino_child_ids.end())
      continue;

    neutrino_child_ids.push_back(pfp->Self());

    std::vector<art::Ptr<recob::Shower>> showers = Assoc_PFParticleShower->at(pfp.key());
    if(showers.size() != 1) continue;
    
    FittedV fitresult = FitVToShower(showers.at(0));

  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

FittedV ShowerReanalyzer::FitVToShower(const art::Ptr<recob::Shower> shower) const {

  std::map<art::Ptr<recob::SpacePoint>,art::Ptr<recob::Hit>> hitspacepointmap = MakeSpacePointHitMap(shower);

  FittedV fitted_v;// = Fitter.DoFit(shower->ShowerStart(),hitspacepointmap);

  return fitted_v; 

} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
