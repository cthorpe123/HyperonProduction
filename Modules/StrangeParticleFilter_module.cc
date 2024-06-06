////////////////////////////////////////////////////////////////////////
// Class:       StrangeParticleFilter
// Plugin Type: filter (art v3_01_02)
// File:        StrangeParticleFilter_module.cc
//
// Purpose: Finds events containing hyperon production (direct and associated)
//
// Generated at Thu Sep 16 09:04:20 2021 by Christopher Thorpe using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"

namespace hyperon {
   class StrangeParticleFilter;
}


class hyperon::StrangeParticleFilter : public art::EDFilter {
   public:
      explicit StrangeParticleFilter(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      StrangeParticleFilter(StrangeParticleFilter const&) = delete;
      StrangeParticleFilter(StrangeParticleFilter&&) = delete;
      StrangeParticleFilter& operator=(StrangeParticleFilter const&) = delete;
      StrangeParticleFilter& operator=(StrangeParticleFilter&&) = delete;

      // Required functions.
      bool filter(art::Event& e) override;

   private:

      fhicl::ParameterSet f_Generator;
};


hyperon::StrangeParticleFilter::StrangeParticleFilter(fhicl::ParameterSet const& p)
   : EDFilter{p},
   f_Generator(p.get<fhicl::ParameterSet>("Generator"))
{
}

bool hyperon::StrangeParticleFilter::filter(art::Event& e)
{
      
      SubModuleGeneratorTruth* Generator_SM = new SubModuleGeneratorTruth(e,f_Generator);
      GeneratorTruth GenT = Generator_SM->GetGeneratorTruth();

      // Return true if event has a hyperon or kaon at a primary vertex inside the TPC
      bool pass = false;
      for(size_t i_h=0;i_h<GenT.HasHyperon.size();i_h++)
        if(GenT.HasHyperon.at(i_h) && GenT.InTPC.at(i_h)) pass = true;
      
      for(size_t i_k=0;i_k<GenT.HasKaon.size();i_k++)
        if(GenT.HasKaon.at(i_k) && GenT.InTPC.at(i_k)) pass = true;
      
      delete Generator_SM;
  
      return pass;

}

DEFINE_ART_MODULE(hyperon::StrangeParticleFilter)
