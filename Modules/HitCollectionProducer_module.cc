////////////////////////////////////////////////////////////////////////
// Class:       HitCollectionProducer
// Plugin Type: producer (art v3_01_02)
// File:        HitCollectionProducer_module.cc
//
// Generated at Wed Aug 21 17:07:38 2019 by Giuseppe Cerati using cetskelgen
// Modified by C Thorpe Sept 2022.
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "ubana/KReco/TrackRebuilder/TrackRebuilder.h"
#include "ubana/HyperonProduction/Tools/LambdaRecoCheat.h"
#include "ubana/HyperonProduction/Tools/ChargedSigmaRecoCheat.h"

#include <memory>

namespace hyperon {
   class HitCollectionProducer;
}

class hyperon::HitCollectionProducer : public art::EDProducer {
   public:
      explicit HitCollectionProducer(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      HitCollectionProducer(HitCollectionProducer const&) = delete;
      HitCollectionProducer(HitCollectionProducer&&) = delete;
      HitCollectionProducer& operator=(HitCollectionProducer const&) = delete;
      HitCollectionProducer& operator=(HitCollectionProducer&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

   private:

   const std::string f_InputHitCollectionLabel;
   kaon_reconstruction::TrackRebuilder trkrebuilder;

   //LambdaRecoCheat hitcollectiontool;
   //ChargedSigmaRecoCheat hitcollectiontool;
   std::unique_ptr<HitCollectionToolBase> hitcollectiontool = nullptr; 

   std::string collectiontooltype;  

};


hyperon::HitCollectionProducer::HitCollectionProducer(fhicl::ParameterSet const& p)
   : EDProducer{p},
   f_InputHitCollectionLabel(p.get<std::string>("InputHitCollectionLabel","gaushit"))
{
   produces<std::vector<recob::Hit>>();
   produces<std::vector<recob::Track>>();

   collectiontooltype = p.get<fhicl::ParameterSet>("HitCollectionTool").get<std::string>("HitCollectionToolType");  

   if(collectiontooltype == "LambdaRecoCheat")
       hitcollectiontool = std::make_unique<LambdaRecoCheat>(p.get<fhicl::ParameterSet>("HitCollectionTool")); 
   else if(collectiontooltype == "ChargedSigmaRecoCheat")
       hitcollectiontool = std::make_unique<ChargedSigmaRecoCheat>(p.get<fhicl::ParameterSet>("HitCollectionTool")); 
  
   if(hitcollectiontool == nullptr)
     throw cet::exception("HitCollectionProducer") << "Not hit collection tool selected" << std::endl;

}

void hyperon::HitCollectionProducer::produce(art::Event& e)
{
   std::unique_ptr<std::vector<recob::Hit>> hitcol(new std::vector<recob::Hit>);
   std::unique_ptr<std::vector<recob::Track>> trackcol(new std::vector<recob::Track>);

   std::vector<std::vector<art::Ptr<recob::Hit>>> hits_v;
   std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>> r_hitspacepointmap_v;
   std::vector<pandora::CartesianVector> r_vertex_v;

   hitcollectiontool->LoadEvent(e);
   hitcollectiontool->MakeHitCollections(hits_v,r_hitspacepointmap_v,r_vertex_v);

   for(size_t i_c=0;i_c<hits_v.size();i_c++){
     if(hits_v.at(i_c).size() < 3) continue;
     std::cout << "Rebuilding track" << std::endl;
     auto status = trkrebuilder.Run(hits_v.at(i_c),r_vertex_v.at(i_c),r_hitspacepointmap_v.at(i_c));
     if(status == STATUS_CODE_SUCCESS){ 
       recob::Track rebuilt_track = trkrebuilder.get_rebuild_reco_track(); 
       trackcol->push_back(rebuilt_track);
       std::cout << "Rebuild track length: " << rebuilt_track.Length() << std::endl;
     }
   }

  // Merge the hit collections together and write to event
  for(size_t i=0;i<hits_v.size();i++)
   for(size_t j=0;j<hits_v.at(i).size();j++)
     hitcol->push_back(*hits_v.at(i).at(j));

  e.put(std::move(hitcol));
  e.put(std::move(trackcol));

  return;

}

DEFINE_ART_MODULE(hyperon::HitCollectionProducer)
