////////////////////////////////////////////////////////////////////////
// Class:       LambdaVertexProducer
// Plugin Type: producer (art v3_01_02)
// File:        LambdaVertexProducer_module.cc
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

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/Utilities/AssociationUtil.h"

//#include "ubana/HyperonProduction/Modules/SubModules/SubModuleReco.h"

#include "ubana/HyperonProduction/Tools/LambdaRecoCheat.h"
#include "ubana/HyperonProduction/Tools/LambdaRecoCheat2.h"
#include "ubana/HyperonProduction/Tools/ChargedSigmaRecoCheat.h"

#include <memory>

namespace hyperon {
  class LambdaVertexProducer;
}

class hyperon::LambdaVertexProducer : public art::EDProducer {
  public:
    explicit LambdaVertexProducer(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    LambdaVertexProducer(LambdaVertexProducer const&) = delete;
    LambdaVertexProducer(LambdaVertexProducer&&) = delete;
    LambdaVertexProducer& operator=(LambdaVertexProducer const&) = delete;
    LambdaVertexProducer& operator=(LambdaVertexProducer&&) = delete;

    // Required functions.
    void produce(art::Event& e) override;

  private:

    // Declare member data here.

    const std::string f_InputHitCollectionLabel;
    std::string collectiontooltype;  
    fhicl::ParameterSet f_HitCollectionTool;

};


hyperon::LambdaVertexProducer::LambdaVertexProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  f_InputHitCollectionLabel(p.get<std::string>("InputHitCollectionLabel","gaushit")),
  f_HitCollectionTool(p.get<fhicl::ParameterSet>("HitCollectionTool"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  produces< std::vector<recob::Vertex> >();
  produces<std::vector<recob::Hit>>();
  produces< art::Assns<recob::Slice,recob::Vertex,void> >();

  collectiontooltype = p.get<fhicl::ParameterSet>("HitCollectionTool").get<std::string>("HitCollectionToolType");  

}

void hyperon::LambdaVertexProducer::produce(art::Event& e)
{
  std::unique_ptr<std::vector<recob::Hit>> hitcol(new std::vector<recob::Hit>);
  std::unique_ptr<std::vector<recob::Vertex>> vertexcol(new std::vector<recob::Vertex>);
  std::unique_ptr<art::Assns<recob::Slice,recob::Vertex>> slicevtxassn(new art::Assns<recob::Slice,recob::Vertex>);

  std::unique_ptr<HitCollectionToolBase> hitcollectiontool = nullptr; 

  if(collectiontooltype == "LambdaRecoCheat")
    hitcollectiontool = std::make_unique<LambdaRecoCheat>(f_HitCollectionTool); 
  else if(collectiontooltype == "LambdaRecoCheat2")
    hitcollectiontool = std::make_unique<LambdaRecoCheat2>(f_HitCollectionTool); 
  else if(collectiontooltype == "ChargedSigmaRecoCheat")
    hitcollectiontool = std::make_unique<ChargedSigmaRecoCheat>(f_HitCollectionTool); 

  if(hitcollectiontool == nullptr)
    throw cet::exception("HitCollectionProducer") << "Not hit collection tool selected" << std::endl;

  std::vector<std::vector<art::Ptr<recob::Hit>>> hits_v;
  std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>> r_hitspacepointmap_v;
  std::vector<pandora::CartesianVector> r_vertex_v;

  hitcollectiontool->LoadEvent(e);
  hitcollectiontool->MakeHitCollections(hits_v,r_hitspacepointmap_v,r_vertex_v);

  double xyz[3] = {r_vertex_v.at(0).GetX(),r_vertex_v.at(0).GetY(),r_vertex_v.at(0).GetZ()};
  recob::Vertex newVtx(xyz, -1);
  vertexcol->push_back(newVtx);   

  // Merge the hit collections together and write to event
  for(size_t i=0;i<hits_v.size();i++)
    for(size_t j=0;j<hits_v.at(i).size();j++)
      hitcol->push_back(*hits_v.at(i).at(j));

  
  
  e.put(std::move(hitcol));
  e.put(std::move(vertexcol));
  e.put(std::move(slicevtxassn));

  return;

}

DEFINE_ART_MODULE(hyperon::LambdaVertexProducer)
