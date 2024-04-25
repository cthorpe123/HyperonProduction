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
   ChargedSigmaRecoCheat hitcollectiontool;


};


hyperon::HitCollectionProducer::HitCollectionProducer(fhicl::ParameterSet const& p)
   : EDProducer{p},
   f_InputHitCollectionLabel(p.get<std::string>("InputHitCollectionLabel","gaushit")),
   hitcollectiontool(p.get<fhicl::ParameterSet>("HitCollectionTool"))
{
   produces<std::vector<recob::Hit>>();
   produces<std::vector<recob::Track>>();
}

void hyperon::HitCollectionProducer::produce(art::Event& e)
{
   std::unique_ptr<std::vector<recob::Hit>> hitcol(new std::vector<recob::Hit>);
   std::unique_ptr<std::vector<recob::Track>> trackcol(new std::vector<recob::Track>);

   std::vector<std::vector<art::Ptr<recob::Hit>>> hits_v;
   std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>> r_hitspacepointmap_v;
   std::vector<pandora::CartesianVector> r_vertex_v;

   hitcollectiontool.LoadEvent(e);
   hitcollectiontool.MakeHitCollections(hits_v,r_hitspacepointmap_v,r_vertex_v);

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

/*
   art::Handle<std::vector<recob::Hit>> Handle_Hit;
   std::vector<art::Ptr<recob::Hit>> Vect_Hit;
   
   if(!e.getByLabel(f_InputHitCollectionLabel,Handle_Hit)) 
      throw cet::exception("HitCollectionProducer") << "Input recob::Hit collection not found" << std::endl;

   art::fill_ptr_vector(Vect_Hit,Handle_Hit);

  // Test - get the first track and generate new hit collection for it
  art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
  std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;
  if(!e.getByLabel("pandora",Handle_PFParticle)) 
    throw cet::exception("HitCollectionProducer") << "No PFParticle Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_PFParticle,Handle_PFParticle);

  art::Handle<std::vector<recob::Track>> Handle_Tracks;
  std::vector<art::Ptr<recob::Track>> Vect_Tracks;
  if(!e.getByLabel("pandora",Handle_Tracks))  
    throw cet::exception("HitCollectionProducer") << "No Track data product!" << std::endl;
  art::fill_ptr_vector(Vect_Tracks,Handle_Tracks);

  art::FindManyP<recob::Track>* Assoc_PFParticleTrack = new art::FindManyP<recob::Track>(Vect_PFParticle,e,"pandora");    
  art::FindManyP<recob::Hit>* Assoc_TrackHit = new art::FindManyP<recob::Hit>(Vect_Tracks,e,"pandora");
  art::FindManyP<recob::Vertex> Assoc_PFParticleVertex = art::FindManyP<recob::Vertex>(Vect_PFParticle,e,"pandora");
  art::FindManyP<recob::SpacePoint> Assoc_HitSpacePoint = art::FindManyP<recob::SpacePoint>(Vect_Hit,e,"pandora");

  size_t neutrinoID = 99999;
  pandora::CartesianVector vertex_pos(-1000,-1000,-1000); 
  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
    if(pfp->IsPrimary() && (abs(pfp->PdgCode()) == 14 || abs(pfp->PdgCode()) == 12)){
      neutrinoID = pfp->Self();
      std::vector<art::Ptr<recob::Vertex>> pfpVertex = Assoc_PFParticleVertex.at(pfp.key());
      for(const art::Ptr<recob::Vertex> &vtx : pfpVertex)
        vertex_pos = pandora::CartesianVector(vtx->position().X() , vtx->position().Y() , vtx->position().Z());
    }
  }

  std::vector<art::Ptr<recob::Hit>> hits;
  std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> hitToSpacePointMap;
  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
    if(pfp->Parent() != neutrinoID) continue; 

    // Grab the associated track
    std::vector<art::Ptr<recob::Track>> pfpTracks = Assoc_PFParticleTrack->at(pfp.key());
    if(pfpTracks.size() != 1) continue;
    art::Ptr<recob::Track> track = pfpTracks.at(0);
    std::cout << "Original track length: " << track->Length() << std::endl; 

    std::cout << "Gettig associated hits" << std::endl;
    std::vector<art::Ptr<recob::Hit>> hits_tmp = Assoc_TrackHit->at(track.key());
    for(art::Ptr<recob::Hit> hit : hits_tmp){
      std::vector<art::Ptr<recob::SpacePoint>> spacepoints = Assoc_HitSpacePoint.at(hit.key());
      if(spacepoints.size() != 1) continue;
      hitcol->push_back(*hit); 
      hits.push_back(hit);
      hitToSpacePointMap[hit] = spacepoints.at(0);
    }
    std::cout << "Done getting associated hits" << std::endl;

    break;
  } 


  // get tracks
  std::cout << "Getting tracks" <<  std::endl;
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > trackList;
  if(e.getByLabel("pandora",trackListHandle)) {
    art::fill_ptr_vector(trackList, trackListHandle);
  }
  art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, e,"pandora");
  std::cout << "Done getting tracks" << std::endl;
 */

/*
  if(hits.size()){
    std::cout << "Rebuilding track" << std::endl;
    auto status = trkrebuilder.Run(hits,vertex_pos,hitToSpacePointMap);
    if(status == STATUS_CODE_SUCCESS){ 
      std::cout << "Rebuilt track" << std::endl;
      recob::Track rebuilt_track = trkrebuilder.get_rebuild_reco_track(); 
      trackcol->push_back(rebuilt_track);
      std::cout << "Rebuild track length: " << rebuilt_track.Length() << std::endl;
    }
  }
*/

  e.put(std::move(hitcol));
  e.put(std::move(trackcol));

  return;

}

DEFINE_ART_MODULE(hyperon::HitCollectionProducer)
