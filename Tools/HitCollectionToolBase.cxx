#ifndef _HitCollectionToolBase_cxx_
#define _HitCollectionToolBase_cxx_

#include "HitCollectionToolBase.h"

using namespace hyperon;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HitCollectionToolBase::LoadEvent(art::Event const& e){

  Vect_G4.clear();
  Vect_Hit.clear();  

  if(!e.getByLabel(params.get<std::string>("HitModuleLabel"),Handle_Hit)) 
    throw cet::exception("HitCollectionToolBase") << "No Hit Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_Hit,Handle_Hit);

  Assoc_HitSpacePoint = new art::FindManyP<recob::SpacePoint>(Vect_Hit,e,params.get<std::string>("HitSpacePointAssnLabel"));

  for(const art::Ptr<recob::Hit>& hit : Vect_Hit){
    std::vector<art::Ptr<recob::SpacePoint>> spacepoints = Assoc_HitSpacePoint->at(hit.key());
    if(spacepoints.size() != 1 || std::isnan(spacepoints.at(0)->XYZ()[0]) || std::isnan(spacepoints.at(0)->XYZ()[1]) || std::isnan(spacepoints.at(0)->XYZ()[2])) continue;
    m_HitsSpacePoints[hit] = spacepoints.at(0); 
    m_SpacePointsHits[spacepoints.at(0)] = hit;
  }


  if(!IsData){

    ParticlesPerHit = new art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>(Handle_Hit,e,params.get<std::string>("HitTruthAssnLabel"));

    if(!e.getByLabel(params.get<std::string>("G4ModuleLabel"),Handle_G4)) 
      throw cet::exception("HitCollectionToolBase") << "No Geant4 Data Products Found!" << std::endl;

    art::fill_ptr_vector(Vect_G4,Handle_G4);

    partByID.clear();

    for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4)
      partByID.insert(std::make_pair(g4p->TrackId(),g4p));

  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool HitCollectionToolBase::WriteHit(std::vector<art::Ptr<recob::Hit>>& r_hits,std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap,art::Ptr<recob::Hit> hit) const {

  if(m_HitsSpacePoints.find(hit) == m_HitsSpacePoints.end()){
    std::cout << "Hit has associated space point" << std::endl;
    return false;
  }
  r_hits.push_back(hit);
  //r_hitspacepointmap[hit] = spacepoints.at(0);
  r_hitspacepointmap[hit] = m_HitsSpacePoints.at(hit);

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HitCollectionToolBase::GetTruthMatchedHits(const unsigned int& trackid,std::vector<art::Ptr<recob::Hit>>& r_hits,std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap) const {

  // Find all of the hits these particles deposit energy in
  std::vector<simb::MCParticle const*> particleVec;
  std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

  for(art::Ptr<recob::Hit> hit : Vect_Hit){

    particleVec.clear();
    matchVec.clear();
    ParticlesPerHit->get(hit.key(),particleVec,matchVec);
    std::unordered_map<int,double>  trkide;

    // Find the MCParticle that deposited the most energy in this hit
    unsigned int id = 0;
    double maxe = 0.0;
    for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){
      if(matchVec[i_particle]->energy > maxe){
        id = particleVec[i_particle]->TrackId(); 
        maxe = matchVec[i_particle]->energy;
      }
    }

    if(trackid == id)
      WriteHit(r_hits,r_hitspacepointmap,hit);

  } 

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HitCollectionToolBase::GetTruthMatchedHits2(const unsigned int& trackid,std::vector<art::Ptr<recob::Hit>>& r_hits,std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap) const {

  // Find all of the hits these particles deposit energy in
  std::vector<simb::MCParticle const*> particleVec;
  std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

  for(art::Ptr<recob::Hit> hit : Vect_Hit){

    particleVec.clear();
    matchVec.clear();
    ParticlesPerHit->get(hit.key(),particleVec,matchVec);
    std::unordered_map<int,double>  trkide;

    double tote = 0.0;
    for(size_t i_particle=0;i_particle<particleVec.size();++i_particle)
      tote += matchVec[i_particle]->energy;

    for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){
      if((unsigned int)particleVec[i_particle]->TrackId() == trackid && matchVec[i_particle]->energy > HitFraction*tote){
        WriteHit(r_hits,r_hitspacepointmap,hit);
        break;
      }
    }

  } 

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HitCollectionToolBase::GetTruthMatchedHits3(const unsigned int& trackid,std::vector<art::Ptr<recob::Hit>>& r_hits,std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>& r_hitspacepointmap) const {

  // Go through all of the trajectory points along mc particle's path inside TPC, space charge correct them
  if(partByID.find(trackid) == partByID.end()) return;
  art::Ptr<simb::MCParticle> mcp = partByID.at(trackid); 
  std::vector<TVector3> true_pos; 

  // For efficiency - calculate min/max xyz values of tajectory - we can ignore any hits outside of this region 
  TVector3 minpos(1e10,1e10,1e10);
  TVector3 maxpos(-1e10,-1e10,-1e10);

  for(unsigned int i_p=0;i_p<mcp->NumberTrajectoryPoints();i_p++){

    TVector3 pos(mcp->Vx(i_p),mcp->Vy(i_p),mcp->Vz(i_p));
    if(!inActiveTPC(pos)) continue;

    // if space between this traj point and the previous is longer than the TrueRecoTolerance, add
    // some extra points via interpolation
    if(true_pos.size() && (true_pos.back()-pos).Mag() > 0.5*TrueRecoTolerance){
      TVector3 step = pos - true_pos.back();
      int steps = floor(step.Mag()/0.5/TrueRecoTolerance);     
      for(int i_s=0;i_s<steps;i_s++)
        true_pos.push_back(true_pos.back() + step*(1.0/steps));
    }

    true_pos.push_back(pos);    

    // Update dimensions of ROI
    if(true_pos.back().X() < minpos.X()) minpos.SetX(true_pos.back().X());
    if(true_pos.back().Y() < minpos.Y()) minpos.SetY(true_pos.back().Y());
    if(true_pos.back().Z() < minpos.Z()) minpos.SetZ(true_pos.back().Z());
    if(true_pos.back().X() > maxpos.X()) maxpos.SetX(true_pos.back().X());
    if(true_pos.back().Y() > maxpos.Y()) maxpos.SetY(true_pos.back().Y());
    if(true_pos.back().Z() > maxpos.Z()) maxpos.SetZ(true_pos.back().Z());
  }

  std::map<art::Ptr<recob::SpacePoint>,art::Ptr<recob::Hit>>::const_iterator it_sp;

  std::vector<simb::MCParticle const*> particleVec;
  std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

  // Iterate through all spacepoints in event - save any that are within TrueRecoTolerance of
  // true trajectory of particle we want to reconstruct
  for(it_sp = m_SpacePointsHits.begin();it_sp != m_SpacePointsHits.end();it_sp++){
    TVector3 sp(it_sp->first->XYZ()[0],it_sp->first->XYZ()[1],it_sp->first->XYZ()[2]);
    if(!inActiveTPC(sp)) continue;

    geo::Point_t point = {sp.X(),sp.Y(),sp.Z()};
    if(std::isnan(point.X()) || std::isnan(point.Y()) || std::isnan(point.Z())) continue;
 
    geo::Vector_t sce_corr = SCE->GetCalPosOffsets(point);

    // corrrect for spacecharge
    sp.SetX(sp.X()-sce_corr.X()-1.76); 
    sp.SetY(sp.Y()+sce_corr.Y()); 
    sp.SetZ(sp.Z()+sce_corr.Z()-0.1); 

    // Ignore any points outside of ROI for speed
    if(sp.X() > maxpos.X()+5 || sp.X() < minpos.X()-5) continue;
    if(sp.Y() > maxpos.Y()+5 || sp.Y() < minpos.Y()-5) continue;
    if(sp.Z() > maxpos.Z()+5 || sp.Z() < minpos.Z()-5) continue;

    // Check if this spacepoint is close to the mc trajectory
    for(TVector3 true_pos : true_pos){
      if((true_pos-sp).Mag() < TrueRecoTolerance){

        WriteHit(r_hits,r_hitspacepointmap,it_sp->second);

        // If spacepoint is close enough to true trajectory, confirm 
        // it has enough energy from partice we're trying to reconstruct
        particleVec.clear();
        matchVec.clear();
        ParticlesPerHit->get(it_sp->second.key(),particleVec,matchVec);
        std::unordered_map<int,double>  trkide;
        double tote = 0.0;
        for(size_t i_particle=0;i_particle<particleVec.size();++i_particle) tote += matchVec[i_particle]->energy;
        for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){
          if((unsigned int)particleVec[i_particle]->TrackId() == trackid && matchVec[i_particle]->energy > HitFraction*tote){
            WriteHit(r_hits,r_hitspacepointmap,it_sp->second);
            break;
          }
        }

        break;
      }
    }

  } 

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
