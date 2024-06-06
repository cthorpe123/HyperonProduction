#ifndef _SubModuleGeneratorTruth_cxx_
#define _SubModuleGeneratorTruth_cxx_

#include "SubModuleGeneratorTruth.h"

using namespace hyperon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleGeneratorTruth::SubModuleGeneratorTruth(art::Event const& e,fhicl::ParameterSet pset,bool particlegunmode) :
ParticleGunMode(particlegunmode)
{

   if(!e.getByLabel(pset.get<std::string>("GeneratorModuleLabel","generator"),Handle_MCTruth))  
      throw cet::exception("SubModuleGeneratorTruth") << "No MC Truth data product!" << std::endl;

   art::fill_ptr_vector(Vect_MCTruth,Handle_MCTruth);  

   HyperonPDGs = pset.get<std::vector<int>>("HyperonPDGs",{3122,3212,3112,3222});
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

GeneratorTruth SubModuleGeneratorTruth::GetGeneratorTruth(){

   if(!Vect_MCTruth.size()){
      std::cout << "MCTruth vector is empty" << std::endl;
      return theTruth;
   }

   // Set sizes of all vector valued outputs
   theTruth.NMCTruths = Vect_MCTruth.size();
   theTruth.Mode.resize(Vect_MCTruth.size()); 
   theTruth.CCNC.resize(Vect_MCTruth.size()); 
   theTruth.TruePrimaryVertex_X.resize(Vect_MCTruth.size()); 
   theTruth.TruePrimaryVertex_Y.resize(Vect_MCTruth.size()); 
   theTruth.TruePrimaryVertex_Z.resize(Vect_MCTruth.size()); 
   theTruth.InTPC.resize(Vect_MCTruth.size(),false); 
   theTruth.HasHyperon.resize(Vect_MCTruth.size(),false); 
   theTruth.HasKaon.resize(Vect_MCTruth.size(),false); 

   int i_truth=0;
   for(const art::Ptr<simb::MCTruth> &theMCTruth : Vect_MCTruth){

      simb::MCNeutrino Nu = theMCTruth->GetNeutrino();

      int mode = Nu.Mode();
      int ccnc = Nu.CCNC();

      if(ccnc == 0) theTruth.CCNC.at(i_truth) = "CC";
      else theTruth.CCNC.at(i_truth) = "NC";

      if(mode == 0) theTruth.Mode.at(i_truth) = "QEL";
      else if(mode == 1) theTruth.Mode.at(i_truth) = "RES";
      else if(mode == 2) theTruth.Mode.at(i_truth) = "DIS";
      else if(mode == 3) theTruth.Mode.at(i_truth) = "COH";
      else if(mode == 5) theTruth.Mode.at(i_truth) = "ElectronScattering";
      else if(mode == 10) theTruth.Mode.at(i_truth) = "MEC";
      else if(mode == 11) theTruth.Mode.at(i_truth) = "Diffractive";
      else if(mode == 1095) theTruth.Mode.at(i_truth) = "HYP";
      else theTruth.Mode.at(i_truth) = "Other";

      for(int k_particles=0;k_particles<theMCTruth->NParticles();k_particles++){

         simb::MCParticle Part = theMCTruth->GetParticle(k_particles);

         //if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1) 
         // theTruth.TruePrimaryVertex.SetXYZ(Part.Vx(),Part.Vy(),Part.Vz());

         if(isNeutrino(Part.PdgCode()) && Part.StatusCode() == 0){
            SimParticle P = MakeSimParticle(Part);
            P.Origin = 0;
            P.MCTruthIndex = i_truth;
            theTruth.Neutrino.push_back(P);
         }

         // If there is a hyperon in the final state in a QEL event, change mode to HYP
         if(isHyperon(Part.PdgCode()) && Part.StatusCode() == 1 && mode == 0) theTruth.Mode.at(i_truth) = "HYP";

         if(Part.StatusCode() == 1 && Part.PdgCode() == 2112) theTruth.EventHasFinalStateNeutron = true;
         if(Part.StatusCode() == 1 && isHyperon(Part.PdgCode()) && std::find(HyperonPDGs.begin(),HyperonPDGs.end(),abs(Part.PdgCode())) != HyperonPDGs.end()){
            theTruth.EventHasHyperon = true;
            theTruth.HasHyperon.at(i_truth) = true;
         }
         if(Part.StatusCode() == 1 && isKaon(Part.PdgCode())){
           theTruth.EventHasKaon = true;        
           theTruth.HasKaon.at(i_truth) = true;
         } 
         if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1) {
           theTruth.TruePrimaryVertex_X.at(i_truth) = Part.Vx();
           theTruth.TruePrimaryVertex_Y.at(i_truth) = Part.Vy();
           theTruth.TruePrimaryVertex_Z.at(i_truth) = Part.Vz();
         }
      }

      if(inActiveTPC(TVector3(theTruth.TruePrimaryVertex_X.at(i_truth),theTruth.TruePrimaryVertex_Y.at(i_truth),theTruth.TruePrimaryVertex_Z.at(i_truth)))){
        theTruth.NMCTruthsInTPC++;
        theTruth.InTPC.at(i_truth) = true;
      }

      i_truth++;
   }

   if(!ParticleGunMode && theTruth.Neutrino.size() != Vect_MCTruth.size())         
      throw cet::exception("SubModuleGeneratorTruth") << "Sim Neutrino/MCTruth vector size mismatch" << std::endl;

   return theTruth;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
