#ifndef __ELEC_TREE_H_
#define __ELEC_TREE_H_


/*

Author: Anil P Singh
 Panjab University

*/

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/eSelector.h"

using namespace std;
using namespace edm;

class patElecTree{


 public:
  patElecTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~patElecTree();
  void Fill(const edm::Event& iEvent);
  double IsoDeposit(const pat::Electron ele,std::string type,double radius,double coneradius,double detastrip);
  void SetBranches();
  void Clear();
  void ConversionRejection(const edm::Event&, const edm::EventSetup&, const pat::Electron&);

 private:
//TTree* tree_;
  //Dont Allow User to Call the Default Constructor.
  patElecTree();


  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<std::string>*, std::string name);
void AddBranch(std::vector<int>*, std::string name);
  edm::InputTag dcsTag_;
  bool isData_;

  edm::InputTag patElecLabel_;
  edm::InputTag beamSpotLabel_;
  edm::InputTag EleRhoLabel_;
  eSelector e2012ID_;
  int Nelecs;
  std::vector<int> patElecType;
 //pt, eta, phi, M : Enough to caluclate
  //px,py,pz etc.==>No need to store later
  TTree* tree_;
  double nele; 
  std::vector<double> fastJetRho_;
  std::vector<double> PFchsJetRho_;
  std::vector<double> lepIsoRho_;
  std::vector<double> patElecEt_;
  std::vector<double> patElecEnergy_;
  std::vector<double> patElecPt_;
  std::vector<double> patElecEta_;
  std::vector<double> patElecPhi_;
  std::vector<double> patElecM_;
  std::vector<double> patElecScEn_;
  std::vector<double> patElecScEt_;
  std::vector<double> patElecScEta_;
  std::vector<double> patElecScPhi_;
  std::vector<double> patElecisEcalDriven_;
  std::vector<double> patElecisTrackerDriven_;
  std::vector<double> patElecSigIhIh_;
  std::vector<double> patElecDelEtaIn_;
  std::vector<double> patElecDelPhiIn_;
  std::vector<double> patElecHoE_;
  std::vector<double> patElecTrkIso_;
  std::vector<double> patElecHcalIso_;
  std::vector<double> patElecEcalIso_;
  std::vector<double> patElecRelIsoComb_;
  std::vector<double> patElecCharge_;
  std::vector<double> patElecChi2Ndoff_;
  std::vector<double> patElecNhits_;
  std::vector<double> patElecQoverP_;
  std::vector<double> patElecDxy_;
  std::vector<double> patElecD0_;
  std::vector<double> patElecDsz_;
  std::vector<double> patElecDz_;
  std::vector<double> patElecDxyBS_;
  std::vector<double> patElecDzBS_;
  std::vector<double> patElecDszBS_;
  std::vector<double> patElecVx_;
  std::vector<double> patElecVy_;
  std::vector<double> patElecVz_;
//  std::vector<double> patElecGenMatchPt_;
 // std::vector<double> patElecGenMatchEta_;
 // std::vector<double> patElecGenMatchPhi_;
 // std::vector<double> patElecGenMotherId_;
 // std::vector<double> patElecGenMatchCharge_;
   std::vector<double> patElecChHadSumPt03_;
  std::vector<double> patElecNeHadSumPt03_;
  std::vector<double> patElecGamSumPt03_;
  std::vector<double> patElecChHadSumPt04_;
  std::vector<double> patElecChHadIso_;
  std::vector<double> patElecNeHadIso_;
  std::vector<double> patElecGamIso_;
  std::vector<double> patElecNeHadSumPt04_;
  std::vector<double> patElecGamSumPt04_;
   std::vector<double> patElecChHadSumPt05_;
  std::vector<double> patElecNeHadSumPt05_;
  std::vector<double> patElecGamSumPt05_;
   

  std::vector<double> patElecMva_;
  //std::vector<double> patElecNumMisHits_;

  // conversion rejection
  std::vector<double> patElecMissingHits_ ;
  std::vector<double> patElecDist_ ;
  std::vector<double> patElecDeltaCotTheta_ ;
  std::vector<double> patElecConvRadius_ ;
  std::vector<double> patElecInBarrel_;
  std::vector<double> patElecInEndcap_;
  
  //add for VBF
  std::vector<double> patElectrackMomentumAtVtxMag2_;
  std::vector<double> patElececalEnergy_;
  std::vector<double> patElechasMatchConv_;
  std::vector<int> patElecID_; 

};

#endif

