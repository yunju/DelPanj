#ifndef __MUON_TREE_H_
#define __MUON_TREE_H_
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
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DelPanj/TreeMaker/interface/muSelector.h"
using namespace std;
using namespace edm;

class patMuonTree{


 public:
  patMuonTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~patMuonTree();
  void Fill(const edm::Event& iEvent);
  double IsoDeposit(const pat::Muon Mu,std::string type,double radius);
  void SetBranches();
  void Clear();

 private:
//TTree* tree_;
  //Dont Allow User to Call the Default Constructor.
  patMuonTree();


  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<std::string>*, std::string name);
  void AddBranch(std::vector<int>*, std::string name);

  edm::InputTag patMuonLabel_;
  edm::InputTag beamSpotLabel_;
muSelector mu2012ID_;

  int Nmuons;
  std::vector<int> patMuonType;
 //pt, eta, phi, M : Enough to caluclate
  //px,py,pz etc.==>No need to store later
  TTree* tree_;
  double nmu; 
  std::vector<std::string> patMuonType_; 
  std::vector<double> patMuonPt_;
  std::vector<double> patMuonEta_;
  std::vector<double> patMuonPhi_;
  std::vector<double> patMuonM_;
  std::vector<double> patMuonSumTrkPt_;
  std::vector<double> patMuonTrkIso_;
  std::vector<double> patMuonHcalIso_;
  std::vector<double> patMuonEcalIso_;
  std::vector<double> patMuonCharge_;
  std::vector<double> patMuonNumChambers_;
  std::vector<double> patMuonNumMatches_;
  std::vector<double> patMuonStationMask_;
  std::vector<double> patMuonNumSegments_;
  
  std::vector<double> patMuonQoverP_;
  std::vector<double> patMuonTheta_;
  std::vector<double> patMuonLambda_;
  std::vector<double> patMuonD0_;
  std::vector<double> patMuonDsz_;
  std::vector<double> patMuonDxyBS_;
  std::vector<double> patMuonDzBS_;
  std::vector<double> patMuonDszBS_;
  std::vector<double> patMuonVx_;
  std::vector<double> patMuonVy_;
  std::vector<double> patMuonVz_;

  std::vector<double> patMuonGenMatchPt_;
  std::vector<double> patMuonGenMatchEta_;
  std::vector<double> patMuonGenMatchPhi_;
  std::vector<double> patMuonGenMotherId_;
  std::vector<double> patMuonGenMatchCharge_;

  //particle based isolation.

  std::vector<double> patMuonChHadSumPt03_;
  std::vector<double> patMuonNeHadSumPt03_;
  std::vector<double> patMuonGamSumPt03_;


  std::vector<double> patMuonChHadSumPt05_;
  std::vector<double> patMuonNeHadSumPt05_;
  std::vector<double> patMuonGamSumPt05_;
 //YunJu's variable 
  std::vector<int> patMuonID_;

  std::vector<int> isGloMu_;
  std::vector<int> isSTDMu_;
  std::vector<int> isTraMu_;
  
  std::vector<double> patMuonChi2Ndoff_;
  std::vector<double> patMuonNhits_;
  std::vector<double> patMuonNMatStat_;
  std::vector<double> patMuondB_;
//  std::vector<double> patMuonDxy_;
  std::vector<double> patMuonDz_;
  std::vector<double> patMuonNumPix_;
  std::vector<double> patMuonNumTlayer_;
//YJ for  muon pfiso 
  std::vector<double> patMuonChHadSumPt04_;
  std::vector<double> patMuonNeHadSumPt04_;
  std::vector<double> patMuonGamSumPt04_;
  std::vector<double> patMuonPUSumPt04_;

};

#endif

