#ifndef __ZZ_TREE_HH_
#define __ZZ_TREE_HH_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include "TTree.h"
#include "TLorentzVector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/muSelector.h"
#include "DelPanj/TreeMaker/interface/eSelector.h"
#include "pandolf/BTagSFUtil/src/BTagSFUtil.h" 

#include "DataFormats/PatCandidates/interface/Jet.h"


using namespace std;
using namespace edm;
using namespace reco;

class ZZTree{


 public:
  ZZTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~ZZTree();
  void Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup);
  void SetBranches();
  void Clear();

 private:
  ZZTree();

  BTagSFUtil* btsfutiljp;
  void AddBranch(double* x, std::string name);
  void AddBranch(int* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<int>*, std::string name);
  void AddBranchArray(const int arraySize, double* x, std::string name);

  bool passLooseJetID(const pat::Jet* recJet);

  edm::InputTag hzzeejj_;
  edm::InputTag hzzmmjj_;
  edm::InputTag eleRhoIsoInputTag_;
  edm::InputTag muoRhoIsoInputTag_;
  eSelector e2012ID_;
  muSelector mu2012ID_;

  TTree* tree_;

  double eleRho_;
  double muoRho_;
  double metSig_;
 
  std::vector<double> higgsPt_;
  std::vector<double> higgsEta_;
  std::vector<double> higgsPhi_;
  std::vector<double> higgsM_;
  std::vector<double> higgsMRefit_;

  std::vector<double> zllPt_;
  std::vector<double> zllEta_;
  std::vector<double> zllPhi_;
  std::vector<double> zllM_;
  std::vector<double> zlldR_; // deltaR between two leptons

  std::vector<double> zjjPt_;
  std::vector<double> zjjEta_;
  std::vector<double> zjjPhi_;
  std::vector<double> zjjM_;
  std::vector<double> zjjMRefit_;
  std::vector<double> zjjdR_; // deltaR between two jets   

  std::vector<int>    HjetIndex_;
  std::vector<int>    HjetHiggsIndex_;
  std::vector<double> HjetE_;
  std::vector<double> HjetPt_;
  std::vector<double> HjetEta_;
  std::vector<double> HjetPhi_;
  
  std::vector<int>    LeptonsIndex_;
  std::vector<int>    LeptonsHiggsIndex_;
  std::vector<double> LeptonsE_;
  std::vector<double> LeptonsPt_;
  std::vector<double> LeptonsEta_;
  std::vector<double> LeptonsPhi_;
  std::vector<double> LeptonsType_;

  // input for the angular likelihood 
  std::vector<double> heliLD_;
  std::vector<double> costhetaNT1_;
  std::vector<double> costhetaNT2_;
  std::vector<double> phiNT_;
  std::vector<double> phiNT1_;
  std::vector<double> costhetastarNT_;

  // refitted
  std::vector<double> heliLDRefit_;
  std::vector<double> costhetaNT1Refit_;
  std::vector<double> costhetaNT2Refit_;
  std::vector<double> phiNTRefit_;
  std::vector<double> phiNT1Refit_;
  std::vector<double> costhetastarNTRefit_;

  std::vector<int>    nBTags_;
  std::vector<int>    lepType_;
  std::vector<int>    passBit_;

  /*


  std::vector<double> eID01;
  std::vector<double> eID02;
  std::vector<double> eID03;
  std::vector<double> eID04;
  std::vector<double> eID05;
  std::vector<double> eID06;
  std::vector<double> eID07;
  std::vector<double> eID08;
  std::vector<double> eID09;
  std::vector<double> eID10;
  std::vector<double> eID11;
  std::vector<double> eID12;
  std::vector<double> eID13;
  std::vector<double> eID14;
  std::vector<double> eID15;
  std::vector<double> eID16;
  std::vector<double> eID17;


  std::vector<double> muID01;
  std::vector<double> muID02;
  std::vector<double> muID03;
  std::vector<double> muID04;
  std::vector<double> muID05;
  std::vector<double> muID06;
  std::vector<double> muID07;
  std::vector<double> muID08;
  std::vector<double> muID09;
  std::vector<double> muID10;
  std::vector<double> muID11;
  std::vector<double> muID12;
  std::vector<double> muID13;
  std::vector<double> muID14;
  std::vector<double> muID15;

  */
  
};

#endif

