#ifndef __MET_TREE_H_
#define __MET_TREE_H_
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
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace std;
using namespace edm;



class patMetTree{

 public:
  patMetTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~patMetTree();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
  
 private:
  patMetTree();
  
  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<std::string>*, std::string name);
  
  edm::InputTag patMetLabel_;

  TTree* tree_;
  std::vector<std::string> patMetType_; 
  std::vector<double> patMetSumEt_;
  std::vector<double> patMetPt_;
  std::vector<double> patMetEta_;
  std::vector<double> patMetPhi_;
  std::vector<double> patMetM_;
};

#endif

