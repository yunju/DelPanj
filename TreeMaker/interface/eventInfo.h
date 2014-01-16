#ifndef __eventInfo__
#define __eventInfo__

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
using namespace std;

class eventInfo{

 public:
  eventInfo(std::string name, TTree* tree);
  ~eventInfo();
  void Fill(const edm::Event& iEvent); 
  void Clear();
 private:
  eventInfo(){};//Don't allow user
  void SetBranches();
  void AddBranch(int* x, std::string name);
  void AddBranch(std::vector<double>* vec, std::string name);
  void AddBranch(std::vector<bool>* vec, std::string name);  
  
  TTree *tree_;
  std::string name_;
  int nEvt_;
  int nRun_;
  int nLumiS_;
  int bunchX_;

  int nVtx_;
  std::vector<double> vertexX_;
  std::vector<double> vertexY_;
  std::vector<double> vertexZ_;		
  int HLTDoubleEle_;
  int HLTDoubleMu_;
};

#endif

