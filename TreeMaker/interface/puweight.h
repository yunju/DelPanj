#ifndef __puweight__
#define __puweight__

#include <memory>
#include <string>
#include <iostream>
#include "TTree.h" 
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


class puweight{

 public:
  puweight(std::string name, TTree* tree);
  ~puweight();
  void Fill(const edm::Event& iEvent); 
  void Clear();
 private:
  puweight(){};//Don't allow user
  void SetBranches();
  void AddBranch(double* x, std::string name);
  void AddBranch(int* x, std::string name);
    
  TTree *tree_;
  std::string name_;

  double nTrueInt_;
  int    nPUVert_;

};

#endif

