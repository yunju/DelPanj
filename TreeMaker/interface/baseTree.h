#ifndef __BASE_TREE_H__
#define __BASE_TREE_H__

/*
Anil Singh
Panjab University

Log: This is a base class to define the common members 
for all the individual object fillers. 

*/

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

typedef edm::ParameterSet pset;

class baseTree{
 public:
  baseTree(std::string identifier, TTree* tree);
 ~baseTree(){delete tree_;};
  //virtual void Fill(const edm::Event& iEvent);
  //virtual void SetBranches();
  //virtual void Clear();
 protected:
  baseTree(){}
  TTree* tree_;
  std::string identifier_; //eg. "pat"/"reco"
  void AddBranch(double* x, std::string name);
  void AddBranch(int* x, std::string name);
  void AddBranch(std::vector<double>* vec, std::string name);
  void AddBranch(std::vector<int>* vec, std::string name);

  };

/*
baseTree::baseTree(std::string identifier, TTree* tree){
  identifier_  = identifier; 
  tree_ = tree;
}

void baseTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),vec);
}


void baseTree::AddBranch(double* x, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}
*/
#endif
