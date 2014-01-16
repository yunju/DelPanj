#include "DelPanj/TreeMaker/interface/baseTree.h"

baseTree::baseTree(std::string identifier, TTree* tree){
  identifier_  = identifier; 
  tree_ = tree;
}

void baseTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),vec);
}


void baseTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),vec);
}


void baseTree::AddBranch(double* x, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}


void baseTree::AddBranch(int* x, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}



