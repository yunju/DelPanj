#include "DelPanj/TreeMaker/interface/patMetTree.h"

patMetTree::patMetTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig){
  tree_=tree; 
  patMetLabel_ = iConfig.getParameter<edm::InputTag>("patMet");
  SetBranches();
}


patMetTree::~patMetTree(){
  delete tree_;
} 


void
patMetTree::Fill(const edm::Event& iEvent){
  Clear();
  edm::Handle<pat::METCollection> patMetHandle;
  if(not iEvent.getByLabel(patMetLabel_,patMetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<patMetLabel_<<std::endl; exit(0);}
  pat::METCollection::const_iterator met=patMetHandle.product()->begin();
  //  patMetSumEt_.push_back(met->sumEt());
  patMetPt_.push_back(met->pt());
  patMetEta_.push_back(met->eta());
  //std::cout<<met->eta()<<std::endl;
  patMetPhi_.push_back(met->phi());
  patMetM_.push_back(met->mass());
} 

void
patMetTree::AddBranch(std::vector<double>* vec, std::string name){
  tree_->Branch(name.c_str(),vec);
}

void 
patMetTree::AddBranch(double* x, std::string name){
  tree_->Branch(name.c_str(),x,(name+"/D").c_str());
}

void 
patMetTree::AddBranch(std::vector<std::string>* vec, std::string name){
  tree_->Branch(name.c_str(),vec);
}

void 
patMetTree::SetBranches(){
  AddBranch(&patMetSumEt_,"patMetSumEt_");
  AddBranch(&patMetPt_, "patMetPt_");
  AddBranch(&patMetEta_, "patMetEta_");
  AddBranch(&patMetPhi_, "patMetPhi_");
  AddBranch(&patMetM_, "patMetM_");
}


void
patMetTree::Clear(){
patMetSumEt_.clear();
patMetPt_.clear();
patMetEta_.clear();
patMetPhi_.clear();
patMetM_.clear();
}
