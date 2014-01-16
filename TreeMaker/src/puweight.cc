/*
  Lovedeep Kaur Saini
  Panjab University, 
  Chandigarh
*/
#include <iostream>
#include "DelPanj/TreeMaker/interface/puweight.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

puweight::puweight(std::string name, TTree* tree){
  name_=name;
  tree_=tree;
  SetBranches();

}


puweight::~puweight(){
  delete tree_;
}


void
puweight::Fill(const edm::Event& iEvent){
  Clear();

  if(iEvent.isRealData()) {
    nTrueInt_ = -999999999.;
    nPUVert_  = -999999999;
  }
  else {
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;


    double npT=-1.;
    int npIT=-1.;

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

      int BX = PVI->getBunchCrossing();
//       std::cout << "BX = " << BX << std::endl;
      if(BX == 0) {
	npT = PVI->getTrueNumInteractions();
	//	std::cout << "npT = " << npT << std::endl;
	npIT = PVI->getPU_NumInteractions();
      }
    }

    nTrueInt_ = npT;
    nPUVert_  = npIT;
  }
}

void
puweight::SetBranches(){
  AddBranch(&nTrueInt_,  "nTrueInt");
  AddBranch(&nPUVert_,   "nPUVert");
}



void
puweight::AddBranch(double* x, std::string name){
  std::string brName="PU_"+name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}

void
puweight::AddBranch(int* x, std::string name){
  std::string brName="PU_"+name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}



void 
puweight::Clear(){
  nTrueInt_=-99999.;
  nPUVert_ =-99999.;
}

