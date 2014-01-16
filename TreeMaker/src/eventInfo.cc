/*
Lovedeep Kaur Saini
Panjab University, 
Chandigarh
*/
#include <iostream>
#include <cstdlib>
#include "DelPanj/TreeMaker/interface/eventInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

eventInfo::eventInfo(std::string name, TTree* tree){
  name_=name;
  tree_=tree;
  SetBranches();

}


eventInfo::~eventInfo(){
  delete tree_;
}


void
eventInfo::Fill(const edm::Event& iEvent){
  Clear();

  nEvt_   = iEvent.id().event();
  nRun_   = iEvent.id().run();
  nLumiS_ = iEvent.luminosityBlock();
  bunchX_ = iEvent.bunchCrossing();
  
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);
  for(unsigned int ind=0;ind<recVtxs->size();ind++) {
    if (!((*recVtxs)[ind].isFake()) && ((*recVtxs)[ind].ndof()>4) 
	&& (fabs((*recVtxs)[ind].z())<=24.0) &&  
	((*recVtxs)[ind].position().Rho()<=2.0) ) {   
      nVtx_++;
      vertexX_.push_back((*recVtxs)[ind].x());
      vertexY_.push_back((*recVtxs)[ind].y());
      vertexZ_.push_back((*recVtxs)[ind].z());
    }
  }
   bool isData = iEvent.isRealData();  
   if(isData)
   {  
     edm::Handle<edm::TriggerResults> hltresults;
     iEvent.getByLabel(edm::InputTag("TriggerResults", "", "HLT"),hltresults);
     edm::TriggerNames TrigNames = iEvent.triggerNames(*hltresults);

     for ( size_t itr = 0; itr < hltresults->size(); ++itr ) {
     std::string passedPathName = TrigNames.triggerName(itr);
     //std::cout<<nEvt_ <<" Trigger names :"<< passedPathName<<std::endl;
     if (passedPathName.find("HLT_Mu17_Mu8_v")!=std::string::npos && hltresults->accept(itr) ) 
     {
      //std::cout<<"fire Mu"<<std::endl;
      HLTDoubleMu_=1;
     }
     else if ((passedPathName.find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL")!=std::string::npos && hltresults->accept(itr)))
     {
        HLTDoubleEle_=1;
       //  std::cout<<"fire Elec"<<std::endl;      
      //    cin.get(); 
     }


 
    }
   }
    else
    {
      //std::cout<<"RUN mc"<<std::endl;
          
    }
    //cin.get();

}

void
eventInfo::SetBranches(){
  AddBranch(&nEvt_,"EventNum");
  AddBranch(&nRun_,  "RunNum");
  AddBranch(&nLumiS_, "LumiSection");
  AddBranch(&bunchX_, "BunchXing");
  AddBranch(&nVtx_, "NumVtx");
  AddBranch(&HLTDoubleEle_, "HLTDoubleEle");
  AddBranch(&HLTDoubleMu_, "HLTDoubleMu");
  AddBranch(&vertexX_, "VertexX_");
  AddBranch(&vertexY_, "VertexY_");
  AddBranch(&vertexZ_, "VertexZ_");
}

void
eventInfo::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName = "EvtInfo_"+name;
  tree_->Branch(brName.c_str(),vec);
}

void
eventInfo::AddBranch(std::vector<bool>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}



void 
eventInfo::AddBranch(int* x, std::string name){
  std::string brName="EvtInfo_"+name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}

void 
eventInfo::Clear(){
  nEvt_   = -99999;
  nRun_   = -99999;
  nLumiS_ = -99999;
  bunchX_ = -99999;
  nVtx_ = 0;
  vertexX_.clear();
  vertexY_.clear();
  vertexZ_.clear();
  HLTDoubleEle_=-999;
  HLTDoubleMu_=-999; 
}

