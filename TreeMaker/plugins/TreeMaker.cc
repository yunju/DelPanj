// -*- C++ -*-
//
// Package:    TreeMaker
// Class:      TreeMaker
// Original Author:  Anil Singh, Ashok Kumar
//                   Panjab University, Delhi University
//         Created:  Tue Jul  6 21:04:59 CEST 2010
//


// system include files
#include <memory>
#include <string>

#include "DelPanj/TreeMaker/interface/TreeMaker.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/plugins/EventCountProducer.cc"
TreeMaker::TreeMaker(const edm::ParameterSet& iConfig)

{
   fillPUweightInfo_=0;
   fillEventInfo_=0;
   fillGenInfo_=0;
   fillTrigInfo_=0;
   fillMuonInfo_=0;
   fillElecInfo_=0;
   fillElecIsoInfo_ = 0;
   fillJetInfo_=0;
   fillMetInfo_=0;
   fillPhotInfo_=0; 
   fillZZInfo_  = 0;
   fillYJHiggInfo_=0; 
   fillYJTestSkim_=0;
   fillJecYJ_=0;   

   

   fillPUweightInfo_ = iConfig.getParameter<bool>("fillPUweightInfo_");
   fillEventInfo_ = iConfig.getParameter<bool>("fillEventInfo_");
   fillGenInfo_   = iConfig.getParameter<bool>("fillGenInfo_");
   fillMuonInfo_  = iConfig.getParameter<bool>("fillMuonInfo_");
   fillElecInfo_  = iConfig.getParameter<bool>("fillElecInfo_");
   fillElecIsoInfo_ = iConfig.getParameter<bool>("fillElecIsoInfo_");
   fillJetInfo_   = iConfig.getParameter<bool>("fillJetInfo_");
   fillMetInfo_   = iConfig.getParameter<bool>("fillMetInfo_");
   fillTrigInfo_  = iConfig.getParameter<bool>("fillTrigInfo_");
   fillPhotInfo_  = iConfig.getParameter<bool>("fillPhotInfo_");

   fillZZInfo_     = iConfig.getParameter<bool>("fillZZInfo_");
   fillYJHiggInfo_ = iConfig.getParameter<bool>("fillYJHiggInfo_");
   fillYJTestSkim_= iConfig.getParameter<bool>("fillYJTestSkim_");
   fillJecYJ_ =iConfig.getParameter<bool>("fillJecYJ_");
   //outFileName_= iConfig.getParameter<std::string>("outFileName");
   edm::Service<TFileService> fs;

   //file = new TFile(outFileName_.c_str(),"recreate");
;
   //tree_ = new TTree("tree","tree");

   tree_ = fs->make<TTree>("tree","tree");
   if( fillPUweightInfo_) puweight_=new puweight("puweight",tree_);
   if( fillEventInfo_) eventInfo_=new eventInfo("eventInfo",tree_); 
   if( fillGenInfo_)   genInfoTree_ = new genInfoTree("gen",tree_,iConfig);
   if( fillMuonInfo_)  patMuTree_= new patMuonTree("mu",tree_,iConfig);
   if( fillElecInfo_)  patElecTree_ = new patElecTree("ele",tree_,iConfig);
   if( fillElecIsoInfo_)  patElecIsoTree_ = new patElecIsoTree("eleIso",tree_,iConfig);
   if( fillMetInfo_)   patMetTree_= new patMetTree("met",tree_,iConfig);
   if( fillJetInfo_)   jetTree_=new jetTree("Jet",tree_,iConfig);
   if( fillTrigInfo_)  patHltTree_ = new patHltTree("Hlt",tree_); 
   if( fillPhotInfo_)  photonTree_ = new photonTree("patPhoton", tree_, iConfig); 
   if( fillZZInfo_)    ZZTree_ = new ZZTree("zz",tree_,iConfig);
   if( fillYJHiggInfo_) YJHiggTree_ = new YJHiggTree("YJHigg",tree_,iConfig); 
   if( fillYJTestSkim_) YJHiggTestskim_ = new YJHiggTestskim("YJTestSkim",tree_,iConfig);
   if( fillJecYJ_) YJJECHiggTree_ = new YJJECHiggTree ("YJJec",tree_,iConfig);
}


TreeMaker::~TreeMaker()
{
  //delete tree_;
}

void
TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  if( fillPUweightInfo_) puweight_  ->Fill(iEvent);
  if( fillEventInfo_) eventInfo_   ->Fill(iEvent);
  if( fillGenInfo_)   genInfoTree_ ->Fill(iEvent);
  if( fillElecInfo_)  patElecTree_ ->Fill(iEvent);
  if( fillElecIsoInfo_) patElecIsoTree_->Fill(iEvent);
  if( fillMuonInfo_)  patMuTree_   ->Fill(iEvent);
  if( fillJetInfo_)   jetTree_     ->Fill(iEvent, iSetup);
  if( fillMetInfo_)   patMetTree_  ->Fill(iEvent);
  if( fillTrigInfo_)  patHltTree_  ->Fill(iEvent);
  if( fillPhotInfo_)  photonTree_  ->Fill(iEvent);
  if( fillZZInfo_)    ZZTree_->Fill(iEvent, iSetup);
  if( fillYJHiggInfo_) 
  {
     YJHiggTree_->Fill(iEvent, iSetup);
          
  }
  if( fillYJTestSkim_)    YJHiggTestskim_->Fill(iEvent, iSetup);
  if(fillJecYJ_)
  {YJJECHiggTree_->Fill(iEvent,iSetup);}



  tree_->Fill();
}


/*
void TreeMaker::endLuminosityBlock(const edm::LuminosityBlock & lumi, const EventSetup & setup) {
// Total number of events is the sum of the events in each of these luminosity blocks
Handle <int>nEventsTotalCounter;
lumi.getByLabel("prePathCounter", nEventsTotalCounter);
cout<<" prepath "<<nEventsTotalCounter<<endl ;

//edm::Handle nEventsFilteredCounter;
//lumi.getByLabel("postPathCounter", nEventsFilteredCounter);
//cout<<" prepath "<<nEventsFilteredCounter<<endl ; 
}
*/









void
TreeMaker::beginJob(){
}




void
TreeMaker::endJob() {
//  file->cd();
//  file->Write();
//  file->Close();
  if( fillYJHiggInfo_){

   cout<<"Report from YJHiggTree: "<<endl;
  cout<<"   - Number of processed events: "<<YJHiggTree_->_nEvents<<endl;
  cout<<"   - Number of rejected events due to trigger: "<<YJHiggTree_->_nFailTrig<<endl;
   cout<<"   - Number of rejected events due to Mets10: "<<YJHiggTree_->_nFailMet<<endl;
  cout<<"   - Number of rejected events due to preselection: "<<YJHiggTree_->_nRejected<<endl;
  cout<<"   - Number of passed events: "<<YJHiggTree_->_nPassed<< "; Num of SB:"<<YJHiggTree_->_nSB <<" Num  total jetevt  save: "<<YJHiggTree_->_nJetEvtSave<<endl;
  cout<<"   - From Dimu trigger: "<<YJHiggTree_->_nPFinalHLTDoubleMu<<" From Diele trigger: "<<YJHiggTree_->_nPFinalHLTDoubleEle<<" From Mu17TkMu8 trigger: "<<YJHiggTree_->_nPFinalHLTMu17TkMu8<<endl;
}

  if( fillJecYJ_){

   cout<<"Report from YJJec: "<<endl;
  cout<<"   - Number of processed events: "<<YJJECHiggTree_->_nEvents<<endl;
  cout<<"   - Number of rejected events due to trigger: "<<YJJECHiggTree_->_nFailTrig<<endl;
   cout<<"   - Number of rejected events due to Mets10: "<<YJJECHiggTree_->_nFailMet<<endl;
  cout<<"   - Number of rejected events due to preselection: "<<YJJECHiggTree_->_nRejected<<endl;
  cout<<"   - Number of passed events: "<<YJJECHiggTree_->_nPassed<< "; Num of SB:"<<YJJECHiggTree_->_nSB <<" Num  total jetevt  save: "<<YJJECHiggTree_->_nJetEvtSave<<endl;
  cout<<"   - Pass final selection(SR) from Dimu trigger: "<<YJJECHiggTree_->_nPFinalHLTDoubleMu<<" From Diele trigger: "<<YJJECHiggTree_->_nPFinalHLTDoubleEle<<" From Mu17TkMu8 trigger: "<<YJJECHiggTree_->_nPFinalHLTMu17TkMu8<<endl;
}







}

DEFINE_FWK_MODULE(TreeMaker);
