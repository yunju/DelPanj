#ifndef  TREE_MAKER_H
#define  TREE_MAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DelPanj/TreeMaker/interface/puweight.h"
#include "DelPanj/TreeMaker/interface/eventInfo.h"
#include "DelPanj/TreeMaker/interface/genInfoTree.h"
#include "DelPanj/TreeMaker/interface/patMuonTree.h"
#include "DelPanj/TreeMaker/interface/patElecTree.h"
#include "DelPanj/TreeMaker/interface/patMetTree.h"
#include "DelPanj/TreeMaker/interface/patHltTree.h"
#include "DelPanj/TreeMaker/interface/jetTree.h"
#include "DelPanj/TreeMaker/interface/photonTree.h"
#include "DelPanj/TreeMaker/interface/patElecIsoTree.h"
#include "DelPanj/TreeMaker/interface/ZZTree.hh"
#include "DelPanj/TreeMaker/interface/YJHiggTree.hh"
#include "DelPanj/TreeMaker/interface/YJHiggTestskim.hh"
#include "DelPanj/TreeMaker/interface/YJJECHiggTree.hh"
#include "TTree.h"
#include "TFile.h"

class TreeMaker : public edm::EDAnalyzer {
   public:
      explicit TreeMaker(const edm::ParameterSet&);
      ~TreeMaker();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup& );
      virtual void endJob() ;
      TFile* file;
      TTree* tree_;
      std::string outFileName_ ;

      bool fillPUweightInfo_; 
      bool fillEventInfo_;

      bool fillGenInfo_;
      bool fillTrigInfo_;
      bool fillMuonInfo_;
      bool fillElecInfo_;
      bool fillElecIsoInfo_;
      bool fillMetInfo_;
      bool fillJetInfo_;
      bool fillPhotInfo_;
      bool fillZZInfo_;
      bool fillYJHiggInfo_;
      bool fillYJTestSkim_;
      bool fillJecYJ_;

      puweight *puweight_;
      eventInfo   *eventInfo_;
      genInfoTree *genInfoTree_;    
      patMuonTree *patMuTree_;
      patElecTree *patElecTree_;
      patElecIsoTree * patElecIsoTree_;
      patMetTree  *patMetTree_;
      jetTree     *jetTree_;
      photonTree  *photonTree_;
      patHltTree  *patHltTree_;      
      ZZTree      *ZZTree_;
      YJHiggTree  *YJHiggTree_;
      YJHiggTestskim *YJHiggTestskim_; 
      YJJECHiggTree *YJJECHiggTree_;    
};


#endif
