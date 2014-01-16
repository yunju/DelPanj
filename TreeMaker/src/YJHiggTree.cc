/// @file
/// File containing the definition of the methods associated to the class.
///
///For VBF analysis YunJu Lu 20130505 


#include "DelPanj/TreeMaker/interface/YJHiggTree.hh"
#include "DelPanj/TreeMaker/interface/cutvalues.h"
//#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "DelPanj/TreeMaker/interface/MuonEffectiveArea.h"

// system include files

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/CandUtils/interface/CenterOfMassBooster.h"
#include "PhysicsTools/CandUtils/interface/Booster.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// ROOT classes
#include <TMath.h>
#include <TVector3.h>
#include <algorithm>
#include <fstream>
#include <Math/VectorUtil.h>
#include <TLegend.h>
#include <TCanvas.h>


// for Btagging scale factors
#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include <time.h>
//for trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//angles
//#include "DelPanj/TreeMaker/plugins/RooBkgd2L2JV2.h"
//#include "DelPanj/TreeMaker/plugins/RooSpinZero5DV2.h"

#include "DelPanj/TreeMaker/interface/Helicity.h"
#include "DelPanj/TreeMaker/interface/HelicityLikelihoodDiscriminant.h"


typedef std::vector< edm::Handle< edm::ValueMap<double> > >             
IsoDepositVals;


void
YJHiggTree::AddBranch(int* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}


void
YJHiggTree::AddBranch(double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}

//---------------------------------------------------
//---------------------------------------------------
void
YJHiggTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

void
YJHiggTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}


//---------------------------------------------------
//---------------------------------------------------
void 
YJHiggTree::AddBranchArray(const int arraySize, double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,Form("%s[%d]/D",brName.data(),arraySize));
}

//---------------------------------------------------------------
//---------------------------------------------------------------

YJHiggTree::YJHiggTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  e2012ID_ ( iConfig.getParameter<edm::ParameterSet>("e2012IDSet")),
  mu2012ID_ ( iConfig.getParameter<edm::ParameterSet>("mu2012IDSet")),
  hzzeejj_(iConfig.getParameter<edm::InputTag>("hzzeejjTag")),
  hzzmmjj_ (iConfig.getParameter<edm::InputTag>("hzzmmjjTag")),
  eleRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("eleRhoIso")),
  muoRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("muoRhoIso"))//,
{
  tree_=tree; 
  SetBranches();

  
  // the second argument is the random seed, any reason to set it 
  // differently or the same for the 3 taggers
  srand ( time(NULL) );
  int seed = rand(); //712687782, 727743360
  std::cout << "seed = " << seed << std::endl;
  btsfutiljp = new BTagSFUtil("JP", seed);
  OnlyMuon=0;
  _nEvents=0;
  _nFailTrig=0;
  _nCandidates=0;
  _nRejected=0;
   _nPassed=0;
   _nSB=0;
   _nJetEvtSave=0;
  _nPFinalHLTDoubleMu=0;
 _nPFinalHLTDoubleEle=0;
 _nPFinalHLTMu17TkMu8=0;
 _nFailMet=0;



/*
   reader->AddVariable( "MaxEtajj", &MaxEtajj );
   reader->AddVariable( "MaxMjj", &MaxMjj );

   reader->AddVariable("Hpt",&Hpt);
   reader->AddVariable("numJets",&numJets);
   reader->AddVariable( "costhetaNT1", &costhetaNT1 );
   reader->AddVariable( "costhetaNT2", &costhetaNT2 );
   reader->AddVariable( "costhetastarNT",&costhetastarNT );
   reader->AddVariable( "phiNT1", &phiNT1 );
   reader->AddVariable( "phiNT2", &phiNT2 );
   reader->BookMVA( "BDT","/home/yunju/HZZ/TMVAClassification_BDT_VBF230_0_A_GGH230_0DY3JDY4JDY2J.weights.xml");
*/

}
void YJHiggTree::endJob (void)
{
   cout<<"Report from YJ: "<<endl;
  cout<<"   - Number of processed events: "<<_nEvents<<endl;
  cout<<"   - Number of rejected events due to trigger: "<<_nFailTrig<<endl;
  cout<<"   - Number of rejected events due to preselection: "<<_nRejected<<endl;
  cout<<"   - Number of passed events: "<<_nPassed<<endl;

}

YJHiggTree::~YJHiggTree()
{


  delete tree_;
  delete btsfutiljp;
}



// ------------ method called to for each event  ------------
void YJHiggTree::Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup)
{
  _nEvents++;
  Clear();
  int eventNum = iEvent.id().event();
  
  //============================================================================
  // 
  //       OBTAIN EVENT-LEVEL VARIABLES
  //  
  //============================================================================
  
  bool isData = iEvent.isRealData();
   HLTDoubleMu_=-999;
   HLTDoubleEle_=-999;
   HLTMu17TkMu8_=-999;
//   if(isData)
//   {  
     edm::Handle<edm::TriggerResults> hltresults;
     iEvent.getByLabel(edm::InputTag("TriggerResults", "", "HLT"),hltresults);
     edm::TriggerNames TrigNames = iEvent.triggerNames(*hltresults);

     for ( size_t itr = 0; itr < hltresults->size(); ++itr ) {
     std::string passedPathName = TrigNames.triggerName(itr);

     //if (passedPathName.find("HLT_Mu17_Mu8_v")!=std::string::npos )std::cout<<iEvent.id().event() <<" Trigger names :"<< passedPathName<<" "<<hltresults->accept(itr)<<std::endl;
    // if ((passedPathName.find("HLT_Mu17_Mu8_v")!=std::string::npos || passedPathName.find("HLT_Mu17_TkMu8")!=std::string::npos) && hltresults->accept(itr) ) 
      if (passedPathName.find("HLT_Mu17_Mu8_v")!=std::string::npos && hltresults->accept(itr) )
      {
     // std::cout<<"fire Mu"<<std::endl;
      HLTDoubleMu_=1;
     }
    else if ((passedPathName.find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")!=std::string::npos && hltresults->accept(itr)))
     {
        HLTDoubleEle_=1;
       //  std::cout<<"fire Elec"<<std::endl;      
      //    cin.get(); 
     }
     else if ((passedPathName.find("HLT_Mu17_TkMu8")!=std::string::npos && hltresults->accept(itr)))
     {
        HLTMu17TkMu8_=1;
       //  std::cout<<"fire Elec"<<std::endl;      
      //    cin.get(); 
     }
    

 
    }//Loop trigger bit
//  }//is data
    
//  cout<<"trigger "<<passDiMuTrig<<endl;


//cin.get();
  bool passTrig=false;
//  if(!isData||HLTDoubleMu_==1||HLTDoubleEle_==1||HLTMu17TkMu8_==1){ passTrig =true;}
  if(HLTDoubleMu_==1||HLTDoubleEle_==1||HLTMu17TkMu8_==1){ passTrig =true;}
  else{_nFailTrig++;}  

genHmass_=-999.0;
//save gen higgs mass
if(!isData) {
     genHmass_=-1.0;
     //get the generated higgs mass
     Handle<std::vector<reco::GenParticle>> gParticles;
     iEvent.getByLabel("genParticles", gParticles);
    
      for( unsigned int k = 0; k < gParticles->size(); k++ )
      {
         const reco::GenParticle & genP = (*gParticles)[ k ];
         if (fabs(genP.pdgId())==25 && genP.status()==3)
         {
          genHmass_=genP.mass();
         }//if it's gen H
      }//genparticle loop
}// if is MC




//START TO CHECK HIGG CAND
/*
   bool hasHMu=false;
  //for YJ to compare
    Handle<std::vector<pat::CompositeCandidate> > hzzlljj;
    iEvent.getByLabel("hzzmmjj","h",hzzlljj); 
    cout<<"################################### Start Event "<<iEvent.id().event()<<endl; 
    cout<<"----[1]Higgs loop "<<hzzmmjj_<<endl;
    for(unsigned i=0; i<hzzlljj->size(); i++)
    {
       const pat::CompositeCandidate & h = (*hzzlljj)[i];
        for(unsigned int imuo=0; imuo < 2; imuo++){
	 
	  const pat::Muon* myMuo
	    = dynamic_cast<const pat::Muon*>(h.daughter(LEPZ)->daughter(imuo)->masterClone().get());
//	  cout<<iEvent.id().event()<<" "<<i<<"th higgs; Mu pt: "<<myMuo->pt()<<endl;
          std::map<std::string, bool> Pass = mu2012ID_.CutRecord(*myMuo);
	  int passOrNot = PassAll(Pass);
	  cout<<iEvent.id().event()<<" "<<i<<"th higgs; Mu pt: "<<myMuo->pt()<<" "<<myMuo->eta()<<" pass:"<<passOrNot<<endl;
         if(passOrNot==0)continue; // 2012 tight muon ID	  
	  hasHMu=true; 

          
          	  	  
	} // end of loop over muon

    }//end for higgs cand loop
    
    //start jet loop
    edm::Handle<std::vector<pat::Jet> > JetHandle;
    iEvent.getByLabel("customPFJetsNoPUSub",JetHandle);
   // iEvent.getByLabel("cleanPatJetsNoPUIsoLept",JetHandle);
    const std::vector<pat::Jet>* jets = JetHandle.product();
    
    int nj1=0;
    int nGoodJets=0;
    for(std::vector<pat::Jet>::const_iterator jet1 =jets->begin();jet1!=jets->end();jet1++)
    {
      nj1++;
      if(jet1->pt()<30) continue;
      if(fabs(jet1->eta())>2.5) continue;
      if(!passLooseJetID(&*jet1)) continue;
      if(jet1->userFloat("puBeta")<0.2) continue;
      nGoodJets++;
      TLorentzVector Vj1;
      Vj1.SetPtEtaPhiE(jet1->pt(),jet1->eta(),jet1->phi(),jet1->energy());
      int nj2=0;
      cout<<"jet index : "<<nj1<<" pt : "<<jet1->pt()<<" eta : "<<jet1->eta()<<endl; 
      for(std::vector<pat::Jet>::const_iterator jet2 =jets->begin();jet2!=jets->end();jet2++)  
      {
        nj2++;
        if(jet2->pt()<30) continue;
        if(fabs(jet2->eta())>2.5) continue;
        if(!passLooseJetID(&*jet2)) continue;
        if(jet2->userFloat("puBeta")<0.2) continue;
        TLorentzVector Vj2;
        Vj2.SetPtEtaPhiE(jet2->pt(),jet2->eta(),jet2->phi(),jet2->energy());
        TLorentzVector Vjj;
        Vjj=Vj1+Vj2;
      //  cout<<"jet1 index1 : "<<nj1<<" jet2 index : "<<nj2<<" mjj:"<<Vjj.M()<<endl;
      }


    }
    if(nGoodJets<2) return;
    
    cout<<"----[2]Userdatamuon Loop"<<endl;

    bool hasMMu=false;
    edm::Handle<pat::MuonCollection> patMuonHandle;
    iEvent.getByLabel("userDataSelectedMuons",patMuonHandle) ;
 
    pat::MuonCollection muColl(*(patMuonHandle.product()));
    std::sort(muColl.begin(),muColl.end(),PtGreater());

    pat::MuonCollection::const_iterator mu1;
    pat::MuonCollection::const_iterator mu2;

    bool passMuWin=false;
    int nmu1=0;
    for(mu1=muColl.begin(); mu1!=muColl.end(); mu1++){
    nmu1++;
   // cout<<nmu1<<" ---1st  "<<endl;
    TLorentzVector Vmu1;
    Vmu1.SetPtEtaPhiM(mu1->pt(),mu1->eta(),mu1->phi(),mu1->mass());
    int nmu2=0;
    for(mu2=muColl.begin(); mu2!=muColl.end(); mu2++){

       //cin.get();
       nmu2++;
     //   cout<<nmu2<<" ---2nd  "<<endl;
       if(nmu1>=nmu2) continue;
     //  cout<<nmu1<<" ---  "<<nmu2<<endl;

      TLorentzVector Vmu2;
       Vmu2.SetPtEtaPhiM(mu2->pt(),mu2->eta(),mu2->phi(),mu2->mass());
       TLorentzVector Vmumu;
       Vmumu=Vmu1+Vmu2;
       cout<<"mu 1 index : "<<nmu1<<" mu 2 index : "<<nmu2<<" mll :"<<Vmumu.M()<<endl;
       if(Vmumu.M()>20) passMuWin=true; 
    }
  }
    
   if(!passMuWin) return;

    int nmu=0;
    
    pat::MuonCollection::const_iterator mu;

     for(mu=muColl.begin(); mu!=muColl.end(); mu++)
    {
     nmu++;
     std::map<std::string, bool> Pass = mu2012ID_.CutRecord(*mu);
     int passOrNot = PassAll(Pass);
     if(fabs(mu->eta())>2.4) continue;
     if(mu->pt()<10) continue;
     //if(passOrNot&&fabs(mu->eta())<2.4) cout<<"UserDataMu Pt :"<<iEvent.id().event()<<" "<<" ;pass:"<<passOrNot <<" ;pt:"<<mu->pt()<<" "<<mu->eta()<<" Charge:"<<mu->charge()<<endl;
     cout<<iEvent.id().event() <<" UserDataMu index : "<<nmu<<" "<<" ;pass:"<<passOrNot <<" ;pt:"<<mu->pt()<<" "<<mu->eta()<<" Charge:"<<mu->charge()<<endl;
     if(!passOrNot) continue;
     hasMMu=true;
    }
       


    if(hasMMu==true&&hasHMu==false) OnlyMuon ++;
    
    cout<<"------------------------------End of event----------------------------- "<<OnlyMuon <<endl;
   // cin.get();
*/



//JET HANDLE

    edm::Handle<std::vector<pat::Jet> > JetHandle;
    iEvent.getByLabel("customPFJetsNoPUSub",JetHandle);
//   iEvent.getByLabel("cleanPatJetsNoPUIsoLept",JetHandle);
    const std::vector<pat::Jet>* jets = JetHandle.product();
//MUON HANDLE    
    edm::Handle<pat::MuonCollection> patMuonHandle;
    iEvent.getByLabel("userDataSelectedMuons",patMuonHandle) ;
    pat::MuonCollection muColl(*(patMuonHandle.product()));
    std::sort(muColl.begin(),muColl.end(),PtGreater());
    
//ELECTRON HANDLE AND SET RHO
    Handle<std::vector<pat::Electron> > patElecHandle;
    iEvent.getByLabel("userDataSelectedElectrons",patElecHandle);
    pat::ElectronCollection eleColl(*(patElecHandle.product()));
    std::sort(eleColl.begin(),eleColl.end(),PtGreater());
 

    // rho for electron
    edm::Handle<double> ele_rho_event;
    iEvent.getByLabel(eleRhoIsoInputTag_,ele_rho_event);
    eleRho_ = *(ele_rho_event.product());
    e2012ID_.SetData(isData);
    e2012ID_.SetRho(eleRho_);

//MET

 Handle<float>  metSignificanceHandle;
  iEvent.getByLabel("metInfoProducer","metSignificance", metSignificanceHandle);
  float metsignificance = *metSignificanceHandle;
  metSig_=metsignificance;

//Pile up

    int pileup =0,pileup_true=0;
    if(!isData){
    edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);
      pileup =-1;
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
        // on time bunch
        if(PVI->getBunchCrossing()==0){
          pileup = PVI->getPU_NumInteractions();
          pileup_true  = PVI->getTrueNumInteractions();
          continue;
        }
      }
    }
    n_pileup_=pileup;
    n_pileup_true_=pileup_true;
   //cout<<n_pileup_<<""<<n_pileup_true_<<endl;



   // check jet pt, id, vtx,lep overlap, if satified is a good jet and put it into LepLapJets 
   // check jet eta <4.7 for VBF jet require
   // note that eta< 2.4 and "isTaggable" and "puBeat"(has eta cut) are applied in dijet loops later after saving VBF information 
   
    int nj1=0;
    int nGoodJets=0;
    bool passJetWin= false;
    bool passJetSB= false;
    std:vector<int> LepLapJets; //a good jet loose id ,pt 30 eta4.7, no lepton overlap
    for(std::vector<pat::Jet>::const_iterator jet1 =jets->begin();jet1!=jets->end();jet1++)
    {
      nj1++;
      if(jet1->pt()<30) continue;
      if(fabs(jet1->eta())>4.7) continue;
      if(!passLooseJetID(&*jet1)) continue;
      //if(jet1->userInt("isTaggable")==0) continue;//use charge tracks has |eta|<2.4
      //if(jet1->userFloat("puBeta")<0.2) continue;
    
      //check overlap with muon
         bool Overlapjets=false;
          for( pat::MuonCollection::const_iterator muLap =muColl.begin();muLap!=muColl.end();muLap++)
          {
             if(muLap->pt()<10) continue;
             if(fabs(muLap->eta())>2.4) continue;
             std::map<std::string, bool> PassMu= mu2012ID_.CutRecord(*muLap);
             int passOrNotMu = PassAll(PassMu);
        //cout<<"event: "<<iEvent.id().event()<<" mu: "<<nmu1<<" pt: "<<mu1->pt() <<"id:"<<passOrN    ot<<endl;
        //cin.get();
             if(passOrNotMu==0) continue;
             if(deltaR (muLap->eta(),muLap->phi(),jet1->eta(), jet1->phi() )<0.5 ){Overlapjets=true;}
          }
      //check overlap with ele
          for(pat::ElectronCollection::const_iterator eleLap =eleColl.begin();eleLap!=eleColl.end();eleLap++) 
          {
             if(eleLap->pt()<10) continue;
             if(fabs(eleLap->eta())>2.5) continue;
             if(eleLap->userFloat("cutIDCode") <= 1) continue;
             if(eleLap->userFloat("passTriggerTight") <= 0) continue;
             if(deltaR (eleLap->eta(),eleLap->phi(),jet1->eta(), jet1->phi() )<0.5 ){Overlapjets=true;} 

          }
          if(Overlapjets) continue;
          (LepLapJets.push_back(nj1));
     // cout<<"Good jets"<<nj1<<" : "<<jet1->pt()<<endl; 
      nGoodJets++;
    }

    //check di jet win and maxmjj maxetajj in 1 event
   //find out the highest btag di-jet pair to save 
//cout<<"################################### Start Event "<<iEvent.id().event()<<endl;
//cout<<<<endl;

//Save inform of central jets: Hjet, select jet pair
//Save inform of Mjj sb jets 
  heliLD_=-999;
  costhetaNT1_= -999;
  costhetaNT2_= -999;
  phiNT_= -999;
  phiNT1_= -999;
  costhetastarNT_= -999;
  Hpt_=-999;
  numJets_=-999;
  EvtType_=-999;
  NBtagYJ_=-999;
int HcentralJet1=0;
int HcentralJet2=0;
int nbtag_discriminant=-999;//for signal region  
float zjj_discriminant=999.0;//for signal region




int SBnbtag_discriminant=-999;//for SB region  
float SBzjj_discriminant=999.0;//for SB region
int SBcentralJet1=0;
int SBcentralJet2=0;



double maxmjj=-99.0;
double maxetajj=-99.0;
int njetM1=0;

for(std::vector<pat::Jet>::const_iterator jetM1 =jets->begin();jetM1!=jets->end();jetM1++)
{
    njetM1++;
    bool isM1GoodJet=false;
    for(unsigned int n=0;n<LepLapJets.size();n++)
    {
        if(njetM1==LepLapJets.at(n))
        isM1GoodJet=true; 
        //cout<<" no over jets: "<<LepLapJets.at(n)<<endl;
    }
    if(!isM1GoodJet) continue;
    //cout<<"M1 no over jets: "<<njetM1<<" : "<<jetM1->pt()<<endl;
    int njetM2=0;
    for(std::vector<pat::Jet>::const_iterator jetM2 =jets->begin();jetM2!=jets->end();jetM2++)
    {
       njetM2++;
       if(njetM1>=njetM2) continue; 
       bool isM2GoodJet=false; 
       for(unsigned int n=0;n<LepLapJets.size();n++) 
       {
           if(njetM2==LepLapJets.at(n))
           isM2GoodJet=true;
       }
       if(!isM2GoodJet) continue;
      // cout<<"M2 no over jets: "<<njetM2<<" : "<<jetM2->pt()<<endl;
       TLorentzVector Vj1;
       Vj1.SetPtEtaPhiE(jetM1->pt(),jetM1->eta(),jetM1->phi(),jetM1->energy());     
       TLorentzVector Vj2;
       Vj2.SetPtEtaPhiE(jetM2->pt(),jetM2->eta(),jetM2->phi(),jetM2->energy());
       TLorentzVector Vjj;
       Vjj=Vj1+Vj2;
       //cout<<"Mjj:"<<Vjj.M()<<" etajj:"<<fabs(jetM2->eta()-jetM1->eta())<<endl;
       if(Vjj.M()>maxmjj) maxmjj= Vjj.M();
       if(fabs(jetM2->eta()-jetM1->eta())> maxetajj) maxetajj= fabs(jetM2->eta()-jetM1->eta());
       //check if there is central jets(and pass taggable and puBeta ) 
       if(fabs(jetM1->eta())>2.4||fabs(jetM2->eta())>2.4) continue;
       if(jetM1->userInt("isTaggable")==0||jetM2->userInt("isTaggable")==0) continue;
       if(jetM1->userFloat("puBeta")<0.2||jetM2->userFloat("puBeta")<0.2) continue;
       //get nbtag 
        int nbtags=0,tagbit=0;
        bool jet1_tagged_medium=false;
        bool jet1_tagged_loose=false;
        bool jet2_tagged_medium=false;
        bool jet2_tagged_loose=false;
        if (jetM1->bDiscriminator("jetProbabilityBJetTags")>0.545)jet1_tagged_medium=true;
        if (jetM1->bDiscriminator("jetProbabilityBJetTags")>0.275)jet1_tagged_loose=true;
        if (jetM2->bDiscriminator("jetProbabilityBJetTags")>0.545)jet2_tagged_medium=true;
        if (jetM2->bDiscriminator("jetProbabilityBJetTags")>0.275)jet2_tagged_loose=true;
          //btag category
          if ((jet1_tagged_medium && jet2_tagged_loose)|| (jet2_tagged_medium && jet1_tagged_loose)) nbtags=2;
          else if (jet1_tagged_loose || jet2_tagged_loose) nbtags =1;
          else nbtags=0;
        //for signal region
      if(Vjj.M()<111.0&&Vjj.M()>71.0){         
       passJetWin=true;
//check btag first then check mjj-z 
        //cout<<njetM1<<" jet1batg: "<<jetM1->bDiscriminator("jetProbabilityBJetTags")<<";" <<njetM2<<" jet2batg:"<<jetM2->bDiscriminator("jetProbabilityBJetTags")<<" nbtags:"<<nbtags<<" fabs(Vjj.M()-91.2):"<<fabs(Vjj.M()-91.2)<<endl;
       //cout<<"nbtags: "<<nbtags<<" nbtag_discriminant:"<<nbtag_discriminant<<" zjj_discriminant:"<<zjj_discriminant<<endl;
       if(nbtags>nbtag_discriminant){             
          nbtag_discriminant=nbtags;
          zjj_discriminant=fabs(Vjj.M()-91.2);
          HcentralJet1=njetM1;
          HcentralJet2=njetM2;   
       }
       else if (nbtags==nbtag_discriminant)
       {
             //cout<<"event3 "<<" zjj_discriminant:"<<zjj_discriminant<<endl;
             if (fabs(Vjj.M()-91.2)<zjj_discriminant)
             {
                 
               zjj_discriminant=fabs(Vjj.M()-91.2);
               HcentralJet1=njetM1;
               HcentralJet2=njetM2;                                        
             }
       }
     }//if in signal region 71-111
     else if ((Vjj.M()<130.0&&Vjj.M()>111.0)||(Vjj.M()>60.0&&Vjj.M()<71.0))//for sideband
     {
       passJetSB =true;

       if(nbtags>SBnbtag_discriminant){             
          SBnbtag_discriminant=nbtags;
          SBzjj_discriminant=fabs(Vjj.M()-91.2);
          SBcentralJet1=njetM1;
          SBcentralJet2=njetM2;   
       }
       else if (nbtags==SBnbtag_discriminant)
       {
             //cout<<"event3 "<<" SBzjj_discriminant:"<<SBzjj_discriminant<<endl;
             if (fabs(Vjj.M()-91.2)<SBzjj_discriminant)
             {
                 
               SBzjj_discriminant=fabs(Vjj.M()-91.2);
               SBcentralJet1=njetM1;
               SBcentralJet2=njetM2;                                        
             }
       }


     }//if in SB region  
 
    }//jetm2 loop
}// jetm1 loop



EvtMaxMjj_=maxmjj;
EvtMaxEtajj_=maxetajj;
if(passJetWin) NBtagYJ_=nbtag_discriminant;
if(passJetSB)  NBtagYJ_=SBnbtag_discriminant;

//cout<<" maxmjj: "<<EvtMaxMjj_<<" maxetajj:"<<EvtMaxEtajj_<<endl;
    //start to check muon
    float zll_discriminant=999.0;
    int HcentralLep1=0;    
    int HcentralLep2=0;
    int HLepType=999;
    bool hasMMu=false;

    pat::MuonCollection::const_iterator mu1;
    pat::MuonCollection::const_iterator mu2;

    bool passMuWin=false;
    bool passMuBase =false;
    int nGoodMuL=0;
    int nGoodMuH=0;
    int nmu1=0;
    for(mu1=muColl.begin(); mu1!=muColl.end(); mu1++){
       nmu1++;
       // cout<<nmu1<<" ---1st  "<<endl;
       if(mu1->pt()<20) continue;
       if(fabs(mu1->eta())>2.4) continue;
       //nGoodMuL++; 
       //if(mu1->pt()<40) continue;
       // nGoodMuH++;
       
       std::map<std::string, bool> Pass = mu2012ID_.CutRecord(*mu1);
       int passOrNot = PassAll(Pass);
       //cout<<"event: "<<iEvent.id().event()<<" mu: "<<nmu1<<" pt: "<<mu1->pt() <<"id:"<<passOrNot<<endl;
       //cin.get(); 
       if(passOrNot==0) continue;
       
       TLorentzVector Vmu1;
       Vmu1.SetPtEtaPhiM(mu1->pt(),mu1->eta(),mu1->phi(),mu1->mass());
       int nmu2=0;
       for(mu2=muColl.begin(); mu2!=muColl.end(); mu2++){
          //cin.get();
          nmu2++;
          //   cout<<nmu2<<" ---2nd  "<<endl;
          if(nmu1>=nmu2) continue;
          if(mu2->pt()<20) continue;
          if(mu1->pt()<40&&mu2->pt()<40) continue;
          if(fabs(mu2->eta())>2.4) continue;
          std::map<std::string, bool> Pass2 = mu2012ID_.CutRecord(*mu2); 
          int passOrNot2 = PassAll(Pass2);
          if(!passOrNot2) continue;
          if(mu1->charge()*mu2->charge()==1) continue;  
          
          //jets which is actually muon is clean at jet loop by 0.5 
          //not necessary to do it again 
           /*
          //delta R
          bool Overlapjets=false;
          for(std::vector<pat::Jet>::const_iterator jetLap =jets->begin();jetLap!=jets->end();jetLap++)
          {
           if(jetLap->pt()<30) continue;
           if(fabs(jetLap->eta())>2.4) continue;  
           if(jetLap->userInt("isTaggable")==0) continue;
           if(jetLap->userFloat("puBeta")<0.2) continue;
          //
          cout<<deltaR (mu1->eta(),mu1->phi(),jetLap->eta(), jetLap->phi() )<<endl;
          if(deltaR (mu1->eta(),mu1->phi(),jetLap->eta(), jetLap->phi() )<0.5 || deltaR (mu2->eta(),mu2->phi(),jetLap->eta(), jetLap->phi() )<0.5){Overlapjets=true;}
          }
          cin.get(); 
          //if(Overlapjets) continue;
          */
          //  cout<<nmu1<<" ---  "<<nmu2<<endl;
          passMuBase=true;
          TLorentzVector Vmu2;
          Vmu2.SetPtEtaPhiM(mu2->pt(),mu2->eta(),mu2->phi(),mu2->mass());
          TLorentzVector Vmumu;
          Vmumu=Vmu1+Vmu2;
         // cout<<"mu 1 index : "<<nmu1<<" mu 2 index : "<<nmu2<<" mll :"<<Vmumu.M()<<endl;

          if(Vmumu.M()<20) continue;
          if(Vmumu.M()>106||Vmumu.M()<76) continue;
          
           passMuWin=true;
       // cout<<"mu "<<nmu1<<" : "<<nmu2<<";" <<fabs(Vmumu.M()-91.2)<<endl;
 

         //check Z ll window 
          if (fabs(Vmumu.M()-91.2)<zll_discriminant)
          {
               zll_discriminant=fabs(Vmumu.M()-91.2);
               HcentralLep1=nmu1;
               HcentralLep2=nmu2;
               HLepType=0;
          }

       }//2nd mu loop
    }//1st muon loop
    

     //start to check electrons

    pat::ElectronCollection::const_iterator ele1;
    pat::ElectronCollection::const_iterator ele2;


    bool passEleBase =false;
    
    int nele1=0;
    for(ele1=eleColl.begin(); ele1!=eleColl.end(); ele1++)
    {
       nele1++;
       if(ele1->pt()<20)continue;
       std::map<std::string, bool> Passele1 =  e2012ID_.CutRecord(*ele1);
       if (!PassAll(Passele1)) continue;


       int nele2=0; 
       for(ele2=eleColl.begin(); ele2!=eleColl.end(); ele2++)
       {
         nele2++;
         if(nele1>=nele2) continue;
         if(ele1->charge()*ele2->charge()==1) continue;
         if(ele2->pt()<20)continue;
         if(ele1->pt()<40&&ele2->pt()<40) continue;
         std::map<std::string, bool> Passele2 =  e2012ID_.CutRecord(*ele2);
         if (!PassAll(Passele2)) continue;  
         
         TLorentzVector Vele1;
         Vele1.SetPtEtaPhiM(ele1->pt(),ele1->eta(),ele1->phi(),ele1->mass());
         TLorentzVector Vele2;
         Vele2.SetPtEtaPhiM(ele2->pt(),ele2->eta(),ele2->phi(),ele2->mass());
         TLorentzVector VEleEle;
         VEleEle=Vele1+Vele2;
         if(VEleEle.M()<20) continue;
         if(VEleEle.M()>106||VEleEle.M()<76) continue;
         passEleBase =true; 
       // cout<<"e "<<nele1<<" : "<<nele2<<";" <<fabs(VEleEle.M()-91.2)<<endl; 
         if (fabs(VEleEle.M()-91.2)<zll_discriminant)
         {
               zll_discriminant=fabs(VEleEle.M()-91.2);
               HcentralLep1=nele1;
               HcentralLep2=nele2;
               if(HLepType!=999) {cout<<"found!!"<<endl;cin.get();}
               HLepType=1;
          }
           




       }//2nd ele loop
    }//1st ele loop


  if(!passTrig) return;

  if(metsignificance>10) _nFailMet++;

  if(metsignificance<=10&&(passEleBase||passMuWin)&&(passJetWin||passJetSB))
  {
   
   if(passJetWin) {_nPassed++;EvtType_=1;}
   else if(passJetSB) {_nSB++;EvtType_=2;}
   
   if(HLTDoubleMu_==1&&passJetWin) {_nPFinalHLTDoubleMu++;}
   if(HLTDoubleEle_==1&&passJetWin){_nPFinalHLTDoubleEle++;}
   if(HLTMu17TkMu8_==1&&passJetWin){_nPFinalHLTMu17TkMu8++;}

//LV for angular VB
   TLorentzVector LVJ1;
   TLorentzVector LVJ2;
   bool hasFillLVJ1=false;

//save jet information
bool hassave=false;
int njetSave=0;
for(std::vector<pat::Jet>::const_iterator jetSave =jets->begin();jetSave!=jets->end();jetSave++)
{
    njetSave++;
    bool isSaveGoodJet=false;
    for(unsigned int n=0;n<LepLapJets.size();n++)
    {
        if(njetSave==LepLapJets.at(n))
        isSaveGoodJet=true; 
        //cout<<" no over jets: "<<LepLapJets.at(n)<<endl;
    }
    if(!isSaveGoodJet) continue;
    int isHjet=0;
    int isSBjet=0; 

    if(passJetWin){    
      if(njetSave==HcentralJet1||njetSave==HcentralJet2)
      {
        isHjet=1;
        // cout<<"## "<<iEvent.id().event()<<" njetSave "<<njetSave<<endl;
        //Fill LV for angular vars
        if(!hasFillLVJ1)
        {
           LVJ1.SetPtEtaPhiE(jetSave->pt(),jetSave->eta(),jetSave->phi(),jetSave->energy());
           hasFillLVJ1=true;
        }
        else
        {
           LVJ2.SetPtEtaPhiE(jetSave->pt(),jetSave->eta(),jetSave->phi(),jetSave->energy());
        }
  
      }
    }//if signal region
    else if(passJetSB)
    {
      if(njetSave==SBcentralJet1||njetSave==SBcentralJet2)
      {
         isSBjet=1;
         if(!hasFillLVJ1)
         {
           LVJ1.SetPtEtaPhiE(jetSave->pt(),jetSave->eta(),jetSave->phi(),jetSave->energy());
           hasFillLVJ1=true;
         }
         else
         {
           LVJ2.SetPtEtaPhiE(jetSave->pt(),jetSave->eta(),jetSave->phi(),jetSave->energy());
         }         
      }
    }//if SB region

    hassave =true; 
    JetPt_.push_back(jetSave->pt());
    JetEta_.push_back(jetSave->eta());
    JetPhi_.push_back(jetSave->phi());
    JetM_.push_back(jetSave->mass());
    JetEn_.push_back(jetSave->energy());
    JetJetProb_.push_back(jetSave->bDiscriminator("jetProbabilityBJetTags"));
    JetFromH_.push_back(isHjet);
    JetFromSB_.push_back(isSBjet);    
    JetpuBetaYJ_.push_back(jetSave->userFloat("puBeta"));
}//jetSave loop
if(hassave==true) _nJetEvtSave++;
  //Set Lep LV 
  TLorentzVector LVLep1;
  TLorentzVector LVLep2;
  bool hasFillLVLep1=false;
  int LVLep1Charge=-99;
  if(HLepType==0)
  {
     pat::MuonCollection::const_iterator muSave;
     int nmuSave =0;
     for(muSave=muColl.begin(); muSave!=muColl.end(); muSave++)
     {
        nmuSave++;
        if(nmuSave==HcentralLep1||nmuSave==HcentralLep2)
        {
        //  cout<<"## "<<iEvent.id().event()<<" nLepSaveMu "<<nmuSave<<endl;
	  LeptonsE_.push_back(muSave->energy());
	  LeptonsPt_.push_back(muSave->pt());
	  LeptonsEta_.push_back(muSave->eta());
	  LeptonsPhi_.push_back(muSave->phi());
          LeptonsType_.push_back(HLepType);	  
          if(!hasFillLVLep1)
          {
             LVLep1.SetPtEtaPhiE(muSave->pt(),muSave->eta(),muSave->phi(),muSave->energy());
             hasFillLVLep1=true;
             LVLep1Charge = muSave->charge();
           
          }
          else
          {
             LVLep2.SetPtEtaPhiE(muSave->pt(),muSave->eta(),muSave->phi(),muSave->energy());  
          }
         }//if is h lep
   }//MUSAVE LOOP
  }//if Hlep ==mu
  else if (HLepType==1)
  {
    pat::ElectronCollection::const_iterator eleSave;
    int neleSave=0;
    for(eleSave=eleColl.begin(); eleSave!=eleColl.end(); eleSave++)  
    {
      neleSave++; 
      if(neleSave==HcentralLep1||neleSave==HcentralLep2)
      {
          //cout<<"## "<<iEvent.id().event()<<" nLepSaveEle "<<neleSave<<endl;
	  LeptonsE_.push_back(eleSave->energy());
	  LeptonsPt_.push_back(eleSave->pt());
	  LeptonsEta_.push_back(eleSave->eta());
	  LeptonsPhi_.push_back(eleSave->phi());
          LeptonsType_.push_back(HLepType);	  
     
          if(!hasFillLVLep1)
          {
             LVLep1.SetPtEtaPhiE(eleSave->pt(),eleSave->eta(),eleSave->phi(),eleSave->energy());
             hasFillLVLep1=true;
             LVLep1Charge = eleSave->charge();
          }
          else
          {
             LVLep2.SetPtEtaPhiE(eleSave->pt(),eleSave->eta(),eleSave->phi(),eleSave->energy());  
          }
        
      }

    }//elesave loop

  }//Hlep ==ele

//add angular variable

numJets_=JetPt_.size();
TLorentzVector LV2l2q;
LV2l2q=LVJ1+LVJ2+LVLep1+LVLep2;
Hpt_=LV2l2q.Pt();
Hm_=LV2l2q.M();

    double CalcosthetaNT1=-100.0;
    double CalcosthetaNT2=-100.0;
    double CalcosthetastarNT=-100.0;
    double CalphiNT=-100.0;
    double CalphiNT1=-100.0;

    Helicity myAngles;
    if(LVLep1Charge==-1)
    {
    myAngles.calculateAngles(LVLep1,LVLep2 ,LVJ1 , LVJ2, CalcosthetaNT1, CalcosthetaNT2, CalcosthetastarNT, CalphiNT, CalphiNT1);
    }
    else if(LVLep1Charge==1)
    {
    myAngles.calculateAngles(LVLep2,LVLep1 ,LVJ1 , LVJ2, CalcosthetaNT1, CalcosthetaNT2, CalcosthetastarNT, CalphiNT, CalphiNT1);
     }
     else
     {
     cout<<"Lep Charge error"<<endl;
     cin.get();
     }

  costhetaNT1_= CalcosthetaNT1;
  costhetaNT2_= CalcosthetaNT2;
  phiNT_= CalphiNT;
  phiNT1_= CalphiNT1;
  costhetastarNT_= CalcosthetastarNT;
//  cout<<"angular "<<costhetaNT1_<<" "<<costhetaNT2_<<" "<<phiNT_<<" "<<phiNT1_<<" "<<costhetastarNT_<<endl;  
   double BDTD=0.0;
    // BDTD=reader->EvaluateMVA( "BDT" );

  float ldSig = -100.0;float ldBkg=-100.0;float helyLD=-100.0;
  HelicityLikelihoodDiscriminant LD_;
  HelicityLikelihoodDiscriminant::HelicityAngles myha;
  myha.helCosTheta1    = costhetaNT1_;
  myha.helCosTheta2    = costhetaNT2_;
  myha.helCosThetaStar = costhetastarNT_;
  myha.helPhi          = phiNT_;
  myha.helPhi1         = phiNT1_;
  myha.mzz             = Hm_;
  LD_.setMeasurables(myha);
  ldSig = LD_.getSignalProbability();
  ldBkg = LD_.getBkgdProbability();
  helyLD = ldSig / (ldSig + ldBkg);

  heliLD_=helyLD;
/*
  cout<<"helLD : "<<heliLD_<<" angular "<<costhetaNT1_<<" "<<costhetaNT2_<<" "<<costhetastarNT_<<" "<<phiNT_<<" "<<phiNT1_<<" "<<endl; 
  cout<<"lep 1  "<<LVLep1.Pt()<<" "<<LVLep1.Eta()<<" "<<LVLep1.Phi()<<" "<<LVLep1.E()<<endl;
  cout<<"lep 2  "<<LVLep2.Pt()<<" "<<LVLep2.Eta()<<" "<<LVLep2.Phi()<<" "<<LVLep2.E()<<endl;
  
  cout<<"J 1  "<<LVJ1.Pt()<<" "<<LVJ1.Eta()<<" "<<LVJ1.Phi()<<" "<<LVJ1.E()<<endl; 
  cout<<"J 2  "<<LVJ2.Pt()<<" "<<LVJ2.Eta()<<" "<<LVJ2.Phi()<<" "<<LVJ2.E()<<endl;
cin.get();
*/
  }//if pass selection
  else
  {
    _nRejected++;
  //  return;
  }
 
LepLapJets.clear();

} // end of Fill()

  //-----------------------------------------------------------------------
void  
YJHiggTree::SetBranches(){

  AddBranch(&eleRho_, "eleRho");
  AddBranch(&muoRho_, "muoRho");
  AddBranch(&metSig_, "metSig");

  AddBranch(&EvtMaxMjj_, "EvtMaxMjj");
  AddBranch(&EvtMaxEtajj_, "EvtMaxEtajj");
  AddBranch(&NBtagYJ_,"NBtagYJ_"); 
  
  AddBranch(&JetpuBetaYJ_,"JetpuBetaYJ_");  
 
  AddBranch(&JetPt_, "JetPt_");
  AddBranch(&JetEta_, "JetEta_");
  AddBranch(&JetPhi_, "JetPhi_");
  AddBranch(&JetM_, "JetM_");
  AddBranch(&JetEn_, "JetEn_");
  AddBranch(&JetJetProb_,"JetJetProb_");
  AddBranch(&JetFromH_, "JetFromH_"); 
  AddBranch(&JetFromSB_, "JetFromSB_"); 
 
  AddBranch(&LeptonsE_,"LeptonsE_");
  AddBranch(&LeptonsPt_,"LeptonsPt_");
  AddBranch(&LeptonsEta_,"LeptonsEta_");
  AddBranch(&LeptonsPhi_,"LeptonsPhi_");
 AddBranch(&LeptonsType_,"LeptonsType_"); 

  AddBranch(&heliLD_,"heliLD_");
  AddBranch(&costhetaNT1_,"costhetaNT1_");
  AddBranch(&costhetaNT2_,"costhetaNT2_");
  AddBranch(&phiNT_,"phiNT_");
  AddBranch(&phiNT1_,"phiNT1_");
  AddBranch(&costhetastarNT_,"costhetastarNT_");
   
  AddBranch(&Hpt_,"Hpt_"); 
  AddBranch(&Hm_,"Hm_");
  AddBranch(&numJets_,"numJets_");

  AddBranch(&EvtType_,"EvtType_");

  AddBranch(&n_pileup_,"n_pileup_");
  AddBranch(&n_pileup_true_,"n_pileup_true_");

  AddBranch(&genHmass_,"genHmass_");
  

   AddBranch(&HLTDoubleMu_,"HLTDoubleMu_");
   AddBranch(&HLTDoubleEle_,"HLTDoubleEle_");
   AddBranch(&HLTMu17TkMu8_,"HLTMu17TkMu8_");

/*
  AddBranch(&higgsPt_,"higgsPt");
  AddBranch(&higgsEta_,"higgsEta");
  AddBranch(&higgsPhi_,"higgsPhi");
  AddBranch(&higgsM_,"higgsM");
  AddBranch(&higgsMRefit_,"higgsMRefit");

  AddBranch(&zllPt_,"zllPt");
  AddBranch(&zllEta_,"zllEta");
  AddBranch(&zllPhi_,"zllPhi");
  AddBranch(&zllM_,"zllM");
  AddBranch(&zlldR_,"zlldR");

  AddBranch(&zjjPt_,"zjjPt");
  AddBranch(&zjjEta_,"zjjEta");
  AddBranch(&zjjPhi_,"zjjPhi");
  AddBranch(&zjjM_,"zjjM");
  AddBranch(&zjjMRefit_,"zjjMRefit");
  AddBranch(&zjjdR_,"zjjdR");

  AddBranch(&HjetIndex_,"HjetIndex");
  AddBranch(&HjetHiggsIndex_,"HjetHiggsIndex_");
  AddBranch(&HjetE_,"HjetE");
  AddBranch(&HjetPt_,"HjetPt");
  AddBranch(&HjetEta_,"HjetEta");
  AddBranch(&HjetPhi_,"HjetPhi");
  
  AddBranch(&LeptonsIndex_,"LeptonsIndex");
  AddBranch(&LeptonsHiggsIndex_,"LeptonsHiggsIndex_");
 
  AddBranch(&heliLD_,"heliLD");
  AddBranch(&costhetaNT1_,"costhetaNT1");
  AddBranch(&costhetaNT2_,"costhetaNT2");
  AddBranch(&phiNT_,"phiNT");
  AddBranch(&phiNT1_,"phiNT1");
  AddBranch(&costhetastarNT_,"costhetastarNT");
  
  AddBranch(&heliLDRefit_,"heliLDRefit");
  AddBranch(&costhetaNT1Refit_,"costhetaNT1Refit");
  AddBranch(&costhetaNT2Refit_,"costhetaNT2Refit");
  AddBranch(&phiNTRefit_,"phiNTRefit");
  AddBranch(&phiNT1Refit_,"phiNT1Refit");
  AddBranch(&costhetastarNTRefit_,"costhetastarNTRefit");

  AddBranch(&nBTags_,"nBTags");
  AddBranch(&lepType_,"lepType");
  AddBranch(&passBit_,"passBit");
*/

}


void  
YJHiggTree::Clear(){

  eleRho_ = DUMMY;
  muoRho_ = DUMMY;

  metSig_ = DUMMY;

  EvtMaxMjj_= DUMMY;
  EvtMaxEtajj_= DUMMY;
  NBtagYJ_= DUMMY;



 
  JetpuBetaYJ_.clear();
  JetPt_.clear();
  JetEta_.clear();
  JetPhi_.clear();
  JetM_.clear();
  JetEn_.clear();
 JetJetProb_.clear();
  JetFromH_.clear();
  JetFromSB_.clear();


  higgsPt_.clear();
  higgsEta_.clear();
  higgsPhi_.clear();
  higgsM_.clear();
  higgsMRefit_.clear();

  zllPt_.clear();
  zllEta_.clear();
  zllPhi_.clear();
  zllM_.clear();
  zlldR_.clear();

  zjjPt_.clear();
  zjjEta_.clear();
  zjjPhi_.clear();
  zjjM_.clear();
  zjjMRefit_.clear();
  zjjdR_.clear();

  HjetIndex_.clear();
  HjetHiggsIndex_.clear();
  HjetE_.clear();
  HjetPt_.clear();
  HjetEta_.clear();
  HjetPhi_.clear();
  
  LeptonsIndex_.clear();
  LeptonsHiggsIndex_.clear();
  LeptonsE_.clear();
  LeptonsPt_.clear();
  LeptonsEta_.clear();
  LeptonsPhi_.clear();
  LeptonsType_.clear();
  


  heliLD_=DUMMY;
  costhetaNT1_= DUMMY;
  costhetaNT2_= DUMMY;
  phiNT_= DUMMY;
  phiNT1_= DUMMY;
  costhetastarNT_= DUMMY;
  
  Hpt_=DUMMY;
  Hm_=DUMMY;
  numJets_= DUMMY;  

  EvtType_=DUMMY;

  n_pileup_=DUMMY;
  n_pileup_true_=DUMMY;


  genHmass_=DUMMY;

   HLTDoubleMu_=DUMMY;
   HLTDoubleEle_=DUMMY;
   HLTMu17TkMu8_=DUMMY;


/* 
  heliLDRefit_.clear();
  costhetaNT1Refit_.clear();
  costhetaNT2Refit_.clear();
  phiNTRefit_.clear();
  phiNT1Refit_.clear();
  costhetastarNTRefit_.clear();

  nBTags_.clear();
  lepType_.clear();
  passBit_.clear();
*/
}


bool YJHiggTree::passLooseJetID(const pat::Jet* recjet)
{
  double eta = recjet->eta();
  if(recjet->getPFConstituents().size() <= 1)return false;                                                                               
  if(recjet->neutralHadronEnergyFraction() >= 0.99)return false;
  if(recjet->neutralEmEnergyFraction() >= 0.99)return false;
  //   // for the tracker region
  if(fabs(eta)<2.4 && recjet->chargedHadronEnergyFraction()<= 0.0)return false;
  if(fabs(eta)<2.4 && recjet->chargedEmEnergyFraction() >= 0.99)return false;
  if(fabs(eta)<2.4 && recjet->chargedMultiplicity() <= 0)return false;
  return true;

}

/*
  eID01.push_back(myEle->pt());
  eID02.push_back(myEle->superCluster()->eta());
  eID03.push_back(myEle->deltaEtaSuperClusterTrackAtVtx());
  eID04.push_back(myEle->deltaPhiSuperClusterTrackAtVtx());
  eID05.push_back(myEle->sigmaIetaIeta());
  eID06.push_back(myEle->hadronicOverEm());
  eID07.push_back(myEle->userFloat("dxy"));
  eID08.push_back(myEle->userFloat("dz"));
  eID09.push_back(fabs(1.0/myEle->ecalEnergy() - 
  myEle->eSuperClusterOverP()/myEle->ecalEnergy()));
  eID10.push_back(myEle->convDcot());
  eID11.push_back(myEle->convDist());
  eID12.push_back(myEle->userFloat("hasMatchConv"));
  eID13.push_back(myEle->gsfTrack().get()->trackerExpectedHitsInner().numberOfHits());
	  
  double iso1 = myEle->chargedHadronIso();
  double iso2 = myEle->neutralHadronIso();
  double iso3 = myEle->photonIso();
	  
  ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = 
  isData? ElectronEffectiveArea::kEleEAData2011:
  ElectronEffectiveArea::kEleEAFall11MC;
      
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaType_ =
  ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03;      
 
  double AEff = ElectronEffectiveArea::GetElectronEffectiveArea
  (effAreaType_, fabs(myEle->superCluster()->eta()), effAreaTarget_);

  double iso4 = iso1 + std::max(0.0, iso2+iso3-ele_rho*AEff);
	  
  eID14.push_back(iso1);
  eID15.push_back(iso2);
  eID16.push_back(iso3);
  eID17.push_back(iso4);
*/


/*
  double pt    =  myMuo->pt();
  double eta   =  myMuo->eta();
  double iso1  =  myMuo->pfIsolationR04().sumChargedHadronPt;
  double iso2  =  myMuo->pfIsolationR04().sumNeutralHadronEt;
  double iso3  =  myMuo->pfIsolationR04().sumPhotonEt;
  double isoPU =  myMuo->pfIsolationR04().sumPUPt;    
  double iso4Beta = iso1 + std::max(iso3+iso2-0.5*isoPU,0.);
  MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = 
  isData? MuonEffectiveArea::kMuEAData2012:
  MuonEffectiveArea::kMuEAData2012;
  // 	MuonEffectiveArea::kMuEAFall11MC;
w
  MuonEffectiveArea::MuonEffectiveAreaType effAreaType_= 
  MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;
  double Area = MuonEffectiveArea::GetMuonEffectiveArea(
  effAreaType_, fabs(eta), effAreaTarget_);

  double iso4Rho =  iso1 + std::max(iso3+iso2-muo_rho*Area,0.);
  double isoBeta = iso4Beta/pt;	  
  double isoRho  = iso4Rho/pt;


  muID01.push_back(isoBeta);
  // 	  muID01.push_back(myMuo->isGlobalMuon());
  muID02.push_back(myMuo->isPFMuon());
  muID03.push_back(myMuo->isTrackerMuon());
  if(myMuo->isTrackerMuon() && myMuo->isGlobalMuon()){
  muID04.push_back(myMuo->globalTrack()->normalizedChi2());
  muID05.push_back(myMuo->globalTrack()->hitPattern().numberOfValidMuonHits());
  muID06.push_back(myMuo->numberOfMatchedStations());
  muID07.push_back(myMuo->dB());
  double dzV =myMuo->userFloat("dzVtx") ;
  muID08.push_back(dzV);
  muID09.push_back(myMuo->innerTrack()->hitPattern().numberOfValidPixelHits());
  muID10.push_back(myMuo->innerTrack()->hitPattern().trackerLayersWithMeasurement());
  }
  muID11.push_back(myMuo->pfIsolationR04().sumChargedHadronPt);
  muID12.push_back(myMuo->pfIsolationR04().sumNeutralHadronEt);
  muID13.push_back(myMuo->pfIsolationR04().sumPhotonEt);
  muID14.push_back(myMuo->pfIsolationR04().sumPUPt);
  muID15.push_back(isoRho);
*/

/*

higgsPt_  = DUMMY;
higgsEta_ = DUMMY;
higgsPhi_ = DUMMY;
higgsM_   = DUMMY;

zllPt_  = DUMMY;
zllEta_ = DUMMY;
zllPhi_ = DUMMY;
zllM_   = DUMMY;

zjjPt_  = DUMMY;
zjjEta_ = DUMMY;
zjjPhi_ = DUMMY;
zjjM_   = DUMMY;


int arraySize = sizeof(jetPt_)/sizeof(jetPt_[0]);

for(int i=0; i<arraySize;i++)
{
jetPt_[i] =DUMMY;
jetEta_[i]=DUMMY;
jetPhi_[i]=DUMMY;
jetE_[i]  =DUMMY;
jetRefitPt_[i] =DUMMY;
jetRefitEta_[i]=DUMMY;
jetRefitPhi_[i]=DUMMY;
jetRefitE_[i]=DUMMY;

}
lepType_ = -1;
  
arraySize = sizeof(lepPt_)/sizeof(lepPt_[0]);

for(int i=0; i<arraySize;i++)
{
lepPt_[i] =DUMMY;
lepEta_[i]=DUMMY;
lepPhi_[i]=DUMMY;
lepE_[i]  =DUMMY;

}

AddBranch(&eID01, "eID01");
AddBranch(&eID02, "eID02");
AddBranch(&eID03, "eID03");
AddBranch(&eID04, "eID04");
AddBranch(&eID05, "eID05");
AddBranch(&eID06, "eID06");
AddBranch(&eID07, "eID07");
AddBranch(&eID08, "eID08");
AddBranch(&eID09, "eID09");
AddBranch(&eID10, "eID10");
AddBranch(&eID11, "eID11");
AddBranch(&eID12, "eID12");
AddBranch(&eID13, "eID13");
AddBranch(&eID14, "eID14");
AddBranch(&eID15, "eID15");
AddBranch(&eID16, "eID16");
AddBranch(&eID17, "eID17");

eID01.clear();
eID02.clear();
eID03.clear();
eID04.clear();
eID05.clear();
eID06.clear();
eID07.clear();
eID08.clear();
eID09.clear();
eID10.clear();
eID11.clear();
eID12.clear();
eID13.clear();
eID14.clear();
eID15.clear();
eID16.clear();
eID17.clear();

AddBranch(&muID01, "muID01");
AddBranch(&muID02, "muID02");
AddBranch(&muID03, "muID03");
AddBranch(&muID04, "muID04");
AddBranch(&muID05, "muID05");
AddBranch(&muID06, "muID06");
AddBranch(&muID07, "muID07");
AddBranch(&muID08, "muID08");
AddBranch(&muID09, "muID09");
AddBranch(&muID10, "muID10");
AddBranch(&muID11, "muID11");
AddBranch(&muID12, "muID12");
AddBranch(&muID13, "muID13");
AddBranch(&muID14, "muID14");
AddBranch(&muID15, "muID15");

muID01.clear();
muID02.clear();
muID03.clear();
muID04.clear();
muID05.clear();
muID06.clear();
muID07.clear();
muID08.clear();
muID09.clear();
muID10.clear();
muID11.clear();
muID12.clear();
muID13.clear();
muID14.clear();
muID15.clear();


*/
