/// @file
/// File containing the definition of the methods associated to the class.
///
///For 539 VBF analysis YunJu Lu 20131225 

#include "DataFormats/Common/interface/View.h"//for refAt
#include "DelPanj/TreeMaker/interface/YJJECHiggTree.hh"
#include "DelPanj/TreeMaker/interface/cutvalues.h"
//#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "DelPanj/TreeMaker/interface/MuonEffectiveArea.h"

// system include files

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/MET.h"


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
//#include "HiggsAnalysis/Higgs2l2b/plugins/RooBkgd2L2JV2.cc"
//#include "HiggsAnalysis/Higgs2l2b/plugins/RooSpinZero5DV2.cc"

//#include "HiggsAnalysis/Higgs2l2b/plugins/Helicity.cc"
//#include "HiggsAnalysis/Higgs2l2b/plugins/HelicityLikelihoodDiscriminant.cc"

//#include "HZZ2l2qAnalysis/Higgs2l2qCode/interface/Helicity.h"
//#include "HZZ2l2qAnalysis/Higgs2l2qCode/interface/HelicityLikelihoodDiscriminant.h"


#include "DelPanj/TreeMaker/interface/Helicity.h"
#include "DelPanj/TreeMaker/interface/HelicityLikelihoodDiscriminant.h"

typedef std::vector< edm::Handle< edm::ValueMap<double> > >             
IsoDepositVals;

typedef std::vector< pat::Jet > PFJetCollectionAB;
void
YJJECHiggTree::AddBranch(int* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}

void
YJJECHiggTree::AddBranch(double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}

//---------------------------------------------------
//---------------------------------------------------
void
YJJECHiggTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}
void
YJJECHiggTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}


//---------------------------------------------------
//---------------------------------------------------
void 
YJJECHiggTree::AddBranchArray(const int arraySize, double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,Form("%s[%d]/D",brName.data(),arraySize));
}

//---------------------------------------------------------------
//---------------------------------------------------------------
YJJECHiggTree::YJJECHiggTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  e2012ID_ ( iConfig.getParameter<edm::ParameterSet>("e2012IDSet")),
  mu2012ID_ ( iConfig.getParameter<edm::ParameterSet>("mu2012IDSet")),
  hzzeejj_(iConfig.getParameter<edm::InputTag>("hzzeejjTag")),
  hzzmmjj_ (iConfig.getParameter<edm::InputTag>("hzzmmjjTag")),
  eleRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("eleRhoIso")),
  muoRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("muoRhoIso")),
  SetJEC_C(iConfig.getParameter<int>("SetJEC_PY")) //,
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
void YJJECHiggTree::endJob (void)
{
   cout<<"Report from YJ: "<<endl;
  cout<<"   - Number of processed events: "<<_nEvents<<endl;
  cout<<"   - Number of rejected events due to trigger: "<<_nFailTrig<<endl;
  cout<<"   - Number of rejected events due to preselection: "<<_nRejected<<endl;
  cout<<"   - Number of passed events: "<<_nPassed<<endl;

}
YJJECHiggTree::~YJJECHiggTree()
{


  delete tree_;
  delete btsfutiljp;
}



// ------------ method called to for each event  ------------
void YJJECHiggTree::Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup)
{
  _nEvents++;
  Clear();
  int eventNum = iEvent.id().event();
  
  //============================================================================
  // 
  //       OBTAIN EVENT-LEVEL VARIABLES
  //  
  //============================================================================

  EvtType_=-1;//1==signal region; 2 = sb region
  EvtLepType_=-1;//0==electron; 1 = muon
  bool hasAPassHCand=false;
  int IndexthetheOneH=-999;// index of  the best H candidate in 1 event




  



  int maxNBTag=-999;
  float Zdiff=999; 
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
    iEvent.getByLabel("cleanPatJetsNoPUIsoLept",JetHandle);
//   iEvent.getByLabel("cleanPatJetsNoPUIsoLept",JetHandle);
    const std::vector<pat::Jet>* patjetsVec = JetHandle.product();


//test pu jet // not working, not save properly in skim?
/*
edm::Handle<edm::View<pat::Jet> > jetsview;
iEvent.getByLabel("cleanPatJetsNoPUIsoLept",jetsview);

Handle<ValueMap<float> > puJetIdMVA;
iEvent.getByLabel("puJetMva","full53xDiscriminant",puJetIdMVA);

Handle<ValueMap<int> > puJetIdFlag;
iEvent.getByLabel("puJetMva","full53xId",puJetIdFlag);

 //pick from the event the input jet collection 
  Handle< PFJetCollectionAB > jetColl;
  //edm::Handle<edm::View<pat::Jet> > jetColl;
  iEvent.getByLabel("cleanPatJetsNoPUIsoLept", jetColl);


    std::auto_ptr< std::vector< pat::Jet > > outputPFJets(new std::vector< pat::Jet >(*jetColl));
  uint njetInColl(0);

  for(PFJetCollectionAB::iterator ijet=outputPFJets->begin(); 
      ijet!=outputPFJets->end();++ijet,++njetInColl ){

edm::RefToBase<pat::Jet> jetRef(edm::Ref<PFJetCollectionAB>(jetColl,njetInColl));

    cout<<"The mva "<<(*puJetIdMVA)[jetRef]<<endl;

    cout<<"The ID "<<(*puJetIdFlag)[jetRef]<<endl;

}
*/

/*
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
*/
//MET
/*
 Handle<float>  metSignificanceHandle;
  iEvent.getByLabel("metInfoProducer","metSignificance", metSignificanceHandle);
  float metsignificance = *metSignificanceHandle;
  metSig_=metsignificance;

*/
///////////////////////////////////JES part by yunju////////////////////////////

      //Yun-Ju for taking jec directlly
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl);
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      //*uncGetter=0;
      uncGetter = new JetCorrectionUncertainty(JetCorPar);

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
   
 /*   int nj1=0;
    int nGoodJets=0;
    bool passJetWin= false;
    bool passJetSB= false;
    std:vector<int> LepLapJets; //a good jet loose id ,pt 30 eta4.7, no lepton overlap

    //cout<<"TD "<<endl;
	    bool UP =true;
    for(std::vector<pat::Jet>::const_iterator jet1Itr =patjetsVec->begin();jet1Itr!=patjetsVec->end();jet1Itr++)
    {
      nj1++;

    //ScaleJet(pat::Jet & dest, const pat::Jet * j,const bool UP);
    ScaleJet(j1,dynamic_cast<const pat::Jet *>(&*jet1Itr),SetJEC_C);
    //cout<<j1.pt()<<" jet1 "<<jet1->pt()<<" nJet "<<nj1<<endl;


      //if(jet1->pt()<30) continue;
      if(j1.pt()<30) continue;//j1 is fater jec
      if(fabs(j1.eta())>4.7) continue;
      if(!passLooseJetID(&*jet1Itr)) continue;

 //     cout<<nj1<<" PUID "<<jet1->userInt("puJetIdFlag")<<" "<<jet1->userFloat("puJetIdMVA")<<endl;     

      //if(jet1->userInt("isTaggable")==0) continue;//use charge tracks has |eta|<2.4
      //if(jet1->userFloat("puBeta")<0.2) continue;
      nGoodJets++;
    }//jet1Itr loop
*/
    //check di jet win and maxmjj maxetajj in 1 event
   //find out the highest btag di-jet pair to save 
//cout<<"################################### Start Event "<<iEvent.id().event()<<endl;
//cout<<<<endl;

const int LEPZ=0;
const int HADZ=1;


const  double MIN_LEPPT1 = 40.0;
const  double MIN_LEPPT2 = 20.0;

const double MIN_JETPT  =30.0;
const double MIN_JETBETA= 0.2;
const double MAX_JETETA = 2.4;
//for VBF study
//const double MAX_JETETA = 4.7;



const double MIN_MZ_LL=76.0;
const double MAX_MZ_LL=106.0;

const double MIN_MZ_JJ=71.0;
const double MAX_MZ_JJ=111.0;


const double LOOSE_MIN_MZ_JJ=60.0;
const double LOOSE_MAX_MZ_JJ=130.0;
const double HELICUT=0.5;


  //initialize variables
  int hcand = -1;
  //LOOP IN 2 KINDS OF LEPTONS: ELECTRONS AND MUONS 
  for (int ilep=0; ilep<2; ilep++) {

    //check cand which pass trigger path  
    if(HLTDoubleEle_!=1&&ilep==0) continue;    
    if( (HLTDoubleMu_!=1&&HLTMu17TkMu8_!=1)&&ilep==1) continue;

    hcand=-1;

    //GET THE HIGGS->ZZ->LLJJ COLLECTION
    Handle<std::vector<pat::CompositeCandidate> > hzzlljj;
    if(ilep==0) iEvent.getByLabel("hzzeejj","h",hzzlljj);
    else iEvent.getByLabel("hzzmmjj","h" ,hzzlljj);





    // LOOP OVER HIGGS CANDIDATES

    for(unsigned iHcand=0; iHcand<hzzlljj->size(); iHcand++){
      const pat::CompositeCandidate & h = (*hzzlljj)[iHcand];	
/*
     //Do jec on h cand
     const pat::Jet * j1 = GetJet1(h);
     const pat::Jet * j2 = GetJet2(h);
*/


      const reco::Candidate*  myLepton[2];  
      for(unsigned int il=0; il < 2; il++)
	myLepton[il]= h.daughter(LEPZ)->daughter(il)->masterClone().get();

       const pat::Jet * myJECJet[2];//JEC jet
       const pat::Jet * mydefaultJet[2];//default jet

       for(unsigned int ijet=0; ijet < 2; ijet++)
       {
	mydefaultJet[ijet] =
	  dynamic_cast<const pat::Jet *>(h.daughter(HADZ)->daughter(ijet)->
					 masterClone().get());


        if(ijet==0) myJECJet[ijet]= GetJet1(h);
        else if (ijet==1) myJECJet[ijet]= GetJet2(h);             

        //cout<<"mydefaultJet:"<<mydefaultJet[ijet]->pt()<<" "<<mydefaultJet[ijet]->eta()<<endl;
       // cout<<"    myJECJet:"<<myJECJet[ijet]->pt()<<" "<<myJECJet[ijet]->eta()<<endl; 
      }
      // LOOK FOR TWO GOOD CHARGED LEPTONS
      int nLepPtHi=0;
      int nLepPtLo=0;

      for(unsigned int il=0; il < 2; il++){
	
	if(deltaR(myLepton[il]->eta(), myLepton[il]->phi(),
		  myJECJet[0]->eta(), myJECJet[0]->phi()) < MIN_DR_JETLEP)continue;

	if(deltaR(myLepton[il]->eta(), myLepton[il]->phi(),
		  myJECJet[1]->eta(), myJECJet[1]->phi()) < MIN_DR_JETLEP)continue;


	double pt = myLepton[il]->pt();	  
	if(pt >MIN_LEPPT1) nLepPtHi++;
	if(pt >MIN_LEPPT2) nLepPtLo++;

      }

      if(nLepPtHi < 1)continue;
      if(nLepPtLo < 2)continue;
      
      bool OppCharge = myLepton[0]->charge()*myLepton[1]->charge() < 0 ? 
	true: false;

      if(!OppCharge)continue;

      // they need to pass ID cuts
      int nPassID=0;

      if(ilep==0){ // electron
	nPassID=0;

	for(unsigned int iele=0; iele < 2; iele++){
	  
           const pat::Electron* myEle
	    = dynamic_cast<const pat::Electron*>(h.daughter(LEPZ)->daughter(iele)->masterClone().get());
	  
	  if(myEle->userInt("cutIDCode")<2) continue; // 539 skim	  
	  nPassID++;
	}
      } // if it's an electron type

      else if(ilep==1){ // muon

	nPassID=0;
	for(unsigned int imuo=0; imuo < 2; imuo++){
	 
	  const pat::Muon* myMuo
	    = dynamic_cast<const pat::Muon*>(h.daughter(LEPZ)->daughter(imuo)->masterClone().get());
	  //std::map<std::string, bool> Pass = mu2012ID_.CutRecord(*myMuo);
	  //int passOrNot = PassAll(Pass);
	  if(myMuo->userInt("isIsolated")!=1) continue; // 539 skim (also pass id)	  
	  nPassID++;	  	  
	} // end of loop over muon
      } // if is a muon type
     
      
      if(nPassID < 2)continue;

      // LOOK FOR 2 JETS PASSING BETA CUTS

      int nGoodJets=0;
      int nLooseBTags=0;
      int nMediumBTags=0;

      for(unsigned int ijet=0; ijet < 2; ijet++){	 


	
	double pt  = myJECJet[ijet]->pt();
//        double pt =HJJEC.pt();  
	if(pt < MIN_JETPT)continue;
	
	double eta = myJECJet[ijet]->eta();
//        double eta=HJJEC.pt();
	if(fabs(eta)> MAX_JETETA)continue;
	
	// to suppress jets from pileups
	double puBeta = myJECJet[ijet]->userFloat("puBeta");
	if(puBeta < MIN_JETBETA)continue;

	
	if(deltaR(myLepton[0]->eta(), myLepton[0]->phi(),
		  myJECJet[ijet]->eta(), myJECJet[ijet]->phi()) < MIN_DR_JETLEP)continue;

	if(deltaR(myLepton[1]->eta(), myLepton[1]->phi(),
		  myJECJet[ijet]->eta(), myJECJet[ijet]->phi()) < MIN_DR_JETLEP)continue;
	  
	
	if( !passLooseJetID(myJECJet[ijet]) )continue;


	bool isLoose  = false;
	bool isMedium = false;
	double jpBTag = myJECJet[ijet]->bDiscriminator("jetProbabilityBJetTags");
	int flavor    = myJECJet[ijet]->partonFlavour();
	
	if(jpBTag > MIN_BTAG_JP_LOOSE)isLoose=true;
	if(jpBTag > MIN_BTAG_JP_MEDIUM)isMedium=true;
/*	
 	if(!isData)
	  {
	    double phi = myJECJet[ijet]->phi();
	    double sin_phi = sin(phi*1000000);
	    // 	    int seed = abs(static_cast<int>(sin_phi*100000));	    
	    // 	    BTagSFUtil* btsfutiljp = new BTagSFUtil("JP", seed);
	
	    float Btageff_SF_ = 1.0;

	    for( size_t iMeasure = 0; iMeasure < measureName.size(); iMeasure++ ){
	      //Setup our measurement
	      iSetup.get<BTagPerformanceRecord>().get( measureName[ iMeasure ],perfH);
	      const BtagPerformance & perf = *(perfH.product());
	      BinningPointByMap measurePoint;
	      measurePoint.reset();
	      ///// pass in the et of the jet
	      double jetEt = myJECJet[ijet]->et();
	      measurePoint.insert(BinningVariables::JetEt, jetEt);  
	      ///// pass in the absolute eta of the jet
	      measurePoint.insert(BinningVariables::JetEta, fabs(eta));       
	      //this is the correction for 2012
 	      float SFL_JPL_2012corr = 1.01744  + (0.000813491*jetEt)-
  		(6.01592e-07)*jetEt*jetEt;
  	      float SFL_JPM_2012corr = 0.964487 + (0.00134038*jetEt)-
  		(1.43995e-06)*jetEt*jetEt;

	      // 2011
	      // 	      float SFL_JPL_2012corr = 1.00;
	      //  	      float SFL_JPM_2012corr = 1.00;

	      // Extract the mistag eff value
	      if ( measureType[ iMeasure ] == "BTAGLEFF") {
		// add a suffix eff so we can distingiush it from other values
		std::string suffix = "eff"; 
		ScaleFactorsEff[ijet][ measureName[ iMeasure ] + suffix ] = perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
	      }
	      else{ 
		// Extract the mistag and btag SF
		// The factor Btageff_SF_ is used for Btagging systematics 
		// and should be set to 1.0 as default
		if(measureName[ iMeasure ] == "MISTAGJPM")
		  ScaleFactors[ijet][ measureName[ iMeasure ] ] = 
		    SFL_JPM_2012corr*Btageff_SF_*perf.getResult( measureMap[ measureType[ iMeasure] ],measurePoint);

		else if(measureName[ iMeasure ] == "MISTAGJPL")
		  ScaleFactors[ijet][ measureName[ iMeasure ] ] = 
		    SFL_JPL_2012corr*Btageff_SF_*perf.getResult( measureMap[ measureType[ iMeasure] ],measurePoint);
		else
		  ScaleFactors[ijet][ measureName[ iMeasure ] ] = Btageff_SF_*perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
	    
	      }
	    }
	    
   	    btsfutiljp->modifyBTagsWithSF( 
   					  isLoose, isMedium, flavor, 
   					  ScaleFactors[ijet]["MUJETSWPBTAGJPL"],
   					  ScaleFactors[ijet]["MUJETSWPBTAGJPM"], 
   					  ScaleFactors[ijet]["MISTAGJPL"],
   					  ScaleFactors[ijet]["MISTAGJPM"],
 					  ScaleFactorsEff[ijet]["MISTAGJPLeff"],
   					  ScaleFactorsEff[ijet]["MISTAGJPMeff"]);
	    
	    // 	    delete btsfutiljp;
	  } // if it's MC
*/	
	if(isLoose)  nLooseBTags++;
	if(isMedium) nMediumBTags++;

	nGoodJets++; 
      } // end of loop over jets
      
      // there must be at least two jets
      if(nGoodJets < 2)continue;
      // number of btags
      int nBTags = 0;
      if(nMediumBTags >= 1 && nLooseBTags == 2)nBTags=2;
      else if(nLooseBTags >=1)nBTags=1;
      else nBTags=0;

      hcand = iHcand;
	
      
      const pat::CompositeCandidate & goodH = (*hzzlljj)[hcand];
      const reco::Candidate * Zll = goodH.daughter(LEPZ);

      const reco::Candidate * Zjj = goodH.daughter(HADZ);
      //Final selection in categories

      double higgsPt_local  = goodH.pt();
      double higgsEta_local = goodH.eta();
      double higgsPhi_local = goodH.phi();
      double higgsM_local   = goodH.mass();

      double zllPt_local  = Zll->pt();
      double zllEta_local = Zll->eta();
      double zllPhi_local = Zll->phi();
      double zllM_local   = Zll->mass();
      double zlldR_local  = deltaR(myLepton[0]->eta(),myLepton[0]->phi(),
				   myLepton[1]->eta(),myLepton[1]->phi());
      
      double zjjPt_local  = Zjj->pt();
      double zjjEta_local = Zjj->eta();
      double zjjPhi_local = Zjj->phi();
      double zjjM_local   = Zjj->mass();
      double zjjdR_local  = deltaR(myJECJet[0]->eta(),myJECJet[0]->phi(),
				   myJECJet[1]->eta(),myJECJet[1]->phi());


      
      double heliLD_local = goodH.userFloat("helyLD");


      TLorentzVector VJecJ1;
      VJecJ1.SetPtEtaPhiE(myJECJet[0]->pt(),myJECJet[0]->eta(),myJECJet[0]->phi(),myJECJet[0]->energy());      
      TLorentzVector VJecJ2; 
      VJecJ2.SetPtEtaPhiE(myJECJet[1]->pt(),myJECJet[1]->eta(),myJECJet[1]->phi(),myJECJet[1]->energy());

      TLorentzVector VJECJJ;
      VJECJJ=VJecJ1+VJecJ2;
      double zjjJECM=0.0;
      zjjJECM = VJECJJ.M();   
     // cout<<" zjjM"<<zjjJECM<<" "<<zjjM_local<<endl;

      // check the dilepton Z mass
             if(zllM_local < MIN_MZ_LL)continue;
             if(zllM_local > MAX_MZ_LL)continue;

      // check the mass of the two jets
             /*
             if(zjjM_local < LOOSE_MIN_MZ_JJ)continue;
             if(zjjM_local > LOOSE_MAX_MZ_JJ)continue;
             */
             
             if(zjjJECM < LOOSE_MIN_MZ_JJ)continue;
             if(zjjJECM > LOOSE_MAX_MZ_JJ)continue; 
             

      //helicity cut
       if(heliLD_local> HELICUT ) continue;


       hasAPassHCand=true;//there is at least 1 higg cand pass selection


       //the part to save theoneH

       if(nBTags>maxNBTag)
       {
          maxNBTag=nBTags;
          Zdiff=fabs(zjjJECM-91.2);
          IndexthetheOneH=hcand;
          EvtLepType_=ilep;
       }
       else if(nBTags==maxNBTag)
       {
          if(fabs(zjjJECM-91.2)<Zdiff)    
          {
            Zdiff=fabs(zjjJECM-91.2);
            IndexthetheOneH=hcand;
            EvtLepType_=ilep; 
          } 

       }
      
      //To print out for testing

      //cout<<"event id (all hig)"<<_nEvents<<" Leptype"<< ilep<<" All higgs "<<higgsPt_local<<" Btag "<<nBTags<<" Mjj "<<zjjM_local<<endl;



      
      /*
      double qgld= 99999.0; // needs to be fixed
      if(qgld > MIN_QUARK_GLUON_LD_0BTAG)
	thisBit |= QG_LD;
     e
      // need to implement cuts on 4-body mass, need to be fixed
      double mH_input = higgsM_refit_local;
      if(higgsM_refit_local > MIN_MH_RATIO*mH_input &&
	 higgsM_refit_local < MAX_MH_RATIO*mH_input)
	thisBit |= MH_SIGNAL;
*/
    } // end of loop over Higgs candidates  
    
 } //end of looping over leptype 

     if(!hasAPassHCand) return;
    
      //SAVE THE ONE HIGGS
     
 
      //GET THE HIGGS->ZZ->LLJJ COLLECTION FOR THE ONE
      Handle<std::vector<pat::CompositeCandidate> > hzzlljj;
      if(EvtLepType_==0) iEvent.getByLabel("hzzeejj","h", hzzlljj);
      else iEvent.getByLabel("hzzmmjj","h", hzzlljj);

         

      const pat::CompositeCandidate & theOneH = (*hzzlljj)[IndexthetheOneH];


      const reco::Candidate * Zll = theOneH.daughter(LEPZ);
      const reco::Candidate * Zjj = theOneH.daughter(HADZ);
 
      const reco::Candidate*  theOneLep[2];  
      for(unsigned int il=0; il < 2; il++)
      {  	theOneLep[il]= theOneH.daughter(LEPZ)->daughter(il)->masterClone().get();
        HLeptonsE_.push_back(theOneLep[il]->energy());
        HLeptonsPt_.push_back(theOneLep[il]->pt());
        HLeptonsEta_.push_back(theOneLep[il]->eta());
        HLeptonsPhi_.push_back(theOneLep[il]->phi());
         
      }

      




 
      const pat::Jet * theOneJet[2];
      for(unsigned int ijet=0; ijet < 2; ijet++)
      {	theOneJet[ijet] =
	  dynamic_cast<const pat::Jet *>(theOneH.daughter(HADZ)->daughter(ijet)->
      					 masterClone().get());
          HJetE_.push_back(theOneJet[ijet]->energy());
          HJetPt_.push_back(theOneJet[ijet]->pt());
          HJetEta_.push_back(theOneJet[ijet]->eta());
          HJetPhi_.push_back(theOneJet[ijet]->phi());
      }
       
       theOneHNBtag_.push_back(maxNBTag);
       theOnehiggsPt_.push_back( theOneH.pt());
       theOnehiggsEta_.push_back( theOneH.eta());
       theOnehiggsPhi_.push_back( theOneH.phi());
       theOnehiggsM_.push_back( theOneH.mass());
       theOnehiggsMRefit_.push_back( theOneH.userFloat("HZZRefitMass"));

       theOneHzllPt_.push_back(  Zll->pt());
       theOneHzllEta_.push_back( Zll->eta());
       theOneHzllPhi_.push_back( Zll->phi());
       theOneHzllM_.push_back(  Zll->mass());
      
       theOneHzjjPt_.push_back( Zjj->pt());
       theOneHzjjEta_.push_back( Zjj->eta());
       theOneHzjjPhi_.push_back( Zjj->phi());
       theOneHzjjM_.push_back( Zjj->mass());
       theOneHzjjMRefit_.push_back( theOneH.userFloat("ZjjRefitMass"));
     

      // CHECK the one HIGGS CANDIDATE is SR or SB

      if(Zjj->mass() > MIN_MZ_JJ && Zjj->mass() < MAX_MZ_JJ)
      {  
         EvtType_=1;
         _nPassed++;
        if(HLTDoubleMu_==1) {_nPFinalHLTDoubleMu++;}
        if(HLTDoubleEle_==1){_nPFinalHLTDoubleEle++;}
        if(HLTMu17TkMu8_==1){_nPFinalHLTMu17TkMu8++;}
      }
      else
      { 
	EvtType_ =2;
        _nSB++;  
      } 

      







      //cout<<"event id (the one)"<<_nEvents<<" Leptype"<<EvtLepType_ <<" All higgs "<<theOnehiggsPt_<<" Btag "<<maxNBTag<<" Mjj "<<Zjj->mass()<<endl;


//cin.get();
    for(std::vector<pat::Jet>::const_iterator jet1Itr =patjetsVec->begin();jet1Itr!=patjetsVec->end();jet1Itr++)
    {
      

      ScaleJet(j1JEC,dynamic_cast<const pat::Jet *>(&*jet1Itr),SetJEC_C);


      if(j1JEC.pt()<30) continue;//j1JEC is afater jec
      if(fabs(j1JEC.eta())>4.7) continue;
      if(!passLooseJetID(&*jet1Itr)) continue;

      
      int idflag=jet1Itr->userInt("puJetIdFlag");
      //cout<<"test PU "<<PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )<<endl; 
      if(!PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )) continue;
      bool isTheOneHjet=false; 
      if(fabs(theOneJet[0]->pt()-jet1Itr->pt())<0.001||fabs(theOneJet[1]->pt()-jet1Itr->pt())<0.001)
      isTheOneHjet=true;

      JetPt_.push_back(j1JEC.pt());
      JetEta_.push_back(j1JEC.eta());
      JetPhi_.push_back(j1JEC.phi());
      JetEn_.push_back(j1JEC.energy());
      JetFromtheOneH_.push_back(isTheOneHjet);
      JetPUMVA_.push_back(jet1Itr->userFloat("puJetIdMVA")); 
      }//jet1Itr loop



 delete uncGetter;
} // end of Fill()

  //-----------------------------------------------------------------------
void  
YJJECHiggTree::SetBranches(){

  AddBranch(&EvtType_,"EvtType");
  AddBranch(&EvtLepType_,"EvtLepType");
  AddBranch(&genHmass_,"genHmass");
  AddBranch(&HLTDoubleMu_,"HLTDoubleMu");
  AddBranch(&HLTDoubleEle_,"HLTDoubleEle");
  AddBranch(&HLTMu17TkMu8_,"HLTMu17TkMu8");
  
  AddBranch(&n_pileup_,"n_pileup");
  AddBranch(&n_pileup_true_,"n_pileup_true");

  //start vector variable
  
  AddBranch(&theOneHNBtag_,"theOneHNBtag"); 

  AddBranch(&theOnehiggsPt_,"theOnehiggsPt");
  AddBranch(&theOnehiggsEta_,"theOnehiggsEta");
  AddBranch(&theOnehiggsPhi_,"theOnehiggsPhi");
  AddBranch(&theOnehiggsM_,"theOnehiggsM");
  AddBranch(&theOnehiggsMRefit_,"theOnehiggsMRefit");

  AddBranch(&theOneHzllPt_,"theOneHzllPt");
  AddBranch(&theOneHzllEta_,"theOneHzllEta");
  AddBranch(&theOneHzllPhi_,"theOneHzllPhi");
  AddBranch(&theOneHzllM_,"theOneHzllM");
  AddBranch(&theOneHzlldR_,"theOneHzlldR");

  AddBranch(&theOneHzjjPt_,"theOneHzjjPt");
  AddBranch(&theOneHzjjEta_,"theOneHzjjEta");
  AddBranch(&theOneHzjjPhi_,"theOneHzjjPhi");
  AddBranch(&theOneHzjjM_,"theOneHzjjM");
  AddBranch(&theOneHzjjMRefit_,"theOneHzjjMRefit");
  AddBranch(&theOneHzjjdR_,"theOneHzjjdR");


  AddBranch(&HJetE_,"HJetE");
  AddBranch(&HJetPt_,"HJetPt");
  AddBranch(&HJetEta_,"HJetEta");
  AddBranch(&HJetPhi_,"HJetPhi");
  
  AddBranch(&HLeptonsE_,"HLeptonsE");
  AddBranch(&HLeptonsPt_,"HLeptonsPt");
  AddBranch(&HLeptonsEta_,"HLeptonsEta");
  AddBranch(&HLeptonsPhi_,"HLeptonsPhi");
  
 
  AddBranch(&JetPt_, "JetPt");
  AddBranch(&JetEta_, "JetEta");
  AddBranch(&JetPhi_, "JetPhi");
  AddBranch(&JetEn_, "JetEn");
  AddBranch(&JetFromtheOneH_, "JetFromtheOneH"); 
  AddBranch(&JetPUMVA_,"JetPUMVA");
 
 


  AddBranch(&heliLD_,"heliLD");
  AddBranch(&costhetaNT1_,"costhetaNT1");
  AddBranch(&costhetaNT2_,"costhetaNT2");
  AddBranch(&phiNT_,"phiNT");
  AddBranch(&phiNT1_,"phiNT1");
  AddBranch(&costhetastarNT_,"costhetastarNT");
   

  






}

void  
YJJECHiggTree::Clear(){


   EvtType_=DUMMY;
  EvtLepType_=DUMMY;
  n_pileup_=DUMMY;
  n_pileup_true_=DUMMY;
   genHmass_=DUMMY;
  

   HLTDoubleMu_=DUMMY;
   HLTDoubleEle_=DUMMY;
   HLTMu17TkMu8_=DUMMY;


   theOneHNBtag_.clear();

   theOnehiggsPt_.clear();
   theOnehiggsEta_.clear();
   theOnehiggsPhi_.clear();
   theOnehiggsM_.clear();
   theOnehiggsMRefit_.clear();

   theOneHzllPt_.clear();
   theOneHzllEta_.clear();
   theOneHzllPhi_.clear();
   theOneHzllM_.clear();
   theOneHzlldR_.clear(); // deltaR between two leptons

   theOneHzjjPt_.clear();
   theOneHzjjEta_.clear();
   theOneHzjjPhi_.clear();
   theOneHzjjM_.clear();
   theOneHzjjMRefit_.clear();
   theOneHzjjdR_.clear(); // deltaR between two jets   

  HJetE_.clear();
  HJetPt_.clear();
  HJetEta_.clear();
  HJetPhi_.clear();

  HLeptonsE_.clear();
  HLeptonsPt_.clear();
  HLeptonsEta_.clear();
  HLeptonsPhi_.clear();


 
  JetPt_.clear();
  JetEta_.clear();
  JetPhi_.clear();
  JetEn_.clear();
  JetFromtheOneH_.clear(); 
  JetPUMVA_.clear();

 
  
  heliLD_=DUMMY;
  costhetaNT1_= DUMMY;
  costhetaNT2_= DUMMY;
  phiNT_= DUMMY;
  phiNT1_= DUMMY;
  costhetastarNT_= DUMMY;

}

bool YJJECHiggTree::passLooseJetID(const pat::Jet* recjet)
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

void YJJECHiggTree::ScaleJet(pat::Jet & dest, const pat::Jet * j,const int SetJEC)
{
     if(SetJEC!=-1&&SetJEC!=0&&SetJEC!=1){
      std::cout<<"FATAL EXCEPTION: "<<"SetJEC:wrong input "
               <<std::endl; exit(0);}


    bool UP=false;  
    if(SetJEC==-1){UP=false;} 
    else if(SetJEC==1){UP=true;}

    //cout<<"scale!!!!!!!!!!!!!"<<endl;  
    dest = pat::Jet(*j); // Copying the jet
    reco::Candidate::LorentzVector p = j->p4();
    double jeta = j->eta();
         if (jeta >  5.19) jeta =  5.19;
    else if (jeta < -5.19) jeta = -5.19;
    try  // getUncertainty crashes if out of bounds -- catching and ignoring 
    {
        uncGetter->setJetPt (j->pt());
        uncGetter->setJetEta(jeta);
        double unc =0;
        if(SetJEC==0){unc = 1.;}
        else if(SetJEC==-1||SetJEC==1){      
        unc   = 1. + (UP?1.:-1.) * uncGetter->getUncertainty(UP);}
        //cout<<"##########################################################Uncertainty :"<<unc<<endl;
        p *= unc;
    }catch(...){}
    dest.setP4(p);
}

