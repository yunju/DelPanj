#include "DelPanj/TreeMaker/interface/patElecTree.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
// for conversion finder
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"




/*
NOTE: For the users of particle flow.

Really a large number of issues have to be 
explained for the case of the particle flow
electrons. This module at the moment treats
the electrons as if they were the non particle flow
objects. We will make all our selections similar to 
the regular objects. Particle flow electrons need 
deeper study.
*/



patElecTree::patElecTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
e2012ID_ ( iConfig.getParameter<edm::ParameterSet>("e2012IDSet"))
{
 tree_=tree; 
  patElecLabel_ = iConfig.getParameter<edm::InputTag>("patElectronsPY");
  beamSpotLabel_ = iConfig.getParameter<edm::InputTag> ("beamSpotLabel");
  EleRhoLabel_ = iConfig.getParameter<edm::InputTag> ("EleRhoPY");
  SetBranches();
}


patElecTree::~patElecTree(){
  delete tree_;

} 

void
patElecTree::ConversionRejection(const edm::Event& iEvent, const edm::EventSetup& iSetup, const pat::Electron& ele)
{
/*  double evt_bField;
  
  edm::Handle<reco::TrackCollection> tracks_h;
  iEvent.getByLabel("generalTracks", tracks_h);
  
  if (isData_) 
  {
    edm::Handle<DcsStatusCollection> dcsHandle;
    iEvent.getByLabel(dcsTag_, dcsHandle);
    
    // scale factor = 3.801/18166.0 which are average values taken over a stable two week period
    float currentToBFieldScaleFactor = 2.09237036221512717e-04;
    float current = (*dcsHandle)[0].magnetCurrent();
    evt_bField = current*currentToBFieldScaleFactor;
  }
  else 
  {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
        
    evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
  }
  
  ConversionFinder convFinder;
  ConversionInfo convInfo = convFinder.getConversionInfo(ele, tracks_h, evt_bField);

  patElecMissingHits_.push_back( ele.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits() );
  
  patElecDist_.push_back(convInfo.dist() );
  patElecDeltaCotTheta_.push_back(convInfo.dcot() ) ;
  patElecConvRadius_.push_back(convInfo.radiusOfConversion() ) ;*/
}


double
patElecTree::IsoDeposit(const pat::Electron lepton,std::string type,double radius,double coneradius, double detastrip){

  //This portion is highly error-prone and need extensive 
  //checking before we can really trust it. 

  //Motto: Make Vetos and then calculate the pt deposits
  typedef reco::isodeposit::AbsVetos AbsVetos;
  AbsVetos Vetos;
  reco::isodeposit::ThresholdVeto *thveto;
  reco::isodeposit::RectangularEtaPhiVeto *rveto;
  reco::isodeposit::ConeVeto *cveto;
  if(std::string::npos!=type.find("chHad")){
    thveto = new reco::isodeposit::ThresholdVeto(0.5);
    //rveto  = new reco::isodeposit::RectangularEtaPhiVeto(reco::isodeposit::Direction() ,0,0,0,0);
    //cveto  = new reco::isodeposit::ConeVeto(reco::isodeposit::Direction() ,0);
    Vetos.push_back(thveto);
    //Vetos.push_back(rveto);
    //Vetos.push_back(cveto);
    double chIso = (lepton.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin( radius,Vetos, false ).first);
    return chIso;
  }

  if(std::string::npos!=type.find("neuHad")){
    thveto = new reco::isodeposit::ThresholdVeto(0.5);
    //rveto  = new reco::isodeposit::RectangularEtaPhiVeto(reco::isodeposit::Direction() ,0,0,0,0);
    cveto  = new reco::isodeposit::ConeVeto(reco::isodeposit::Direction() ,coneradius);
    Vetos.push_back(thveto);
    //Vetos.push_back(rveto);
    Vetos.push_back(cveto);
    double nhIso = (lepton.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(radius,Vetos, false ).first);
    return nhIso;
  }
  
  if(std::string::npos!=type.find("gam")){
    thveto = new reco::isodeposit::ThresholdVeto(0.5);
    rveto  = new reco::isodeposit::RectangularEtaPhiVeto(reco::isodeposit::Direction() ,-detastrip,detastrip,-100.0,100.0);
    //    cveto  = new reco::isodeposit::ConeVeto(reco::isodeposit::Direction() ,0);
    Vetos.push_back(thveto);
    Vetos.push_back(rveto);
    //    Vetos.push_back(cveto);
    double phIso = (lepton.isoDeposit(pat::PfGammaIso)->depositAndCountWithin( radius,Vetos ,false ).first);
    return phIso;
  }
  
  return -9999;
}


void
patElecTree::Fill(const edm::Event& iEvent){
Clear();
  
//cout<<"START"<<endl;

Handle<std::vector<pat::Electron> > patElecHandle;
//cout<<"STARTv"<<endl;
  if(not iEvent.getByLabel(patElecLabel_,patElecHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<patElecLabel_<<std::endl; exit(0);}
//cout<<"START2"<<endl;

/*
   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel(beamSpotLabel_, beamSpotHandle);
   reco::BeamSpot beamSpot = (*beamSpotHandle);
*/
   double fastJetRho=-999999.9;
   double PFchsJetRho=-999999.9;
   double  lepIsoRho =  -999999.9;
 
  /////// Pileup density "rho" in the event from fastJet pileup calculation /////
  edm::Handle<double> rho;
//  const edm::InputTag eventrho("kt6PFJets", "rho");
//  iEvent.getByLabel(eventrho,rho);
//  if( *rho == *rho) fastJetRho = *rho;
//  else  fastJetRho =  -999999.9;
  fastJetRho_.push_back(fastJetRho);

  /////// Pileup density "rho" in the event from PFchs pileup calculation /////
  edm::Handle<double> rhochs;
//cout<<"pwwwwwwwwwwdddddww3"<<endl;
//  const edm::InputTag eventrhochs(JetForElecTree_, "rho");
//  const edm::InputTag eventrhochs("kt6PFJets", "rho");                    
  /*
  if(not iEvent.getByLabel(JetForElecTree_,rhochs)){
  std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
                <<JetForElecTree_<<std::endl; exit(0);}
*/
//not use rho for now 12/20
/*
   iEvent.getByLabel(eventrhochs,rhochs);


  if( *rhochs == *rhochs) PFchsJetRho = *rhochs;
  else  PFchsJetRho =  -999999.9;


  PFchsJetRho_.push_back(PFchsJetRho);
*/
PFchsJetRho_.push_back(-999);

  /////// Pileup density "rho" for lepton isolation subtraction /////
//not use rho for now 12/20
/*
  edm::Handle<double> rhoLepIso;
  const edm::InputTag eventrhoLepIso("kt6PFJets", "rho");
  iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
  if( *rhoLepIso == *rhoLepIso) lepIsoRho = *rhoLepIso;
  else  lepIsoRho =  -999999.9;
 
 lepIsoRho_.push_back(lepIsoRho);
*/

edm::Handle<double> ele_rho_event;
iEvent.getByLabel(EleRhoLabel_,ele_rho_event);

lepIsoRho_.push_back(*(ele_rho_event.product()));
//cout<<"lepIsoRho_"<<*(ele_rho_event.product())<<endl; 
  pat::ElectronCollection eleColl(*(patElecHandle.product()));
  std::sort(eleColl.begin(),eleColl.end(),PtGreater());

  nele=eleColl.size();
  pat::ElectronCollection::const_iterator ele;
  // cout <<"El: "<<nele<< endl;
  
  for(ele=eleColl.begin(); ele!=eleColl.end(); ele++){
  // cout <<"El333"<<endl;
  
   //eleid
   std::map<std::string, bool> Pass =  e2012ID_.CutRecord(*ele);
    int passOrNot = PassAll(Pass);
  // cout <<"Elec passOrNot: "<<passOrNot<<endl;
     patElecID_.push_back(passOrNot);
    patElecEt_.push_back(ele->et());
    patElecEnergy_.push_back(ele->energy());
    patElecPt_.push_back(ele->p4().pt());
    patElecEta_.push_back(ele->eta());
    patElecPhi_.push_back(ele->phi());
    patElecM_.push_back(ele->mass());
    patElecTrkIso_.push_back(ele->dr03TkSumPt());
    patElecHcalIso_.push_back(ele->dr03HcalTowerSumEt());
    patElecEcalIso_.push_back(ele->dr03EcalRecHitSumEt());
    patElecCharge_.push_back(ele->threeCharge());

//     if ( fabs(ele->eta()) < 1.479 ) // barrel
//       patElecRelIsoComb_.push_back( ( ele->dr03TkSumPt() + std::max(0., ele->dr03EcalRecHitSumEt() - 1.) + ele->dr03HcalTowerSumEt() ) / ele->p4().Pt() ) ;
//     else
//       patElecRelIsoComb_.push_back( ( ele->dr03TkSumPt() + ele->dr03EcalRecHitSumEt() + ele->dr03HcalTowerSumEt() ) / ele->p4().Pt() ) ;
     
     
     double supercluster_e =-999;
    double sigihih = -999;
    double deletain = -999;
    double delphiin = -999;
    double hoe = -999;
     double supercluster_et =-999;
     double supercluster_eta =-999;
     double supercluster_phi =-999;
    
     if(ele->superCluster().isNonnull()){
       supercluster_et  = ele->superCluster()->energy()/cosh(ele->superCluster()->eta());
       supercluster_eta = ele->superCluster()->eta();
       supercluster_phi = ele->superCluster()->phi();
       supercluster_e   = ele->superCluster()->energy();
    sigihih          = ele->sigmaIetaIeta();
    delphiin         = ele->deltaPhiSuperClusterTrackAtVtx();
    deletain         = ele->deltaEtaSuperClusterTrackAtVtx();
    hoe              = ele->hadronicOverEm();
     }
    
//     //Don't know how this gonna fare in case of track-only PFElectrons.
//     patElecScEn_.push_back(supercluster_e);
    patElecSigIhIh_.push_back(sigihih);
    patElecDelEtaIn_.push_back(deletain);
    patElecDelPhiIn_.push_back(delphiin);
    patElecHoE_.push_back(hoe);
     patElecScEt_.push_back(supercluster_et);
     patElecScEta_.push_back(supercluster_eta);
     patElecScPhi_.push_back(supercluster_phi);
//     patElecisEcalDriven_.push_back(ele->ecalDrivenSeed());
//     patElecisTrackerDriven_.push_back(ele->trackerDrivenSeed());
    patElecMissingHits_.push_back( ele->gsfTrack().get()->trackerExpectedHitsInner().numberOfHits() );
    patElecDist_.push_back(ele->convDist() );
    patElecDeltaCotTheta_.push_back(ele->convDcot() ) ;
//     patElecConvRadius_.push_back(ele->convRadius() ) ;
     patElecDxy_.push_back(ele->userFloat("dxy"));
//     patElecDsz_.push_back(ele->gsfTrack()->dsz());
     patElecDz_.push_back(ele->userFloat("dz"));
//     patElecD0_.push_back(ele.userFloat("dz"));
//     patElecDxyBS_.push_back(ele->gsfTrack()->dxy(beamSpot.position()));
//     patElecDszBS_.push_back(ele->gsfTrack()->dsz(beamSpot.position()));
//     patElecDzBS_.push_back(ele->gsfTrack()->dz(beamSpot.position()));
//     patElecMva_.push_back(ele->pfCandidateRef()->mva_e_pi());
 //  cout <<"Elec p2assOrNot: "<<passOrNot<<endl;

/*
    patElecChHadSumPt03_.push_back(IsoDeposit(*ele,"chHad",0.3,0.07,0.025));
    patElecNeHadSumPt03_.push_back(IsoDeposit(*ele,"neuHad",0.3,0.07,0.025));
    patElecGamSumPt03_.push_back(IsoDeposit(*ele,"gam",0.3,0.07,0.025));
    patElecChHadSumPt04_.push_back(IsoDeposit(*ele,"chHad",0.4,0.07,0.025));
    patElecNeHadSumPt04_.push_back(IsoDeposit(*ele,"neuHad",0.4,0.07,0.025));
    patElecGamSumPt04_.push_back(IsoDeposit(*ele,"gam",0.4,0.07,0.025));
    patElecChHadSumPt05_.push_back(IsoDeposit(*ele,"chHad",0.5,0.07,0.025));
    patElecNeHadSumPt05_.push_back(IsoDeposit(*ele,"neuHad",0.5,0.07,0.025));
    patElecGamSumPt05_.push_back(IsoDeposit(*ele,"gam",0.5,0.07,0.025));
*/
    patElecChHadSumPt03_.push_back(-999);
    patElecNeHadSumPt03_.push_back(-999);
    patElecGamSumPt03_.push_back(-999);
    patElecChHadSumPt04_.push_back(-999);
    patElecNeHadSumPt04_.push_back(-999);
    patElecGamSumPt04_.push_back(-999);
    patElecChHadSumPt05_.push_back(-999);
    patElecNeHadSumPt05_.push_back(-999);
    patElecGamSumPt05_.push_back(-999);


    patElecChHadIso_.push_back(ele->chargedHadronIso());
    patElecNeHadIso_.push_back(ele->neutralHadronIso());
    patElecGamIso_.push_back(ele->photonIso());
    patElecInBarrel_.push_back(ele->isEB());
    patElecInEndcap_.push_back(ele->isEE());
  
   
//cout<<"pss3"<<endl;
  patElectrackMomentumAtVtxMag2_.push_back(ele->trackMomentumAtVtx().mag2());
  patElececalEnergy_.push_back(ele->ecalEnergy());
  patElechasMatchConv_.push_back(ele->userFloat("hasMatchConv"));

//  cout<<"p3"<<endl;
  }//electron loop
}
void
patElecTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}


void
patElecTree::AddBranch(std::vector<double>* vec, std::string name){
  //std::string brName="patElec"+name;
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

void 
patElecTree::AddBranch(double* x, std::string name){
  std::string brName="Electron"+name;
  tree_->Branch(name.c_str(),x,(brName+"/D").c_str());
}

void 
patElecTree::AddBranch(std::vector<std::string>* vec, std::string name){
  std::string brName="patElec"+name;
  tree_->Branch(name.c_str(),vec);
}





void 
patElecTree::SetBranches(){
  AddBranch(&fastJetRho_, "fastJetRho_");
  AddBranch(&PFchsJetRho_, "PFchsJetRho_");
  AddBranch(&lepIsoRho_, "lepIsoRho_");
  AddBranch(&patElecEt_, "patElecEt_");
  AddBranch(&patElecEnergy_, "patElecEnergy_");
  AddBranch(&patElecPt_, "patElecPt_");
  AddBranch(&patElecEta_, "patElecEta_");
  AddBranch(&patElecPhi_, "patElecPhi_");
  AddBranch(&patElecM_, "patElecM_");
  AddBranch(&patElecScEn_, "patElecScEn_");
  AddBranch(&patElecScEt_, "patElecScEt_");
  AddBranch(&patElecScEta_, "patElecScEta_");
  AddBranch(&patElecScPhi_, "patElecScPhi_");
  AddBranch(&patElecisEcalDriven_, "patElecisEcalDriven_");
  AddBranch(&patElecisTrackerDriven_, "patElecisTrackerDriven_");
  AddBranch(&patElecSigIhIh_, "patElecSigIhIh_");
  AddBranch(&patElecDelEtaIn_, "patElecDelEtaIn_");
  AddBranch(&patElecDelPhiIn_, "patElecDelPhiIn_");
  AddBranch(&patElecHoE_, "patElecHoE_");
  AddBranch(&patElecTrkIso_, "patElecTrkIso_");
  AddBranch(&patElecHcalIso_, "patElecHcalIso_");
  AddBranch(&patElecEcalIso_, "patElecEcalIso_");
  AddBranch(&patElecCharge_, "patElecCharge_");
  AddBranch(&patElecRelIsoComb_, "patElecRelIsoComb_");
  AddBranch(&patElecDxy_, "patElecDxy_");
  AddBranch(&patElecD0_, "patElecD0_");
  AddBranch(&patElecDsz_, "patElecDsz_");
  AddBranch(&patElecDz_, "patElecDz_");
  AddBranch(&patElecDxyBS_,"patElecDxyBS_");
  AddBranch(&patElecDszBS_,"patElecDszBS_");
  AddBranch(&patElecDzBS_,"patElecDzBS_");
  //Some variables for particle flow electron selection
  AddBranch(&patElecMva_ ,"patElecMva_");
  //  AddBranch(&patElecNumMisHits_ ,"patElecNumMisHits_");
  AddBranch(&patElecChHadSumPt03_, "patElecChHadSumPt03_");
  AddBranch(&patElecNeHadSumPt03_, "patElecNeHadSumPt03_");
  AddBranch(&patElecGamSumPt03_, "patElecGamSumPt03_");
  AddBranch(&patElecChHadSumPt04_, "patElecChHadSumPt04_");
  AddBranch(&patElecChHadIso_, "patElecChHadIso_");
  AddBranch(&patElecNeHadIso_, "patElecNeHadIso_");
  AddBranch(&patElecGamIso_, "patElecGamIso_");
  AddBranch(&patElecNeHadSumPt04_, "patElecNeHadSumPt04_");
  AddBranch(&patElecGamSumPt04_, "patElecGamSumPt04_");
  AddBranch(&patElecChHadSumPt05_, "patElecChHadSumPt05_");
  AddBranch(&patElecNeHadSumPt05_, "patElecNeHadSumPt05_");
  AddBranch(&patElecGamSumPt05_, "patElecGamSumPt05_");
  
  // conversion rejection
  AddBranch(&patElecMissingHits_, "patElecMissingHits_");
  AddBranch(&patElecDist_, "patElecDist_" );
  AddBranch(&patElecDeltaCotTheta_, "patElecDeltaCotTheta_");
  AddBranch(&patElecConvRadius_, "patElecConvRadius_");
  AddBranch(&patElecInBarrel_,"patElecInBarrel_");
  AddBranch(&patElecInEndcap_,"patElecInEndcap_");
  //addfor VBF
  AddBranch(&patElectrackMomentumAtVtxMag2_,"patElectrackMomentumAtVtxMag2_"); 
  AddBranch(&patElececalEnergy_,"patElececalEnergy_");
  AddBranch(&patElechasMatchConv_,"patElechasMatchConv_");
  AddBranch(&patElecID_,"patElecID_");
}


void
patElecTree::Clear(){
  fastJetRho_.clear();
  PFchsJetRho_.clear();
  lepIsoRho_.clear();
  patElecEt_.clear();
  patElecEnergy_.clear();
  patElecPt_.clear();
  patElecEta_.clear();
  patElecPhi_.clear();
  patElecM_.clear();
  patElecScEn_.clear();
  patElecScEt_.clear();
  patElecScEta_.clear();
  patElecScPhi_.clear();
  patElecisEcalDriven_.clear();
  patElecisTrackerDriven_.clear();
  patElecSigIhIh_.clear();
  patElecDelEtaIn_.clear();
  patElecDelPhiIn_.clear();
  patElecHoE_.clear();
  patElecTrkIso_.clear();
  patElecHcalIso_.clear();
  patElecEcalIso_.clear();
  patElecRelIsoComb_.clear();
  patElecCharge_.clear();
  patElecChi2Ndoff_.clear();
  patElecNhits_.clear();
  patElecQoverP_.clear();
  patElecDxy_.clear();
  patElecD0_.clear();
  patElecDsz_.clear();
  patElecDz_.clear();
  patElecDxyBS_.clear();
  patElecDzBS_.clear();
  patElecDszBS_.clear();
  patElecVx_.clear();
  patElecVy_.clear();
  patElecVz_.clear();
//  patElecGenMatchPt_.clear();
  //patElecGenMatchEta_.clear();
  //patElecGenMatchPhi_.clear();
  //patElecGenMotherId_.clear();
  //patElecGenMatchCharge_.clear();
   patElecChHadSumPt03_.clear();
  patElecNeHadSumPt03_.clear();
  patElecGamSumPt03_.clear();
  patElecChHadSumPt04_.clear();
  patElecChHadIso_.clear();
  patElecNeHadIso_.clear();
  patElecGamIso_.clear();
  patElecNeHadSumPt04_.clear();
  patElecGamSumPt04_.clear();
  patElecChHadSumPt05_.clear();
  patElecNeHadSumPt05_.clear();
  patElecGamSumPt05_.clear();
  patElecMva_.clear();
//  patElecNumMisHits_.clear();
  // conversion rejection
  patElecMissingHits_.clear();
  patElecDist_.clear();
  patElecDeltaCotTheta_.clear();
  patElecConvRadius_.clear();
 
  patElecInBarrel_.clear();
  patElecInEndcap_.clear();

  patElectrackMomentumAtVtxMag2_.clear();
  patElececalEnergy_.clear();
  patElechasMatchConv_.clear();
patElecID_.clear();
}


















  //AddBranch(&patElecScEt_, "patElecScEt_");
  //AddBranch(&patElecCalo_, "patElecCalo_");
  //AddBranch(&patElecChi2Ndoff_, "patElecChi2Ndoff_");
  //AddBranch(&patElecNhits_, "patElecNhits_");
  //AddBranch(&patElecQoverP_, "patElecQoverP_");
  //AddBranch(&patElecDxy_, "patElecDxy_");
  //AddBranch(&patElecD0_, "patElecD0_");
  //AddBranch(&patElecDsz_, "patElecDsz_");
  //AddBranch(&patElecDs_, "patElecDs_");
  //AddBranch(&patElecDxyBS_, "patElecDxyBS_");
  //AddBranch(&patElecDzBS_, "patElecDzBS_");
  //AddBranch(&patElecDszBS_, "patElecDszBS_");
  //AddBranch(&patElecVx_, "patElecVx_");
  //AddBranch(&patElecVy_, "patElecVy_");
  //AddBranch(&patElecVz_, "patElecVz_");
