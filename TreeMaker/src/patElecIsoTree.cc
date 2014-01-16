#include "DelPanj/TreeMaker/interface/patElecIsoTree.h"

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

double GetCombinedPFIso(pat::Electron& ele, double rad){
  reco::isodeposit::AbsVetos Vetos;
  reco::isodeposit::ThresholdVeto *thveto = new reco::isodeposit::ThresholdVeto(0.5);
  Vetos.push_back(thveto);
  float chDep =(ele.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(rad,Vetos, false ).first);
  float neDep = (ele.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(rad,Vetos, false ).first);
  float gamDep= (ele.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(rad,Vetos, false ).first);
  return (chDep+neDep+gamDep);
}

double GetCombinedPFIsoBetaCorr(pat::Electron& ele, double rad){
  reco::isodeposit::AbsVetos Vetos;
  reco::isodeposit::ThresholdVeto *thveto = new reco::isodeposit::ThresholdVeto(0.5);
  Vetos.push_back(thveto);
  float chDep =(ele.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(rad,Vetos, false ).first);
  float neDep = (ele.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(rad,Vetos, false ).first);
  float gamDep= (ele.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(rad,Vetos, false ).first);
  reco::isodeposit::AbsVetos Vetos1;
  float puDep =(ele.isoDeposit(pat::PfPUChargedHadronIso)->depositAndCountWithin(rad,Vetos1, false ).first);
       
  float two = neDep+gamDep-0.5*puDep;
  if(two<0)two = 0; 
 
  return (chDep+two);
}

patElecIsoTree::patElecIsoTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  ewp1_(iConfig.getParameter<edm::ParameterSet>("leadElecPset_"))
{
 tree_=tree; 
  patElecLabel_ = iConfig.getParameter<edm::InputTag>("patElectrons");
  beamSpotLabel_ = iConfig.getParameter<edm::InputTag> ("beamSpotLabel");
  JetForElecTree_ = iConfig.getParameter<edm::InputTag> ("JetForElecTreePY");
SetBranches();
}


patElecIsoTree::~patElecIsoTree(){
  delete tree_;
  
} 
void
patElecIsoTree::Fill(const edm::Event& iEvent){
  Clear();
  edm::Handle<pat::ElectronCollection> patElecHandle;
  if(not iEvent.getByLabel(patElecLabel_,patElecHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<patElecLabel_<<std::endl; exit(0);}
  
//  edm::Handle<reco::BeamSpot> beamSpotHandle;
//  iEvent.getByLabel(beamSpotLabel_, beamSpotHandle);
//  reco::BeamSpot beamSpot = (*beamSpotHandle);
  
   double fastJetRho=-999999.9;
   double pfChsJetRho=-999999.9;
   double  lepIsoRho =  -999999.9;
 
  
  /////// Pileup density "rho" in the event from fastJet pileup calculation /////
  edm::Handle<double> rho;
  const edm::InputTag eventrho("kt6PFJets", "rho");
  iEvent.getByLabel(eventrho,rho);
  if( *rho == *rho) fastJetRho = *rho;
  else  fastJetRho =  -999999.9;
  fastJetRho_.push_back(fastJetRho);
  
  /////// Pileup density "rho" in the event from PFchs pileup calculation /////
  edm::Handle<double> rhochs;
  const edm::InputTag eventrhochs("kt6PFJetsAK5","rho");
  /*
  const edm::InputTag eventrhochs(,"rho");
  
 if(not iEvent.getByLabel(JetForElecTree_,rhochs)){
  std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
               <<JetForElecTree_<<std::endl; exit(0);}
*/


  iEvent.getByLabel(eventrhochs,rhochs);
  if( *rhochs == *rhochs) pfChsJetRho = *rhochs;
  else  pfChsJetRho =  -999999.9;
  pfChsJetRho_.push_back(pfChsJetRho);
  
  /////// Pileup density "rho" for lepton isolation subtraction /////
  edm::Handle<double> rhoLepIso;
  const edm::InputTag eventrhoLepIso("kt6PFJetsForIsolation", "rho");
  iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
  if( *rhoLepIso == *rhoLepIso) lepIsoRho = *rhoLepIso;
  else  lepIsoRho =  -999999.9;
  lepIsoRho_.push_back(lepIsoRho);
  
  pat::ElectronCollection eleColl(*(patElecHandle.product()));
  std::sort(eleColl.begin(),eleColl.end(),PtGreater());
  
  nele=eleColl.size();
  pat::ElectronCollection::const_iterator ele;
  bool index = 0;
  for(ele=eleColl.begin(); ele!=eleColl.end(); ele++){
    //Apply The Electron Selection here...
    pat::Electron e(*ele);
    std::map<std::string, bool> leadPass = ewp1_.CutRecord(e);
    if(!(PassAll(leadPass)))continue; //supply void isolation cuts from outside.
    patElecPt_.push_back(ele->p4().pt());
    patElecEta_.push_back(ele->eta());
    patElecTrkIso_.push_back(ele->dr03TkSumPt());
    patElecHcalIso_.push_back(ele->dr03HcalTowerSumEt());
    patElecEcalIso_.push_back(ele->dr03EcalRecHitSumEt());
    patElecCombPFIso02_.push_back(GetCombinedPFIso(e,0.2));
    patElecCombPFIso03_.push_back(GetCombinedPFIso(e,0.3));
    patElecCombPFIso04_.push_back(GetCombinedPFIso(e,0.4));
    patElecCombPFIso05_.push_back(GetCombinedPFIso(e,0.5));
    patElecCombPFIso06_.push_back(GetCombinedPFIso(e,0.6));
    patElecCombPFIso07_.push_back(GetCombinedPFIso(e,0.7));
    patElecCombPFIso08_.push_back(GetCombinedPFIso(e,0.8));
    patElecCombPFIsoDBeta02_.push_back(GetCombinedPFIsoBetaCorr(e,0.2));
    patElecCombPFIsoDBeta03_.push_back(GetCombinedPFIsoBetaCorr(e,0.3));
    patElecCombPFIsoDBeta04_.push_back(GetCombinedPFIsoBetaCorr(e,0.4));
    patElecCombPFIsoDBeta05_.push_back(GetCombinedPFIsoBetaCorr(e,0.5));
    patElecCombPFIsoDBeta06_.push_back(GetCombinedPFIsoBetaCorr(e,0.6));
    patElecCombPFIsoDBeta07_.push_back(GetCombinedPFIsoBetaCorr(e,0.7));
    patElecCombPFIsoDBeta08_.push_back(GetCombinedPFIsoBetaCorr(e,0.8));
  }
}

void
patElecIsoTree::AddBranch(std::vector<double>* vec, std::string name){
  //std::string brName="patElec"+name;
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

void 
patElecIsoTree::AddBranch(double* x, std::string name){
  std::string brName="Electron"+name;
  tree_->Branch(name.c_str(),x,(brName+"/D").c_str());
}

void 
patElecIsoTree::AddBranch(std::vector<std::string>* vec, std::string name){
  std::string brName="patElec"+name;
  tree_->Branch(name.c_str(),vec);
}



void 
patElecIsoTree::SetBranches(){

  AddBranch(&fastJetRho_, "fastJetRho_");
  AddBranch(&pfChsJetRho_, "pfChsJetRho_");
  AddBranch(&lepIsoRho_, "lepIsoRho_");
  AddBranch(&patElecPt_, "patElecPt_");
  AddBranch(&patElecEta_, "patElecEta_");
  AddBranch(&patElecTrkIso_, "patElecTrkIso_");
  AddBranch(&patElecHcalIso_, "patElecHcalIso_");
  AddBranch(&patElecEcalIso_, "patElecEcalIso_");
  AddBranch(&patElecCombPFIso02_,"patElecCombPFIso02_");
  AddBranch(&patElecCombPFIsoDBeta02_,"patElecCombPFIsoDBeta02_");
  AddBranch(&patElecCombPFIso03_,"patElecCombPFIso03_");
  AddBranch(&patElecCombPFIsoDBeta03_,"patElecCombPFIsoDBeta03_");
  AddBranch(&patElecCombPFIso04_,"patElecCombPFIso04_");
  AddBranch(&patElecCombPFIsoDBeta04_,"patElecCombPFIsoDBeta04_");
  AddBranch(&patElecCombPFIso05_,"patElecCombPFIso05_");
  AddBranch(&patElecCombPFIsoDBeta05_,"patElecCombPFIsoDBeta05_");
  AddBranch(&patElecCombPFIso06_,"patElecCombPFIso06_");
  AddBranch(&patElecCombPFIsoDBeta06_,"patElecCombPFIsoDBeta06_");
  AddBranch(&patElecCombPFIso07_,"patElecCombPFIso07_");
  AddBranch(&patElecCombPFIsoDBeta07_,"patElecCombPFIsoDBeta07_");
  AddBranch(&patElecCombPFIso08_,"patElecCombPFIso08_");
  AddBranch(&patElecCombPFIsoDBeta08_,"patElecCombPFIsoDBeta08_");

}


void
patElecIsoTree::Clear(){
  fastJetRho_.clear();
  pfChsJetRho_.clear();
  lepIsoRho_.clear();
  patElecPt_.clear();
  patElecEta_.clear();
  patElecTrkIso_.clear();
  patElecHcalIso_.clear();
  patElecEcalIso_.clear();
  patElecCombPFIso02_.clear();
  patElecCombPFIsoDBeta02_.clear();
  patElecCombPFIso03_.clear();
  patElecCombPFIsoDBeta03_.clear();
  patElecCombPFIso04_.clear();
  patElecCombPFIsoDBeta04_.clear();
  patElecCombPFIso05_.clear();
  patElecCombPFIsoDBeta05_.clear();
  patElecCombPFIso06_.clear();
  patElecCombPFIsoDBeta06_.clear();
  patElecCombPFIso07_.clear();
  patElecCombPFIsoDBeta07_.clear();
  patElecCombPFIso08_.clear();
  patElecCombPFIsoDBeta08_.clear();
  /*
  patElecCharHadIso03_.clear();
  patElecNeuHadIso03_.clear();
  patElecGammaIso03_.clear();
  patElecPfPUCharHadIso03_.clear();
  patElecCharHadVetoIso03_.clear();
  patElecNeuHadVetoIso03_.clear();
  patElecGammaVetoIso03_.clear();
  patElecPfPUCharHadVetoIso03_.clear();
  patElecCharHadIso04_.clear();
  patElecNeuHadIso04_.clear();
  patElecGammaIso04_.clear();
  patElecPfPUCharHadIso04_.clear();
  patElecCharHadVetoIso04_.clear();
  patElecNeuHadVetoIso04_.clear();
  patElecGammaVetoIso04_.clear();
  patElecPfPUCharHadVetoIso04_.clear();
  patElecCharHadIso04_.clear();
  patElecNeuHadIso04_.clear();
  patElecGammaIso04_.clear();
  patElecPfPUCharHadIso04_.clear();
  patElecCharHadVetoIso04_.clear();
  patElecNeuHadVetoIso04_.clear();
  patElecGammaVetoIso04_.clear();
  patElecPfPUCharHadVetoIso04_.clear();
  patElecCharHadIso05_.clear();
  patElecNeuHadIso05_.clear();
  patElecGammaIso05_.clear();
  patElecPfPUCharHadIso05_.clear();
  patElecCharHadVetoIso05_.clear();
  patElecNeuHadVetoIso05_.clear();
  patElecGammaVetoIso05_.clear();
  patElecPfPUCharHadVetoIso05_.clear();
  patElecCharHadIso06_.clear();
  patElecNeuHadIso06_.clear();
  patElecGammaIso06_.clear();
  patElecPfPUCharHadIso06_.clear();
  patElecCharHadVetoIso06_.clear();
  patElecNeuHadVetoIso06_.clear();
  patElecGammaVetoIso06_.clear();
  patElecPfPUCharHadVetoIso06_.clear();
  patElecCharHadIso07_.clear();
  patElecNeuHadIso07_.clear();
  patElecGammaIso07_.clear();
  patElecPfPUCharHadIso07_.clear();
  patElecCharHadVetoIso07_.clear();
  patElecNeuHadVetoIso07_.clear();
  patElecGammaVetoIso07_.clear();
  patElecPfPUCharHadVetoIso07_.clear();
  patElecCharHadIso08_.clear();
  patElecNeuHadIso08_.clear();
  patElecGammaIso08_.clear();
  patElecPfPUCharHadIso08_.clear();
  patElecCharHadVetoIso08_.clear();
  patElecNeuHadVetoIso08_.clear();
  patElecGammaVetoIso08_.clear();
  patElecPfPUCharHadVetoIso08_.clear();
  */
}



