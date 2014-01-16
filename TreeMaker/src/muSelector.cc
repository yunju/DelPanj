#include "DelPanj/TreeMaker/interface/muSelector.h"
#include "DelPanj/TreeMaker/interface/MuonEffectiveArea.h"
/*
  This works according to the official muon recipe in the following page:
  https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon_selection/
*/


muSelector::muSelector(const edm::ParameterSet ps):      
  ptX_  (ps.getParameter<double>("pt")),
  etaX_ (ps.getParameter<double>("eta")),
  reqTrigMatch_ (ps.getParameter<bool>("requireTrigMatch")),
  idPar_(ps.getParameter<edm::ParameterSet> ("idPar")),
  isoPar_(ps.getParameter<edm::ParameterSet> ("isoPar")),
  
  //isolation parameters 
  usePfIso_(ps.getParameter<bool>("usePfIso")),  
  useRelIso_(ps.getParameter<bool>("useRelIso")),
  useCombIso_(ps.getParameter<bool>("useCombIso")),
  coneRadius_(ps.getParameter<double>("coneRad")),
  doDelBetaCorr_(ps.getParameter<bool>("delBCorr"))  
{
  normalizedChi2X_=(idPar_.getParameter<double>("normalizedChi2"));
  muonHitsX_=(idPar_.getParameter<double>("muonHits"));
  nMatchesX_=(idPar_.getParameter<double>("nMatches"));
  dxyX_ =(idPar_.getParameter<double>("dxy"));
  dzX_ =(idPar_.getParameter<double>("dz"));

  trackerHitsX_=(idPar_.getParameter<double>("trackerHits"));
  pixelHitsX_=(idPar_.getParameter<double>("pixelHits"));
  
  
  if(!useCombIso_){
    iso1X_ =(isoPar_.getParameter<double>("iso1"));
    iso2X_ =(isoPar_.getParameter<double>("iso2"));
    iso3X_ =(isoPar_.getParameter<double>("iso3"));
  }  
  else 
    iso4X_ =(isoPar_.getParameter<double>("iso4"));

  isData_ = true;
  rho_    = 0;
}


std::map<std::string, bool>  muSelector::CutRecord(const pat::Muon& mu){
  
  std::map<std::string, bool> cuts;
  double pt = mu.pt();
  double eta = mu.eta();
  
   
  //Kinematic and Fiducial Acceptance
  cuts["ptx"]         = pt>ptX_;
  cuts["etax"]        = fabs(eta)<etaX_;
  
  //MuonIdentification
  cuts["dxy"]         = fabs(mu.dB())<dxyX_;

  //
  if(dzX_<0.0)
    cuts["dz"]        = 1;
  else // userfloat exists
    cuts["dz"]        = fabs(mu.userFloat("dzVtx")) < dzX_;

  bool isPFMuon      = mu.isPFMuon();
  bool isTrackMuon   = mu.isTrackerMuon();
  bool isGlobalMuon  = mu.isGlobalMuon();

  if(isTrackMuon && isGlobalMuon && isPFMuon){//If global as well as tracker muon...
    cuts["trkLayers"]     = mu.innerTrack()->hitPattern().trackerLayersWithMeasurement()>trackerHitsX_;
    cuts["pixelHits"]     = mu.innerTrack()->hitPattern().numberOfValidPixelHits()>pixelHitsX_;
    cuts["muonHits"]      = mu.globalTrack()->hitPattern().numberOfValidMuonHits()>muonHitsX_;
    cuts["nMatches"]      = mu.numberOfMatchedStations()>nMatchesX_;
    cuts["normalizedChi2"]= mu.globalTrack()->normalizedChi2()<normalizedChi2X_; 
  }
  else{
    cuts["trkLayers"]      = 0;
    cuts["pixelHits"]      = 0;
    cuts["muonHits"]       = 0;
    cuts["nMatches"]       = 0;    
    cuts["normalizedChi2"] = 0;
  }
  

  //Muon Isolation... configurable.
  double iso1 = 999.;
  double iso2 = 999.;
  double iso3 = 999.;
  double iso4 = 999.;
  double isoPU = -999.; 
  
  //If user want pfIso
  if(usePfIso_){   
    iso1  =  mu.pfIsolationR04().sumChargedHadronPt;
    iso2  =  mu.pfIsolationR04().sumNeutralHadronEt;
    iso3  =  mu.pfIsolationR04().sumPhotonEt;
    isoPU =  mu.pfIsolationR04().sumPUPt;    
  }
  //If user does not want pfIso
  else {
    iso1   = mu.trackIso();
    iso2   = mu.ecalIso(); 
    iso3   = mu.hcalIso(); 
  }
  
  //Calculate combined Iso
  iso4  = (iso1+iso2+iso3);
  if(doDelBetaCorr_)
    {  //If user wants betacorrected
      iso4 = iso1 + std::max(iso3+iso2-0.5*isoPU,0.);
    }
  else  // do the rho correction
    {
      MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = 
	MuonEffectiveArea::kMuEAData2012;

      MuonEffectiveArea::MuonEffectiveAreaType effAreaType_= 
	MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;
      double Area = MuonEffectiveArea::GetMuonEffectiveArea
	(effAreaType_, fabs(eta), effAreaTarget_);
      iso4 = iso1 + std::max(0.0, iso2+iso3-rho_*Area);
    }
  

  //If user wants relative iso
  if(useRelIso_){
    iso1 = iso1/std::max(0.1,pt);
    iso2 = iso2/std::max(0.1,pt);
    iso3 = iso3/std::max(0.1,pt);
    iso4 = iso4/std::max(0.1,pt);
  }

  if(useCombIso_){
    cuts["iso1"]   = 1;
    cuts["iso2"]   = 1;
    cuts["iso3"]   = 1;
    cuts["iso4"]   = iso4 < iso4X_;
  }
  else{
    cuts["iso1"]   = iso1 < iso1X_;
    cuts["iso2"]   = iso2 < iso2X_;
    cuts["iso3"]   = iso3 < iso3X_;
    cuts["iso4"]   = 1;
  }
  
  //Do the trigger Match stuff.
  if(reqTrigMatch_){
    cuts["isTrigMatch"]= mu.triggerObjectMatches().size() > 0 ;
  }
  else 
    cuts["isTrigMatch"]= 1;
  
  return cuts;
}
