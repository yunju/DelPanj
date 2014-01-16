#include "DelPanj/TreeMaker/interface/eSelector.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

eSelector::eSelector(const edm::ParameterSet ps):      
  ptX_  (ps.getParameter<double>("ptx")),
  etaX_ (ps.getParameter<double>("etax")),
  //identification parameters
  idBrl_(ps.getParameter<edm::ParameterSet> ("idBrl")),   
  idEcp_(ps.getParameter<edm::ParameterSet> ("idEcp")),
  isoBrl_(ps.getParameter<edm::ParameterSet> ("isoBrl")),
  isoEcp_(ps.getParameter<edm::ParameterSet> ("isoEcp")),
  //isolation parameters 
  usePfIso_(ps.getParameter<bool>("usePfIso")),  
  useRelIso_(ps.getParameter<bool>("useRelIso")),
  useCombIso_(ps.getParameter<bool>("useCombIso")),
  coneRadius_(ps.getParameter<double>("coneRad")),
  doRhoCorr_(ps.getParameter<bool>("rhoCorr"))
{     
  detainBrlX_= (idBrl_.getParameter<double>("detain"));
  delphiBrlX_= (idBrl_.getParameter<double>("delphi"));  
  sieieBrlX_ = (idBrl_.getParameter<double>("sieie"));
  hoeBrlX_   = (idBrl_.getParameter<double>("hoe"));
  d0vtxBrlX_ = (idBrl_.getParameter<double>("d0vtx"));
  dzvtxBrlX_ = (idBrl_.getParameter<double>("dzvtx"));
  ooemoopBrlX_= (idBrl_.getParameter<double>("ooemoop"));
  distBrlX_   = (idBrl_.getParameter<double>("dist"));
  dcotBrlX_   = (idBrl_.getParameter<double>("dcot"));
  hasConvBrlX_= (idBrl_.getParameter<double>("hasConv"));
  nmisHitBrlX_= (idBrl_.getParameter<double>("nmisHit"));


  detainEcpX_= (idEcp_.getParameter<double>("detain"));
  delphiEcpX_= (idEcp_.getParameter<double>("delphi"));  
  sieieEcpX_ = (idEcp_.getParameter<double>("sieie"));
  hoeEcpX_   = (idEcp_.getParameter<double>("hoe"));
  d0vtxEcpX_ = (idEcp_.getParameter<double>("d0vtx"));
  dzvtxEcpX_ = (idEcp_.getParameter<double>("dzvtx"));
  ooemoopEcpX_= (idEcp_.getParameter<double>("ooemoop"));
  distEcpX_   = (idEcp_.getParameter<double>("dist"));
  dcotEcpX_   = (idEcp_.getParameter<double>("dcot"));
  hasConvEcpX_= (idEcp_.getParameter<double>("hasConv"));
  nmisHitEcpX_= (idEcp_.getParameter<double>("nmisHit"));

  
  if(!useCombIso_){
    iso1BrlX_ =(isoBrl_.getParameter<double>("iso1"));
    iso2BrlX_ =(isoBrl_.getParameter<double>("iso2"));
    iso3BrlX_ =(isoBrl_.getParameter<double>("iso3"));
  }else 
    iso4BrlX_ =(isoBrl_.getParameter<double>("iso4"));
  if(!useCombIso_){
    iso1EcpX_ =(isoEcp_.getParameter<double>("iso1"));
    iso2EcpX_ =(isoEcp_.getParameter<double>("iso2"));
    iso3EcpX_ =(isoEcp_.getParameter<double>("iso3"));
  }else
    iso4EcpX_ =(isoEcp_.getParameter<double>("iso4"));

  isData_ = true;
  rho_    = 0;

}


std::map<std::string, bool> 
eSelector::CutRecord(const pat::Electron& e){

  std::map<std::string, bool> cuts;
// cout<<"cu :"<<cuts["tx"]<<endl;
       
  //good idea to call the eta and pt methods for once.
  double pt = e.pt();
  double eta = e.superCluster()->eta();
  bool ingap= fabs(eta)>=1.4442 && fabs(eta)<= 1.566;
  //       std::cout<<"isoBrlX: "<< dcotBrlX_<<std::endl;
  if(ingap)
    {
      cuts["ptx"]    = 0;
      cuts["etax"]   = 0;
      cuts["detain"] = 0;
      cuts["delphi"] = 0;
      cuts["sieie"]  = 0;
      cuts["hoe"]    = 0;
      cuts["d0vtx"]  = 0;
      cuts["dzvtx"]  = 0;
      cuts["ooemoop"] = 0;
      cuts["dist"] = 0;	
      cuts["dcot"] = 0;
      cuts["hasConv"] = 0;
      cuts["nmshits"] = 0;
      cuts["iso1"]   = 0;
      cuts["iso2"]   = 0;
      cuts["iso3"]   = 0;
      cuts["iso4"]   = 0;
      cuts["trig"]   = 0;     
       return cuts;
    }
  
  //Default values must fail any sensible IsoCut
  double iso1 = 999.;
  double iso2 = 999.;
  double iso3 = 999.;
  double iso4 = 999.;
  
  bool passTriggerTight = 
    (e.dr03TkSumPt()/std::max(0.1,pt) < 0.2) &&
    (e.dr03EcalRecHitSumEt()/std::max(0.1,pt) < 0.2) &&
    (e.dr03HcalTowerSumEt()/std::max(0.1,pt) < 0.2);
  

  if(usePfIso_){  
    // will be fixed later, Eiko
    //std::cout<<"Using Particle Flow Iso"<<std::endl;
    iso1  =  e.chargedHadronIso(); 
    iso2  =  e.neutralHadronIso();
    iso3  =  e.photonIso();    
  }
  else {
    iso1   = e.dr03TkSumPt();
    iso2   = e.dr03EcalRecHitSumEt();
    iso3   = e.dr03HcalTowerSumEt();
  }
      
  //Calculate combined Iso
  iso4  = (iso1+iso2+iso3);

  if(doRhoCorr_)
    {
     //old effective area at 2011
     /*
      ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = 
  	isData_? ElectronEffectiveArea::kEleEAData2011:
  	ElectronEffectiveArea::kEleEAFall11MC;

      
      ElectronEffectiveArea::ElectronEffectiveAreaType effAreaType_ =
	ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03;      
 
      double AEff = ElectronEffectiveArea::GetElectronEffectiveArea
	(effAreaType_, fabs(eta), effAreaTarget_);
      iso4 = iso1 + std::max(0.0, iso2+iso3-rho_*AEff);
      
      */
      // for2013 effective area pf iso 
      float EffectiveArea;
      if(abs(e.superCluster()->eta())<1.0 )                                               EffectiveArea = 0.13 ;
      if(abs(e.superCluster()->eta())>1.0     && abs(e.superCluster()->eta())<1.479 )     EffectiveArea = 0.14 ;
      if(abs(e.superCluster()->eta())>1.479   && abs(e.superCluster()->eta())<2.0 )       EffectiveArea = 0.07 ;
      if(abs(e.superCluster()->eta())>2.0     && abs(e.superCluster()->eta())<2.2 )       EffectiveArea = 0.09 ;
      if(abs(e.superCluster()->eta())>2.2     && abs(e.superCluster()->eta())<2.3 )       EffectiveArea = 0.11 ;
      if(abs(e.superCluster()->eta())>2.3     && abs(e.superCluster()->eta())<2.4 )       EffectiveArea = 0.11 ;
      if(abs(e.superCluster()->eta())>2.4 )                                               EffectiveArea = 0.14 ; 



      iso4 = iso1 + std::max(0.0, iso2+iso3-max(rho_,0.0)*EffectiveArea);
      //cout<<iso4<<" "<<iso1<<" "<<iso2<<" "<<iso3<<" "<<rho_<<" "<<EffectiveArea<<endl;
    }
  else 
    iso4 = (iso1+iso2+iso3);
  
  /*@anil: Add the choice for doing rho correction as well*/

  //If user wants relative iso
  if(useRelIso_){
    iso1 = iso1/std::max(0.1,pt);
    iso2 = iso2/std::max(0.1,pt);
    iso3 = iso3/std::max(0.1,pt);
    iso4 = iso4/std::max(0.1,pt);
  }

  //  cout << " p = " << sqrt(e.trackMomentumAtVtx().mag2()) << "\t" << 
  //    e.ecalEnergy() / e.eSuperClusterOverP() << endl;
  double ooemoop = fabs(
			1.0/std::max((double)1e-3,(double)e.ecalEnergy()) - 
			1.0/std::max((double)1e-3,(double)sqrt(e.trackMomentumAtVtx().mag2()))
			);


  cuts["ptx"]    = pt>ptX_;
  cuts["etax"]   = !(ingap)&&fabs(eta)<etaX_;
  cuts["trig"]   =passTriggerTight;
     
  if(e.isEB()){

    cuts["detain"] = fabs(e.deltaEtaSuperClusterTrackAtVtx())<detainBrlX_;
    cuts["delphi"] = fabs(e.deltaPhiSuperClusterTrackAtVtx())<delphiBrlX_;
    cuts["sieie"]  = fabs(e.sigmaIetaIeta())<sieieBrlX_;
    cuts["hoe"]    = e.hadronicOverEm() <hoeBrlX_;

    cuts["d0vtx"]  = d0vtxBrlX_ < 0 ? 1: 
      (fabs(e.userFloat("dxy")) < d0vtxBrlX_);
    cuts["dzvtx"]  = dzvtxBrlX_ < 0 ? 1: 
      (fabs(e.userFloat("dz")) < dzvtxBrlX_);

    cuts["ooemoop"] = ooemoop < ooemoopBrlX_;
// We're using userFloat("hasMatchConv") now and the hasConvBrlX is set to be 1e-6.
// Only e.userFloat("hasMatchConv"))<hasConvBrlX_ is matter now.	 
    cuts["dcot"]   = hasConvBrlX_ < 0 ? 
      !(fabs(e.convDcot())< dcotBrlX_ && fabs(e.convDist())< distBrlX_):1;
    cuts["dist"]   = hasConvBrlX_ < 0 ? 
      !(fabs(e.convDcot())< dcotBrlX_ && fabs(e.convDist())< distBrlX_):1;
    cuts["hasConv"] = hasConvBrlX_ > 0? 
      (e.userFloat("hasMatchConv"))<hasConvBrlX_:1;
    cuts["nmshits"]= e.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits()<=nmisHitBrlX_;
	 
    if(useCombIso_){
      cuts["iso1"]   = 1;
      cuts["iso2"]   = 1;
      cuts["iso3"]   = 1;
  //    cuts["iso4"]   = iso4 < iso4BrlX_ && passTriggerTight;
      cuts["iso4"]   = iso4 < iso4BrlX_ ;//the trigger is using only for barrel before, is it correct? 
    }
    else{
      cuts["iso1"]   = iso1 < iso1BrlX_;
      cuts["iso2"]   = iso2 < iso2BrlX_;
      cuts["iso3"]   = iso3 < iso3BrlX_;
      cuts["iso4"]   = 1;

    }
	 
  }
       
  else if(e.isEE()){

    cuts["detain"] = fabs(e.deltaEtaSuperClusterTrackAtVtx())<detainEcpX_;
    cuts["delphi"] = fabs(e.deltaPhiSuperClusterTrackAtVtx())<delphiEcpX_;
    cuts["sieie"]  = fabs(e.sigmaIetaIeta())<sieieEcpX_;
    cuts["hoe"]    = e.hadronicOverEm() < hoeEcpX_;

    cuts["d0vtx"]  = d0vtxEcpX_ < 0 ? 1: 
      (fabs(e.userFloat("dxy")) < d0vtxEcpX_);
    cuts["dzvtx"]  = dzvtxEcpX_ < 0 ? 1: 
      (fabs(e.userFloat("dz")) < dzvtxEcpX_);

    cuts["ooemoop"] = ooemoop < ooemoopEcpX_;

    cuts["dcot"]   = hasConvEcpX_ < 0 ? 
      !(fabs(e.convDcot())< dcotEcpX_ && fabs(e.convDist())< distEcpX_):1;
    cuts["dist"]   = hasConvEcpX_ < 0 ? 
      !(fabs(e.convDcot())< dcotEcpX_ && fabs(e.convDist())< distEcpX_):1;
    cuts["hasConv"] = hasConvEcpX_ > 0? 
      (e.userFloat("hasMatchConv"))<hasConvEcpX_:1;
    cuts["nmshits"]= e.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits()<=nmisHitEcpX_;

    if(useCombIso_){
      cuts["iso1"]   = 1;
      cuts["iso2"]   = 1;
      cuts["iso3"]   = 1;
      cuts["iso4"]   = iso4 < iso4EcpX_;
    }
    else{
      cuts["iso1"]   = iso1 < iso1EcpX_;
      cuts["iso2"]   = iso2 < iso2EcpX_;
      cuts["iso3"]   = iso3 < iso3EcpX_;
      cuts["iso4"]   = 1;
    }

  }
    
     
  return cuts;
}
