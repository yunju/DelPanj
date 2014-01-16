#include "PFJetSmearFactory.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include <iostream>

using namespace std;


cmg::PFJetSmearFactory::event_ptr cmg::PFJetSmearFactory::create(const edm::Event& iEvent, 
						       const edm::EventSetup& iES) {
	
  
  cmg::PFJetSmearFactory::event_ptr result =PFJetFactory_.create(iEvent,iES);
  JetCorrectionUncertainty *jecUnc_=0;  
  if(applyScale_&&jecUnc_==0){
    //if applyJEC unc but name of txt file is empty, read them from the db
    
    if(!applyScaleDB_){
      jecUnc_=new JetCorrectionUncertainty(scaleFile_.c_str());
    } 
    else{
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      iES.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      jecUnc_ = new JetCorrectionUncertainty(JetCorPar);
    }
  }//end if applyScale and 1st event
      


  for(cmg::PFJetFactory::collection::iterator mi = result->begin();
      mi != result->end(); ++mi){
    
    
    if(applyResolution_){
      bool  doGaussian = true;
      JetResolution ptResol(resolutionFile_, doGaussian);
      
      float geneta = mi->eta() ;
      float genpt = mi->pt();  
      float jetpt = genpt;
      if(genpt > 10){ // no reliable resolution numbers that low
	//std::cout << "eta/pt " << geneta <<" " << genpt << std::endl;
	TF1* fPtResol = ptResol.resolutionEtaPt(geneta,genpt);
	fPtResol->Print("v");
	float rndm = fPtResol->GetRandom();
	float etafactor = 0;
	if(fabs(geneta) <1.1)//numbers from JME-10-010
	  etafactor = 0.04;
	else if(fabs(geneta) <1.7)
	  etafactor = 0.02;
	else if(fabs(geneta) <2.3)
	  etafactor = 0.09;
	else
	  etafactor = 0.04;
	
	if(resolutionOverride_ > 0.)
	  etafactor = resolutionOverride_;
	float smearfactor = fabs(1.- (1.-rndm)*etafactor) ;
	//std::cout << rndm  << " "<<smearfactor << std::endl;
	jetpt = smearfactor*genpt;
	const math::PtEtaPhiMLorentzVector jetp4(jetpt, mi->eta(), mi->phi(), mi->mass());
	//const math::PtEtaPhiMLorentzVector jetp4(mi->pt(), mi->eta(), mi->phi(), mi->mass());
	mi->setP4(jetp4);
      }
    }
    
    if(applyScale_){   
      double jetPt = mi->pt();
      double jetEta = mi->eta();
      jecUnc_->setJetEta(jetEta);
      jecUnc_->setJetPt(jetPt);
      double unc = jecUnc_->getUncertainty(true);
      float smearfactor = unc*nSigmaScale_ ;
      //      std::cout<<"Smearing factor = "<<smearfactor<<std::endl;
      jetPt = (1+smearfactor)*jetPt;
      
      const math::PtEtaPhiMLorentzVector jetp4 (jetPt, mi->eta(), mi->phi(), mi->mass());
      mi->setP4(jetp4);  
    }//end if applyScale
    
  }
  delete jecUnc_;//too bad we cannot declare it as data member...
  return result;
}
