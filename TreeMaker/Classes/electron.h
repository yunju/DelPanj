#ifndef __ELECTRON_HH_
#define __ELECTRON_HH_

#include "lepton.h"

template <class T>
class electron: public lepton<T>{
 public:
  electron();
  electron(T* const tree, unsigned int pos);
  electron(const electron& mu);
  electron& operator =(electron mu);
  bool isBarrel(){return ((fabs(this->eta_)<1.4442));}
  bool isEndcap(){return ((fabs(this->eta_)>1.566 && fabs(this->eta_)<2.5));}
  
  double Mva(){return mva_;}
  double RelTrkIso(){return (this->tIso_/this->pt_);}
  double RelEcalIso(){return (this->eIso_/this->pt_);}
  double RelHcalIso(){return (this->hIso_/this->pt_);}
  double Sihih(){return sihih_;}
  double DetaIn(){return  detain_;}
  double DphiIn() {return dphiin_;}
  double HoE(){return dhoe_;}
  double MissHits() {return missingHits_;}
  double Dist(){return dist_;}
  double Dcot(){return dcot_;}

double Q(){
return q_;
}
  bool MvaX(double mva){return mva_>mva;}
  bool PassConvRej(){
    bool passConvRej=0;
    if(missingHits_<= 0 && fabs(dist_)>=0.02 && fabs(dcot_)>=0.02)
      passConvRej = 1.0; 
    return passConvRej;
  };
  bool PassWP80(){
    bool passWP80=0;
    if((isBarrel() &&
       (
	missingHits_<= 0
	&& fabs(dist_)>0.02
	&& fabs(dcot_)>0.02
 &&	 (this->tIso_/this->pt_) < 0.09
	&& (this->eIso_/this->pt_) < 0.07
	&& (this->hIso_/this->pt_) < 0.10
	&& sihih_<0.01
	 && fabs(dphiin_)<0.06
	 && fabs(detain_)<0.004
	 && dhoe_<0.04
	 ))
       ||(isEndcap() &&
	  (
	   missingHits_<= 0
	   && fabs(dist_)>0.02
	   && fabs(dcot_)>0.02
&&	    (this->tIso_/this->pt_) < 0.04
	   && (this->eIso_/this->pt_) < 0.05
	   && (this->hIso_/this->pt_) < 0.025
	   && sihih_<0.03
	        && fabs(dphiin_)<0.03
	    && fabs(detain_)<0.007
	     && dhoe_<0.15  //0.025
	    )
	  ) )   
       passWP80 = 1.0; 
    return passWP80;
  };

  bool PassWP95(){
    bool passWP95=0;
    if((isBarrel() &&
       (
	missingHits_<= 1
	//&& fabs(dist_)>=0.02
	//&& fabs(dcot_)>=0.02
	&& (this->tIso_/this->pt_) < 0.15
	&& (this->eIso_/this->pt_) < 2.0
	&& (this->hIso_/this->pt_) < 0.12
	&&sihih_<0.01
	//&& fabs(dphiin_)<=0.06
	&& fabs(detain_)<=0.007
	&& dhoe_<0.15
	))
       ||(isEndcap() &&
	  (
	   missingHits_<= 1
	   //&& fabs(dist_)>=0.02
	   //&& fabs(dcot_)>=0.02
	   && (this->tIso_/this->pt_) < 0.08
	   && (this->eIso_/this->pt_) < 0.06
	   && (this->hIso_/this->pt_) < 0.05
	   &&sihih_<0.03
	   //&& fabs(dphiin_)<=0.01
	      && fabs(detain_)<=0.01
	   && dhoe_<0.07
	  )
	  )  )     
       passWP95 = 1.0; 
    return passWP95;
  };


  private:
  void createFromTree(T* const tree, unsigned int pos);
  double  scen_;
  double  sihih_;
  double  detain_;
  double  dphiin_;
  double  dhoe_;
  double  mva_;
  double  missingHits_;
  double  dist_;
  double  dcot_;
  double eIso_;
double tIso_;
double hIso_;
double q_;
};


template<class T>
electron<T>::electron():lepton<T>(){}

template<class T>
electron<T>::electron(T* const tree, unsigned int pos):lepton<T>(){createFromTree(tree,pos);}

template <class T>
electron<T>::electron(const electron<T>& el):lepton<T>(el){}

template <class T>
electron<T>& electron<T>::operator =(electron<T>  el){
  if (this != &el) {//guard against self assignment.
    this->lepton<T>::operator=(el);
  }
  return *this;
}
  
template <class T>
void electron<T>::createFromTree(T* const tree, unsigned int pos){

  //  pos_=pos;
  double pt=tree->patElecPt_->at(pos);
  double et=tree->patElecEt_->at(pos);
  double eta= tree->patElecEta_->at(pos) ;
  double phi= tree->patElecPhi_->at(pos) ;
  double m= tree->patElecM_->at(pos) ;

  this->particle::init(pt,eta,phi,m,0);
  
  this->q_=tree->patElecCharge_->at(pos);
/*  
  this->lepton<T>::normChi2_= 99;//tree->patElecChi2Ndoff_->at(pos)  ;
  this->lepton<T>::nHits_= 99;//tree->patElectronNumExpHits_->at(pos) ;
  this->lepton<T>::corrDxy_= 99;//tree->patElecDxyBS_->at(pos) ;
  this->lepton<T>::dxy_= 99;//tree->patElecDxy_->at(pos)  ;
 */ 
  this->tIso_= tree->patElecTrkIso_->at(pos)  ;
  this->eIso_= tree->patElecEcalIso_->at(pos)  ;
  this->hIso_= tree->patElecHcalIso_->at(pos)  ;
  /*
  this->lepton<T>::charHadSumPt03_ = tree->patElecChHadSumPt03_->at(pos);
  this->lepton<T>::neutHadSumPt03_ = tree->patElecNeHadSumPt03_->at(pos);
  this->lepton<T>::gammaSumPt03_   = tree->patElecGamSumPt03_->at(pos);
  
  this->lepton<T>::charHadSumPt04_ = tree->patElecChHadSumPt04_->at(pos);
  this->lepton<T>::neutHadSumPt04_ = tree->patElecNeHadSumPt04_->at(pos);
  this->lepton<T>::gammaSumPt04_   = tree->patElecGamSumPt04_->at(pos);
  
  this->lepton<T>::charHadSumPt05_ = tree->patElecChHadSumPt05_->at(pos);
  this->lepton<T>::neutHadSumPt05_ = tree->patElecNeHadSumPt05_->at(pos);
  this->lepton<T>::gammaSumPt05_   = tree->patElecGamSumPt05_->at(pos);
*/

  scen_   = tree->patElecScEn_->at(pos);
  sihih_  = tree->patElecSigIhIh_->at(pos);
  detain_ = tree->patElecDelEtaIn_->at(pos);
  dphiin_ = tree->patElecDelPhiIn_->at(pos);
  dhoe_   = tree->patElecHoE_->at(pos);
  mva_    = tree->patElecMva_->at(pos);
  missingHits_ = tree->patElecMissingHits_->at(pos);
  dist_ = tree->patElecDist_->at(pos);
  dcot_ = tree->patElecDeltaCotTheta_->at(pos);
}
#endif
