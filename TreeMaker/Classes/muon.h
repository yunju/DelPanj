#ifndef __MUON_HH_
#define __MUON_HH_

#include "lepton.h"

//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
template <class T>
class muon: public lepton<T>{
 public:
  muon();
  muon(T* const tree, unsigned int pos);
  muon(const muon& mu);
  muon& operator =(muon mu);  
  double Q(){return q_;}
 private:
  void createFromTree(T* const tree, unsigned int pos);
  double q_;
};


//===============================================================
//===============================================================
template<class T>
muon<T>::muon():lepton<T>(){}

//===============================================================
//===============================================================
template<class T>
muon<T>::muon(T* const tree, unsigned int pos):lepton<T>(){
  createFromTree(tree,pos);
}

//===============================================================
//===============================================================
template <class T>
muon<T>::muon(const muon<T>& el):lepton<T>(el){}

//===============================================================
//===============================================================
template <class T>
muon<T>& muon<T>::operator =(muon<T>  el){
  if (this != &el) {//guard against self assignment.
    this->lepton<T>::operator=(el);
  }
  return *this;
}
  
//===============================================================
//===============================================================
template <class T>
void muon<T>::createFromTree(T* const tree, unsigned int pos){
  double pt=tree->patMuonPt_->at(pos);
  double eta= tree->patMuonEta_->at(pos) ;
  double phi= tree->patMuonPhi_->at(pos) ;
  double m= tree->patMuonM_->at(pos) ;
  q_ = tree->patMuonCharge_->at(pos);
  //this->lepton<T>::pos_=pos;  
  this->lepton<T>::init(pt,eta,phi,m,0);
/*
  this->lepton<T>::q_=tree->patMuonCharge_->at(pos);
  
  this->lepton<T>::normChi2_= tree->patMuonChi2Ndoff_->at(pos)  ;
  this->lepton<T>::nHits_= tree->patMuonNhits_->at(pos) ;
  this->lepton<T>::corrDxy_= tree->patMuonDxyBS_->at(pos) ;
  this->lepton<T>::dxy_= tree->patMuonDxy_->at(pos)  ;
  
  this->lepton<T>::tIso_= tree->patMuonTrkIso_->at(pos)  ;
  this->lepton<T>::eIso_= tree->patMuonEcalIso_->at(pos)  ;
  this->lepton<T>::hIso_= tree->patMuonHcalIso_->at(pos)  ;
  
  this->lepton<T>::charHadSumPt03_= tree->patMuonChHadSumPt03_->at(pos);
  this->lepton<T>::neutHadSumPt03_= tree->patMuonNeHadSumPt03_->at(pos);
  this->lepton<T>::gammaSumPt03_= tree->patMuonGamSumPt03_->at(pos);
  
  this->lepton<T>::charHadSumPt04_= tree->patMuonChHadSumPt04_->at(pos);
  this->lepton<T>::neutHadSumPt04_= tree->patMuonNeHadSumPt04_->at(pos);
  this->lepton<T>::gammaSumPt04_= tree->patMuonGamSumPt04_->at(pos);
  
  this->lepton<T>::charHadSumPt05_= tree->patMuonChHadSumPt05_->at(pos);
  this->lepton<T>::neutHadSumPt05_= tree->patMuonNeHadSumPt05_->at(pos);
  this->lepton<T>::gammaSumPt05_= tree->patMuonGamSumPt05_->at(pos);
*/
}

#endif
