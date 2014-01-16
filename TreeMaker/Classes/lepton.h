#ifndef __LEPTON_HH_
#define __LEPTON_HH_

//Class Lepton is common for both the
//electron and muon objects. The dynamic
//binding has to be used for the constructor
//that extract information from the tree.

#include "particle.h"

template <class T>
class lepton: public particle{

 public:
  lepton(){}
/*lepton(T* const tree, unsigned int pos){createFromTree(tree,pos);}
  lepton(const lepton& lep);
  lepton& operator=(lepton lep);
  virtual ~lepton(){};
    

  double Q(){return q_;}
  double Nhits(){return nHits_;}
  double NChi2(){return normChi2_;}
  double Dxy(){return dxy_;}
  double DxyCorr(){return corrDxy_;}
  double Charge(){return q_;}
  double IsoSumPt(std::string type, double radius);
  double IsoDeposit(std::string type, double radius);

  bool NHitsX(unsigned int hits){return nHits_>hits;}  
  bool NormChi2X(unsigned int chix){return normChi2_<chix;}  
  bool DxyCorrX(double dxycorrx)
  {return (abs(corrDxy_)<dxycorrx||abs(corrDxy_)==dxycorrx);}  
  bool IsoSumPtX(std::string type, double radius, double);
  bool IsoDepositX(std::string type, double radius, double);
  
  // protected:
  
  virtual void createFromTree(T* const tree, unsigned int pos){};
 
  double q_;//charge

  //Tracking and It's Quality
  double normChi2_;//normalized chi2 for global fit.
  double nHits_;//number of hits in the inner track.
  double corrDxy_;//transverse impact parameter
  double dxy_;//transverse impact paramete wrt nominal IP

  //Isolation Variables.
  double tIso_;
  double eIso_;
  double hIso_;
  
  double charHadSumPt03_;
  double charHadSumPt04_;
  double charHadSumPt05_;
  
  double neutHadSumPt03_;
  double neutHadSumPt04_;
  double neutHadSumPt05_;
  

  double gammaSumPt03_;
  double gammaSumPt04_;
  double gammaSumPt05_;
*/
};
/*
template<class T>
lepton<T>::lepton():
  particle(),
  q_(0),
  normChi2_(999),
  nHits_(0),
  corrDxy_(999),
  dxy_(999),
  
  tIso_(999),eIso_(999),hIso_(999),	     
  charHadSumPt03_(999),charHadSumPt04_(999),charHadSumPt05_(999),
  neutHadSumPt03_(999),neutHadSumPt04_(999),neutHadSumPt05_(999),
  gammaSumPt03_(999),gammaSumPt04_(999),gammaSumPt05_(999){}

template <class T>
lepton<T>& lepton<T>::operator =(lepton<T>  lep){
  this-> pos_= lep.pos_;
  this-> pt_= lep.pt_;
  this-> eta_= lep.eta_;
  this-> phi_= lep.phi_;
  this-> m_= lep.m_;
  this-> q_= lep.q_;
  
  this-> normChi2_= lep.normChi2_;
  this-> nHits_= lep.nHits_;
  this-> corrDxy_= lep.corrDxy_;
  this-> dxy_= lep.dxy_;
  
  this-> tIso_= lep.tIso_;
  this-> eIso_= lep.eIso_;
  this-> hIso_= lep.hIso_;
  
  this-> charHadSumPt03_= lep.charHadSumPt03_;
  this-> charHadSumPt04_= lep.charHadSumPt04_;
  this-> charHadSumPt05_= lep.charHadSumPt05_;
  
  this-> neutHadSumPt03_= lep.neutHadSumPt03_;
  this-> neutHadSumPt04_= lep.neutHadSumPt04_;
  this-> neutHadSumPt05_= lep.neutHadSumPt05_;
  
  
  this-> gammaSumPt03_= lep.gammaSumPt03_;
  this-> gammaSumPt04_= lep.gammaSumPt04_;
  this-> gammaSumPt05_= lep.gammaSumPt05_;

  return *this;
}

template <class T>
double lepton<T>::IsoSumPt(std::string type, double radius){

  if(type=="neutHad"){
    if(radius==0.3)return neutHadSumPt03_;
    else if(radius==0.4)return neutHadSumPt04_;
    else if(radius==0.5)return neutHadSumPt05_;
    else {
      std::cout<<"NeutHadSumPt("<<radius<<") not stored."<<std::endl;
      std::cout<<"Returning NeutHadSumPt(0.4) instead."<<std::endl;
      return neutHadSumPt03_;
    }
  }
  
  else if(type=="charHad"){
    if(radius==0.3)return charHadSumPt03_;
    else if(radius==0.4)return charHadSumPt04_;
    else if(radius==0.5)return charHadSumPt05_;
    else {
      std::cout<<"CharHadSumPt("<<radius<<") not stored."<<std::endl;
      std::cout<<"Returning CharHadSumPt(0.4) instead."<<std::endl;
      return charHadSumPt04_;
      
    }
  }
  
  else if(type=="gamma"){
    if(radius==0.3)return gammaSumPt03_;
    else if(radius==0.4)return gammaSumPt04_;
    else if(radius==0.5)return gammaSumPt05_;
    else {
      std::cout<<"GammaSumPt("<<radius<<") not stored."<<std::endl;
      std::cout<<"Returning GammaSumPt(0.4) instead."<<std::endl;
      return gammaSumPt03_;
    }
  }
  
  else  if(type=="combRelSumPt"){
    if(radius==0.3)
      return (gammaSumPt03_+charHadSumPt03_+neutHadSumPt03_)/pt_;
    else if(radius==0.4)
      return (gammaSumPt04_+charHadSumPt04_+neutHadSumPt04_)/pt_;
    else if(radius==0.5)
      return (gammaSumPt05_+charHadSumPt05_+neutHadSumPt05_)/pt_;
    else {
      std::cout<<"combRelSumPt("<<radius<<") not stored."<<std::endl;
      std::cout<<"Returning combRelSumPt(0.4) instead."<<std::endl;
      return (gammaSumPt04_+charHadSumPt04_+neutHadSumPt04_)/pt_;
    }
  }
  else return 999999;
}


template <class T>
bool lepton<T>::IsoSumPtX(std::string type, double radius, double iso){
  double isolation = IsoSumPt(type,radius);
  return isolation < iso;
}*/
#endif

