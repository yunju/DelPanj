#ifndef __PFJET_HH_
#define __PFJET_HH_



#include <cmath>
#include <string>
#include <sstream>
#include "TLorentzVector.h"
#include "particle.h"

template<class T>
class pfjet : public particle{

 public:
  pfjet();
  pfjet(const pfjet<T>& j);
  pfjet(T* const tree, unsigned int pos);
  pfjet& operator =(pfjet<T>  j);


  double chf(){return chf_;}
  double nhf(){return nhf_;}
  double cmult(){return cmult_;}
  double cemf(){return cemf_;}
  double nemf(){return nemf_;}

  bool  loosePfJetIdX(){return ((pt_>10)&&(AbsEta()<2.4)
				&&(chf_>0)&&(nhf_<0.99)&&(cmult_>0.0)
				&&(cemf_<0.99)&&(nemf_<0.99));}
  bool  tightPfJetIdX(){return 1;}
  void ShiftPt(double percent){pt_=pt_+pt_*(percent/100.);}


 private:
  double chf_;
  double nhf_;
  double cmult_;
  double cemf_;
  double nemf_;
};

template<class T>
pfjet<T>::pfjet():
  particle(),
  chf_(0),
  nhf_(0),
  cmult_(0),
  cemf_(0),
  nemf_(0){}

//===============================================================

template<class T>
pfjet<T>::pfjet(const pfjet<T>& j)
  :particle(j),
   chf_(j.chf_),
   nhf_(j.nhf_),
   cmult_(j.cmult_),
   cemf_(j.cemf_),
   nemf_(j.nemf_){}


template <class T>
pfjet<T>& pfjet<T>::operator =(pfjet<T>  j){
  if (this != &j) {//guard against self assignment.
    this->particle::operator=(j);
    chf_=j.chf_;
    nhf_=j.nhf_;
    cmult_=j.cmult_;
    cemf_=j.cemf_;
    nemf_=j.nemf_;
  }
  return *this;
}
  
  
  template<class T>
  pfjet<T>::pfjet(T* const tree, unsigned int pos):particle(){
    //  pos_=pos;
    pt_= tree->patJetPfAk05Pt_->at(pos);
    eta_= tree->patJetPfAk05Eta_->at(pos) ;
    phi_= tree->patJetPfAk05Phi_->at(pos) ;
    m_= tree->patJetPfAk05M_->at(pos) ;
    chf_=tree->patJetPfAk05CharHadEFr_->at(pos);
    nhf_=tree->patJetPfAk05NeutHadEFr_->at(pos);
    cmult_=tree->patJetPfAk05CharMulti_->at(pos);
    cemf_=tree->patJetPfAk05CharEmEFr_->at(pos);
    nemf_=tree->patJetPfAk05NeutEmEFr_->at(pos);
  }

#endif
