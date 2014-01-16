#ifndef __PARTICLE_HH_
#define __PARTICLE_HH_

//This is a general class that contain 
//kinematics of the particles.
//We can apply the selection cuts on it.
//We can add two particles to obtain 
//the invariant mass.

#include<iostream>
#include<vector>
//include following to provide "fabs(double)"
#include <math.h>
#include <cmath>
#include "TLorentzVector.h"
#include "JetUtilMC.h"

class particle{

 public:

  //Constructors
  particle();
  particle(double pt, double eta, double phi,
	   double m, double e);
  particle(const particle& p);

  //Accessors
  double  Pt(){return pt_;}
  double  Eta(){return eta_;}
  double  AbsEta(){return fabs(eta_);}
  double  Rapidity(){return rap_;}
  double  Energy(){return e_;}
  double  Phi(){return phi_;}
  double  Px(){return px_;}
  double  Py(){return py_;}
  double  Pz(){return pz_;}
  double  M(){return m_;}
  double  Theta(){return theta_;}

  //bools to make cuts
  bool  PtX(double pt){return pt_>pt;}
  bool  EtaX(double eta){return ((fabs(eta_)<eta));}
  bool  elecEtaX(double eta){return ((fabs(eta_)<eta) 
				   && !(fabs(eta_)>1.4442 && fabs(eta_)<1.566));}
  bool  EtaEBX(double eta){return ((fabs(eta_)<=eta));}
  bool  EtaEEX(double eta){return ((fabs(eta_)>=1.566 && fabs(eta_)<eta));}
  bool  PhiX(double phi){return fabs(phi_)<phi;}
  bool  MX(double m1, double m2){return (m_>m1&&m_<m2);}

  //SomeCalculations
  double deltaR(const particle p);
  double deltaRMin(std::vector<particle>& pvec);
  //double deltaTheta(const particle p);
  void speak();
  


  //Operators
  particle& operator=(const particle&);
  particle operator+(const particle&);



 protected:
  void init(double pt, double eta, double phi,
	    double m, double e);
  double pt_;
  double eta_;
  double phi_;
  double px_;
  double py_;
  double pz_;
  double m_;
  double e_;
  double rap_;
  double theta_;

};

//---------------------------------------------
//---------------------------------------------
particle::particle():
  pt_(0.),eta_(0.),phi_(0.),
  px_(0.),py_(0.),pz_(0.),m_(0.),e_(0.),rap_(0),
  theta_(0)
{

}

//---------------------------------------------
//---------------------------------------------
particle::particle(double pt, double eta, double phi,
		   double m, double e){
  init( pt,eta,phi,m,e);
}

//---------------------------------------------
//---------------------------------------------
particle::particle(const particle& p):
  pt_(p.pt_),eta_(p.eta_),phi_(p.phi_),
  px_(p.px_),py_(p.py_),pz_(p.pz_),m_(p.m_),
  e_(p.e_),rap_(p.rap_),theta_(p.theta_)
{
}

//---------------------------------------------
//---------------------------------------------
particle& 
particle::operator=(const particle& p){
  if(this==&p)return *this;
  init( p.pt_, p.eta_,p.phi_,
        p.m_,p.e_);
  return *this;
}


//---------------------------------------------
//---------------------------------------------
particle 
particle::operator+(const particle& p){

  TLorentzVector* lw = new TLorentzVector();
  lw->SetPtEtaPhiM(p.pt_,p.eta_,p.phi_,p.m_);
  
  TLorentzVector* nw = new TLorentzVector();
  nw->SetPtEtaPhiM(pt_,eta_,phi_,m_);
  
  TLorentzVector* rl = new TLorentzVector();
  (*rl) = (*lw)+(*nw);
  
  
  particle par;
  double e = rl->E();
  double phi = rl->Phi();
  double pt = rl->Pt();  
  double eta = rl->Eta();
  double m = rl->M();


  par.init(pt,eta,phi,m,e);

    return par;
}


//---------------------------------------------
//---------------------------------------------
double 
particle::deltaR(const particle p){
  // double deleta = eta_-p.eta_;
  // double delphi = phi_-p.phi_;
  double dR = radius(eta_,phi_,p.eta_,p.phi_);
  return dR;
}

//---------------------------------------------
//---------------------------------------------
double 
particle::deltaRMin(std::vector<particle>& pvec){
  double delRMin=9999;
  for(unsigned int i=0; i!=pvec.size();i++){
    double delR = deltaR(pvec[i]);
    if(delRMin>delR)delRMin=delR;
  }
  return delRMin;
}


//---------------------------------------------
//---------------------------------------------
void 
particle::init(double pt, double eta, double phi,
	       double m, double e){

  TLorentzVector* p = new TLorentzVector();

  p->SetPtEtaPhiM(pt,eta,phi,m);

  pt_=pt;
  eta_=eta;
  phi_=phi;
  px_=p->Px();
  py_=p->Py();
  pz_=p->Pz();
  m_=m;
  e_=p->E();
  rap_=p->Rapidity();
  theta_=p->Theta();
  phi_=p->Phi();
  e=p->Energy();//just to avoid warning
}


//---------------------------------------------
//---------------------------------------------
void 
particle::speak(){
  std::cout
    <<"pt "<<pt_<<"  "
    <<"eta "<<eta_<<"  "
    <<"phi "<<phi_<<"  "
    <<"px "<<px_<<"  "
    <<"py "<<py_<<"  "
    <<"pz "<<pz_<<"  "
    <<"m " <<m_<<"  "
    <<std::endl;
}
#endif
