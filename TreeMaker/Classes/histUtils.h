#ifndef _HISTO_UTILS_H_
#define _HISTO_UTILS_H_


//Note: Keep these classes simple.
//No need to provide =, + or even the
//copy constructors.

#include "particle.h"
#include "lepton.h"
#include "electron.h"
#include "TH1D.h"
#include <map>
#include "prettyHistograms.h"
template<class T>
class partHisto{

public:
  partHisto(std::string token,std::string level,double wt=1);
  virtual ~partHisto();
  virtual void Fill(particle& p);
/*  void SetNBins(std::string key, int nbins, double min, double max){
    histMap_[key]->SetBins(nbins,min,max);
  }*/

protected:
  partHisto();
  double weight;
  TH1D* hPt_;
  //TH1D* hEta_;
  //TH1D* hPhi_;
  //TH1D* hRap_;
  TH1D* hM_;
};


template<class T>
partHisto<T>::partHisto(std::string token, std::string level,double wt){

  weight=wt;
  hPt_=prettyHistograms(hPt_ ,token,"p_{T}",level,"p_{T} [GeV]",40,0,200);
  //hEta_=prettyHistograms(hEta_,token,"PseudoRapidity",level,"#eta",40,-4,4);
  //hPhi_=prettyHistograms(hPhi_,token,"Phi",level,"#Phi",50,-5,+5);
  hM_=prettyHistograms(hM_,token,"Mass",level,"Mass [GeV]",60,60,120);
  //hRap_=prettyHistograms(hRap_,token,"Rapidity",level,"Y",50,-5,5);
}

template<class T>
void partHisto<T>::Fill(particle& p){
  hPt_->Fill(p.Pt(),weight);
  //hEta_->Fill(p.Eta(),weight);
  //hPhi_->Fill(p.Phi(),weight);
  hM_->Fill(p.M(),weight);
  //hRap_->Fill(p.Rapidity(),weight);
}

template<class T>
partHisto<T>::~partHisto(){
   delete hPt_;
   //delete hEta_;
   //delete hPhi_;
   //delete hRap_;
   delete hM_;
}

/*
template<class T>
class elecHisto{

public:
  elecHisto(std::string token,std::string level,double wt=1);
  ~elecHisto();
  void Fill(electron<T>& p)  ;

protected:
  elecHisto(){};
  double weight;
 
  TH1D* hMissHits_;
  TH1D* hDcot_;
  TH1D* hDist_;
  TH1D* hRelTrkIso_;
  TH1D* hRelEcalIso_;
  TH1D* hRelHcalIso_;
  TH1D* hSihih_;
  TH1D* hDphiIn_;
  TH1D* hDetaIn_;
  TH1D* hHoE_;

};

template<class T>
elecHisto<T>::elecHisto(std::string token,std::string level,double wt){
  weight=wt;
  hMissHits_   =prettyHistograms(hMissHits_,token,"Missing-Hits",level,"N_{MisHits}",5,0,5);
  hDcot_       =prettyHistograms(hDcot_,token,"dCot",level,"dCot",100,-0.05,0.05);
  hDist_       =prettyHistograms(hDist_,token,"Dist",level,"Dist",100,-0.05,0.05);
  hSihih_      =prettyHistograms(hSihih_,token,"sigma_{ieta ieta}",level,"#sigma_{i#eta i#eta}",100,0.0,0.1);
  hDphiIn_     =prettyHistograms(hDphiIn_,token,"Delta(Phi)_{in}",level,"#Delta(#Phi)_{in}",100,0.0,0.1);
  hDetaIn_     =prettyHistograms(hDetaIn_,token,"Delta(eta)_{in}",level,"#Delta(#eta)_{in}",100,0.0,0.1);
  hHoE_        =prettyHistograms(hHoE_,token,"HoE",level,"H/E",100,0.0,0.2);
  hRelTrkIso_  =prettyHistograms(hRelTrkIso_ ,token,"RelTkIso",level,"Relative Track Isolation",100,0.0,3.0);
  hRelEcalIso_ =prettyHistograms(hRelEcalIso_,token,"RelEcalIso",level,"Relative ECal Isolation",100,0.0,3.0);
  hRelHcalIso_ =prettyHistograms(hRelHcalIso_,token,"RelHcalIso",level,"Relative HCal Isolation",100,0.0,3.0);
  
}

template<class T>
elecHisto<T>::~elecHisto(){
  delete hMissHits_;
  delete hDcot_;
  delete hDist_;
  delete hRelTrkIso_;
  delete hRelEcalIso_;
  delete hRelHcalIso_;
  delete hSihih_;
  delete hDphiIn_;
  delete hDetaIn_;
  delete hHoE_;
}


template<class T>
void elecHisto<T>::Fill(electron<T>& p) {
  hMissHits_->Fill(p.MissHits(),weight);
  hDcot_->Fill(p.Dcot(),weight);
  hDist_->Fill(p.Dist(),weight);
  hRelTrkIso_->Fill(p.RelTrkIso(),weight);
  hRelEcalIso_->Fill(p.RelEcalIso(),weight);
  hRelHcalIso_->Fill(p.RelHcalIso(),weight);
  hSihih_->Fill(p.Sihih(),weight);
  hDphiIn_->Fill(p.DphiIn(),weight);
  hDetaIn_->Fill(p.DetaIn(),weight);
  hHoE_->Fill(p.HoE(),weight);
}

*/

/*


template<class T>
class jetHisto{

public:
  jetHisto(std::string token,std::string level,double wt=1);
  ~jetHisto();
  void Fill(pfjet<T>& p)  ;

protected:
  jetHisto(){};
  double weight;
  TH1D*  hchf_;
  TH1D*  hnhf_;
  TH1D*  hcmult_;
  TH1D*  hcemf_;
  TH1D*  hnemf_; 
  
};

template<class T>
jetHisto<T>::jetHisto(std::string token,std::string level,double wt){
  weight=wt;
  std::string title = "h_"+token+"chf_"+level;
  hchf_ = new TH1D(title.c_str(),title.c_str(),100,-1,1);

  title = "h_"+token+"nhf_"+level;
  hnhf_ = new TH1D(title.c_str(),title.c_str(),100,-1,1);

  title = "h_"+token+"cmult_"+level;
  hcmult_ = new TH1D(title.c_str(),title.c_str(),100,-1,10);

  title = "h_"+token+"cemf_"+level;
  hcemf_ = new TH1D(title.c_str(),title.c_str(),100,-1,1);

  title = "h_"+token+"nemf_"+level;
  hnemf_ = new TH1D(title.c_str(),title.c_str(),100,-1,1);

}

template<class T>
jetHisto<T>::~jetHisto(){
  delete hchf_;
  delete hnhf_;
  delete hcmult_;
  delete hcemf_;
  delete hnemf_; 

}


template<class T>
void jetHisto<T>::Fill(pfjet<T>& p) {
  hchf_->Fill(p.chf(),weight);
  hnhf_->Fill(p.nhf(),weight);
  hcmult_->Fill(p.cmult(),weight);
  hcemf_->Fill(p.cemf(),weight);
  hnemf_->Fill(p.nemf(),weight);
}

*/

#endif
