#ifndef __ESELECTOR_H__
#define __ESELECTOR_H__

/*
  ElectronSelector
  Optimal Usage: pat Electrons.

  Anil Singh,
  Panjab University.
  Shin-Shan Eiko Yu,
  National Central University
*/



#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class eSelector{
   
 public:
  eSelector(const edm::ParameterSet ps);  
  std::map<std::string, bool> CutRecord(const pat::Electron& e);
  
  void SetRho(double rho){rho_ = rho;}
  void SetData(bool isData){isData_ = isData;}
  
  ~eSelector(){}
  
  double ptX_;
  double etaX_;

  double detainBrlX_;
  double delphiBrlX_;  
  double sieieBrlX_;
  double hoeBrlX_;
  double d0vtxBrlX_;
  double dzvtxBrlX_;
  double ooemoopBrlX_;
  double distBrlX_;
  double dcotBrlX_;
  double hasConvBrlX_;
  double nmisHitBrlX_;

  double detainEcpX_;
  double delphiEcpX_;  
  double sieieEcpX_;
  double hoeEcpX_;
  double d0vtxEcpX_;
  double dzvtxEcpX_;
  double ooemoopEcpX_;
  double distEcpX_;
  double dcotEcpX_;
  double hasConvEcpX_;
  double nmisHitEcpX_;

  
  double iso1BrlX_;
  double iso2BrlX_;
  double iso3BrlX_;
  double iso4BrlX_;
  
  double iso1EcpX_;
  double iso2EcpX_;
  double iso3EcpX_;
  double iso4EcpX_;
  
  
   
  //some flags.
  bool usePfIso_;
  bool useRelIso_;
  bool useCombIso_;
  bool coneRadius_;
  bool doRhoCorr_;   
  
  //some pset declaration
  edm::ParameterSet idBrl_;
  edm::ParameterSet idEcp_;
  edm::ParameterSet isoBrl_;
  edm::ParameterSet isoEcp_; 

  double rho_;
  double isData_;

};

#endif
