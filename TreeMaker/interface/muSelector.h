#ifndef __MUSELECTOR_H__
#define __MUSELECTOR_H__

/*
MuonSelector
Optimal Usage: patMuons.

Lovedeep Kaur,
Panjab University.

Shin-Shan Eiko Yu,
National Central University

*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
// for tracks
#include "DataFormats/TrackReco/interface/Track.h"

// for vertexing
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// for making histograms
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"


class muSelector{
   
 public:
  muSelector(const edm::ParameterSet ps);
  std::map<std::string, bool> CutRecord(const pat::Muon& mu);

  void SetRho(double rho){rho_ = rho;}
  void SetData(bool isData){isData_ = isData;}

  
  ~muSelector(){}
  
  double ptX_;
  double etaX_;
  double dxyX_;
  double dzX_;
  double normalizedChi2X_;
  double trackerHitsX_;
  double pixelHitsX_;
  double muonHitsX_;
  double nMatchesX_;

  double iso1X_;
  double iso2X_;
  double iso3X_;
  double iso4X_;
  
  //some flags.
  bool usePfIso_;
  bool useRelIso_;
  bool useCombIso_;
  bool coneRadius_;
  bool applyGlbTrk_;    
  bool doDelBetaCorr_;   
  bool reqTrigMatch_;
  edm::ParameterSet idPar_;
  edm::ParameterSet isoPar_;

  double rho_;
  double isData_;


};

#endif
