/*
The Isolation information is quite large and voluminous
so we store it here as an independent tree
*/

#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/eSelector.h"
using namespace std;
using namespace edm;

class patElecIsoTree{

 public:
  patElecIsoTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~patElecIsoTree();
  void Fill(const edm::Event& iEvent);
  double IsoDeposit(const pat::Electron ele,std::string type,double radius,double coneradius,double detastrip);
  void SetBranches();
  void Clear();
  //  void ConversionRejection(const edm::Event&, const edm::EventSetup&, const pat::Electron&);

 private:
//TTree* tree_;
  //Dont Allow User to Call the Default Constructor.
  patElecIsoTree();


  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<std::string>*, std::string name);
  edm::InputTag dcsTag_;
  edm::InputTag patElecLabel_;
  edm::InputTag beamSpotLabel_;
  edm::InputTag JetForElecTree_;

  bool isData_;
  int Nelecs;
  eSelector ewp1_;//for lead ele.
 
 

 //pt, eta, phi, M : Enough to caluclate
  //px,py,pz etc.==>No need to store later
  TTree* tree_;
  double nele; 

  std::vector<double> fastJetRho_;
  std::vector<double> pfChsJetRho_;
  std::vector<double> lepIsoRho_;


  std::vector<double> patElecPt_;
  std::vector<double> patElecEta_;
  std::vector<double> patElecTrkIso_;
  std::vector<double> patElecHcalIso_;
  std::vector<double> patElecEcalIso_;

  std::vector<double> patElecCombPFIso02_;
  std::vector<double> patElecCombPFIsoDBeta02_;
  std::vector<double> patElecCombPFIso03_;
  std::vector<double> patElecCombPFIsoDBeta03_;
  std::vector<double> patElecCombPFIso04_;
  std::vector<double> patElecCombPFIsoDBeta04_;
  std::vector<double> patElecCombPFIso05_;
  std::vector<double> patElecCombPFIsoDBeta05_;
  std::vector<double> patElecCombPFIso06_;
  std::vector<double> patElecCombPFIsoDBeta06_;
  std::vector<double> patElecCombPFIso07_;
  std::vector<double> patElecCombPFIsoDBeta07_;
  std::vector<double> patElecCombPFIso08_;
  std::vector<double> patElecCombPFIsoDBeta08_;






  /*
  std::vector<double> patElecCharHadIso03_;
  std::vector<double> patElecNeuHadIso03_;
  std::vector<double> patElecGammaIso03_;
  std::vector<double> patElecPfPUCharHadIso03_;
  std::vector<double> patElecCharHadVetoIso03_;
  std::vector<double> patElecNeuHadVetoIso03_;
  std::vector<double> patElecGammaVetoIso03_;
  std::vector<double> patElecPfPUCharHadVetoIso03_;

  std::vector<double> patElecCharHadIso04_;
  std::vector<double> patElecNeuHadIso04_;
  std::vector<double> patElecGammaIso04_;
  std::vector<double> patElecPfPUCharHadIso04_;
  std::vector<double> patElecCharHadVetoIso04_;
  std::vector<double> patElecNeuHadVetoIso04_;
  std::vector<double> patElecGammaVetoIso04_;
  std::vector<double> patElecPfPUCharHadVetoIso04_;

  std::vector<double> patElecCharHadIso05_;
  std::vector<double> patElecNeuHadIso05_;
  std::vector<double> patElecGammaIso05_;
  std::vector<double> patElecPfPUCharHadIso05_;
  std::vector<double> patElecCharHadVetoIso05_;
  std::vector<double> patElecNeuHadVetoIso05_;
  std::vector<double> patElecGammaVetoIso05_;
  std::vector<double> patElecPfPUCharHadVetoIso05_;

  std::vector<double> patElecCharHadIso06_;
  std::vector<double> patElecNeuHadIso06_;
  std::vector<double> patElecGammaIso06_;
  std::vector<double> patElecPfPUCharHadIso06_;
  std::vector<double> patElecCharHadVetoIso06_;
  std::vector<double> patElecNeuHadVetoIso06_;
  std::vector<double> patElecGammaVetoIso06_;
  std::vector<double> patElecPfPUCharHadVetoIso06_;

  std::vector<double> patElecCharHadIso07_;
  std::vector<double> patElecNeuHadIso07_;
  std::vector<double> patElecGammaIso07_;
  std::vector<double> patElecPfPUCharHadIso07_;
  std::vector<double> patElecCharHadVetoIso07_;
  std::vector<double> patElecNeuHadVetoIso07_;
  std::vector<double> patElecGammaVetoIso07_;
  std::vector<double> patElecPfPUCharHadVetoIso07_;


  std::vector<double> patElecCharHadIso08_;
  std::vector<double> patElecNeuHadIso08_;
  std::vector<double> patElecGammaIso08_;
  std::vector<double> patElecPfPUCharHadIso08_;
  std::vector<double> patElecCharHadVetoIso08_;
  std::vector<double> patElecNeuHadVetoIso08_;
  std::vector<double> patElecGammaVetoIso08_;
  std::vector<double> patElecPfPUCharHadVetoIso08_;
  */
};
