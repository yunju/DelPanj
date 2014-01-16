/// @file
/// File containing the definition of the methods associated to the class.
///

#include "DelPanj/TreeMaker/interface/ZZTree.hh"
#include "DelPanj/TreeMaker/interface/cutvalues.h"
//#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "DelPanj/TreeMaker/interface/MuonEffectiveArea.h"

// system include files

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/CandUtils/interface/CenterOfMassBooster.h"
#include "PhysicsTools/CandUtils/interface/Booster.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// ROOT classes
#include <TMath.h>
#include <TVector3.h>
#include <algorithm>
#include <fstream>
#include <Math/VectorUtil.h>
#include <TLegend.h>
#include <TCanvas.h>


// for Btagging scale factors
#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include <time.h>

typedef std::vector< edm::Handle< edm::ValueMap<double> > >             
IsoDepositVals;


void
ZZTree::AddBranch(int* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}


void
ZZTree::AddBranch(double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}

//---------------------------------------------------
//---------------------------------------------------
void
ZZTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

void
ZZTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}


//---------------------------------------------------
//---------------------------------------------------
void 
ZZTree::AddBranchArray(const int arraySize, double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,Form("%s[%d]/D",brName.data(),arraySize));
}

//---------------------------------------------------------------
//---------------------------------------------------------------

ZZTree::ZZTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  e2012ID_ ( iConfig.getParameter<edm::ParameterSet>("e2012IDSet")),
  mu2012ID_ ( iConfig.getParameter<edm::ParameterSet>("mu2012IDSet")),
  hzzeejj_(iConfig.getParameter<edm::InputTag>("hzzeejjTag")),
  hzzmmjj_ (iConfig.getParameter<edm::InputTag>("hzzmmjjTag")),
  eleRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("eleRhoIso")),
  muoRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("muoRhoIso"))//,
{
  tree_=tree; 
  SetBranches();

  
  // the second argument is the random seed, any reason to set it 
  // differently or the same for the 3 taggers
  srand ( time(NULL) );
  int seed = rand(); //712687782, 727743360
  std::cout << "seed = " << seed << std::endl;
  btsfutiljp = new BTagSFUtil("JP", seed);

}


ZZTree::~ZZTree()
{
  delete tree_;
  delete btsfutiljp;
}


// ------------ method called to for each event  ------------
void ZZTree::Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup)
{
  Clear();
  int eventNum = iEvent.id().event();
  
  //============================================================================
  // 
  //       OBTAIN EVENT-LEVEL VARIABLES
  //  
  //============================================================================
  
  bool isData = iEvent.isRealData();

  // missing Et significance
  edm::Handle<pat::METCollection> met_H;
  iEvent.getByLabel("patMETsAK5", met_H);
  metSig_ = (&(met_H->front()))->significance();

  // rho for electron
  edm::Handle<double> ele_rho_event;
  iEvent.getByLabel(eleRhoIsoInputTag_,ele_rho_event);
  eleRho_ = *(ele_rho_event.product());
  e2012ID_.SetData(isData);
  e2012ID_.SetRho(eleRho_);

  // rho for muon
  edm::Handle<double> muo_rho_event;
  iEvent.getByLabel(muoRhoIsoInputTag_,muo_rho_event);
  muoRho_ = *(muo_rho_event.product());

  mu2012ID_.SetData(isData);
  mu2012ID_.SetRho(muoRho_);


  //============================================================================
  // 
  //   B-TAGGING Scale Factors for DB setup
  //  
  //============================================================================


  ///////////////////////////////
  ///   Begin DB setup
  ///////////////////////////////

  //// This is needed for the DB

  std::map<std::string,PerformanceResult::ResultType> measureMap;
  measureMap["BTAGBEFF"]=PerformanceResult::BTAGBEFF;
  measureMap["BTAGBERR"]=PerformanceResult::BTAGBERR;
  measureMap["BTAGCEFF"]=PerformanceResult::BTAGCEFF;
  measureMap["BTAGCERR"]=PerformanceResult::BTAGCERR;
  measureMap["BTAGLEFF"]=PerformanceResult::BTAGLEFF;
  measureMap["BTAGLERR"]=PerformanceResult::BTAGLERR;
  measureMap["BTAGNBEFF"]=PerformanceResult::BTAGNBEFF;
  measureMap["BTAGNBERR"]=PerformanceResult::BTAGNBERR;
  measureMap["BTAGBEFFCORR"]=PerformanceResult::BTAGBEFFCORR;
  measureMap["BTAGBERRCORR"]=PerformanceResult::BTAGBERRCORR;
  measureMap["BTAGCEFFCORR"]=PerformanceResult::BTAGCEFFCORR;
  measureMap["BTAGCERRCORR"]=PerformanceResult::BTAGCERRCORR;
  measureMap["BTAGLEFFCORR"]=PerformanceResult::BTAGLEFFCORR;
  measureMap["BTAGLERRCORR"]=PerformanceResult::BTAGLERRCORR;
  measureMap["BTAGNBEFFCORR"]=PerformanceResult::BTAGNBEFFCORR;
  measureMap["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  measureMap["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  measureMap["MUEFF"]=PerformanceResult::MUEFF;
  measureMap["MUERR"]=PerformanceResult::MUERR;
  measureMap["MUFAKE"]=PerformanceResult::MUFAKE;
  measureMap["MUEFAKE"]=PerformanceResult::MUEFAKE;

  edm::ESHandle<BtagPerformance> perfH;
  std::vector<std::string> measureName;
  std::vector<std::string> measureType;

  // Define which Btag and Mistag algorithm you want to use. These are not user defined and need to be exact
  measureName.push_back("MISTAGJPM");
  measureName.push_back("MUJETSWPBTAGJPM");
  measureName.push_back("MISTAGJPL");
  measureName.push_back("MUJETSWPBTAGJPL");
  measureName.push_back("MISTAGJPM");
  measureName.push_back("MISTAGJPL");

  // Tell DB you want the SF. These are not user defined and need to be exact
  measureType.push_back("BTAGLEFFCORR");
  measureType.push_back("BTAGBEFFCORR");
  measureType.push_back("BTAGLEFFCORR");
  measureType.push_back("BTAGBEFFCORR");
  measureType.push_back("BTAGLEFF");
  measureType.push_back("BTAGLEFF");

  // These are user defined maps that we will use to store the SF
  std::map<std::string,float> ScaleFactors[2];   //store the Btag and Mistag SF for jet0 and jet 1
  std::map<std::string,float> ScaleFactorsEff[2];   //store the Mistag eff for jet0 and jet 1

  ///////////////////////////////
  ///   End DB setup
  ///////////////////////////////


  //============================================================================
  // 
  //    NOW WE START LOOKING FOR HIGGS CANDIDATES
  //  
  //============================================================================

  //initialize variables
  int hcand = -1;
  int bestHCandIndex = -1;

  //LOOP IN 2 KINDS OF LEPTONS: ELECTRONS AND MUONS 
  for (int ilep=0; ilep<2; ilep++) {
    hcand=-1;

    //GET THE HIGGS->ZZ->LLJJ COLLECTION
    Handle<std::vector<pat::CompositeCandidate> > hzzlljj;
    if(ilep==0) iEvent.getByLabel(hzzeejj_, hzzlljj);
    else iEvent.getByLabel(hzzmmjj_, hzzlljj);

    // LOOP OVER HIGGS CANDIDATES

    for(unsigned i=0; i<hzzlljj->size(); i++){
      const pat::CompositeCandidate & h = (*hzzlljj)[i];	

      const reco::Candidate*  myLepton[2];  
      for(unsigned int il=0; il < 2; il++)
	myLepton[il]= h.daughter(LEPZ)->daughter(il)->masterClone().get();

      const pat::Jet * myJet[2];
      for(unsigned int ijet=0; ijet < 2; ijet++)
	myJet[ijet] =
	  dynamic_cast<const pat::Jet *>(h.daughter(HADZ)->daughter(ijet)->
					 masterClone().get());
      
      // LOOK FOR TWO GOOD CHARGED LEPTONS
      int nLepPtHi=0;
      int nLepPtLo=0;

      for(unsigned int il=0; il < 2; il++){
	
	if(deltaR(myLepton[il]->eta(), myLepton[il]->phi(),
		  myJet[0]->eta(), myJet[0]->phi()) < MIN_DR_JETLEP)continue;

	if(deltaR(myLepton[il]->eta(), myLepton[il]->phi(),
		  myJet[1]->eta(), myJet[1]->phi()) < MIN_DR_JETLEP)continue;


	double pt = myLepton[il]->pt();	  
	if(pt >MIN_LEPPT1) nLepPtHi++;
	if(pt >MIN_LEPPT2) nLepPtLo++;

      }

      if(nLepPtHi < 1)continue;
      if(nLepPtLo < 2)continue;
      
      bool OppCharge = myLepton[0]->charge()*myLepton[1]->charge() < 0 ? 
	true: false;

      if(!OppCharge)continue;


      // they need to pass ID cuts
      int nPassID=0;

      if(ilep==0){ // electron
	nPassID=0;

	for(unsigned int iele=0; iele < 2; iele++){
	  
	  const pat::Electron* myEle
	    = dynamic_cast<const pat::Electron*>(h.daughter(LEPZ)->daughter(iele)->masterClone().get());

	  std::map<std::string, bool> Pass    = e2012ID_.CutRecord(*myEle); 
	  int passOrNot = PassAll(Pass);

	  // 	      std::map<std::string, bool>::iterator iterCut= Pass.begin();
	  // 	      for(;iterCut!=Pass.end();iterCut++)
	  // 		std::cout<< myEle->eta() << "-->"<<iterCut->first<<"\t"
	  // 			 <<iterCut->second<<std::endl;            


	  if(passOrNot==0)continue; // 2012 loose electron ID	  
	  nPassID++;
	}
      } // if it's an electron type

      else if(ilep==1){ // muon

	nPassID=0;
	for(unsigned int imuo=0; imuo < 2; imuo++){
	 
	  const pat::Muon* myMuo
	    = dynamic_cast<const pat::Muon*>(h.daughter(LEPZ)->daughter(imuo)->masterClone().get());
	  std::map<std::string, bool> Pass = mu2012ID_.CutRecord(*myMuo);
	  int passOrNot = PassAll(Pass);
	  if(passOrNot==0)continue; // 2012 tight muon ID	  
	  nPassID++;	  	  
	} // end of loop over muon
      } // if is a muon type
     
      
      if(nPassID < 2)continue;


      // LOOK FOR 2 JETS PASSING BETA CUTS

      int nGoodJets=0;
      int nLooseBTags=0;
      int nMediumBTags=0;

      for(unsigned int ijet=0; ijet < 2; ijet++){	 
	
	double pt  = myJet[ijet]->pt();
	if(pt < MIN_JETPT)continue;
	
	double eta = myJet[ijet]->eta();
	if(fabs(eta)> MAX_JETETA)continue;
	
	// to suppress jets from pileups
	double puBeta = myJet[ijet]->userFloat("puBeta");
	if(puBeta < MIN_JETBETA)continue;

	
	if(deltaR(myLepton[0]->eta(), myLepton[0]->phi(),
		  myJet[ijet]->eta(), myJet[ijet]->phi()) < MIN_DR_JETLEP)continue;

	if(deltaR(myLepton[1]->eta(), myLepton[1]->phi(),
		  myJet[ijet]->eta(), myJet[ijet]->phi()) < MIN_DR_JETLEP)continue;
	  
	
	if( !passLooseJetID(myJet[ijet]) )continue;


	bool isLoose  = false;
	bool isMedium = false;
	double jpBTag = myJet[ijet]->bDiscriminator("jetProbabilityBJetTags");
	int flavor    = myJet[ijet]->partonFlavour();
	
	if(jpBTag > MIN_BTAG_JP_LOOSE)isLoose=true;
	if(jpBTag > MIN_BTAG_JP_MEDIUM)isMedium=true;
	
 	if(!isData)
	  {
	    double phi = myJet[ijet]->phi();
	    double sin_phi = sin(phi*1000000);
	    // 	    int seed = abs(static_cast<int>(sin_phi*100000));	    
	    // 	    BTagSFUtil* btsfutiljp = new BTagSFUtil("JP", seed);
	
	    float Btageff_SF_ = 1.0;

	    for( size_t iMeasure = 0; iMeasure < measureName.size(); iMeasure++ ){
	      //Setup our measurement
	      iSetup.get<BTagPerformanceRecord>().get( measureName[ iMeasure ],perfH);
	      const BtagPerformance & perf = *(perfH.product());
	      BinningPointByMap measurePoint;
	      measurePoint.reset();
	      ///// pass in the et of the jet
	      double jetEt = myJet[ijet]->et();
	      measurePoint.insert(BinningVariables::JetEt, jetEt);  
	      ///// pass in the absolute eta of the jet
	      measurePoint.insert(BinningVariables::JetEta, fabs(eta));       
	      //this is the correction for 2012
 	      float SFL_JPL_2012corr = 1.01744  + (0.000813491*jetEt)-
  		(6.01592e-07)*jetEt*jetEt;
  	      float SFL_JPM_2012corr = 0.964487 + (0.00134038*jetEt)-
  		(1.43995e-06)*jetEt*jetEt;

	      // 2011
	      // 	      float SFL_JPL_2012corr = 1.00;
	      //  	      float SFL_JPM_2012corr = 1.00;

	      // Extract the mistag eff value
	      if ( measureType[ iMeasure ] == "BTAGLEFF") {
		// add a suffix eff so we can distingiush it from other values
		std::string suffix = "eff"; 
		ScaleFactorsEff[ijet][ measureName[ iMeasure ] + suffix ] = perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
	      }
	      else{ 
		// Extract the mistag and btag SF
		// The factor Btageff_SF_ is used for Btagging systematics 
		// and should be set to 1.0 as default
		if(measureName[ iMeasure ] == "MISTAGJPM")
		  ScaleFactors[ijet][ measureName[ iMeasure ] ] = 
		    SFL_JPM_2012corr*Btageff_SF_*perf.getResult( measureMap[ measureType[ iMeasure] ],measurePoint);

		else if(measureName[ iMeasure ] == "MISTAGJPL")
		  ScaleFactors[ijet][ measureName[ iMeasure ] ] = 
		    SFL_JPL_2012corr*Btageff_SF_*perf.getResult( measureMap[ measureType[ iMeasure] ],measurePoint);
		else
		  ScaleFactors[ijet][ measureName[ iMeasure ] ] = Btageff_SF_*perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
	    
	      }
	    }
	    
   	    btsfutiljp->modifyBTagsWithSF( 
   					  isLoose, isMedium, flavor, 
   					  ScaleFactors[ijet]["MUJETSWPBTAGJPL"],
   					  ScaleFactors[ijet]["MUJETSWPBTAGJPM"], 
   					  ScaleFactors[ijet]["MISTAGJPL"],
   					  ScaleFactors[ijet]["MISTAGJPM"],
 					  ScaleFactorsEff[ijet]["MISTAGJPLeff"],
   					  ScaleFactorsEff[ijet]["MISTAGJPMeff"]);
	    
	    // 	    delete btsfutiljp;
	  } // if it's MC
	
	if(isLoose)  nLooseBTags++;
	if(isMedium) nMediumBTags++;

	nGoodJets++; 
      } // end of loop over jets
      
      // there must be at least two jets
    //  if(nGoodJets < 2)continue;

      // number of btags
      int nBTags = 0;
      if(nMediumBTags >= 1 && nLooseBTags == 2)nBTags=2;
      else if(nLooseBTags >=1)nBTags=1;
      else nBTags=0;

      hcand = i;
	
      
      const pat::CompositeCandidate & goodH = (*hzzlljj)[hcand];
      const reco::Candidate * Zll = goodH.daughter(LEPZ);

      const reco::Candidate * Zjj = goodH.daughter(HADZ);
      //Final selection in categories

      double higgsPt_local  = goodH.pt();
      double higgsEta_local = goodH.eta();
      double higgsPhi_local = goodH.phi();
      double higgsM_local   = goodH.mass();
      double higgsM_refit_local = goodH.userFloat("HZZRefitMass");

      double zllPt_local  = Zll->pt();
      double zllEta_local = Zll->eta();
      double zllPhi_local = Zll->phi();
      double zllM_local   = Zll->mass();
      double zlldR_local  = deltaR(myLepton[0]->eta(),myLepton[0]->phi(),
				   myLepton[1]->eta(),myLepton[1]->phi());
      
      double zjjPt_local  = Zjj->pt();
      double zjjEta_local = Zjj->eta();
      double zjjPhi_local = Zjj->phi();
      double zjjM_local   = Zjj->mass();
      double zjjM_refit_local = goodH.userFloat("ZjjRefitMass");
      double zjjdR_local  = deltaR(myJet[0]->eta(),myJet[0]->phi(),
				   myJet[1]->eta(),myJet[1]->phi());


      
      double heliLD_local = goodH.userFloat("helyLD");
      double heliLD_refit_local = goodH.userFloat("helyLDRefit");

      // check the dilepton Z mass
      //       if(zllM_local < MIN_MZ_LL)continue;
      //       if(zllM_local > MAX_MZ_LL)continue;

      // check the mass of the two jets
      //       if(zjjM_local < LOOSE_MIN_MZ_JJ)continue;
      //       if(zjjM_local > LOOSE_MAX_MZ_JJ)continue;

      bestHCandIndex++;
    
      higgsPt_.push_back(higgsPt_local);
      higgsEta_.push_back(higgsEta_local);
      higgsPhi_.push_back(higgsPhi_local);
      higgsM_.push_back(higgsM_local);
      higgsMRefit_.push_back(higgsM_refit_local);

      zllPt_.push_back(zllPt_local);
      zllEta_.push_back(zllEta_local);
      zllPhi_.push_back(zllPhi_local);
      zllM_.push_back(zllM_local);
      zlldR_.push_back(zlldR_local);

      zjjPt_.push_back(zjjPt_local);
      zjjEta_.push_back(zjjEta_local);
      zjjPhi_.push_back(zjjPhi_local);
      zjjM_.push_back(zjjM_local);
      zjjMRefit_.push_back(zjjM_refit_local);
      zjjdR_.push_back(zjjdR_local);

      for(int isjet=0; isjet<2; isjet++)
	{

	  HjetIndex_.push_back(isjet);
	  HjetHiggsIndex_.push_back(bestHCandIndex);
	  HjetE_.push_back(myJet[isjet]->energy());
	  HjetPt_.push_back(myJet[isjet]->pt());
	  HjetEta_.push_back(myJet[isjet]->eta());
	  HjetPhi_.push_back(myJet[isjet]->phi());
	  
	  
	}

      for(int isLep=0; isLep<2; isLep++)
	{

	  LeptonsIndex_.push_back(isLep);
	  LeptonsHiggsIndex_.push_back(bestHCandIndex);
	  LeptonsE_.push_back(myLepton[isLep]->energy());
	  LeptonsPt_.push_back(myLepton[isLep]->pt());
	  LeptonsEta_.push_back(myLepton[isLep]->eta());
	  LeptonsPhi_.push_back(myLepton[isLep]->phi());
          LeptonsType_.push_back(ilep);	  
	  
	}


      heliLD_.push_back(heliLD_local);
      heliLDRefit_.push_back(heliLD_refit_local);

      costhetaNT1_.push_back(goodH.userFloat("costhetaNT1"));
      costhetaNT2_.push_back(goodH.userFloat("costhetaNT2"));
      phiNT_.push_back(goodH.userFloat("phiNT"));
      phiNT1_.push_back(goodH.userFloat("phiNT1"));
      costhetastarNT_.push_back(goodH.userFloat("costhetastarNT"));
      
      costhetaNT1Refit_.push_back(goodH.userFloat("costhetaNT1Refit"));
      costhetaNT2Refit_.push_back(goodH.userFloat("costhetaNT2Refit"));
      phiNTRefit_.push_back(goodH.userFloat("phiNTRefit"));
      phiNT1Refit_.push_back(goodH.userFloat("phiNT1Refit"));
      costhetastarNTRefit_.push_back(goodH.userFloat("costhetastarNTRefit"));

      nBTags_.push_back(nBTags);
      lepType_.push_back(ilep);

      // CHECK WHICH CUTS THIS HIGGS CANDIDATE PASSES
      int thisBit = 0;

      if(zjjM_local > MIN_MZ_JJ && zjjM_local < MAX_MZ_JJ)
	thisBit |= MZJJ_SIGNAL;
      else
	thisBit |= MZJJ_SIDEBAND;

      if(metSig_ < MAX_MET_SIG[nBTags])
	thisBit |= PFMET_SIG;

      double heliLD_cutValue = -9999;
      heliLD_cutValue = MIN_HELD_CONSTANT[nBTags] + 
	MIN_HELD_SLOPE[nBTags]* higgsM_refit_local;
      
      if(heliLD_refit_local > heliLD_cutValue)
	thisBit |= HELI_LD;

      double qgld= 99999.0; // needs to be fixed
      if(qgld > MIN_QUARK_GLUON_LD_0BTAG)
	thisBit |= QG_LD;
      
      // need to implement cuts on 4-body mass, need to be fixed
      double mH_input = higgsM_refit_local;
      if(higgsM_refit_local > MIN_MH_RATIO*mH_input &&
	 higgsM_refit_local < MAX_MH_RATIO*mH_input)
	thisBit |= MH_SIGNAL;


      // now check all cuts

      if( (thisBit & MZJJ_SIGNAL) &&
	  (thisBit & MH_SIGNAL)   &&
	  (thisBit & PFMET_SIG)   &&	      
	  (thisBit & HELI_LD)     
//	 && (thisBit & QG_LD)
        )
	thisBit |= ALL_SIGNAL;
      
      else if( (thisBit & MZJJ_SIDEBAND) &&
	       (thisBit & MH_SIGNAL)   &&
	       (thisBit & PFMET_SIG)   &&	      
	       (thisBit & HELI_LD)     &&
	       (thisBit & QG_LD))
	thisBit |= ALL_SIDEBAND;
      
      passBit_.push_back(thisBit);



    } // end of loop over Higgs candidates  
    
  } // end of looping over lepton types

} // end of Fill()

  //-----------------------------------------------------------------------
void  
ZZTree::SetBranches(){

  AddBranch(&eleRho_, "eleRho");
  AddBranch(&muoRho_, "muoRho");
  AddBranch(&metSig_, "metSig");

  AddBranch(&higgsPt_,"higgsPt");
  AddBranch(&higgsEta_,"higgsEta");
  AddBranch(&higgsPhi_,"higgsPhi");
  AddBranch(&higgsM_,"higgsM");
  AddBranch(&higgsMRefit_,"higgsMRefit");

  AddBranch(&zllPt_,"zllPt");
  AddBranch(&zllEta_,"zllEta");
  AddBranch(&zllPhi_,"zllPhi");
  AddBranch(&zllM_,"zllM");
  AddBranch(&zlldR_,"zlldR");

  AddBranch(&zjjPt_,"zjjPt");
  AddBranch(&zjjEta_,"zjjEta");
  AddBranch(&zjjPhi_,"zjjPhi");
  AddBranch(&zjjM_,"zjjM");
  AddBranch(&zjjMRefit_,"zjjMRefit");
  AddBranch(&zjjdR_,"zjjdR");

  AddBranch(&HjetIndex_,"HjetIndex");
  AddBranch(&HjetHiggsIndex_,"HjetHiggsIndex_");
  AddBranch(&HjetE_,"HjetE");
  AddBranch(&HjetPt_,"HjetPt");
  AddBranch(&HjetEta_,"HjetEta");
  AddBranch(&HjetPhi_,"HjetPhi");
  
  AddBranch(&LeptonsIndex_,"LeptonsIndex");
  AddBranch(&LeptonsHiggsIndex_,"LeptonsHiggsIndex_");
  AddBranch(&LeptonsE_,"LeptonsE");
  AddBranch(&LeptonsPt_,"LeptonsPt");
  AddBranch(&LeptonsEta_,"LeptonsEta");
  AddBranch(&LeptonsPhi_,"LeptonsPhi");
 AddBranch(&LeptonsType_,"LeptonsType"); 
 
  AddBranch(&heliLD_,"heliLD");
  AddBranch(&costhetaNT1_,"costhetaNT1");
  AddBranch(&costhetaNT2_,"costhetaNT2");
  AddBranch(&phiNT_,"phiNT");
  AddBranch(&phiNT1_,"phiNT1");
  AddBranch(&costhetastarNT_,"costhetastarNT");
  
  AddBranch(&heliLDRefit_,"heliLDRefit");
  AddBranch(&costhetaNT1Refit_,"costhetaNT1Refit");
  AddBranch(&costhetaNT2Refit_,"costhetaNT2Refit");
  AddBranch(&phiNTRefit_,"phiNTRefit");
  AddBranch(&phiNT1Refit_,"phiNT1Refit");
  AddBranch(&costhetastarNTRefit_,"costhetastarNTRefit");

  AddBranch(&nBTags_,"nBTags");
  AddBranch(&lepType_,"lepType");
  AddBranch(&passBit_,"passBit");


}


void  
ZZTree::Clear(){

  eleRho_ = DUMMY;
  muoRho_ = DUMMY;

  metSig_ = DUMMY;

  higgsPt_.clear();
  higgsEta_.clear();
  higgsPhi_.clear();
  higgsM_.clear();
  higgsMRefit_.clear();

  zllPt_.clear();
  zllEta_.clear();
  zllPhi_.clear();
  zllM_.clear();
  zlldR_.clear();

  zjjPt_.clear();
  zjjEta_.clear();
  zjjPhi_.clear();
  zjjM_.clear();
  zjjMRefit_.clear();
  zjjdR_.clear();

  HjetIndex_.clear();
  HjetHiggsIndex_.clear();
  HjetE_.clear();
  HjetPt_.clear();
  HjetEta_.clear();
  HjetPhi_.clear();
  
  LeptonsIndex_.clear();
  LeptonsHiggsIndex_.clear();
  LeptonsE_.clear();
  LeptonsPt_.clear();
  LeptonsEta_.clear();
  LeptonsPhi_.clear();
  LeptonsType_.clear();
  


  heliLD_.clear();
  costhetaNT1_.clear();
  costhetaNT2_.clear();
  phiNT_.clear();
  phiNT1_.clear();
  costhetastarNT_.clear();
 
 
  heliLDRefit_.clear();
  costhetaNT1Refit_.clear();
  costhetaNT2Refit_.clear();
  phiNTRefit_.clear();
  phiNT1Refit_.clear();
  costhetastarNTRefit_.clear();

  nBTags_.clear();
  lepType_.clear();
  passBit_.clear();


}


bool ZZTree::passLooseJetID(const pat::Jet* recjet)
{
  double eta = recjet->eta();
  if(recjet->getPFConstituents().size() <= 1)return false;                                                                               
  if(recjet->neutralHadronEnergyFraction() >= 0.99)return false;
  if(recjet->neutralEmEnergyFraction() >= 0.99)return false;
  //   // for the tracker region
  if(fabs(eta)<2.4 && recjet->chargedHadronEnergyFraction()<= 0.0)return false;
  if(fabs(eta)<2.4 && recjet->chargedEmEnergyFraction() >= 0.99)return false;
  if(fabs(eta)<2.4 && recjet->chargedMultiplicity() <= 0)return false;
  return true;

}

/*
  eID01.push_back(myEle->pt());
  eID02.push_back(myEle->superCluster()->eta());
  eID03.push_back(myEle->deltaEtaSuperClusterTrackAtVtx());
  eID04.push_back(myEle->deltaPhiSuperClusterTrackAtVtx());
  eID05.push_back(myEle->sigmaIetaIeta());
  eID06.push_back(myEle->hadronicOverEm());
  eID07.push_back(myEle->userFloat("dxy"));
  eID08.push_back(myEle->userFloat("dz"));
  eID09.push_back(fabs(1.0/myEle->ecalEnergy() - 
  myEle->eSuperClusterOverP()/myEle->ecalEnergy()));
  eID10.push_back(myEle->convDcot());
  eID11.push_back(myEle->convDist());
  eID12.push_back(myEle->userFloat("hasMatchConv"));
  eID13.push_back(myEle->gsfTrack().get()->trackerExpectedHitsInner().numberOfHits());
	  
  double iso1 = myEle->chargedHadronIso();
  double iso2 = myEle->neutralHadronIso();
  double iso3 = myEle->photonIso();
	  
  ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = 
  isData? ElectronEffectiveArea::kEleEAData2011:
  ElectronEffectiveArea::kEleEAFall11MC;
      
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaType_ =
  ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03;      
 
  double AEff = ElectronEffectiveArea::GetElectronEffectiveArea
  (effAreaType_, fabs(myEle->superCluster()->eta()), effAreaTarget_);

  double iso4 = iso1 + std::max(0.0, iso2+iso3-ele_rho*AEff);
	  
  eID14.push_back(iso1);
  eID15.push_back(iso2);
  eID16.push_back(iso3);
  eID17.push_back(iso4);
*/


/*
  double pt    =  myMuo->pt();
  double eta   =  myMuo->eta();
  double iso1  =  myMuo->pfIsolationR04().sumChargedHadronPt;
  double iso2  =  myMuo->pfIsolationR04().sumNeutralHadronEt;
  double iso3  =  myMuo->pfIsolationR04().sumPhotonEt;
  double isoPU =  myMuo->pfIsolationR04().sumPUPt;    
  double iso4Beta = iso1 + std::max(iso3+iso2-0.5*isoPU,0.);
  MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = 
  isData? MuonEffectiveArea::kMuEAData2012:
  MuonEffectiveArea::kMuEAData2012;
  // 	MuonEffectiveArea::kMuEAFall11MC;
  MuonEffectiveArea::MuonEffectiveAreaType effAreaType_= 
  MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;
  double Area = MuonEffectiveArea::GetMuonEffectiveArea(
  effAreaType_, fabs(eta), effAreaTarget_);

  double iso4Rho =  iso1 + std::max(iso3+iso2-muo_rho*Area,0.);
  double isoBeta = iso4Beta/pt;	  
  double isoRho  = iso4Rho/pt;


  muID01.push_back(isoBeta);
  // 	  muID01.push_back(myMuo->isGlobalMuon());
  muID02.push_back(myMuo->isPFMuon());
  muID03.push_back(myMuo->isTrackerMuon());
  if(myMuo->isTrackerMuon() && myMuo->isGlobalMuon()){
  muID04.push_back(myMuo->globalTrack()->normalizedChi2());
  muID05.push_back(myMuo->globalTrack()->hitPattern().numberOfValidMuonHits());
  muID06.push_back(myMuo->numberOfMatchedStations());
  muID07.push_back(myMuo->dB());
  double dzV =myMuo->userFloat("dzVtx") ;
  muID08.push_back(dzV);
  muID09.push_back(myMuo->innerTrack()->hitPattern().numberOfValidPixelHits());
  muID10.push_back(myMuo->innerTrack()->hitPattern().trackerLayersWithMeasurement());
  }
  muID11.push_back(myMuo->pfIsolationR04().sumChargedHadronPt);
  muID12.push_back(myMuo->pfIsolationR04().sumNeutralHadronEt);
  muID13.push_back(myMuo->pfIsolationR04().sumPhotonEt);
  muID14.push_back(myMuo->pfIsolationR04().sumPUPt);
  muID15.push_back(isoRho);
*/

/*

higgsPt_  = DUMMY;
higgsEta_ = DUMMY;
higgsPhi_ = DUMMY;
higgsM_   = DUMMY;

zllPt_  = DUMMY;
zllEta_ = DUMMY;
zllPhi_ = DUMMY;
zllM_   = DUMMY;

zjjPt_  = DUMMY;
zjjEta_ = DUMMY;
zjjPhi_ = DUMMY;
zjjM_   = DUMMY;


int arraySize = sizeof(jetPt_)/sizeof(jetPt_[0]);

for(int i=0; i<arraySize;i++)
{
jetPt_[i] =DUMMY;
jetEta_[i]=DUMMY;
jetPhi_[i]=DUMMY;
jetE_[i]  =DUMMY;
jetRefitPt_[i] =DUMMY;
jetRefitEta_[i]=DUMMY;
jetRefitPhi_[i]=DUMMY;
jetRefitE_[i]=DUMMY;

}
lepType_ = -1;
  
arraySize = sizeof(lepPt_)/sizeof(lepPt_[0]);

for(int i=0; i<arraySize;i++)
{
lepPt_[i] =DUMMY;
lepEta_[i]=DUMMY;
lepPhi_[i]=DUMMY;
lepE_[i]  =DUMMY;

}

AddBranch(&eID01, "eID01");
AddBranch(&eID02, "eID02");
AddBranch(&eID03, "eID03");
AddBranch(&eID04, "eID04");
AddBranch(&eID05, "eID05");
AddBranch(&eID06, "eID06");
AddBranch(&eID07, "eID07");
AddBranch(&eID08, "eID08");
AddBranch(&eID09, "eID09");
AddBranch(&eID10, "eID10");
AddBranch(&eID11, "eID11");
AddBranch(&eID12, "eID12");
AddBranch(&eID13, "eID13");
AddBranch(&eID14, "eID14");
AddBranch(&eID15, "eID15");
AddBranch(&eID16, "eID16");
AddBranch(&eID17, "eID17");

eID01.clear();
eID02.clear();
eID03.clear();
eID04.clear();
eID05.clear();
eID06.clear();
eID07.clear();
eID08.clear();
eID09.clear();
eID10.clear();
eID11.clear();
eID12.clear();
eID13.clear();
eID14.clear();
eID15.clear();
eID16.clear();
eID17.clear();

AddBranch(&muID01, "muID01");
AddBranch(&muID02, "muID02");
AddBranch(&muID03, "muID03");
AddBranch(&muID04, "muID04");
AddBranch(&muID05, "muID05");
AddBranch(&muID06, "muID06");
AddBranch(&muID07, "muID07");
AddBranch(&muID08, "muID08");
AddBranch(&muID09, "muID09");
AddBranch(&muID10, "muID10");
AddBranch(&muID11, "muID11");
AddBranch(&muID12, "muID12");
AddBranch(&muID13, "muID13");
AddBranch(&muID14, "muID14");
AddBranch(&muID15, "muID15");

muID01.clear();
muID02.clear();
muID03.clear();
muID04.clear();
muID05.clear();
muID06.clear();
muID07.clear();
muID08.clear();
muID09.clear();
muID10.clear();
muID11.clear();
muID12.clear();
muID13.clear();
muID14.clear();
muID15.clear();


*/
