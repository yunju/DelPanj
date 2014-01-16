// -*- C++ -*-
//
// Package:    ZZ2l2qFilter
// Class:      ZZ2l2qFilter
// 
/**\class ZZ2l2qFilter ZZ2l2qFilter.cc testFilter/ZZ2l2qFilter/src/ZZ2l2qFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yun-Ju Lu,32 3-C06,+41227679377,
//         Created:  Tue Aug 27 11:57:24 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//for trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include <vector>
 using namespace edm;
 using namespace std;
 using namespace reco;



//
// class declaration
//

class ZZ2l2qFilter : public edm::EDFilter {
   public:
      explicit ZZ2l2qFilter(const edm::ParameterSet&);
      ~ZZ2l2qFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      int Npcrosess;
      int NpassTrig;
      int NpassMetS;
      int NfailTrig;
      int NfailMetS;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZZ2l2qFilter::ZZ2l2qFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  Npcrosess =0;
  NpassTrig =0;
  NpassMetS =0;
  NfailTrig =0;
  NfailMetS =0; 

}


ZZ2l2qFilter::~ZZ2l2qFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZZ2l2qFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

    Npcrosess++;

  
int    HLTDoubleMu_=-999;
int   HLTDoubleEle_=-999;
int   HLTMu17TkMu8_=-999;
//   if(isData)
//   {


     edm::Handle<edm::TriggerResults> hltresults;
     iEvent.getByLabel(edm::InputTag("TriggerResults", "", "HLT"),hltresults);
     edm::TriggerNames TrigNames = iEvent.triggerNames(*hltresults);

     for ( size_t itr = 0; itr < hltresults->size(); ++itr ) {
     std::string passedPathName = TrigNames.triggerName(itr);

     //if (passedPathName.find("HLT_Mu17_Mu8_v")!=std::string::npos )std::cout<<iEvent.id().event() <<" Trigger names :"<< passedPathName<<" "<<hltresults->accept(itr)<<std::endl;
    // if ((passedPathName.find("HLT_Mu17_Mu8_v")!=std::string::npos || passedPathName.find("HLT_Mu17_TkMu8")!=std::string::npos) && hltresults->accept(itr) )
      if (passedPathName.find("HLT_Mu17_Mu8_v")!=std::string::npos && hltresults->accept(itr) )
      {
     // std::cout<<"fire Mu"<<std::endl;
      HLTDoubleMu_=1;
     }
    else if ((passedPathName.find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")!=std::string::npos && hltresults->accept(itr)))
     {
        HLTDoubleEle_=1;
       //  std::cout<<"fire Elec"<<std::endl;
      //    cin.get();
     }
     else if ((passedPathName.find("HLT_Mu17_TkMu8")!=std::string::npos && hltresults->accept(itr)))
     {
     HLTMu17TkMu8_=1;
       //  std::cout<<"fire Elec"<<std::endl;
      //    cin.get();
     }

  }//Loop trigger bit  

  //Met Cut 

//MET

 //edm::Handle<float>  metSignificanceHandle;
  Handle<pat::METCollection> metHandle;
  iEvent.getByLabel("patMETsPFJetsAK5", metHandle);
  
  pat::METCollection met_h = *metHandle;
  
   auto_ptr<float> metSignificance_( new float );
   *metSignificance_ = -999.;
  TMatrixD metmat = met_h.front().getSignificanceMatrix();

  if( (metmat < 1.0e10) && (metmat > -1.0e10) ) { 
  *metSignificance_ = met_h.front().significance();
   }

   if(HLTDoubleMu_==1||HLTDoubleEle_==1||HLTMu17TkMu8_==1)
   { 
     NpassTrig++;      
     if((met_h.front().significance()<=10))
     {
       NpassMetS++; 
       return true;
     }
     else
     {
       NfailMetS++; 
       return false;
     } 
    
   }
   else {
   NfailTrig++;
   return false;}
   
    
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZZ2l2qFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZZ2l2qFilter::endJob() {
std::cout<<"##ZZ2l2qFilter(Trig&Met) " <<"NpassTrig: "<< NpassTrig <<" NpassTrig and MetS 10: "<<NpassMetS<<" Npcrosess: "<<Npcrosess<<std::endl;
std::cout<<"##ZZ2l2qFilter(Trig&Met) " <<"NfailTrig: "<< NfailTrig <<" Nfail MetS 10: "<<NfailMetS<<" Npcrosess: "<<Npcrosess<<std::endl;

}

// ------------ method called when starting to processes a run  ------------
bool 
ZZ2l2qFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ZZ2l2qFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ZZ2l2qFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ZZ2l2qFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZZ2l2qFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZZ2l2qFilter);
