// -*- C++ -*-
//
// Package:    EventCounter
// Class:      EventCounter
// 
/**\class EventCounter EventCounter.cc DelPanj/EventCounter/src/EventCounter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Anil Pratap Singh,,,
//         Created:  Mon Sep 19 05:31:51 CEST 2011
// $Id: EventCounter.cc,v 1.2 2013/03/18 13:49:23 yunju Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// for SHERPA weighting
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//
// class declaration
//

class YJEventCounter : public edm::EDAnalyzer {
   public:
      explicit YJEventCounter(const edm::ParameterSet&);
      ~YJEventCounter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      int instance_;
      TH1D* hCount_;
      TH1D* hWeight_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
//int YJEventCounter::instance_ = 1;
//
// constructors and destructor
//
YJEventCounter::YJEventCounter(const edm::ParameterSet& iConfig)
:instance_(iConfig.getParameter<int>("instance"))
{
//now do what ever initialization is needed
  std::stringstream ss;
  ss<<instance_;
  std::string name = "hCounter_"+ss.str();
  edm::Service<TFileService> fs;
  hCount_=fs->make<TH1D>(name.c_str(), name.c_str(),10.,0.,10.);

  name = "hWeight_"+ss.str();
  hWeight_=fs->make<TH1D>(name.c_str(), name.c_str(),10,0,10.0);


}


YJEventCounter::~YJEventCounter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
YJEventCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   hCount_->Fill(instance_);

   edm::Handle<GenEventInfoProduct>    genEventScale;

   double w = 1.0;
   if(!iEvent.isRealData()){
     if (iEvent.getByLabel("generator", genEventScale)) 
       w = genEventScale->weight();
   }
   hWeight_->Fill(instance_,w);

}


// ------------ method called once each job just before starting event loop  ------------
void 
YJEventCounter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
YJEventCounter::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
YJEventCounter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
YJEventCounter::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
YJEventCounter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
YJEventCounter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
YJEventCounter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(YJEventCounter);
