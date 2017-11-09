//ACESS DATA WITH TRIGGER
// -*- C++ -*-
//
// Package:    JpsiAnalyzerOpen2011
// Class:      JpsiAnalyzerOpen2011
// 
/**\class JpsiAnalyzerOpen2011 JpsiAnalyzerOpen2011.cc CMSOpenData2011/JpsiAnalyzerOpen2011/src/JpsiAnalyzerOpen2011.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sandro Fonseca De Souza,32 4-C14,+41227674949,
//         Created:  Thu Apr 13 21:58:44 CEST 2017
// $Id$
//
//


// system include files
#include <memory>
#include <cmath>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Muon Reco
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
/*
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
*/
//includes ROOT
#include <TROOT.h>
#include <TTree.h>

//
// class declaration
//

class JpsiAnalyzerOpen2011 : public edm::EDAnalyzer {
   public:
      explicit JpsiAnalyzerOpen2011(const edm::ParameterSet&);
      ~JpsiAnalyzerOpen2011();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      bool triggerfired(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname);
      bool triggerfound(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname);
 
     unsigned int triggerIndex(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname);
 
  	// ----------member data ---------------------------

	TTree* AnalysisTree;
	//TBranch* branch;

    //Histogrms
   	TH1F* leadingMuon_Pt ;
  	TH1F* leadingMuon_Eta ;
    TH1F* leadingMuon_Phi ;
    TH1F* leadingMuon_Charge ;
	TH1F* leadingMuon_Mass ;

    TH1F* trailingMuon_Pt ;
    TH1F* trailingMuon_Eta ;
    TH1F* trailingMuon_Phi ;
    TH1F* trailingMuon_Charge ;
	TH1F* trailingMuon_Mass ;

	TH1F* HistoMuon_Pt;
	TH1F* HistoMuon_Eta;
	TH1F* HistoMuon_Phi;
	TH1F* HistoMuon_Charge;
	TH1F* HistoMuon_Mass;

	TH1F* HistoMuonTight_Pt;
	TH1F* HistoMuonTight_Eta;
	TH1F* HistoMuonTight_Phi;
	TH1F* HistoMuonTight_Charge;
	TH1F* HistoMuonTight_Mass;

	TH1F* HistoMuonTightValidHits_Pt;
	TH1F* HistoMuonTightValidHits_Eta;
	TH1F* HistoMuonTightValidHits_Phi;
	TH1F* HistoMuonTightValidHits_Charge;
	TH1F* HistoMuonTightValidHits_Mass;

	TH1F* HistoMuonTightValidHitsPixelLayer_Pt;
	TH1F* HistoMuonTightValidHitsPixelLayer_Eta;
	TH1F* HistoMuonTightValidHitsPixelLayer_Phi;
	TH1F* HistoMuonTightValidHitsPixelLayer_Charge;
	TH1F* HistoMuonTightValidHitsPixelLayer_Mass;

	TH1F* HistoMuonTightValidHitsPixelLayerChi2_Pt;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2_Eta;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2_Phi;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2_Charge;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2_Mass;

	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDz_Pt;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDz_Eta;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDz_Phi;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDz_Charge;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDz_Mass;

	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Pt;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Eta;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Phi;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Charge;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Mass;

	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge;
	TH1F* HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass;

	TH1F* histo_Mll;
	TH1F* histo_MllPt;
	TH1F* histo_MllEta;
	TH1F* histo_MllPhi;

    bool verbose_;
    bool triggerflag_;	
    edm::InputTag primaryVertexProducer_;
    edm::InputTag recoMuons_;
// Reco configs
  double minMuPt_;
  double maxMuEta_;
  double muonLeadPt_, muonTrailPt_;
  double minJPsiMass_ ;
// Trigger
  edm::InputTag hlTriggerEvent_;      // Input tag for TriggerEvent
  std::string triggerName_;
  edm::InputTag hlTriggerResults_;    // Input tag for TriggerResults
  edm::InputTag collectionName_;
    bool debug;
//Counters
    int irun, ievt,lumiBlock;
	int Total_Events = 0;
	int CounterMuons = 0;
	int TrackerMuon	= 0;
	int GlobalMuon = 0;
	int TMOneStationTight = 0;
	int	NumberOfValidMuonHits = 0;
	int pixelLayersWithMeasurement = 0;
	int	normalizedChi2 = 0;
	int db_dz = 0;
	int PFMuon = 0; 
	int TrackerGlobalMuon = 0;
	int nDimuon = 0;

//Trigger
    int countInAccepted = 0;
    int countInTriggered = 0;
//Creating Vectors

	std::vector<int> VectorEvent;
	std::vector<int> VectorRun;
	std::vector<int> VectorlumiBlock;

	std::vector<double> VectorMuon_Pt;
	std::vector<double>	VectorMuon_Eta;
	std::vector<double>	VectorMuon_Phi;
	std::vector<int>	VectorMuon_Charge;
	std::vector<double>	VectorMuon_Mass;

	std::vector<double> TrackerMuonPt;
	std::vector<double> TrackerMuonEta;
	std::vector<double> TrackerMuonPhi;
	std::vector<int> TrackerMuonCharge;

	std::vector<double> GlobalMuonPt;
	std::vector<double> GlobalMuonEta;
	std::vector<double> GlobalMuonPhi;
	std::vector<int> GlobalMuonCharge;

	std::vector<double> VectorMuonTight_Pt;
	std::vector<double> VectorMuonTight_Eta;
	std::vector<double> VectorMuonTight_Phi;
	std::vector<int> VectorMuonTight_Charge;
	std::vector<double> VectorMuonTight_Mass;

	std::vector<double> VectorMuonTightValidHits_Pt;
	std::vector<double> VectorMuonTightValidHits_Eta;
	std::vector<double> VectorMuonTightValidHits_Phi;
	std::vector<int> VectorMuonTightValidHits_Charge;
	std::vector<double> VectorMuonTightValidHits_Mass;

	std::vector<double> VectorMuonTightValidHitsPixelLayer_Pt;
	std::vector<double> VectorMuonTightValidHitsPixelLayer_Eta;
	std::vector<double> VectorMuonTightValidHitsPixelLayer_Phi;
	std::vector<int> VectorMuonTightValidHitsPixelLayer_Charge;
	std::vector<double> VectorMuonTightValidHitsPixelLayer_Mass;

	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2_Pt;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2_Eta;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2_Phi;
	std::vector<int> VectorMuonTightValidHitsPixelLayerChi2_Charge;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2_Mass;

	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi;
	std::vector<int> VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass;

	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi;
	std::vector<int> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass;

	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi;
	std::vector<int> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge;
	std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass;
		
	std::vector<double> VectorleadingMuon_Pt;
	std::vector<double> VectorleadingMuon_Eta;
	std::vector<double> VectorleadingMuon_Phi;
	std::vector<int> VectorleadingMuon_Charge;
	std::vector<double> VectorleadingMuon_Mass;

	std::vector<double> VectortrailingMuon_Pt;
	std::vector<double> VectortrailingMuon_Eta;
	std::vector<double> VectortrailingMuon_Phi;
	std::vector<int> VectortrailingMuon_Charge;
	std::vector<double> VectortrailingMuon_Mass;

	std::vector<double> VectorMll;
	std::vector<double> VectorMllpT;
	std::vector<double> VectorMlleta;
	std::vector<double> VectorMllphi;

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
JpsiAnalyzerOpen2011::JpsiAnalyzerOpen2011(const edm::ParameterSet& iConfig):
verbose_ (iConfig.getParameter< bool > ("verbose")),
triggerflag_ (iConfig.getParameter< bool > ("triggerflag")),
primaryVertexProducer_ (iConfig.getParameter<edm::InputTag>("primaryVertexProducer")),
recoMuons_(iConfig.getParameter<edm::InputTag>("recoMuonsLabel")),
// Reco config with the trigger
minMuPt_ (iConfig.getParameter<double>("minMuPt")),
maxMuEta_ (iConfig.getParameter<double>("maxMuEta")),
muonLeadPt_ (iConfig.getParameter<double>("minMuonLeadPt")),
muonTrailPt_ (iConfig.getParameter<double>("minMuonTrailPt")),
minJPsiMass_ (iConfig.getParameter<double>("minJPsiMass")),
hlTriggerEvent_ (iConfig.getUntrackedParameter<edm::InputTag>("TriggerEventTag", edm::InputTag("hltTriggerSummaryAOD", "", "HLT"))),
triggerName_ (iConfig.getUntrackedParameter<std::string>("PathName","HLT_Dimuon0_Jpsi_v")),
hlTriggerResults_  (iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultsTag", edm::InputTag("TriggerResults", "", "HLT")))
{
// Histos File
edm::Service<TFileService> fs;

// Define Histos
TH1D::SetDefaultSumw2();

	HistoMuon_Pt = fs->make<TH1F>("HistoMuon_Pt"  , "Muon_Pt"  ,   100,   0., 200.);
	HistoMuon_Eta = fs->make<TH1F>("HistoMuon_Eta" , "Muon_Eta" ,   100,  -5.,   5.);
	HistoMuon_Phi = fs->make<TH1F>("HistoMuon_Phi" , "Muon_Phi" ,   100,  -4.,   4.);
	HistoMuon_Charge = fs->make<TH1F>("HistoMuon_Charge" , "Muon_Charge", 100 ,-2., 2.);
	HistoMuon_Mass = fs->make<TH1F>("HistoMuon_Mass" , "Muon_Mass", 100 ,0., 1.);

	HistoMuonTight_Pt = fs->make<TH1F>("HistoMuonTight_Pt"  , "Muon_Pt"  ,   100,   0., 200.);
	HistoMuonTight_Eta = fs->make<TH1F>("HistoMuonTight_Eta" , "Muon_Eta" ,   100,  -5.,   5.);
	HistoMuonTight_Phi = fs->make<TH1F>("HistoMuonTight_Phi" , "Muon_Phi" ,   100,  -4.,   4.);
	HistoMuonTight_Charge = fs->make<TH1F>("HistoMuonTight_Charge" , "Muon_Charge", 100 ,-2., 2.);
	HistoMuonTight_Mass = fs->make<TH1F>("HistoMuonTight_Mass" , "Muon_Mass", 100 ,0., 1.);
	
	HistoMuonTightValidHits_Pt = fs->make<TH1F>("HistoMuonTightValidHits_Pt"  , "Muon_Pt"  ,   100,   0., 200.);
	HistoMuonTightValidHits_Eta = fs->make<TH1F>("HistoMuonTightValidHits_Eta" , "Muon_Eta" ,   100,  -5.,   5.);
	HistoMuonTightValidHits_Phi = fs->make<TH1F>("HistoMuonTightValidHits_Phi" , "Muon_Phi" ,   100,  -4.,   4.);
	HistoMuonTightValidHits_Charge = fs->make<TH1F>("HistoMuonTightValidHits_Charge" , "Muon_Charge", 100 ,-2., 2.);
	HistoMuonTightValidHits_Mass = fs->make<TH1F>("HistoMuonTightValidHits_Mass" , "Muon_Mass", 100 ,0., 1.);

	HistoMuonTightValidHitsPixelLayer_Pt = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayer_Pt"  , "Muon_Pt"  ,   100,   0., 200.);
	HistoMuonTightValidHitsPixelLayer_Eta = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayer_Eta" , "Muon_Eta" ,   100,  -5.,   5.);
	HistoMuonTightValidHitsPixelLayer_Phi = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayer_Phi" , "Muon_Phi" ,   100,  -4.,   4.);
	HistoMuonTightValidHitsPixelLayer_Charge = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayer_Charge" , "Muon_Charge", 100 ,-2., 2.);
	HistoMuonTightValidHitsPixelLayer_Mass = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayer_Mass" , "Muon_Mass", 100 ,0., 1.);

	HistoMuonTightValidHitsPixelLayerChi2_Pt = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2_Pt"  , "Muon_Pt"  ,   100,   0., 200.);
	HistoMuonTightValidHitsPixelLayerChi2_Eta = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2_Eta" , "Muon_Eta" ,   100,  -5.,   5.);
	HistoMuonTightValidHitsPixelLayerChi2_Phi = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2_Phi" , "Muon_Phi" ,   100,  -4.,   4.);
	HistoMuonTightValidHitsPixelLayerChi2_Charge = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2_Charge" , "Muon_Charge", 100 ,-2., 2.);
	HistoMuonTightValidHitsPixelLayerChi2_Mass = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2_Mass" , "Muon_Mass", 100 ,0., 1.);

	HistoMuonTightValidHitsPixelLayerChi2DbDz_Pt = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDz_Pt"  , "Muon_Pt"  ,   100,   0., 200.);
	HistoMuonTightValidHitsPixelLayerChi2DbDz_Eta = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDz_Eta" , "Muon_Eta" ,   100,  -5.,   5.);
	HistoMuonTightValidHitsPixelLayerChi2DbDz_Phi = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDz_Phi" , "Muon_Phi" ,   100,  -4.,   4.);
	HistoMuonTightValidHitsPixelLayerChi2DbDz_Charge = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDz_Charge" , "Muon_Charge", 100 ,-2., 2.);
	HistoMuonTightValidHitsPixelLayerChi2DbDz_Mass = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDz_Mass" , "Muon_Mass", 100 ,0., 1.);

	HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Pt = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Pt"  , "Muon_Pt"  ,   100,   0., 200.);
	HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Eta = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Eta" , "Muon_Eta" ,   100,  -5.,   5.);
	HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Phi = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Phi" , "Muon_Phi" ,   100,  -4.,   4.);
	HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Charge = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Charge" , "Muon_Charge", 100 ,-2., 2.);
	HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Mass = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Mass" , "Muon_Mass", 100 ,0., 1.);

	HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt"  , "Muon_Pt"  ,   100,   0., 200.);
	HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta" , "Muon_Eta" ,   100,  -5.,   5.);
	HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi" , "Muon_Phi" ,   100,  -4.,   4.);
	HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge" , "Muon_Charge", 100 ,-2., 2.);
	HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass = fs->make<TH1F>("HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass" , "Muon_Mass", 100 ,0., 1.);

	leadingMuon_Pt = fs->make<TH1F>("leadingMuon_Pt"  , "mu1_pt"  ,   100,   0., 200.);
	leadingMuon_Eta = fs->make<TH1F>("leadingMuon_Eta" , "mu1_eta" ,   100,  -5.,   5.);
	leadingMuon_Phi = fs->make<TH1F>("leadingMuon_Phi" , "mu1_phi" ,   100,  -4.,   4.);
	leadingMuon_Charge = fs->make<TH1F>("leadingMuon_Charge" , "mu1_charge", 100 ,-2., 2.);
	leadingMuon_Mass = fs->make<TH1F>("leadingMuon_Phi" , "mu1_mass" ,   100,  0.,   1.);
	
	trailingMuon_Pt = fs->make<TH1F>("trailingMuon_Pt"  , "mu2_pt"  ,   100,   0., 200.);
	trailingMuon_Eta = fs->make<TH1F>("trailingMuon_Eta" , "mu2_eta" ,   100,  -5.,   5.);
	trailingMuon_Phi = fs->make<TH1F>("trailingMuon_Phi" , "mu2_phi" ,   100,  -4.,   4.);
	trailingMuon_Charge = fs->make<TH1F>("trailingMuon_Charge" , "mu2_charge", 100 ,-2., 2.);
	trailingMuon_Mass = fs->make<TH1F>("trailingMuon_Phi" , "mu2_mass" ,   100,  -0,   1.);

	histo_Mll = fs->make<TH1F>("histo_Mll"  , "histo_Mll"  ,   100,   0., 100.);
	histo_MllPt = fs->make<TH1F>("histo_MllPt"  , "histo_MllPt"  ,   100,   0., 200.);
	histo_MllEta = fs->make<TH1F>("leadingMuonEta" , "Muon_Eta" ,   100,  -5.,   5.);
	histo_MllPhi = fs->make<TH1F>("histo_MllPhi" , "histo_MllPhi" ,   100,  -4.,   10.);

	//Define Trees
	AnalysisTree = fs->make<TTree>("AnalysisTree","Muon Analysis Tree");
	
}


JpsiAnalyzerOpen2011::~JpsiAnalyzerOpen2011()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
	std::vector<double>().swap(VectorMuon_Pt);
	std::vector<double>().swap(VectorMuon_Eta);
	std::vector<double>().swap(VectorMuon_Phi);
	std::vector<int>().swap(VectorMuon_Charge);
	std::vector<double>().swap(VectorMuon_Mass);

	std::vector<double>().swap(TrackerMuonPt);
	std::vector<double>().swap(TrackerMuonEta);
	std::vector<double>().swap(TrackerMuonPhi);
	std::vector<int>().swap(TrackerMuonCharge);

	std::vector<double>().swap(GlobalMuonPt);
	std::vector<double>().swap(GlobalMuonEta);
	std::vector<double>().swap(GlobalMuonPhi);
	std::vector<int>().swap(GlobalMuonCharge);
		
	std::vector<double>().swap(VectorMuonTight_Pt);
	std::vector<double>().swap(VectorMuonTight_Eta);
	std::vector<double>().swap(VectorMuonTight_Phi);
	std::vector<int>().swap(VectorMuonTight_Charge);
	std::vector<double>().swap(VectorMuonTight_Mass);

	std::vector<double>().swap(VectorMuonTightValidHits_Pt);
	std::vector<double>().swap(VectorMuonTightValidHits_Eta);
	std::vector<double>().swap(VectorMuonTightValidHits_Phi);
	std::vector<int>().swap(VectorMuonTightValidHits_Charge);
	std::vector<double>().swap(VectorMuonTightValidHits_Mass);

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayer_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayer_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayer_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayer_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayer_Mass);

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayerChi2_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2_Mass);	

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass);	

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass);

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass);

	std::vector<double>().swap(VectorleadingMuon_Pt);
	std::vector<double>().swap(VectorleadingMuon_Eta);
	std::vector<double>().swap(VectorleadingMuon_Phi);
	std::vector<int>().swap(VectorleadingMuon_Charge);
	std::vector<double>().swap(VectorleadingMuon_Mass);

	std::vector<double>().swap(VectortrailingMuon_Pt);
	std::vector<double>().swap(VectortrailingMuon_Eta);
	std::vector<double>().swap(VectortrailingMuon_Phi);
	std::vector<int>().swap(VectortrailingMuon_Charge);
	std::vector<double>().swap(VectortrailingMuon_Mass);

	std::vector<double>().swap(VectorMll);
	std::vector<double>().swap(VectorMllpT);
	std::vector<double>().swap(VectorMlleta);
	std::vector<double>().swap(VectorMllphi);
	
}


//TRIGGERS
bool JpsiAnalyzerOpen2011::triggerfired(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname){
  const edm::TriggerNames TrigNames_ = ev.triggerNames(*triggerResultsHandle_);
  const unsigned int ntrigs = triggerResultsHandle_->size();
  for (unsigned int itr=0; itr<ntrigs; itr++){
    TString trigName=TrigNames_.triggerName(itr);
    if (!triggerResultsHandle_->accept(itr)) continue;
    if(trigName.Contains(trigname))      return true;
  }
  return false;
}


bool JpsiAnalyzerOpen2011::triggerfound(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname){
  const edm::TriggerNames TrigNames_ = ev.triggerNames(*triggerResultsHandle_);
  const unsigned int ntrigs = triggerResultsHandle_->size();
  for (unsigned int itr=0; itr<ntrigs; itr++){
    TString trigName=TrigNames_.triggerName(itr);
    if(trigName.Contains(trigname) )     return true;
     if(verbose_) std::cout << "HLT trigger Name avaliable : " << trigName << std::endl;
 
 }
  return false;
}

unsigned int JpsiAnalyzerOpen2011::triggerIndex(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname){
  const edm::TriggerNames TrigNames_ = ev.triggerNames(*triggerResultsHandle_);
  const unsigned int ntrigs = triggerResultsHandle_->size();
  unsigned int itr;
  for (itr=0; itr<ntrigs; itr++){
    TString trigName=TrigNames_.triggerName(itr);
    if(trigName.Contains(trigname))      return itr;
 
 }
  return itr;
}

//
// member functions
//

// ------------ method called for each event  ------------
void
JpsiAnalyzerOpen2011::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ //int countInAccepted = 0;
 // int countInTriggered = 0;

  ievt = iEvent.id().event();
  irun = iEvent.run();
  lumiBlock = iEvent.luminosityBlock();

   if (verbose_) std::cout<<" Run # "<<irun<<" Evt # "<<ievt<<" lumiBlock # "<<lumiBlock<<std::endl;

   using namespace edm;
   const reco::Vertex* vertex = 0;
   edm::Handle< reco::VertexCollection > recoVertices;
   iEvent.getByLabel(primaryVertexProducer_, recoVertices);
   int nVertices = recoVertices->size();
   if (verbose_) std::cout << "nVertices "<< nVertices << std::endl;
   if(nVertices>0) vertex = & ((*recoVertices)[0]);

   edm::Handle< reco::MuonCollection > recoMuons;
   iEvent.getByLabel(recoMuons_, recoMuons);
 
  // *** get Handle to the TriggerEvent
   edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
   iEvent.getByLabel(hlTriggerEvent_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    std::cout << "Error in getting TriggerEvent product from Event!" << std::endl;
    return;
  }

// *** get Handle to the TriggerResults
  edm::Handle<TriggerResults> HLTR; 
  iEvent.getByLabel(hlTriggerResults_, HLTR);
  if (!HLTR.isValid()) {
    std::cout << "HLT TriggerResults with label " << hlTriggerResults_ << " not found!";
    return;
  }


// Only events in which the path actually fired had stored the filter results and products:	  

  bool triggerFound = triggerfound(iEvent,HLTR,triggerName_);
  if(verbose_) std::cout<< "HLT trigger is fired: " << triggerFound << std::endl; 
  if (triggerFound) countInAccepted++;
 bool triggerFired = triggerfired(iEvent,HLTR,triggerName_);
   if(verbose_) std::cout<< "TriggerFired : "<< triggerFired <<" Trigger Path :"<< triggerName_ << std::endl; 
  const unsigned int numberOfHltPaths = HLTR->size();
  //const unsigned int numberOfHltFilters = triggerEventHandle_->sizeFilters();



unsigned int pathIndex = triggerIndex(iEvent,HLTR,triggerName_);
  if (pathIndex>=numberOfHltPaths) {
    std::cout << " WARNING: path " << triggerName_ << " out of range in TriggerResults" << std::endl;
    return;
  }

  if (HLTR->wasrun(pathIndex)) {
    if (!triggerFound) std::cout << " WARNING: path exists in HLTR but not found in TriggerResults" << std::endl;
  }
  else {
    if (triggerFound) std::cout << " WARNING: path found in TriggerResults but it does not exist in HLTR" << std::endl;
  }


	
  //  reco::Vertex vertex;
   std::vector<reco::Muon> myLeptons;

	Total_Events++;
    //   if (triggerFired){
         if(triggerflag_) countInTriggered++;
  // Reco Muons
  for (reco::MuonCollection::const_iterator muon = recoMuons->begin(); muon != recoMuons->end(); muon++) {
                CounterMuons++;               
		VectorMuon_Pt.push_back(muon->pt());
		VectorMuon_Eta.push_back(muon->eta());
		VectorMuon_Phi.push_back(muon->phi());
		VectorMuon_Charge.push_back(muon->charge());
		VectorMuon_Mass.push_back(muon->mass());

		HistoMuon_Pt->Fill(muon->pt()) ;
		HistoMuon_Eta->Fill(muon->eta()) ;
		HistoMuon_Phi->Fill(muon->phi()) ;
		HistoMuon_Charge->Fill(muon->charge()) ;
		HistoMuon_Mass->Fill(muon->mass());
    
		if (muon->isTrackerMuon()){
		TrackerMuon++;
		TrackerMuonPt.push_back(muon->pt());
		TrackerMuonEta.push_back(muon->eta());
		TrackerMuonPhi.push_back(muon->phi());
		TrackerMuonCharge.push_back(muon->charge());	
		}

	 	if (muon->isGlobalMuon()){
		GlobalMuon++;	
		GlobalMuonPt.push_back(muon->pt());
		GlobalMuonEta.push_back(muon->eta());
		GlobalMuonPhi.push_back(muon->phi());
		GlobalMuonCharge.push_back(muon->charge());
	    }


      bool isTMOneStationTight = muon::isGoodMuon(*muon,muon::TMOneStationTight);	
          
		if (isTMOneStationTight) 
		{
			TMOneStationTight++;

			VectorMuonTight_Pt.push_back(muon->pt());
			VectorMuonTight_Eta.push_back(muon->eta());
			VectorMuonTight_Phi.push_back(muon->phi());
			VectorMuonTight_Charge.push_back(muon->charge());
			VectorMuonTight_Mass.push_back(muon->mass());

			HistoMuonTight_Pt->Fill(muon->pt()) ;
			HistoMuonTight_Eta->Fill(muon->eta()) ;
			HistoMuonTight_Phi->Fill(muon->phi()) ;
			HistoMuonTight_Charge->Fill(muon->charge()) ;
			HistoMuonTight_Mass->Fill(muon->mass());

			

			if (verbose_) std::cout<< " dxy: "<< fabs(muon->innerTrack()->dxy(vertex->position()))  << std::endl; 
  			if (verbose_) std::cout<< " dz: "<< fabs(muon->innerTrack()->dz(vertex->position()))  << std::endl;
        	// if (verbose_) std::cout << muon->charge() << std::endl;
       			   
			if(muon->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10)
			{
				NumberOfValidMuonHits++;

				VectorMuonTightValidHits_Pt.push_back(muon->pt());
				VectorMuonTightValidHits_Eta.push_back(muon->eta());
				VectorMuonTightValidHits_Phi.push_back(muon->phi());
				VectorMuonTightValidHits_Charge.push_back(muon->charge());
				VectorMuonTightValidHits_Mass.push_back(muon->mass());

				HistoMuonTightValidHits_Pt->Fill(muon->pt()) ;
				HistoMuonTightValidHits_Eta->Fill(muon->eta()) ;
				HistoMuonTightValidHits_Phi->Fill(muon->phi()) ;
				HistoMuonTightValidHits_Charge->Fill(muon->charge()) ;
				HistoMuonTightValidHits_Mass->Fill(muon->mass());


				if(muon->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 1)
				{
					pixelLayersWithMeasurement++;	

					VectorMuonTightValidHitsPixelLayer_Pt.push_back(muon->pt());
					VectorMuonTightValidHitsPixelLayer_Eta.push_back(muon->eta());
					VectorMuonTightValidHitsPixelLayer_Phi.push_back(muon->phi());
					VectorMuonTightValidHitsPixelLayer_Charge.push_back(muon->charge());
					VectorMuonTightValidHitsPixelLayer_Mass.push_back(muon->mass());

					HistoMuonTightValidHitsPixelLayer_Pt->Fill(muon->pt()) ;
					HistoMuonTightValidHitsPixelLayer_Eta->Fill(muon->eta()) ;
					HistoMuonTightValidHitsPixelLayer_Phi->Fill(muon->phi()) ;
					HistoMuonTightValidHitsPixelLayer_Charge->Fill(muon->charge()) ;
					HistoMuonTightValidHitsPixelLayer_Mass->Fill(muon->mass());

			
					if(muon->innerTrack()->normalizedChi2() < 1.8)
					{
						normalizedChi2++;

						VectorMuonTightValidHitsPixelLayerChi2_Pt.push_back(muon->pt());
						VectorMuonTightValidHitsPixelLayerChi2_Eta.push_back(muon->eta());
						VectorMuonTightValidHitsPixelLayerChi2_Phi.push_back(muon->phi());
						VectorMuonTightValidHitsPixelLayerChi2_Charge.push_back(muon->charge());
						VectorMuonTightValidHitsPixelLayerChi2_Mass.push_back(muon->mass());

						HistoMuonTightValidHitsPixelLayerChi2_Pt->Fill(muon->pt()) ;
						HistoMuonTightValidHitsPixelLayerChi2_Eta->Fill(muon->eta()) ;
						HistoMuonTightValidHitsPixelLayerChi2_Phi->Fill(muon->phi()) ;
						HistoMuonTightValidHitsPixelLayerChi2_Charge->Fill(muon->charge()) ;
						HistoMuonTightValidHitsPixelLayerChi2_Mass->Fill(muon->mass());

	
						
						if( (fabs(muon->innerTrack()->dxy(vertex->position())) < 3./*cm*/) && (fabs(muon->innerTrack()->dz(vertex->position())) < 15./*cm*/) )
						{
							db_dz++;

							VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt.push_back(muon->pt());
							VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta.push_back(muon->eta());
							VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi.push_back(muon->phi());
							VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge.push_back(muon->charge());
							VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass.push_back(muon->mass());

							HistoMuonTightValidHitsPixelLayerChi2DbDz_Pt->Fill(muon->pt()) ;
							HistoMuonTightValidHitsPixelLayerChi2DbDz_Eta->Fill(muon->eta()) ;
							HistoMuonTightValidHitsPixelLayerChi2DbDz_Phi->Fill(muon->phi()) ;
							HistoMuonTightValidHitsPixelLayerChi2DbDz_Charge->Fill(muon->charge()) ;
							HistoMuonTightValidHitsPixelLayerChi2DbDz_Mass->Fill(muon->mass());

							if (muon->isPFMuon())
							{ //Particle-Flow muon id - Can be complemented by muon quality cuts similar to those used in the Tight Muon selection.
								PFMuon++;
				
								VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt.push_back(muon->pt());
								VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta.push_back(muon->eta());
								VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi.push_back(muon->phi());
								VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge.push_back(muon->charge());
								VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass.push_back(muon->mass());

								HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Pt->Fill(muon->pt()) ;
								HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Eta->Fill(muon->eta()) ;
								HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Phi->Fill(muon->phi()) ;
								HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Charge->Fill(muon->charge()) ;
								HistoMuonTightValidHitsPixelLayerChi2DbDzPf_Mass->Fill(muon->mass());
							
								if (muon->isTrackerMuon() || muon->isGlobalMuon())//muon type selection
								{
									TrackerGlobalMuon++;
		
									VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt.push_back(muon->pt());
									VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta.push_back(muon->eta());
									VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi.push_back(muon->phi());
									VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge.push_back(muon->charge());
									VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass.push_back(muon->mass());

									HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->Fill(muon->pt()) ;
									HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->Fill(muon->eta()) ;
									HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->Fill(muon->phi()) ;
									HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->Fill(muon->charge()) ;
									HistoMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass->Fill(muon->mass());

									myLeptons.push_back(*muon); //fill

								}  //end muon type selection					
							} //PF muon
						} //end db_dz
					} //end normalizedChi2
				}//pixelLayersWithMeasurement
			} //end numberOfValidTrackerHits         
    	}//end TMOneStationTight 
   } //Muon Loop
  //}// Trigger selection

std::sort(myLeptons.begin(),myLeptons.end(), [](const reco::Muon &a, const reco::Muon &b){
    return a.pt() > b.pt();
  });

  // Print Reco Muons
  for (std::vector<reco::Muon>::const_iterator muon = myLeptons.begin(); muon != myLeptons.end(); muon++) {
    if(verbose_) std::cout<<"RECO myLeptons->pt():  " << muon->pt() << std::endl;
  }// Muon loop

  if(verbose_) std::cout<<"RECO myLeptons.size() all  " << myLeptons.size() << std::endl;


  // dimuon selection
  if (myLeptons.size() >= 2) { //loop ndimuon
    //recoFilter_ = true;
    nDimuon++;
    if(verbose_) std::cout<<"RECO  Muons Multiplicity:  " << myLeptons.size() << std::endl;
    // if(verbose_) std::cout<<"RECO Dimuons Multiplicity:  " << nDimuon << std::endl;
    reco::Muon leadingMuon = myLeptons[0];
    reco::Muon trailingMuon = myLeptons[1];
    //Dimuons  selection
    // if ((leadingMuon.charge() != trailingMuon.charge())) {
    // } else {return false;}

    if(verbose_) std::cout<< "Leading Muon pt, eta, phi, charge = " << leadingMuon.pt() << " "<< leadingMuon.eta() << " "<< leadingMuon.phi() << " " << leadingMuon.charge() << std::endl;
    if(verbose_) std::cout<< "Trailing Muon  pt, eta, phi,charge = " << trailingMuon.pt() << " " << trailingMuon.eta() << " " << trailingMuon.phi() << " " << trailingMuon.charge()<< std::endl;

	//Invariant mass of dimuons
    double Mll = (leadingMuon.p4() + trailingMuon.p4()).mass();
    double MllpT = (leadingMuon.p4() + trailingMuon.p4()).pt();
    double Mlleta = (leadingMuon.p4() + trailingMuon.p4()).eta();
    double Mllphi = (leadingMuon.p4() + trailingMuon.p4()).phi();
    if(verbose_) std::cout<< "Dimuons Invariant Mass Mll, pT, eta, phi: " << Mll << " " << MllpT << " " << Mlleta << " " << Mllphi << std::endl;

	//Fill Histograms
    leadingMuon_Pt->Fill(leadingMuon.pt());
    leadingMuon_Eta->Fill(leadingMuon.eta());
    leadingMuon_Phi->Fill(leadingMuon.phi());
    leadingMuon_Charge->Fill(leadingMuon.charge());
	leadingMuon_Mass->Fill(leadingMuon.mass());

    trailingMuon_Pt->Fill(trailingMuon.pt());
    trailingMuon_Eta->Fill(trailingMuon.eta());
    trailingMuon_Phi->Fill(trailingMuon.phi());
    trailingMuon_Charge->Fill(trailingMuon.charge());
	trailingMuon_Mass->Fill(trailingMuon.mass());

	histo_Mll -> Fill(Mll);
	histo_MllPt -> Fill(MllpT);
	histo_MllEta -> Fill(Mlleta);
	histo_MllPhi -> Fill(Mllphi);

	//Fill Vectors
	VectorleadingMuon_Pt.push_back(leadingMuon.pt());
	VectorleadingMuon_Eta.push_back(leadingMuon.eta());
	VectorleadingMuon_Phi.push_back(leadingMuon.phi());
	VectorleadingMuon_Charge.push_back(leadingMuon.charge());
	VectorleadingMuon_Mass.push_back(leadingMuon.mass());

	VectortrailingMuon_Pt.push_back(trailingMuon.pt());
	VectortrailingMuon_Eta.push_back(trailingMuon.eta());
	VectortrailingMuon_Phi.push_back(trailingMuon.phi());
	VectortrailingMuon_Charge.push_back(trailingMuon.charge());
	VectortrailingMuon_Mass.push_back(trailingMuon.mass());

	VectorMll.push_back(Mll);
	VectorMllpT.push_back(MllpT);
	VectorMlleta.push_back(Mlleta);
	VectorMllphi.push_back(Mllphi);
	
    }//end loop ndimuons
	


/*	

    if (leadingMuon.pt() >= muonLeadPt_ || trailingMuon.pt() >= muonTrailPt_ ) {

      // ***
      //   // jpsi peak
      //     // ***
      //
      if (Mll > minJPsiMass_ && Mll < maxJPsiMass_){
        nJpsi++;                           
        if(verbose_) std::cout<<" Invariant Mass in JPsi peak, pT, eta, phi " << Mll << " " << MllpT << " " << Mlleta << " " << Mllphi << std::endl;
      } else {return false;}// jpsi selection
    } else {return false;}//lead and trail muon pT cut

*/


}//end analize


// ------------ method called once each job just before starting event loop  ------------
void
JpsiAnalyzerOpen2011::beginJob()
{
	//Create and Fill Branchs
	AnalysisTree->Branch("ievt", &ievt, "ievt/I");
	AnalysisTree->Branch("irun", &irun, "irun/I");
	AnalysisTree->Branch("lumiblock", &lumiBlock, "lumiblock/I" );
	AnalysisTree->Branch("Total_Events", &Total_Events, "Total_Events/I");
	AnalysisTree->Branch("Muons", &CounterMuons, "CounterMuons/I");
	AnalysisTree->Branch("TrackerMuon", &TrackerMuon, "TrackerMuon/I");
	AnalysisTree->Branch("GlobalMuon", &GlobalMuon, "GlobalMuon/I");

	AnalysisTree->Branch("VectorMuon_Pt","std::vector<double>", &VectorMuon_Pt);
	AnalysisTree->Branch("VectorMuon_Eta","std::vector<double>", &VectorMuon_Eta);
	AnalysisTree->Branch("VectorMuon_Phi","std::vector<double>", &VectorMuon_Phi);
	AnalysisTree->Branch("VectorMuon_Charge","std::vector<int>", &VectorMuon_Charge);
	AnalysisTree->Branch("VectorMuon_Mass","std::vector<double>", &VectorMuon_Mass);

	AnalysisTree->Branch("TMOneStationTight", &TMOneStationTight, "TMOneStationTight/I");

	AnalysisTree->Branch("VectorMuonTight_Pt","std::vector<double>", &VectorMuonTight_Pt);
	AnalysisTree->Branch("VectorMuonTight_Eta","std::vector<double>", &VectorMuonTight_Eta);
	AnalysisTree->Branch("VectorMuonTight_Phi","std::vector<double>", &VectorMuonTight_Phi);
	AnalysisTree->Branch("VectorMuonTight_Charge","std::vector<int>", &VectorMuonTight_Charge);
	AnalysisTree->Branch("VectorMuonTight_Mass","std::vector<double>", &VectorMuonTight_Mass);

	AnalysisTree->Branch("NumberOfValidMuonHits", &NumberOfValidMuonHits, "NumberOfValidMuonHits/I");

	AnalysisTree->Branch("VectorMuonTightValidHits_Pt","std::vector<double>", &VectorMuonTightValidHits_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHits_Eta","std::vector<double>", &VectorMuonTightValidHits_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHits_Phi","std::vector<double>", &VectorMuonTightValidHits_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHits_Charge","std::vector<int>", &VectorMuonTightValidHits_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHits_Mass","std::vector<double>", &VectorMuonTightValidHits_Mass);

	AnalysisTree->Branch("pixelLayersWithMeasurement", &pixelLayersWithMeasurement, "pixelLayersWithMeasurement/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayer_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayer_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayer_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayer_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayer_Mass);

	AnalysisTree->Branch("normalizedChi2", &normalizedChi2, "normalizedChi2/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayerChi2_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2_Mass);

	AnalysisTree->Branch("db_dz", &db_dz, "db_dz/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass);

	AnalysisTree->Branch("PFMuon", &PFMuon, "PFMuon/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass);

	AnalysisTree->Branch("TrackerGlobalMuon", &TrackerGlobalMuon, "TrackerGlobalMuon/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass);

	AnalysisTree->Branch("nDimuon", &nDimuon, "nDimuon/I");

	AnalysisTree->Branch("VectorEvent","std::vector<int>", &VectorEvent);
	AnalysisTree->Branch("VectorRun","std::vector<int>", &VectorRun);
	AnalysisTree->Branch("VectorlumiBlock","std::vector<int>", &VectorlumiBlock);


	AnalysisTree->Branch("leadingMuon_Pt","std::vector<double>", &VectorleadingMuon_Pt);
	AnalysisTree->Branch("leadingMuon_Eta","std::vector<double>", &VectorleadingMuon_Eta);
	AnalysisTree->Branch("leadingMuon_Phi","std::vector<double>", &VectorleadingMuon_Phi);
	AnalysisTree->Branch("leadingMuon_Charge","std::vector<int>", &VectorleadingMuon_Charge);
	AnalysisTree->Branch("leadingMuon_Mass","std::vector<double>", &VectorleadingMuon_Mass);

	AnalysisTree->Branch("trailingMuon_Pt","std::vector<double>", &VectortrailingMuon_Pt);
	AnalysisTree->Branch("trailingMuon_Eta","std::vector<double>", &VectortrailingMuon_Eta);
	AnalysisTree->Branch("trailingMuon_Phi","std::vector<double>", &VectortrailingMuon_Phi);
	AnalysisTree->Branch("trailingMuon_Charge","std::vector<int>", &VectortrailingMuon_Charge);
	AnalysisTree->Branch("trailingMuon_Mass","std::vector<double>", &VectortrailingMuon_Mass);

	AnalysisTree->Branch("Mll","std::vector<double>", &VectorMll);
	AnalysisTree->Branch("MllpT","std::vector<double>", &VectorMllpT);
	AnalysisTree->Branch("Mlleta","std::vector<double>", &VectorMlleta);
	AnalysisTree->Branch("Mllphi","std::vector<double>", &VectorMllphi);
        
}

// ------------ method called once each job just after ending the event loop  ------------nalysisTree->Branch("VectorMll","std::vector<double>", &VectorMll);
//
void
JpsiAnalyzerOpen2011::endJob()
{
	//Fill the Trees
	AnalysisTree->Fill();

	std::cout << "Total_Events: " << Total_Events << std::endl; 
	std::cout<<" N of Evts using Trigger Fired :"<< countInTriggered <<" Path: " <<  triggerName_ << std::endl; 
	std::cout << "Muons after Trigger: " << CounterMuons << std::endl;
	std::cout << "TMOneStationTight: "<< TMOneStationTight<<std::endl;
	std::cout << "NumberOfValidMuonHits: " << NumberOfValidMuonHits << std::endl;
	std::cout << "pixelLayersWithMeasurement > 1: " << pixelLayersWithMeasurement << std::endl;
	std::cout << "normalizedChi2: "<< normalizedChi2<<std::endl;
	std::cout << "db < 3cm and dz < 15cm: "<< db_dz<<std::endl;
	std::cout << "PFMuon: " << PFMuon << std::endl;
	std::cout << "TrackerGlobalMuon: " << TrackerGlobalMuon << std::endl;
	std::cout << "nDimuon: " << nDimuon << std::endl;
	 	
}

// ------------ method called when starting to processes a run  ------------
void
JpsiAnalyzerOpen2011::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
JpsiAnalyzerOpen2011::endRun(edm::Run const&, edm::EventSetup const&)
{
}
//define this as a plug-in
DEFINE_FWK_MODULE(JpsiAnalyzerOpen2011);
