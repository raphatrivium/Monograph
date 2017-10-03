//Program to extract data (variables and vector format) of the root file and create histograms. This root file is made of others root files, therefore, has multiple entries which we access with a diferent method than a normal root file.

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atof */
#include <math.h>       /* sin */
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <vector>
#ifndef ROOT_TLatex
#define ROOT_TLatex
#ifndef ROOTiosfwd
#include "Riosfwd.h"
#endif
#ifndef ROOT_TText
#include "TText.h"
#endif
#ifndef ROOT_TAttLine
#include "TAttLine.h"
#endif

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main()
{
	//call a file for a histogram style (Optional)
	gROOT->LoadMacro("styleTDR.C"); 
	setTDRStyle();

	//counters for file f1 (pythia)
	int count_Total_Events_pythia = 0;
	int counter_Muon_pythia = 0; 
	int counter_PtTrigger7_pythia = 0;
	int counter_TMOneStationTight_pythia = 0;
	int counter_NumberOfValidMuonHits_pythia = 0;
	int counter_pixelLayersWithMeasurement_pythia = 0;
	int counter_normalizedChi2_pythia = 0;
	int counter_db_dz_pythia = 0;
	int counter_PFMuon_pythia = 0;
	int counter_TrackerGlobalMuon_pythia = 0;
	int counter_nDimuon_pythia = 0;
	int	count_OppositeCharge_pythia = 0; 
	int count_Eta_pythia = 0;
	int count_Jpsi_pythia = 0; 
	int count_region1_pythia = 0;
	int count_region2_pythia = 0; 
	int count_region3_pythia = 0;
	int count_region4_pythia = 0;
	int count_region5_pythia = 0;
	int count_pico_massa = 0;

	//counters for file f2 (data)
	int count_Total_Events_data = 0;
	int counter_Muon_data = 0;
	int counter_TMOneStationTight_data = 0;
	int counter_NumberOfValidMuonHits_data = 0;
	int counter_pixelLayersWithMeasurement_data = 0;
	int counter_normalizedChi2_data = 0;
	int counter_db_dz_data = 0;
	int counter_PFMuon_data = 0;
	int counter_TrackerGlobalMuon_data = 0;
	int counter_nDimuon_data = 0;
	int	count_OppositeCharge_data = 0; 
	int count_Eta_data = 0;
	int count_Jpsi_data = 0; 
	int count_region1_data = 0;
	int count_region2_data = 0; 
	int count_region3_data = 0;
	int count_region4_data = 0;
	int count_region5_data = 0;    
		
	//Lorentz Vector
	TLorentzVector mu_1;
	TLorentzVector mu_2;	
	
	//Variaveis	and Vectors
	int Total_Events = 0;
	int Muons = 0;
	int PtTrigger7 = 0;
	int TMOneStationTight = 0;
	int NumberOfValidMuonHits = 0;
	int pixelLayersWithMeasurement = 0;
	int normalizedChi2 = 0;
	int db_dz = 0;
	int PFMuon = 0;
	int TrackerGlobalMuon = 0;
	int nDimuon = 0;

	std::vector<double>* VectorMuon_Pt = 0.;
	std::vector<double>* VectorMuon_Eta = 0.;
	std::vector<double>* VectorMuon_Phi = 0.;
	std::vector<int>* VectorMuon_Charge = 0.;
	std::vector<double>* VectorMuon_Mass = 0.;

	std::vector<double>* TrackerMuonPt = 0.;
	std::vector<double>* TrackerMuonEta = 0.;
	std::vector<double>* TrackerMuonPhi = 0.;
	std::vector<int>* TrackerMuonCharge = 0.;

	std::vector<double>* GlobalMuonPt = 0.;
	std::vector<double>* GlobalMuonEta = 0.;
	std::vector<double>* GlobalMuonPhi = 0.;
	std::vector<int>* GlobalMuonCharge = 0.;

	std::vector<double>* VectorMuonTight_Pt = 0.;
	std::vector<double>* VectorMuonTight_Eta = 0.;
	std::vector<double>* VectorMuonTight_Phi = 0.;
	std::vector<int>* VectorMuonTight_Charge = 0.;
	std::vector<double>* VectorMuonTight_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHits_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHits_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHits_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHits_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHits_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayer_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayer_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayer_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayer_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayer_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayerChi2_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass = 0.;
		
	std::vector<double>* leadingMuon_Pt = 0.;
	std::vector<double>* leadingMuon_Eta = 0.;
	std::vector<double>* leadingMuon_Phi = 0.;
	std::vector<int>* leadingMuon_Charge = 0.;
	std::vector<double>* leadingMuon_Mass = 0.;

	std::vector<double>* trailingMuon_Pt = 0.;
	std::vector<double>* trailingMuon_Eta = 0.;
	std::vector<double>* trailingMuon_Phi = 0.;
	std::vector<int>* trailingMuon_Charge = 0.;
	std::vector<double>* trailingMuon_Mass = 0.;

	int Total_Events2 = 0;
	int Muons2 = 0;
	int TMOneStationTight2 = 0;
	int NumberOfValidMuonHits2 = 0;
	int pixelLayersWithMeasurement2 = 0;
	int normalizedChi22 = 0;
	int db_dz2 = 0;
	int PFMuon2 = 0;
	int TrackerGlobalMuon2 = 0;
	int nDimuon2 = 0;

	std::vector<double>* leadingMuon_Pt2 = 0.;
	std::vector<double>* leadingMuon_Eta2 = 0.;
	std::vector<double>* leadingMuon_Phi2 = 0.;
	std::vector<int>* leadingMuon_Charge2 = 0.;
	std::vector<double>* leadingMuon_Mass2 = 0.;

	std::vector<double>* trailingMuon_Pt2 = 0.;
	std::vector<double>* trailingMuon_Eta2 = 0.;
	std::vector<double>* trailingMuon_Phi2 = 0.;
	std::vector<int>* trailingMuon_Charge2 = 0.;
	std::vector<double>* trailingMuon_Mass2 = 0.;
	
	std::vector<double>* VectorMll = 0.;
	std::vector<double>* VectorMllpT = 0.;
	std::vector<double>* VectorMlleta = 0.;
	std::vector<double>* VectorMllphi = 0.;
	
	double M = 0.;
	double Pt = 0.;
	double Eta = 0.;
	double Rapidity = 0.;
	
	//***********************************************************	
	//Criacao dos histogramas

	//Histogramas cinematics quantities of the muons
	TH1F *h_Muon_Pt = new TH1F("h_Muon_Pt","h_Muon_Pt",100,0,100);
	h_Muon_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_Muon_Pt->SetName("h_Muon_Pt");
	TH1F *h_Muon_Eta = new TH1F("h_Muon_Eta","h_Muon_Eta",100,-4,4);
	h_Muon_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_Muon_Eta->SetName("h_Muon_Eta");
	TH1F *h_Muon_Phi = new TH1F("h_Muon_Phi","h_Muon_Phi",100,-4,4);
	h_Muon_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_Muon_Phi->SetName("h_Muon_Phi");
	TH1F *h_Muon_Charge = new TH1F("h_Muon_Charge","h_Muon_Charge",100,-2,2);
	h_Muon_Charge->SetTitle("Distribuicao das cargas Muons; Charge  ; Eventos ");
	h_Muon_Charge->SetName("h_Muon_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTight_Pt = new TH1F("h_MuonTight_Pt","h_MuonTight_Pt",100,0,100);
	h_MuonTight_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTight_Pt->SetName("h_MuonTight_Pt");
	TH1F *h_MuonTight_Eta = new TH1F("h_MuonTight_Eta","h_MuonTight_Eta",100,-4,4);
	h_MuonTight_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTight_Eta->SetName("h_MuonTight_Eta");
	TH1F *h_MuonTight_Phi = new TH1F("h_MuonTight_Phi","h_MuonTight_Phi",100,-4,4);
	h_MuonTight_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTight_Phi->SetName("h_MuonTight_Phi");
	TH1F *h_MuonTight_Charge = new TH1F("h_MuonTight_Charge","h_MuonTight_Charge",100,-2,2);
	h_MuonTight_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTight_Charge->SetName("h_MuonTight_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHits_Pt = new TH1F("h_MuonTightValidHits_Pt","h_MuonTightValidHits_Pt",100,0,100);
	h_MuonTightValidHits_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHits_Pt->SetName("h_MuonTightValidHits_Pt");
	TH1F *h_MuonTightValidHits_Eta = new TH1F("h_MuonTightValidHits_Eta","h_MuonTightValidHits_Eta",100,-4,4);
	h_MuonTightValidHits_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHits_Eta->SetName("h_MuonTightValidHits_Eta");
	TH1F *h_MuonTightValidHits_Phi = new TH1F("h_MuonTightValidHits_Phi","h_MuonTightValidHits_Phi",100,-4,4);
	h_MuonTightValidHits_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHits_Phi->SetName("h_MuonTightValidHits_Phi");
	TH1F *h_MuonTightValidHits_Charge = new TH1F("h_MuonTightValidHits_Charge","h_MuonTightValidHits_Charge",100,-2,2);
	h_MuonTightValidHits_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHits_Charge->SetName("h_MuonTightValidHits_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayer_Pt = new TH1F("h_MuonTightValidHitsPixelLayer_Pt","h_MuonTightValidHitsPixelLayer_Pt",100,0,100);
	h_MuonTightValidHitsPixelLayer_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayer_Pt->SetName("h_MuonTightValidHitsPixelLayer_Pt");
	TH1F *h_MuonTightValidHitsPixelLayer_Eta = new TH1F("h_MuonTightValidHitsPixelLayer_Eta","h_MuonTightValidHitsPixelLayer_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayer_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayer_Eta->SetName("h_MuonTightValidHitsPixelLayer_Eta");
	TH1F *h_MuonTightValidHitsPixelLayer_Phi = new TH1F("h_MuonTightValidHitsPixelLayer_Phi","h_MuonTightValidHitsPixelLayer_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayer_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayer_Phi->SetName("h_MuonTightValidHitsPixelLayer_Phi");
	TH1F *h_MuonTightValidHitsPixelLayer_Charge = new TH1F("h_MuonTightValidHitsPixelLayer_Charge","h_MuonTightValidHitsPixelLayer_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayer_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayer_Charge->SetName("h_MuonTightValidHitsPixelLayer_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayerChi2_Pt = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Pt","h_MuonTightValidHitsPixelLayerChi2_Pt",100,0,100);
	h_MuonTightValidHitsPixelLayerChi2_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2_Pt->SetName("h_MuonTightValidHitsPixelLayerChi2_Pt");
	TH1F *h_MuonTightValidHitsPixelLayerChi2_Eta = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Eta","h_MuonTightValidHitsPixelLayerChi2_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2_Eta->SetName("h_MuonTightValidHitsPixelLayerChi2_Eta");
	TH1F *h_MuonTightValidHitsPixelLayerChi2_Phi = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Phi","h_MuonTightValidHitsPixelLayerChi2_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2_Phi->SetName("h_MuonTightValidHitsPixelLayerChi2_Phi");
	TH1F *h_MuonTightValidHitsPixelLayerChi2_Charge = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Charge","h_MuonTightValidHitsPixelLayerChi2_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayerChi2_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2_Charge->SetName("h_MuonTightValidHitsPixelLayerChi2_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDz_Pt = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDz_Pt","h_MuonTightValidHitsPixelLayerChi2DbDz_Pt",100,0,100);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->SetName("h_MuonTightValidHitsPixelLayerChi2DbDz_Pt");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDz_Eta = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDz_Eta","h_MuonTightValidHitsPixelLayerChi2DbDz_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDz_Eta->SetName("h_MuonTightValidHitsPixelLayerChi2DbDz_Eta");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDz_Phi = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDz_Phi","h_MuonTightValidHitsPixelLayerChi2DbDz_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDz_Phi->SetName("h_MuonTightValidHitsPixelLayerChi2DbDz_Phi");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDz_Charge = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Charge","h_MuonTightValidHitsPixelLayerChi2DbDz_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Charge->SetTitle("Distribuicao das cargas dos muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDz_Charge->SetName("h_MuonTightValidHitsPixelLayerChi2DbDz_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt","h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt",100,0,100);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta","h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi","h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge","h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt","h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt",100,0,100);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta","h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi","h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge","h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge");

//-----------------------------------------------------------------------------
//Cinematics quantities of the Dimuons
	TH1F *h_Dimuons_M = new TH1F("h_Dimuons_Pt","h_Dimuons_M",100,0,10);
	h_Dimuons_M->SetTitle("Distribuicao Massa Invariante dos Dimuons; #mu#mu [GeV] ; Eventos ");
	h_Dimuons_M->SetName("h_Dimuons_M");

	TH1F *h_Dimuons_Pt = new TH1F("h_Dimuons_Pt","h_Dimuons_Pt",100,0,50);
	h_Dimuons_Pt->SetTitle("Distribuicao Pt dos Dimuons; #mu#mu p_{T} [GeV] ; Eventos ");
	h_Dimuons_Pt->SetName("h_Dimuons_Pt");

	TH1F *h_Dimuons_Eta = new TH1F("h_Dimuons_Pt","h_Dimuons_Eta",100,-4,4);
	h_Dimuons_Eta->SetTitle("Distribuicao Pseudo-Rapidez dos Dimuons; #eta ; Eventos ");
	h_Dimuons_Eta->SetName("h_Dimuons_Eta");

	TH1F *h_Dimuons_Rapidity = new TH1F("h_Dimuons_Rapidity","h_Dimuons_Rapidity",100,-4,4);
	h_Dimuons_Rapidity->SetTitle("Distribuicao Rapidez dos Muons; y ; Eventos ");
	h_Dimuons_Rapidity->SetName("h_Dimuons_Rapidity");

//-----------------------------------------------------------------------------

	TH1F *h_DimuonsOppositeCharge_M = new TH1F("h_DimuonsOppositeCharge_Pt","h_DimuonsOppositeCharge_M",100,0,10);
	h_DimuonsOppositeCharge_M->SetTitle("Distribuicao Massa Invariante dos muons; #mu#mu [GeV] ; Eventos ");
	h_DimuonsOppositeCharge_M->SetName("h_DimuonsOppositeCharge_M");

	TH1F *h_DimuonsOppositeCharge_Pt = new TH1F("h_DimuonsOppositeCharge_Pt","h_DimuonsOppositeCharge_Pt",100,0,50);
	h_DimuonsOppositeCharge_Pt->SetTitle("Distribuicao Pt dos muons; #mu#mu p_{T} [GeV] ; Eventos ");
	h_DimuonsOppositeCharge_Pt->SetName("h_DimuonsOppositeCharge_Pt");

	TH1F *h_DimuonsOppositeCharge_Eta = new TH1F("h_DimuonsOppositeCharge_Eta","h_DimuonsOppositeCharge_Eta",100,-4,4);
	h_DimuonsOppositeCharge_Eta->SetTitle("Distribuicao Pseudo-Rapidez dos Dimuons; #eta ; Eventos ");
	h_DimuonsOppositeCharge_Eta->SetName("h_DimuonsOppositeCharge_Eta");

	TH1F *h_DimuonsOppositeCharge_Rapidity = new TH1F("h_DimuonsOppositeCharge_Rapidity","h_DimuonsOppositeCharge_Rapidity",100,-4,4);
	h_DimuonsOppositeCharge_Rapidity->SetTitle("Distribuicao Rapidez dos muons; y ; Eventos ");
	h_DimuonsOppositeCharge_Rapidity->SetName("h_DimuonsOppositeCharge_Rapidity");

//-----------------------------------------------------------------------------

	TH1F *h_DimuonsOppositeChargeEta_M = new TH1F("h_DimuonsOppositeChargeEta_Pt","h_DimuonsOppositeChargeEta_M",100,0,10);
	h_DimuonsOppositeChargeEta_M->SetTitle("Distribuicao Massa Invariante dos muons; #mu#mu [GeV] ; Eventos ");
	h_DimuonsOppositeChargeEta_M->SetName("h_DimuonsOppositeChargeEta_M");

	TH1F *h_DimuonsOppositeChargeEta_Pt = new TH1F("h_DimuonsOppositeChargeEta_Pt","h_DimuonsOppositeChargeEta_Pt",100,0,50);
	h_DimuonsOppositeChargeEta_Pt->SetTitle("Distribuicao Pt dos muons; #mu#mu p_{T} [GeV] ; Eventos ");
	h_DimuonsOppositeChargeEta_Pt->SetName("h_DimuonsOppositeChargeEta_Pt");

	TH1F *h_DimuonsOppositeChargeEta_Eta = new TH1F("h_DimuonsOppositeChargeEta_Eta","h_DimuonsOppositeChargeEta_Eta",100,-4,4);
	h_DimuonsOppositeChargeEta_Eta->SetTitle("Distribuicao Pseudo-Rapidez dos Dimuons; #eta ; Eventos ");
	h_DimuonsOppositeChargeEta_Eta->SetName("h_DimuonsOppositeChargeEta_Eta");

	TH1F *h_DimuonsOppositeChargeEta_Rapidity = new TH1F("h_DimuonsOppositeChargeEta_Rapidity","h_DimuonsOppositeChargeEta_Rapidity",100,-4,4);
	h_DimuonsOppositeChargeEta_Rapidity->SetTitle("Distribuicao Rapidez dos muons; y ; Eventos ");
	h_DimuonsOppositeChargeEta_Rapidity->SetName("h_DimuonsOppositeChargeEta_Rapidity");

//-----------------------------------------------------------------------------
	
	TH1F *h_DimuonsOppositeChargeEtaJpsi_M = new TH1F("h_DimuonsOppositeChargeEtaJpsi_Pt","h_DimuonsOppositeChargeEtaJpsi_M",100,0,10);
	h_DimuonsOppositeChargeEtaJpsi_M->SetTitle("Distribuicao Massa Invariante dos muons; #mu#mu [GeV] ; Eventos ");
	h_DimuonsOppositeChargeEtaJpsi_M->SetName("h_DimuonsOppositeChargeEtaJpsi_M");

	TH1F *h_DimuonsOppositeChargeEtaJpsi_Pt = new TH1F("h_DimuonsOppositeChargeEtaJpsi_Pt","h_DimuonsOppositeChargeEtaJpsi_Pt",100,0,50);
	h_DimuonsOppositeChargeEtaJpsi_Pt->SetTitle("Distribuicao Pt dos muons; #mu#mu p_{T} [GeV] ; Eventos ");
	h_DimuonsOppositeChargeEtaJpsi_Pt->SetName("h_DimuonsOppositeChargeEtaJpsi_Pt");

	TH1F *h_DimuonsOppositeChargeEtaJpsi_Eta = new TH1F("h_DimuonsOppositeChargeEtaJpsi_Eta","h_DimuonsOppositeChargeEtaJpsi_Eta",100,-4,4);
	h_DimuonsOppositeChargeEtaJpsi_Eta->SetTitle("Distribuicao Pseudo-Rapidez dos Dimuons; #eta ; Eventos ");
	h_DimuonsOppositeChargeEtaJpsi_Eta->SetName("h_DimuonsOppositeChargeEtaJpsi_Eta");

	TH1F *h_DimuonsOppositeChargeEtaJpsi_Rapidity = new TH1F("h_DimuonsOppositeChargeEtaJpsi_Rapidity","h_DimuonsOppositeChargeEtaJpsi_Rapidity",100,-4,4);
	h_DimuonsOppositeChargeEtaJpsi_Rapidity->SetTitle("Distribuicao Rapidez dos muons; y  ; Eventos ");
	h_DimuonsOppositeChargeEtaJpsi_Rapidity->SetName("h_DimuonsOppositeChargeEtaJpsi_Rapidity");

//-----------------------------------------------------------------------------
	//TH2F (rapidity x  Transverse Momentum) for J/Psi candidates
	TH2F *h2_Jpsi = new TH2F("h2_Jpsi", "h2_Jpsi", 20, 0, 2.0, 20, 8, 30);
	h2_Jpsi->SetTitle("Distribuicao |y| x p_{T} na Janela J/#psi; |Y| ; p_{T} [GeV] ");
	h2_Jpsi->SetName("p_{T} #times y");

//-----------------------------------------------------------------------------
	//rapidity regions for pythia
	TH1F *h_dimuons_M_y1 = new TH1F("h_dimuons_M_y1","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y1->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y1->SetName("h_dimuons_M_y1");

	TH1F *h_dimuons_M_y2 = new TH1F("h_dimuons_M_y2","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y2->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y2->SetName("h_dimuons_M_y2");
	
	TH1F *h_dimuons_M_y3 = new TH1F("h_dimuons_M_y3","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y3->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y3->SetName("h_dimuons_M_y3");

	TH1F *h_dimuons_M_y4 = new TH1F("h_dimuons_M_y4","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y4->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y4->SetName("h_dimuons_M_y4");

	TH1F *h_dimuons_M_y5 = new TH1F("h_dimuons_M_y5","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y5->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y5->SetName("h_dimuons_M_y5");

	TH1F *h_pico_massa = new TH1F("h_pico_massa","h_dimuons_phi",100,2.8,3.4);
	h_pico_massa->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_pico_massa->SetName("h_pico_massa");

//-----------------------------------------------------------------------------
	//rapidity regions for data
	TH1F *h_dimuons_M_y1_trigger = new TH1F("h_dimuons_M_y1_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y1_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y1_trigger->SetName("h_dimuons_M_y1_trigger");

	TH1F *h_dimuons_M_y2_trigger = new TH1F("h_dimuons_M_y2_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y2_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y2_trigger->SetName("h_dimuons_M_y2_trigger");
	
	TH1F *h_dimuons_M_y3_trigger = new TH1F("h_dimuons_M_y3_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y3_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y3_trigger->SetName("h_dimuons_M_y3_trigger");

	TH1F *h_dimuons_M_y4_trigger = new TH1F("h_dimuons_M_y4_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y4_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y4_trigger->SetName("h_dimuons_M_y4_trigger");

	TH1F *h_dimuons_M_y5_trigger = new TH1F("h_dimuons_M_y5_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y5_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y5_trigger->SetName("h_dimuons_M_y5_trigger");

//-----------------------------------------------------------------------------
	//Behavior of the J/psi candidates
	TH1F *h_Jpsi_Pt = new TH1F("h_Jpsi_Pt","h_Jpsi_Pt",100,10,40);
	h_Jpsi_Pt->SetTitle("Distribuicao p_{T} dos Candidatos a J/#psi; p_{T} [GeV] ; Eventos ");
	h_Jpsi_Pt->SetName("h_Jpsi_Pt");

	TH1F *h_Jpsi_Eta = new TH1F("h_Jpsi_Eta","h_Jpsi_Eta",100,-3,3);
	h_Jpsi_Eta->SetTitle("Distribuicao #eta dos Candidatos a J/#psi; #eta ; Eventos ");
	h_Jpsi_Eta->SetName("h_Jpsi_Eta");
	
	TH1F *h_Jpsi_Rapidity = new TH1F("h_Jpsi_Rapidity","h_Jpsi_Rapidity",100,-3,3);
	h_Jpsi_Rapidity->SetTitle("Distribuicao y dos Candidatos a J/#psi; y ; Eventos ");
	h_Jpsi_Rapidity->SetName("h_Jpsi_Rapidity");

//End Histograms
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------	
	//Reading the root file and the tree
	TFile *f1 = new TFile("histo_pythia_trigger.root");
	TTree *t1 = (TTree*)f1->Get("demo/AnalysisTree");

	//for file2
	TFile *f2 = new TFile("histo_data_1_621.root");
	TTree *t2 = (TTree*)f2->Get("demo/AnalysisTree");

//---------------------------------------------------------------------------------
	// addressing the memory to vector and variables for file f1

	//For Variables
	TBranch *b_Total_Events = t1->GetBranch("Total_Events");
	b_Total_Events->SetAddress(&Total_Events);
	TBranch *b_Muons = t1->GetBranch("Muons");
	b_Muons->SetAddress(&Muons);
	TBranch *b_PtTrigger7 = t1->GetBranch("PtTrigger7");
	b_PtTrigger7->SetAddress(&PtTrigger7);
	TBranch *b_TMOneStationTight = t1->GetBranch("TMOneStationTight");
	b_TMOneStationTight->SetAddress(&TMOneStationTight);
	TBranch *b_NumberOfValidMuonHits = t1->GetBranch("NumberOfValidMuonHits");
	b_NumberOfValidMuonHits->SetAddress(&NumberOfValidMuonHits);
	TBranch *b_pixelLayersWithMeasurement = t1->GetBranch("pixelLayersWithMeasurement");
	b_pixelLayersWithMeasurement->SetAddress(&pixelLayersWithMeasurement);
	TBranch *b_normalizedChi2 = t1->GetBranch("normalizedChi2");
	b_normalizedChi2->SetAddress(&normalizedChi2);
	TBranch *b_db_dz = t1->GetBranch("db_dz");
	b_db_dz->SetAddress(&db_dz);
	TBranch *b_PFMuon = t1->GetBranch("PFMuon");
	b_PFMuon->SetAddress(&PFMuon);
	TBranch *b_TrackerGlobalMuon = t1->GetBranch("TrackerGlobalMuon");
	b_TrackerGlobalMuon->SetAddress(&TrackerGlobalMuon);
	TBranch *b_nDimuon = t1->GetBranch("nDimuon");
	b_nDimuon->SetAddress(&nDimuon);

	//For Vectors	
	TBranch *b_VectorMuon_Pt = t1->GetBranch("VectorMuon_Pt");
	b_VectorMuon_Pt->SetAddress(&VectorMuon_Pt);
	TBranch *b_VectorMuon_Eta = t1->GetBranch("VectorMuon_Eta");
	b_VectorMuon_Eta->SetAddress(&VectorMuon_Eta);
	TBranch *b_VectorMuon_Phi = t1->GetBranch("VectorMuon_Phi");
	b_VectorMuon_Phi->SetAddress(&VectorMuon_Phi);
	TBranch *b_VectorMuon_Charge = t1->GetBranch("VectorMuon_Charge");
	b_VectorMuon_Charge->SetAddress(&VectorMuon_Charge);
	TBranch *b_VectorMuon_Mass = t1->GetBranch("VectorMuon_Mass");
	b_VectorMuon_Mass->SetAddress(&VectorMuon_Mass);

	TBranch *b_VectorMuonTight_Pt = t1->GetBranch("VectorMuonTight_Pt");
	b_VectorMuonTight_Pt->SetAddress(&VectorMuonTight_Pt);
	TBranch *b_VectorMuonTight_Eta = t1->GetBranch("VectorMuonTight_Eta");
	b_VectorMuonTight_Eta->SetAddress(&VectorMuonTight_Eta);
	TBranch *b_VectorMuonTight_Phi = t1->GetBranch("VectorMuonTight_Phi");
	b_VectorMuonTight_Phi->SetAddress(&VectorMuonTight_Phi);
	TBranch *b_VectorMuonTight_Charge = t1->GetBranch("VectorMuonTight_Charge");
	b_VectorMuonTight_Charge->SetAddress(&VectorMuonTight_Charge);
	TBranch *b_VectorMuonTight_Mass = t1->GetBranch("VectorMuonTight_Mass");
	b_VectorMuonTight_Mass->SetAddress(&VectorMuonTight_Mass);

	TBranch *b_VectorMuonTightValidHits_Pt = t1->GetBranch("VectorMuonTightValidHits_Pt");
	b_VectorMuonTightValidHits_Pt->SetAddress(&VectorMuonTightValidHits_Pt);
	TBranch *b_VectorMuonTightValidHits_Eta = t1->GetBranch("VectorMuonTightValidHits_Eta");
	b_VectorMuonTightValidHits_Eta->SetAddress(&VectorMuonTightValidHits_Eta);
	TBranch *b_VectorMuonTightValidHits_Phi = t1->GetBranch("VectorMuonTightValidHits_Phi");
	b_VectorMuonTightValidHits_Phi->SetAddress(&VectorMuonTightValidHits_Phi);
	TBranch *b_VectorMuonTightValidHits_Charge = t1->GetBranch("VectorMuonTightValidHits_Charge");
	b_VectorMuonTightValidHits_Charge->SetAddress(&VectorMuonTightValidHits_Charge);
	TBranch *b_VectorMuonTightValidHits_Mass = t1->GetBranch("VectorMuonTightValidHits_Mass");
	b_VectorMuonTightValidHits_Mass->SetAddress(&VectorMuonTightValidHits_Mass);

	TBranch *b_VectorMuonTightValidHitsPixelLayer_Pt = t1->GetBranch("VectorMuonTightValidHitsPixelLayer_Pt");
	b_VectorMuonTightValidHitsPixelLayer_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayer_Pt);
	TBranch *b_VectorMuonTightValidHitsPixelLayer_Eta = t1->GetBranch("VectorMuonTightValidHitsPixelLayer_Eta");
	b_VectorMuonTightValidHitsPixelLayer_Eta->SetAddress(&VectorMuonTightValidHitsPixelLayer_Eta);
	TBranch *b_VectorMuonTightValidHitsPixelLayer_Phi = t1->GetBranch("VectorMuonTightValidHitsPixelLayer_Phi");
	b_VectorMuonTightValidHitsPixelLayer_Phi->SetAddress(&VectorMuonTightValidHitsPixelLayer_Phi);
	TBranch *b_VectorMuonTightValidHitsPixelLayer_Charge = t1->GetBranch("VectorMuonTightValidHitsPixelLayer_Charge");
	b_VectorMuonTightValidHitsPixelLayer_Charge->SetAddress(&VectorMuonTightValidHitsPixelLayer_Charge);
	TBranch *b_VectorMuonTightValidHitsPixelLayer_Mass = t1->GetBranch("VectorMuonTightValidHitsPixelLayer_Mass");
	b_VectorMuonTightValidHitsPixelLayer_Mass->SetAddress(&VectorMuonTightValidHitsPixelLayer_Mass);

	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2_Pt = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2_Pt");
	b_VectorMuonTightValidHitsPixelLayerChi2_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2_Pt);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2_Eta = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2_Eta");
	b_VectorMuonTightValidHitsPixelLayerChi2_Eta->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2_Eta);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2_Phi = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2_Phi");
	b_VectorMuonTightValidHitsPixelLayerChi2_Phi->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2_Phi);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2_Charge = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2_Charge");
	b_VectorMuonTightValidHitsPixelLayerChi2_Charge->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2_Charge);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2_Mass = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2_Mass");
	b_VectorMuonTightValidHitsPixelLayerChi2_Mass->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2_Mass);

	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass);

	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass);

	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge);
	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass = t1->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass);

	TBranch *b_leadingMuon_Pt = t1->GetBranch("leadingMuon_Pt");
	b_leadingMuon_Pt->SetAddress(&leadingMuon_Pt);
	TBranch *b_leadingMuon_Eta = t1->GetBranch("leadingMuon_Eta");
	b_leadingMuon_Eta->SetAddress(&leadingMuon_Eta);
	TBranch *b_leadingMuon_Phi = t1->GetBranch("leadingMuon_Phi");
	b_leadingMuon_Phi->SetAddress(&leadingMuon_Phi);
	TBranch *b_leadingMuon_Charge = t1->GetBranch("leadingMuon_Charge");
	b_leadingMuon_Charge->SetAddress(&leadingMuon_Charge);
	TBranch *b_leadingMuon_Mass = t1->GetBranch("leadingMuon_Mass");
	b_leadingMuon_Mass->SetAddress(&leadingMuon_Mass);

	TBranch *b_trailingMuon_Pt = t1->GetBranch("trailingMuon_Pt");
	b_trailingMuon_Pt->SetAddress(&trailingMuon_Pt);
	TBranch *b_trailingMuon_Eta = t1->GetBranch("trailingMuon_Eta");
	b_trailingMuon_Eta->SetAddress(&trailingMuon_Eta);
	TBranch *b_trailingMuon_Phi = t1->GetBranch("trailingMuon_Phi");
	b_trailingMuon_Phi->SetAddress(&trailingMuon_Phi);
	TBranch *b_trailingMuon_Charge = t1->GetBranch("trailingMuon_Charge");
	b_trailingMuon_Charge->SetAddress(&trailingMuon_Charge);
	TBranch *b_trailingMuon_Mass = t1->GetBranch("trailingMuon_Mass");
	b_trailingMuon_Mass->SetAddress(&trailingMuon_Mass);

	//=======================================================================================
	//addressing the memory to vector and variables for file f2
	//For Variables
	TBranch *b_Total_Events2 = t2->GetBranch("Total_Events");
	b_Total_Events2->SetAddress(&Total_Events2);
	TBranch *b_Muons2 = t2->GetBranch("Muons");
	b_Muons2->SetAddress(&Muons2);
	TBranch *b_TMOneStationTight2 = t2->GetBranch("TMOneStationTight");
	b_TMOneStationTight2->SetAddress(&TMOneStationTight2);
	TBranch *b_NumberOfValidMuonHits2 = t2->GetBranch("NumberOfValidMuonHits");
	b_NumberOfValidMuonHits2->SetAddress(&NumberOfValidMuonHits2);
	TBranch *b_pixelLayersWithMeasurement2 = t2->GetBranch("pixelLayersWithMeasurement");
	b_pixelLayersWithMeasurement2->SetAddress(&pixelLayersWithMeasurement2);
	TBranch *b_normalizedChi22 = t2->GetBranch("normalizedChi2");
	b_normalizedChi22->SetAddress(&normalizedChi22);
	TBranch *b_db_dz2 = t2->GetBranch("db_dz");
	b_db_dz2->SetAddress(&db_dz2);
	TBranch *b_PFMuon2 = t2->GetBranch("PFMuon");
	b_PFMuon2->SetAddress(&PFMuon2);
	TBranch *b_TrackerGlobalMuon2 = t2->GetBranch("TrackerGlobalMuon");
	b_TrackerGlobalMuon2->SetAddress(&TrackerGlobalMuon2);
	TBranch *b_nDimuon2 = t2->GetBranch("nDimuon");
	b_nDimuon2->SetAddress(&nDimuon2);

	//For Vectors
	TBranch *b_leadingMuon_Pt2 = t2->GetBranch("leadingMuon_Pt");
	b_leadingMuon_Pt2->SetAddress(&leadingMuon_Pt2);
	TBranch *b_leadingMuon_Eta2 = t2->GetBranch("leadingMuon_Eta");
	b_leadingMuon_Eta2->SetAddress(&leadingMuon_Eta2);
	TBranch *b_leadingMuon_Phi2 = t2->GetBranch("leadingMuon_Phi");
	b_leadingMuon_Phi2->SetAddress(&leadingMuon_Phi2);
	TBranch *b_leadingMuon_Charge2 = t2->GetBranch("leadingMuon_Charge");
	b_leadingMuon_Charge2->SetAddress(&leadingMuon_Charge2);
	TBranch *b_leadingMuon_Mass2 = t2->GetBranch("leadingMuon_Mass");
	b_leadingMuon_Mass2->SetAddress(&leadingMuon_Mass2);

	TBranch *b_trailingMuon_Pt2 = t2->GetBranch("trailingMuon_Pt");
	b_trailingMuon_Pt2->SetAddress(&trailingMuon_Pt2);
	TBranch *b_trailingMuon_Eta2 = t2->GetBranch("trailingMuon_Eta");
	b_trailingMuon_Eta2->SetAddress(&trailingMuon_Eta2);
	TBranch *b_trailingMuon_Phi2 = t2->GetBranch("trailingMuon_Phi");
	b_trailingMuon_Phi2->SetAddress(&trailingMuon_Phi2);
	TBranch *b_trailingMuon_Charge2 = t2->GetBranch("trailingMuon_Charge");
	b_trailingMuon_Charge2->SetAddress(&trailingMuon_Charge2);
	TBranch *b_trailingMuon_Mass2 = t2->GetBranch("trailingMuon_Mass");
	b_trailingMuon_Mass2->SetAddress(&trailingMuon_Mass2);
	
	//**********************************************************		
	//Reading Number of tree entries for file f1
	Long64_t nentries = t1->GetEntries();
	cout<< "Number of tree entries: "<< nentries <<std::endl;

	Long64_t GetEntriesFast = t1->GetEntriesFast();
	cout<< "GetEntriesFast: "<< GetEntriesFast <<std::endl;

	Long64_t nbytes = 0, nb = 0, i=0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) //loop tree entries for file f1
	{
      	Long64_t ientry = t1->LoadTree(jentry);
      	std::cout << "nentries " << nentries << " ientry " << jentry << " jentry " <<jentry <<std::endl;
   
      	if (ientry < 0) break;
		b_Total_Events->GetEntry(ientry);
		b_Muons->GetEntry(ientry); 
		b_PtTrigger7->GetEntry(ientry);
		b_TMOneStationTight->GetEntry(ientry);
		b_NumberOfValidMuonHits->GetEntry(ientry);
		b_pixelLayersWithMeasurement->GetEntry(ientry);
		b_normalizedChi2->GetEntry(ientry);
		b_db_dz->GetEntry(ientry);
		b_PFMuon->GetEntry(ientry);
		b_TrackerGlobalMuon->GetEntry(ientry);
		b_nDimuon->GetEntry(ientry);

		counter_Muon_pythia += Muons;
		counter_PtTrigger7_pythia += PtTrigger7;
		counter_TMOneStationTight_pythia += TMOneStationTight;
		counter_NumberOfValidMuonHits_pythia += NumberOfValidMuonHits;
		counter_pixelLayersWithMeasurement_pythia += pixelLayersWithMeasurement;
		counter_normalizedChi2_pythia += normalizedChi2;
		counter_db_dz_pythia += db_dz;
		counter_PFMuon_pythia += PFMuon;
		counter_TrackerGlobalMuon_pythia += TrackerGlobalMuon;
		counter_nDimuon_pythia += nDimuon;

		cout << "Total_Events: "<< Total_Events << endl;

		count_Total_Events_pythia += Total_Events;

       	b_VectorMuon_Pt->GetEntry(ientry);
		b_VectorMuon_Eta->GetEntry(ientry);
		b_VectorMuon_Phi->GetEntry(ientry);
		b_VectorMuon_Charge->GetEntry(ientry);

		b_VectorMuonTight_Pt->GetEntry(ientry);
		b_VectorMuonTight_Eta->GetEntry(ientry);
		b_VectorMuonTight_Phi->GetEntry(ientry);
		b_VectorMuonTight_Charge->GetEntry(ientry);

		b_VectorMuonTightValidHits_Pt->GetEntry(ientry);
		b_VectorMuonTightValidHits_Eta->GetEntry(ientry);
		b_VectorMuonTightValidHits_Phi->GetEntry(ientry);
		b_VectorMuonTightValidHits_Charge->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayer_Pt->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayer_Eta->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayer_Phi->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayer_Charge->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayerChi2_Pt->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2_Eta->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2_Phi->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2_Charge->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->GetEntry(ientry);
		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->GetEntry(ientry);

		b_leadingMuon_Pt->GetEntry(ientry);
		b_leadingMuon_Eta->GetEntry(ientry);
		b_leadingMuon_Phi->GetEntry(ientry);
		b_leadingMuon_Charge->GetEntry(ientry);
		b_leadingMuon_Mass->GetEntry(ientry);

		b_trailingMuon_Pt->GetEntry(ientry);
		b_trailingMuon_Eta->GetEntry(ientry);
		b_trailingMuon_Phi->GetEntry(ientry);
		b_trailingMuon_Charge->GetEntry(ientry);
		b_trailingMuon_Mass->GetEntry(ientry);

		for(Long64_t i=0; i<VectorMuon_Pt->size();i++)
		{  
      		//std::cout<< "VectorMuon_Pt->at("<< i <<") "<< VectorMuon_Pt->at(i) << std::endl;  
       		//double muon_pT = VectorMuon_Pt->at(i);

			h_Muon_Pt->Fill(VectorMuon_Pt->at(i));
			h_Muon_Eta->Fill(VectorMuon_Eta->at(i));
			h_Muon_Phi->Fill(VectorMuon_Phi->at(i));
			h_Muon_Charge->Fill(VectorMuon_Charge->at(i));

       		//myHisto->Fill(VectorMuon_Pt->at(i)); 
   			//nb = fChain->GetEntry(jentry);   nbytes += nb;
      		// if (Cut(ientry) < 0) continue;

  		} 

		for(Long64_t i=0; i<VectorMuonTight_Pt->size();i++)
		{  
			h_MuonTight_Pt->Fill(VectorMuonTight_Pt->at(i));
			h_MuonTight_Eta->Fill(VectorMuonTight_Eta->at(i));
			h_MuonTight_Phi->Fill(VectorMuonTight_Phi->at(i));
			h_MuonTight_Charge->Fill(VectorMuonTight_Charge->at(i));
  		} 

		for(Long64_t i=0; i<VectorMuonTightValidHits_Pt->size();i++)
		{   		
			h_MuonTightValidHits_Pt->Fill(VectorMuonTightValidHits_Pt->at(i));
			h_MuonTightValidHits_Eta->Fill(VectorMuonTightValidHits_Eta->at(i));
			h_MuonTightValidHits_Phi->Fill(VectorMuonTightValidHits_Phi->at(i));
			h_MuonTightValidHits_Charge->Fill(VectorMuonTightValidHits_Charge->at(i));
  		} 

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayer_Pt->size();i++)
		{        		
			h_MuonTightValidHitsPixelLayer_Pt->Fill(VectorMuonTightValidHitsPixelLayer_Pt->at(i));
			h_MuonTightValidHitsPixelLayer_Eta->Fill(VectorMuonTightValidHitsPixelLayer_Eta->at(i));
			h_MuonTightValidHitsPixelLayer_Phi->Fill(VectorMuonTightValidHitsPixelLayer_Phi->at(i));
			h_MuonTightValidHitsPixelLayer_Charge->Fill(VectorMuonTightValidHitsPixelLayer_Charge->at(i));
  		} 

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayerChi2_Pt->size();i++)
		{      		
			h_MuonTightValidHitsPixelLayerChi2_Pt->Fill(VectorMuonTightValidHitsPixelLayer_Pt->at(i));
			h_MuonTightValidHitsPixelLayerChi2_Eta->Fill(VectorMuonTightValidHitsPixelLayer_Eta->at(i));
			h_MuonTightValidHitsPixelLayerChi2_Phi->Fill(VectorMuonTightValidHitsPixelLayer_Phi->at(i));
			h_MuonTightValidHitsPixelLayerChi2_Charge->Fill(VectorMuonTightValidHitsPixelLayer_Charge->at(i));
  		} 

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt->size();i++)
		{       		
			h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt->at(i));
			h_MuonTightValidHitsPixelLayerChi2DbDz_Eta->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta->at(i));
			h_MuonTightValidHitsPixelLayerChi2DbDz_Phi->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi->at(i));
			h_MuonTightValidHitsPixelLayerChi2DbDz_Charge->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge->at(i));
  		} 

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt->size();i++)
		{      		
			h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt->at(i));
			h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta->at(i));
			h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi->at(i));
			h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge->at(i));
  		} 

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->size();i++)
		{       		
			h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->at(i));
			h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->at(i));
			h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->at(i));
			h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->at(i));
  		}

		for(Long64_t i=0; i<trailingMuon_Pt->size();i++) //loop nDimuon
		{  
			//Lorentz Vector   		
			mu_1.SetPtEtaPhiM(leadingMuon_Pt->at(i), leadingMuon_Eta->at(i), leadingMuon_Phi->at(i), leadingMuon_Mass->at(i));
			mu_2.SetPtEtaPhiM(trailingMuon_Pt->at(i), trailingMuon_Eta->at(i), trailingMuon_Phi->at(i), trailingMuon_Mass->at(i));
			//Calculation of the cinematics quantities
			M = (mu_1+mu_2).Mag();		//Invariant Mass  #mu#mu of two Particles
			Pt = (mu_1+mu_2).Pt();      //transverse momentum muon pair
			Eta = (mu_1+mu_2).Eta();      //Pseudo-Rapidity muon pair
			Rapidity = (mu_1+mu_2).Rapidity(); //Rapidity muon pair
			//std::cout<< "Invariante Mass " << i << ": "<< M << std::endl;
			h_Dimuons_M->Fill(M);
			h_Dimuons_Pt->Fill(Pt);
			h_Dimuons_Eta->Fill(Eta);
			h_Dimuons_Rapidity->Fill(Rapidity);

		
			if ( leadingMuon_Charge->at(i) != trailingMuon_Charge->at(i) ) // loop charge
			{
				count_OppositeCharge_pythia++;
				h_DimuonsOppositeCharge_M->Fill(M);
				h_DimuonsOppositeCharge_Pt->Fill(Pt);
				h_DimuonsOppositeCharge_Eta->Fill(Eta);
				h_DimuonsOppositeCharge_Rapidity->Fill(Rapidity);

						
						if ( fabs(leadingMuon_Eta->at(i)) < 2.4 && fabs(trailingMuon_Eta->at(i)) < 2.4 ) // loop eta
						{
							count_Eta_pythia++;
							h_DimuonsOppositeChargeEta_M->Fill(M);
							h_DimuonsOppositeChargeEta_Pt->Fill(Pt);
							h_DimuonsOppositeChargeEta_Eta->Fill(Eta);
							h_DimuonsOppositeChargeEta_Rapidity->Fill(Rapidity);

							if ( (M > 2.8) && (M < 3.4) ) // loop Jpsi
							{
								count_Jpsi_pythia++;
								h_DimuonsOppositeChargeEtaJpsi_M->Fill(M);
								h_DimuonsOppositeChargeEtaJpsi_Pt->Fill(Pt);
								h_DimuonsOppositeChargeEtaJpsi_Eta->Fill(Eta);
								h_DimuonsOppositeChargeEtaJpsi_Rapidity->Fill(Rapidity);

								h2_Jpsi->Fill( fabs(Rapidity), Pt );

								h_Jpsi_Pt->Fill(Pt);
								h_Jpsi_Eta->Fill(Eta);
								h_Jpsi_Rapidity->Fill(Rapidity);
								
								if ( fabs(Rapidity) < 0.3 ) // loop Jpsi
								{	count_region1_pythia++;
									h_dimuons_M_y1->Fill(M);
								}
								if ( fabs(Rapidity) > 0.3 && fabs(Rapidity) < 0.6) // loop Jpsi
								{   count_region2_pythia++;
									h_dimuons_M_y2->Fill(M);
								}
								if ( fabs(Rapidity) > 0.6 && fabs(Rapidity) < 0.9) // loop Jpsi
								{	count_region3_pythia++;
									h_dimuons_M_y3->Fill(M);
								}
								if ( fabs(Rapidity) > 0.9 && fabs(Rapidity) < 1.2) // loop Jpsi
								{	count_region4_pythia++;
									h_dimuons_M_y4->Fill(M);
								}
								if ( fabs(Rapidity) > 1.2 && fabs(Rapidity) < 1.5) // loop Jpsi
								{	count_region5_pythia++;
									h_dimuons_M_y5->Fill(M);
								}

															
								
									
							} //end loop jpsi
						}//end loop eta
			}//end loop charge
  		} //end loop nDimuon

}//End loop tree entries for file f1

//**********************************************************		
	//Reading Number of tree entries for file f2
	Long64_t nentries2 = t2->GetEntries();
	cout<< "Numero de Entradas: "<< nentries2 <<std::endl;

	Long64_t nbytes = 0, nb = 0, i=0;
	for (Long64_t kentry=0; kentry<nentries2;kentry++) // loop tree entries for file f2
	{
      	Long64_t ientry = t2->LoadTree(kentry);
      	std::cout << "nentries " << nentries2 << " ientry " << kentry << " kentry " <<kentry <<std::endl;
   
      	if (ientry < 0) break;

		b_Total_Events2->GetEntry(ientry);
		b_Muons2->GetEntry(ientry);
		b_TMOneStationTight2->GetEntry(ientry);
		b_NumberOfValidMuonHits2->GetEntry(ientry);
		b_pixelLayersWithMeasurement2->GetEntry(ientry);
		b_normalizedChi22->GetEntry(ientry);
		b_db_dz2->GetEntry(ientry);
		b_PFMuon2->GetEntry(ientry);
		b_TrackerGlobalMuon2->GetEntry(ientry);
		b_nDimuon2->GetEntry(ientry);

		count_Total_Events_data += Total_Events2;
		counter_Muon_data += Muons2;
		counter_TMOneStationTight_data += TMOneStationTight2;
		counter_NumberOfValidMuonHits_data += NumberOfValidMuonHits2;
		counter_pixelLayersWithMeasurement_data += pixelLayersWithMeasurement2;
		counter_normalizedChi2_data += normalizedChi22;
		counter_db_dz_data += db_dz2;
		counter_PFMuon_data += PFMuon2;
		counter_TrackerGlobalMuon_data += TrackerGlobalMuon2;
		counter_nDimuon_data += nDimuon2;

		b_leadingMuon_Pt2->GetEntry(ientry);
		b_leadingMuon_Eta2->GetEntry(ientry);
		b_leadingMuon_Phi2->GetEntry(ientry);
		b_leadingMuon_Charge2->GetEntry(ientry);
		b_leadingMuon_Mass2->GetEntry(ientry);

		b_trailingMuon_Pt2->GetEntry(ientry);
		b_trailingMuon_Eta2->GetEntry(ientry);
		b_trailingMuon_Phi2->GetEntry(ientry);
		b_trailingMuon_Charge2->GetEntry(ientry);
		b_trailingMuon_Mass2->GetEntry(ientry);
		
		for(Long64_t i=0; i<trailingMuon_Pt2->size();i++) //loop nDimuon
		{  
			//T		
			mu_1.SetPtEtaPhiM(leadingMuon_Pt2->at(i), leadingMuon_Eta2->at(i), leadingMuon_Phi2->at(i), leadingMuon_Mass2->at(i));
			mu_2.SetPtEtaPhiM(trailingMuon_Pt2->at(i), trailingMuon_Eta2->at(i), trailingMuon_Phi2->at(i), trailingMuon_Mass2->at(i));
			//calculo das grandezas cinematicas
			M = (mu_1+mu_2).Mag();		//Massa Invariante #mu#mu of two Particles
			Pt = (mu_1+mu_2).Pt();      //transverse momentum muon pair
			Eta = (mu_1+mu_2).Eta();      //Pseudo-Rapidity muon pair
			Rapidity = (mu_1+mu_2).Rapidity(); //Rapidity muon pair
			//std::cout<< "Invariant Mass  " << i << ": "<< M << std::endl;
				
			if ( leadingMuon_Charge2->at(i) != trailingMuon_Charge2->at(i) ) // loop charge
			{
				count_OppositeCharge_data++;
					
						
						if ( fabs(leadingMuon_Eta2->at(i)) < 2.4 && fabs(trailingMuon_Eta2->at(i)) < 2.4 ) // loop eta
						{
							count_Eta_data++;
							
							if ( (M > 2.8) && (M < 3.4) ) // loop Jpsi
							{
								count_Jpsi_data++;
																
								if ( fabs(Rapidity) < 0.3 ) // loop Jpsi
								{	count_region1_data++;
									h_dimuons_M_y1_trigger->Fill(M);
								}
								if ( fabs(Rapidity) > 0.3 && fabs(Rapidity) < 0.6) // loop Jpsi
								{   count_region2_data++;
									h_dimuons_M_y2_trigger->Fill(M);
								}
								if ( fabs(Rapidity) > 0.6 && fabs(Rapidity) < 0.9) // loop Jpsi
								{	count_region3_data++;
									h_dimuons_M_y3_trigger->Fill(M);
								}
								if ( fabs(Rapidity) > 0.9 && fabs(Rapidity) < 1.2) // loop Jpsi
								{	count_region4_data++;
									h_dimuons_M_y4_trigger->Fill(M);
								}
								if ( fabs(Rapidity) > 1.2 && fabs(Rapidity) < 1.5) // loop Jpsi
								{	count_region5_data++;
									h_dimuons_M_y5_trigger->Fill(M);
								}
			
								if ( (M > 3.04) && (M < 3.14) ) // loop Jpsi
								{	count_pico_massa++;
									h_pico_massa->Fill(M);	
								}
										
							} //end loop jpsi
						}//end loop eta

			}//end loop charge
  		} //end loop nDimuon

}//End loop tree entries for file f2

	cout << "        " << endl;	
	cout << "======================================================== " << endl;
	cout << "count_Total_Events_pythia: "<< count_Total_Events_pythia << endl;	
	cout << "counter_Muon_pythia: "<< counter_Muon_pythia << endl;
	cout << "counter_PtTrigger7_pythia: "<< counter_PtTrigger7_pythia << endl;
	cout << "counter_TMOneStationTight_pythia: "<< counter_TMOneStationTight_pythia << endl;	
	cout << "counter_NumberOfValidMuonHits_pythia: "<< counter_NumberOfValidMuonHits_pythia << endl;	
	cout << "counter_pixelLayersWithMeasurement_pythia: "<< counter_pixelLayersWithMeasurement_pythia << endl;	
	cout << "counter_normalizedChi2_pythia: "<< counter_normalizedChi2_pythia << endl;	
	cout << "counter_db_dz_pythia: "<< counter_db_dz_pythia << endl;	
	cout << "counter_PFMuon_pythia: "<< counter_PFMuon_pythia << endl;	
	cout << "counter_TrackerGlobalMuon_pythia: "<< counter_TrackerGlobalMuon_pythia << endl;	
	cout << "counter_nDimuon_pythia: "<< counter_nDimuon_pythia << endl;
	cout << "count_OppositeCharge_pythia: "<< count_OppositeCharge_pythia << endl;
	cout << "count_Eta_pythia: "<< count_Eta_pythia << endl;
	cout << "count_Jpsi_pythia: "<< count_Jpsi_pythia << endl;
	cout << "count_region1_pythia: "<< count_region1_pythia << endl;
	cout << "count_region2_pythia: "<< count_region2_pythia << endl;
	cout << "count_region3_pythia: "<< count_region3_pythia << endl;
	cout << "count_region4_pythia: "<< count_region4_pythia << endl;
	cout << "count_region5_pythia: "<< count_region5_pythia << endl;
	cout << "======================================================== " << endl;
	cout << "        " << endl;

	cout << "        " << endl;	
	cout << "======================================================== " << endl;
	cout << "count_Total_Events_data: "<< count_Total_Events_data << endl;	
	cout << "counter_Muon_data: "<< counter_Muon_data << endl;	
	cout << "counter_TMOneStationTight_data: "<< counter_TMOneStationTight_data << endl;	
	cout << "counter_NumberOfValidMuonHits_data: "<< counter_NumberOfValidMuonHits_data << endl;	
	cout << "counter_pixelLayersWithMeasurement_data: "<< counter_pixelLayersWithMeasurement_data << endl;	
	cout << "counter_normalizedChi2_data: "<< counter_normalizedChi2_data << endl;	
	cout << "counter_db_dz_data: "<< counter_db_dz_data << endl;	
	cout << "counter_PFMuon_data: "<< counter_PFMuon_data << endl;	
	cout << "counter_TrackerGlobalMuon_data: "<< counter_TrackerGlobalMuon_data << endl;	
	cout << "counter_nDimuon_data: "<< counter_nDimuon_data << endl;
	cout << "count_OppositeCharge_data: "<< count_OppositeCharge_data << endl;
	cout << "count_Eta_data: "<< count_Eta_data << endl;
	cout << "count_Jpsi_data: "<< count_Jpsi_data << endl;
	cout << "count_region1_data: "<< count_region1_data << endl;
	cout << "count_region2_data: "<< count_region2_data << endl;
	cout << "count_region3_data: "<< count_region3_data << endl;
	cout << "count_region4_data: "<< count_region4_data << endl;
	cout << "count_region5_data: "<< count_region5_data << endl;
	cout << "======================================================== " << endl;
	cout << "        " << endl; 

    //=========================================================================	
	//Creating Canvas
	/*TCanvas* c1 = new TCanvas("c1","Canvas 1 - behavior of the muons during quality selection",1200,600);
	c1->Divide(2,2);
	c1->cd(1);

	h_Muon_Pt->SetLineColor(kRed);
	h_MuonTight_Pt->SetLineColor(kBlue);
	h_MuonTightValidHits_Pt->SetLineColor(kGreen);
	h_MuonTightValidHitsPixelLayer_Pt->SetLineColor(kBlack);
	h_MuonTightValidHitsPixelLayerChi2_Pt->SetLineColor(kYellow);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->SetLineColor(+4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->SetLineColor(kViolet);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->SetLineColor(kBlue);

	TLegend* leg_dimuons_Pt = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Pt->SetFillColor(kWhite);
	leg_dimuons_Pt->SetFillStyle(1001);
	leg_dimuons_Pt->AddEntry(h_Muon_Pt,"#mu","L");
	//leg_dimuons_Pt->AddEntry(h_MuonTight_Pt,"#mu Tight","L");
	//leg_dimuons_Pt->AddEntry(h_MuonTightValidHits_Pt,"#mu numberOfValidTrackerHits > 10","L");
	//leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayer_Pt,"#mu pixelLayersWithMeasurement > 1","L");
	//leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayerChi2_Pt,"#mu Chi2","L");
	//leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDz_Pt,"#mu db < 3cm and dz < 15cm","L");
	leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt,"Soft #mu","L");
	//leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt,"#mu TrackerGlobal","L");
		
	h_Muon_Pt->Draw();
	//h_MuonTight_Pt->Draw("sames");
	//h_MuonTightValidHits_Pt->Draw("sames");
	//h_MuonTightValidHitsPixelLayer_Pt->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2_Pt->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->Draw("sames");
	leg_dimuons_Pt->Draw("same");
	//-------------------------------------------------------------------------
	c1->cd(2);
	h_Muon_Eta->SetLineColor(kRed);
	h_MuonTight_Eta->SetLineColor(kBlue);
	h_MuonTightValidHits_Eta->SetLineColor(kGreen);
	h_MuonTightValidHitsPixelLayer_Eta->SetLineColor(kBlack);
	h_MuonTightValidHitsPixelLayerChi2_Eta->SetLineColor(kYellow);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Eta->SetLineColor(+4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta->SetLineColor(kViolet);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->SetLineColor(kBlue);

	TLegend* leg_dimuons_Eta = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Eta->SetFillColor(kWhite);
	leg_dimuons_Eta->SetFillStyle(1001);
	leg_dimuons_Eta->AddEntry(h_Muon_Eta,"#mu","L");
	//leg_dimuons_Eta->AddEntry(h_MuonTight_Eta,"#mu Tight","L");
	//leg_dimuons_Eta->AddEntry(h_MuonTightValidHits_Eta,"#mu numberOfValidTrackerHits > 10","L");
	//leg_dimuons_Eta->AddEntry(h_MuonTightValidHitsPixelLayer_Eta,"#mu pixelLayersWithMeasurement > 1","L");
	//leg_dimuons_Eta->AddEntry(h_MuonTightValidHitsPixelLayerChi2_Eta,"#mu Chi2","L");
	//leg_dimuons_Eta->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDz_Eta,"#mu db < 3cm and dz < 15cm","L");
	leg_dimuons_Eta->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta,"Soft #mu","L");
	//leg_dimuons_Eta->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta,"#mu TrackerGlobal","L");
	
	
	h_Muon_Eta->Draw();
	//h_MuonTight_Eta->Draw("sames");
	//h_MuonTightValidHits_Eta->Draw("sames");
	//h_MuonTightValidHitsPixelLayer_Eta->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2_Eta->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2DbDz_Eta->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->Draw("sames");
	leg_dimuons_Eta->Draw("same");
	//-------------------------------------------------------------------------
	c1->cd(3);
	h_Muon_Phi->SetLineColor(kRed);
	//h_MuonTight_Phi->SetLineColor(kBlue);
	//h_MuonTightValidHits_Phi->SetLineColor(kGreen);
	//h_MuonTightValidHitsPixelLayer_Phi->SetLineColor(kBlack);
	//h_MuonTightValidHitsPixelLayerChi2_Phi->SetLineColor(kYellow);
	//h_MuonTightValidHitsPixelLayerChi2DbDz_Phi->SetLineColor(+4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi->SetLineColor(kViolet);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->SetLineColor(kBlue);

	TLegend* leg_dimuons_Phi = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Phi->SetFillColor(kWhite);
	leg_dimuons_Phi->SetFillStyle(1001);
	leg_dimuons_Phi->AddEntry(h_Muon_Phi,"#mu","L");
	//leg_dimuons_Phi->AddEntry(h_MuonTight_Phi,"#mu Tight","L");
	//leg_dimuons_Phi->AddEntry(h_MuonTightValidHits_Phi,"#mu numberOfValidTrackerHits > 10","L");
	//leg_dimuons_Phi->AddEntry(h_MuonTightValidHitsPixelLayer_Phi,"#mu pixelLayersWithMeasurement > 1","L");
	//leg_dimuons_Phi->AddEntry(h_MuonTightValidHitsPixelLayerChi2_Phi,"#mu Chi2","L");
	//leg_dimuons_Phi->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDz_Phi,"#mu db < 3cm and dz < 15cm","L");
	leg_dimuons_Phi->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi,"Soft #mu","L");
	//leg_dimuons_Phi->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi,"#mu TrackerGlobal","L");
	
	h_Muon_Phi->Draw();
	//h_MuonTight_Phi->Draw("sames");
	//h_MuonTightValidHits_Phi->Draw("sames");
	//h_MuonTightValidHitsPixelLayer_Phi->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2_Phi->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2DbDz_Phi->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->Draw("sames");
	leg_dimuons_Phi->Draw("same");
	//-------------------------------------------------------------------------
	c1->cd(4);		
	h_Muon_Charge->SetLineColor(kRed);
	h_MuonTight_Charge->SetLineColor(kBlue);
	h_MuonTightValidHits_Charge->SetLineColor(kGreen);
	h_MuonTightValidHitsPixelLayer_Charge->SetLineColor(kBlack);
	h_MuonTightValidHitsPixelLayerChi2_Charge->SetLineColor(kYellow);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Charge->SetLineColor(+4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge->SetLineColor(kViolet);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->SetLineColor(kBlue);

	TLegend* leg_dimuons_Charge = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Charge->SetFillColor(kWhite);
	leg_dimuons_Charge->SetFillStyle(1001);
	leg_dimuons_Charge->AddEntry(h_Muon_Charge,"#mu","L");
	//leg_dimuons_Charge->AddEntry(h_MuonTight_Charge,"#mu Tight","L");
	//leg_dimuons_Charge->AddEntry(h_MuonTightValidHits_Charge,"#mu numberOfValidTrackerHits > 10","L");
	//leg_dimuons_Charge->AddEntry(h_MuonTightValidHitsPixelLayer_Charge,"#mu pixelLayersWithMeasurement > 1","L");
	//leg_dimuons_Charge->AddEntry(h_MuonTightValidHitsPixelLayerChi2_Charge,"#mu Chi2","L");
	//leg_dimuons_Charge->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDz_Charge,"#mu db < 3cm and dz < 15cm","L");
	leg_dimuons_Charge->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge,"Soft #mu","L");
	//leg_dimuons_Charge->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge,"#mu TrackerGlobal","L");
	
	h_Muon_Charge->Draw();
	//h_MuonTight_Charge->Draw("sames");
	//h_MuonTightValidHits_Charge->Draw("sames");
	//h_MuonTightValidHitsPixelLayer_Charge->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2_Charge->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2DbDz_Charge->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->Draw("sames");
	leg_dimuons_Charge->Draw("same");*/

//*******************************************************************************************
	//Creating Canvas
	TCanvas* c2 = new TCanvas("c2","Canvas 2 - behavior of the dimuons after quality selection",1200,600);
	c2->Divide(2,2);
	c2->cd(1);

	h_Dimuons_M->SetLineColor(kRed);
	h_DimuonsOppositeCharge_M->SetLineColor(kBlue);
	h_DimuonsOppositeChargeEta_M->SetLineColor(kBlack);
	h_DimuonsOppositeChargeEtaJpsi_M->SetLineColor(kViolet);
	
	TLegend* leg_dimuons_M = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_M->SetFillColor(kWhite);
	leg_dimuons_M->SetFillStyle(1001);
	leg_dimuons_M->AddEntry(h_Dimuons_M,"#mu#mu","L");
	leg_dimuons_M->AddEntry(h_DimuonsOppositeCharge_M,"#mu^{+}#mu^{-}","L");
	leg_dimuons_M->AddEntry(h_DimuonsOppositeChargeEta_M,"#mu^{+}#mu^{-} Eta<2.4","L");
	leg_dimuons_M->AddEntry(h_DimuonsOppositeChargeEtaJpsi_M,"Candidatos a J/#psi","L");
	
	h_Dimuons_M->Draw();
	h_DimuonsOppositeCharge_M->Draw("sames");
	h_DimuonsOppositeChargeEta_M->Draw("sames");
	h_DimuonsOppositeChargeEtaJpsi_M->Draw("sames");
	leg_dimuons_M->Draw();
	
	//-------------------------------------------------------------------------
	c2->cd(2);

	h_Dimuons_Pt->SetLineColor(kRed);
	h_DimuonsOppositeCharge_Pt->SetLineColor(kBlue);
	h_DimuonsOppositeChargeEta_Pt->SetLineColor(kBlack);
	h_DimuonsOppositeChargeEtaJpsi_Pt->SetLineColor(kViolet);

	TLegend* leg_dimuons_Pt = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Pt->SetFillColor(kWhite);
	leg_dimuons_Pt->SetFillStyle(1001);
	leg_dimuons_Pt->AddEntry(h_Dimuons_Pt,"#mu#mu","L");
	leg_dimuons_Pt->AddEntry(h_DimuonsOppositeCharge_Pt,"#mu^{+}#mu^{-} ","L");
	leg_dimuons_Pt->AddEntry(h_DimuonsOppositeChargeEta_Pt,"#mu^{+}#mu^{-} Eta<2.4","L");
	leg_dimuons_Pt->AddEntry(h_DimuonsOppositeChargeEtaJpsi_Pt,"Candidatos a J/#psi","L");

	h_Dimuons_Pt->Draw();
	h_DimuonsOppositeCharge_Pt->Draw("sames");
	h_DimuonsOppositeChargeEta_Pt->Draw("sames");
	h_DimuonsOppositeChargeEtaJpsi_Pt->Draw("sames");
	leg_dimuons_Pt->Draw();
	//-------------------------------------------------------------------------
	c2->cd(3);

	h_Dimuons_Eta->SetLineColor(kRed);
	h_DimuonsOppositeCharge_Eta->SetLineColor(kBlue);
	h_DimuonsOppositeChargeEta_Eta->SetLineColor(kBlack);
	h_DimuonsOppositeChargeEtaJpsi_Eta->SetLineColor(kViolet);

	TLegend* leg_dimuons_Eta = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Eta->SetFillColor(kWhite);
	leg_dimuons_Eta->SetFillStyle(1001);
	leg_dimuons_Eta->AddEntry(h_Dimuons_Eta,"#mu#mu","L");
	leg_dimuons_Eta->AddEntry(h_DimuonsOppositeCharge_Eta,"#mu^{+}#mu^{-} ","L");
	leg_dimuons_Eta->AddEntry(h_DimuonsOppositeChargeEta_Eta,"#mu^{+}#mu^{-} Eta<2.4","L");
	leg_dimuons_Eta->AddEntry(h_DimuonsOppositeChargeEtaJpsi_Eta,"#Candidatos a J/#psi","L");

	h_Dimuons_Eta->Draw();
	h_DimuonsOppositeCharge_Eta->Draw("sames");
	h_DimuonsOppositeChargeEta_Eta->Draw("sames");
	h_DimuonsOppositeChargeEtaJpsi_Eta->Draw("sames");
	leg_dimuons_Eta->Draw();
	//-------------------------------------------------------------------------
	c2->cd(4);

	h_Dimuons_Rapidity->SetLineColor(kRed);
	h_DimuonsOppositeCharge_Rapidity->SetLineColor(kBlue);
	h_DimuonsOppositeChargeEta_Rapidity->SetLineColor(kBlack);
	h_DimuonsOppositeChargeEtaJpsi_Rapidity->SetLineColor(kViolet);

	TLegend* leg_dimuons_Rapidity = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Rapidity->SetFillColor(kWhite);
	leg_dimuons_Rapidity->SetFillStyle(1001);
	leg_dimuons_Rapidity->AddEntry(h_Dimuons_Rapidity,"#mu#mu","L");
	leg_dimuons_Rapidity->AddEntry(h_DimuonsOppositeCharge_Rapidity,"#mu^{+}#mu^{-} ","L");
	leg_dimuons_Rapidity->AddEntry(h_DimuonsOppositeChargeEta_Rapidity,"#mu^{+}#mu^{-} Eta<2.4","L");
	leg_dimuons_Rapidity->AddEntry(h_DimuonsOppositeChargeEtaJpsi_Rapidity,"#Candidatos a J/#psi","L");

	h_Dimuons_Rapidity->Draw();
	h_DimuonsOppositeCharge_Rapidity->Draw("sames");
	h_DimuonsOppositeChargeEta_Rapidity->Draw("sames");
	h_DimuonsOppositeChargeEtaJpsi_Rapidity->Draw("sames");
	leg_dimuons_Rapidity->Draw();

	//=====================================================================
	TCanvas* c3 = new TCanvas("c3","Canvas 3 - behavior of the J/psi candidates in the THF2 (y x pt) ",1200,600);
	//c3->cd(2);
	gStyle->SetPalette(55);
	h2_Jpsi->SetStats(0); //Remove the statistic box
	h2_Jpsi->Draw("CONT4COLZ");
		
	//=====================================================================================
	TCanvas* c6 = new TCanvas("c6","Canvas 6 - behavior of the J/psi candidates",1200,600);
	c6->Divide(2,2);

	c6->cd(1);
	h_Jpsi_Pt->SetLineColor(kViolet);

	TLegend* leg_Jpsi_Pt = new TLegend(0.75,0.81,0.97,0.97);
	leg_Jpsi_Pt->SetFillColor(kWhite);
	leg_Jpsi_Pt->SetFillStyle(1001);
	leg_Jpsi_Pt->AddEntry(h_Jpsi_Pt,"p_{T} Pythia 6","L");

	h_Jpsi_Pt->Draw();
	leg_Jpsi_Pt->Draw();

	c6->cd(2);
	h_Jpsi_Eta->SetLineColor(kViolet);

	TLegend* leg_Jpsi_Eta = new TLegend(0.75,0.81,0.97,0.97);
	leg_Jpsi_Eta->SetFillColor(kWhite);
	leg_Jpsi_Eta->SetFillStyle(1001);
	leg_Jpsi_Eta->AddEntry(h_Jpsi_Eta,"#eta Pythia 6","L");
	
	h_Jpsi_Eta->Draw();
	leg_Jpsi_Eta->Draw();
	
	c6->cd(3);
	h_Jpsi_Rapidity->SetLineColor(kViolet);

	TLegend* leg_Jpsi_Rapidity = new TLegend(0.75,0.81,0.97,0.97);
	leg_Jpsi_Rapidity->SetFillColor(kWhite);
	leg_Jpsi_Rapidity->SetFillStyle(1001);
	leg_Jpsi_Rapidity->AddEntry(h_Jpsi_Rapidity,"y Pythia 6","L");
	
	h_Jpsi_Rapidity->Draw();
	leg_Jpsi_Rapidity->Draw();

//=====================================================================================
//Canvas para para regioes de rapidez
	TCanvas* c7 = new TCanvas("c7","Canvas 7 - rapidity regions 1 and 2",1200,600);
	c7->Divide(2);
	
	c7->cd(1);
	h_dimuons_M_y1->SetLineColor(kViolet);
	h_dimuons_M_y1->SetMarkerStyle(7);
	h_dimuons_M_y1->SetStats(0);

	h_dimuons_M_y1_trigger->SetMarkerStyle(7);
	h_dimuons_M_y1_trigger->SetStats(0);

	//Arbitrary normalization
	Double_t norm = 1.;
    h_dimuons_M_y1->Scale(norm/h_dimuons_M_y1->Integral(), "width");
    h_dimuons_M_y1->Scale(120.0);
	
	h_dimuons_M_y1_trigger->Draw("E");
	h_dimuons_M_y1->Draw("LSAME");

	TLegend* leg_dimuons_M_y1 = new TLegend(0.60,0.63,0.82,0.8);
   	leg_dimuons_M_y1->SetFillColor(kWhite);
    leg_dimuons_M_y1->SetFillStyle(1001);
    leg_dimuons_M_y1->AddEntry(leg_dimuons_M_y1,"|y| < 0.3","");
	leg_dimuons_M_y1->SetBorderSize(0);

	//TLegend* leg_dimuons_Eta = new TLegend(0.75,0.81,0.97,0.97);	
	leg_dimuons_M_y1->AddEntry(h_dimuons_M_y1_trigger,"Data 2011","E");
	leg_dimuons_M_y1->AddEntry(h_dimuons_M_y1,"MC Pythia 6","L");
	leg_dimuons_M_y1->Draw();
	//----------------------------------------------------------------------
	c7->cd(2);

	h_dimuons_M_y2->SetLineColor(kViolet);
	h_dimuons_M_y2->SetMarkerStyle(7);
	h_dimuons_M_y2->SetStats(0);

	h_dimuons_M_y2_trigger->SetMarkerStyle(7);
	h_dimuons_M_y2_trigger->SetStats(0);

	Double_t norm = 1.;
    h_dimuons_M_y2->Scale(norm/h_dimuons_M_y2->Integral(), "width");
    h_dimuons_M_y2->Scale(120.0);

	TLegend* leg_dimuons_M_y2 = new TLegend(0.60,0.63,0.82,0.8);
   	leg_dimuons_M_y2->SetFillColor(kWhite);
    leg_dimuons_M_y2->SetFillStyle(1001);
    leg_dimuons_M_y2->AddEntry(leg_dimuons_M_y2,"0.3 < |y| < 0.6","");
	leg_dimuons_M_y2->SetBorderSize(0);
	leg_dimuons_M_y2->AddEntry(h_dimuons_M_y2_trigger,"Data 2011","E");
	leg_dimuons_M_y2->AddEntry(h_dimuons_M_y2,"MC Pythia 6","L");

	h_dimuons_M_y2_trigger->Draw("E");
	h_dimuons_M_y2->Draw("LSAME");
	leg_dimuons_M_y2->Draw();
	
//-------------------------------------------------------------------------
	TCanvas* c8 = new TCanvas("c8","Canvas 8 - rapidity regions 3 and 4",1200,600);
	c8->Divide(2);
	c8->cd(1);

	h_dimuons_M_y3->SetLineColor(kBlue);
	h_dimuons_M_y3->SetMarkerStyle(7);
	h_dimuons_M_y3->SetStats(0);

	h_dimuons_M_y3_trigger->SetMarkerStyle(7);
	h_dimuons_M_y3_trigger->SetStats(0);

	//Double_t norm = 1.;
    //h_dimuons_M_y3->Scale(norm/h_dimuons_M_y3->Integral(), "width");
    //h_dimuons_M_y3->Scale(120.0);

	//Double_t scaleMC = 1/h_dimuons_M_y3_trigger->Integral();
	//h_dimuons_M_y3->Scale(scaleMC);
	//h_dimuons_M_y3_trigger->Scale(scaleMC);

	//Double_t scaleMC = 1/h_dimuons_M_y3->Integral();
	//h_dimuons_M_y3->Scale(scaleMC);

	//Double_t scaleDados = 1/h_dimuons_M_y3_trigger->Integral();
	//h_dimuons_M_y3_trigger->Scale(scaleDados);

	//h_dimuons_M_y3_trigger->GetYaxis()->SetRange(0,3);
	//h_dimuons_M_y3_trigger->GetYaxis()->SetRangeUser(0,1)

	//cout << "scaleMC: "<< scaleMC << endl;
	//cout << "scaleDados: "<< scaleDados << endl;

	//h_dimuons_M_y3_trigger->Draw("P");
	//h_dimuons_M_y3->Draw("LSAME");

	//h_dimuons_M_y3->Draw("L");
	//h_dimuons_M_y3_trigger->Draw("PSAME");

	Double_t normD = h_dimuons_M_y3_trigger->GetEntries();
	Double_t normM = h_dimuons_M_y3->GetEntries();
	h_dimuons_M_y3->Scale(normD/normM);

	//Double_t normD = h_pico_massa->GetEntries();
	//Double_t normM = h_dimuons_M_y3->GetEntries();
	//h_dimuons_M_y3->Scale(normD/normM);

	//h_dimuons_M_y3->Draw("L");
	//h_dimuons_M_y3_trigger->Draw("PSAME");

	h_dimuons_M_y3_trigger->SetMarkerStyle(21);
	
	//e para MC
	h_dimuons_M_y3->SetMarkerStyle(21);
	h_dimuons_M_y3->SetMarkerColor(kBlue);

	h_dimuons_M_y3->SetFillColor(kBlue);

	h_dimuons_M_y3->Draw("HIST");
	h_dimuons_M_y3_trigger->Draw("e1pSAME");
	

	TLegend* leg_dimuons_M_y3 = new TLegend(0.60,0.63,0.82,0.8);
   	leg_dimuons_M_y3->SetFillColor(kWhite);
    leg_dimuons_M_y3->SetFillStyle(1001);
    leg_dimuons_M_y3->AddEntry(leg_dimuons_M_y3,"0.6 < |y| < 0.9","");
	leg_dimuons_M_y3->SetBorderSize(0);

	
	leg_dimuons_M_y3->AddEntry(h_dimuons_M_y3_trigger,"Data 2011","e1pSAME");
	leg_dimuons_M_y3->AddEntry(h_dimuons_M_y3,"MC Pythia 6","L");
	
	leg_dimuons_M_y3->Draw();
	
	//-------------------------------------------------------------------
	c8->cd(2);

	h_dimuons_M_y4->SetLineColor(kBlue);
	h_dimuons_M_y4->SetMarkerStyle(7);
	h_dimuons_M_y4->SetStats(0);

	h_dimuons_M_y4_trigger->SetMarkerStyle(7);
	h_dimuons_M_y4_trigger->SetStats(0);

	//Double_t norm = 1.;
   // h_dimuons_M_y4->Scale(norm/h_dimuons_M_y4->Integral(), "width");
    //h_dimuons_M_y4->Scale(110.0);

	//Double_t scaleMC = 1/h_dimuons_M_y4->Integral();
	//h_dimuons_M_y4_trigger->Scale(scaleMC);
	//h_dimuons_M_y4->Scale(scaleMC);

	

	//h_dimuons_M_y4_trigger->Draw("E");
	//h_dimuons_M_y4->Draw("LSAME");

	Double_t normD = h_dimuons_M_y4_trigger->GetEntries();
	Double_t normM = h_dimuons_M_y4->GetEntries();
	h_dimuons_M_y4->Scale(normD/normM);

	h_dimuons_M_y4_trigger->SetMarkerStyle(21);
	
	//e para MC
	h_dimuons_M_y4->SetMarkerStyle(21);
	h_dimuons_M_y4->SetMarkerColor(kBlue);

	h_dimuons_M_y4->SetFillColor(kBlue);

	h_dimuons_M_y4->Draw("HIST");
	h_dimuons_M_y4_trigger->Draw("e1pSAME");


	TLegend* leg_dimuons_M_y4 = new TLegend(0.60,0.63,0.82,0.8);
   	leg_dimuons_M_y4->SetFillColor(kWhite);
    leg_dimuons_M_y4->SetFillStyle(1001);
    leg_dimuons_M_y4->AddEntry(leg_dimuons_M_y4,"0.9 < |y| < 1.2","");
	leg_dimuons_M_y4->SetBorderSize(0);

	//TLegend* leg_dimuons_Eta = new TLegend(0.75,0.81,0.97,0.97);
	
	leg_dimuons_M_y4->AddEntry(h_dimuons_M_y4_trigger,"Data 2011","e1pSAME");
	leg_dimuons_M_y4->AddEntry(h_dimuons_M_y4,"MC Pythia 6","L");
	
	leg_dimuons_M_y4->Draw();
	
	//-------------------------------------------------------------------------

}//end program
