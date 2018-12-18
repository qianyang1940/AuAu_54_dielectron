#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TPythia6.h"
#include "TMCParticle.h"
#include <iostream>
#include <vector>
#include <string>
#include "histogram.h"
//#include "TPythia6/TPythia6.h"
//#include "TPythia6/TMCParticle6.h"

#define eMass 0.000511
#define charmMass 1.275
#define bottomMass 4.18

void runhfevent(Int_t irun=1, Int_t Nevt=1000000, Int_t iseed=789456, Bool_t debug=1){
	TStopwatch*   stopWatch = new TStopwatch();
	stopWatch->Start();
	gSystem->Load("libPhysics");
	//gSystem->Load("libEG.so");
	gSystem->Load("libEGPythia6.so");
	gSystem->Load("libPythia6.so");
	//gSystem->Load("TPythia6/TPythia6.so");

	//========================= initial trees ==========================
	eTree = new TTree("meTree","eTree");
	eTree->SetAutoSave(1000000);

	cout << "Initialize the eTree ... " << endl;

	eTree->Branch("EventId",&meTree.EventId,"EventId/I");

	eTree->Branch("Z0Pt",&meTree.Z0Pt,"Z0Pt/D");
	eTree->Branch("Z0Eta",&meTree.Z0Eta,"Z0Eta/D");
	eTree->Branch("Z0Phi",&meTree.Z0Phi,"Z0Phi/D");
	eTree->Branch("Z0Y",&meTree.Z0Y,"Z0Y/D");

	eTree->Branch("eePt",&meTree.eePt,"eePt/D");
	eTree->Branch("eeEta",&meTree.eeEta,"eeEta/D");
	eTree->Branch("eePhi",&meTree.eePhi,"eePhi/D");
	eTree->Branch("eeRapidity",&meTree.eeRapidity,"eeRapidity/D");
	eTree->Branch("eeM",&meTree.eeM,"eeM/D");
	eTree->Branch("ePosPt",&meTree.ePosPt,"ePosPt/D");
	eTree->Branch("ePosEta",&meTree.ePosEta,"ePosEta/D");
	eTree->Branch("ePosPhi",&meTree.ePosPhi,"ePosPhi/D");
	eTree->Branch("eNegPt",&meTree.eNegPt,"eNegPt/D");
	eTree->Branch("eNegEta",&meTree.eNegEta,"eNegEta/D");
	eTree->Branch("eNegPhi",&meTree.eNegPhi,"eNegPhi/D");

	//======================== initial parameters ==========================

	// Initialize Pythia
	TPythia6 *myPythia6 = new TPythia6();

	//myPythia6->SetMSEL(1); //call PYGIVE('msel = 1')   ! pp min. bias. 
	//myPythia6->SetMSEL(4);//c trigger
	//myPythia6->SetMSEL(5);//b trigger
	myPythia6->SetMSEL(11);//Z0 or gamma* trigger same as below
	//myPythia6->SetMSEL(0);//Z0 or gamma* trigger
	//myPythia6->SetMSUB(1,1);
	myPythia6->SetMSTP(43,1);//qqbar -> gamma* -> l+l-
	//myPythia6->SetMSTP(81,1);//call PYGIVE('mstp(81) = 1') ! multiple interaction on.
	//myPythia6->SetMSTP(82,4);//call PYGIVE('mstp(82) = 4') ! impact parameter choice.
	//myPythia6->SetMSTP(2, 2);//call PYGIVE('mstp(2) = 2') ! 2nd order running of alpha_s.
	//myPythia6->SetMSTP(33,3);//call PYGIVE('mstp(33) = 3') ! K-factors.

	//Double_t guass distribution
	//PARP(83)=0.2
	//PARP(84)=0.5
	//myPythia6->SetPARP(83,0.2);
	//myPythia6->SetPARP(84,0.5);

	//myPythia6->SetPARP(82,3.2);//call PYGIVE('parp(82) = 3.2') ! PT0 multiple distribution tail.
	myPythia6->SetCKIN(1,1.); //M
	//myPythia6->SetCKIN(3,1.0);//Set the high pT trigger 2GeV/c    
	//myPythia6->SetCKIN(7,-1.0);//y min
	//myPythia6->SetCKIN(8,1.0);//y max

	//tuned parameters:
	//myPythia6->SetMSTJ(11,3);//Peterson FF
	//myPythia6->SetPARJ(54,-1e-05);//mstj11 epsilon, make spectrum hard D=-0.05
	//myPythia6->SetPARJ(55,-1e-05);
	myPythia6->SetMSTP(33,1);//common kfactor stored in parp31
	myPythia6->SetMSTP(32,4);//Q^2 scale D:8
	myPythia6->SetMSTP(51,7);//CTEQ5L PDF
	//myPythia6->SetMSTP(51,4);//GRV-LO PDF
	myPythia6->SetPARP(31,1.8); //D:1.5
	myPythia6->SetPARP(91,1.5);//<kt>
	//myPythia6->SetPARP(93,10);//maximum pt broadening D:5
	//myPythia6->SetPMAS(4,1,charmMass);//c-quark mass D=1.5
	//myPythia6->SetPMAS(5,1,bottomMass);//b-quark mass  D=4.8
	//myPythia6->SetPARP(67,4);//mstp32*4 high pT tuned parameter

	//====================== particle decay mode ==========================

	//switch off non-electron decays
	for(Int_t i=162; i<=189; i++) myPythia6->SetMDME(i,1,0); 

	myPythia6->SetMDME(170,1,1);
	myPythia6->SetMDME(182,1,1);


	//============================ initial run =============================

	myPythia6->SetMRPY(1,iseed);
	myPythia6->Initialize("CMS","p","p",54);

	myPythia6->Pystat(2);

	//============================ run events ==============================

	Int_t nevent = 0;
	for(Int_t i=1;i<=Nevt;i++){
		myPythia6->GenerateEvent();

		if(debug){
			cout<<endl;
			cout<<"Working on "<<i<<"-th event ..."<<endl;
			myPythia6->Pylist(1);
		}
		else if(i<10) myPythia6->Pylist(1);

		TObjArray *particles = myPythia6->GetListOfParticles();
		Int_t nParticles = particles->GetEntries();
		if(i%10000==0) {
			cout<<"Working on "<<i<<"-th event ..."<<endl;
		}

		Bool_t isZ0   = kFALSE;
		Bool_t isePos = kFALSE;
		Bool_t iseNeg = kFALSE;

		memset(&meTree,0,sizeof(meTree));

		TMCParticle *mZ0 = (TMCParticle*)particles->At(9);

		if(mZ0->GetKF() == 23) isZ0 = kTRUE;
		Int_t Z0id  = mZ0->GetKF();
		Int_t Z0_ks = mZ0->GetKS();
		if(Z0id==23 && Z0_ks==11){
			Double_t Z0px = mZ0->GetPx();
			Double_t Z0py = mZ0->GetPy();
			Double_t Z0pz = mZ0->GetPz();
			Double_t Z0mass = mZ0->GetMass();
			TLorentzVector Z0FourMom(0,0,0,0);
			Z0FourMom.SetXYZM(Z0px,Z0py,Z0pz,Z0mass);
			meTree.Z0Pt  = Z0FourMom.Pt();
			meTree.Z0Eta = Z0FourMom.PseudoRapidity();
			meTree.Z0Phi = Z0FourMom.Phi();
			meTree.Z0Y   = Z0FourMom.Rapidity();
		}

		for(Int_t j=10; j<14; j++){
			TMCParticle *meTmp = (TMCParticle*)particles->At(j);

			if(meTmp->GetKF() == -11){
				isePos = kTRUE;
				Double_t epospx = meTmp->GetPx();
				Double_t epospy = meTmp->GetPy();
				Double_t epospz = meTmp->GetPz();
				TVector3 ePosMom(epospx,epospy,epospz);
				meTree.ePosPt  = ePosMom.Pt();
				meTree.ePosEta = ePosMom.Eta();
				meTree.ePosPhi = ePosMom.Phi();
			}

			if(meTmp->GetKF() == 11){
				iseNeg = kTRUE;
				Double_t enegpx = meTmp->GetPx();
				Double_t enegpy = meTmp->GetPy();
				Double_t enegpz = meTmp->GetPz();
				TVector3 eNegMom(enegpx,enegpy,enegpz);
				meTree.eNegPt  = eNegMom.Pt();
				meTree.eNegEta = eNegMom.Eta();
				meTree.eNegPhi = eNegMom.Phi();
			}
		}

		if(!isZ0 || !isePos || !iseNeg) continue;

		TLorentzVector ePosFourMom(0,0,0,0);
		TLorentzVector eNegFourMom(0,0,0,0);
		TLorentzVector eeFourMom(0,0,0,0);

		ePosFourMom.SetPtEtaPhiM(meTree.ePosPt,meTree.ePosEta,meTree.ePosPhi,eMass);
		eNegFourMom.SetPtEtaPhiM(meTree.eNegPt,meTree.eNegEta,meTree.eNegPhi,eMass);

		eeFourMom = ePosFourMom + eNegFourMom;

		meTree.EventId       = nevent;
		
		meTree.eePt          = eeFourMom.Pt();
		meTree.eeEta         = eeFourMom.Eta();
		meTree.eePhi         = eeFourMom.Phi();
		meTree.eeRapidity    = eeFourMom.Rapidity();
		meTree.eeM           = eeFourMom.M();

		eTree->Fill();

		nevent ++;

	}//event

	//============================= output ==============================

	char rootfilename[100];
	if(debug) sprintf(rootfilename,"test/pythiaevent_test.root");
	else      sprintf(rootfilename,"output/pythiaevent%d%d%d.root",irun/100,(irun/10)%10,irun%10);
	TFile* file = new TFile(rootfilename,"RECREATE");
	file->cd();
	eTree->Write();
	file->Close();

	myPythia6->Pystat(1);

	stopWatch->Stop();
	stopWatch->Print();
}
