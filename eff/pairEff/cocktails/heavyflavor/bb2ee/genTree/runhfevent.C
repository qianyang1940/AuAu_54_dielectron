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
using std::cout;
using std::endl;

#define eMass 0.000511
#define charmMass 1.275
#define bottomMass 4.18

void runhfevent(Int_t irun=1, Int_t Nevt=1000, Int_t iseed=789456, Bool_t debug=1){
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
	eTree->Branch("bPt",&meTree.bPt,"bPt/D");
	eTree->Branch("bEta",&meTree.bEta,"bEta/D");
	eTree->Branch("bPhi",&meTree.bPhi,"bPhi/D");
	eTree->Branch("bY",&meTree.bY,"bY/D");

	eTree->Branch("bbarPt",&meTree.bbarPt,"bbarPt/D");
	eTree->Branch("bbarEta",&meTree.bbarEta,"bbarEta/D");
	eTree->Branch("bbarPhi",&meTree.bbarPhi,"bbarPhi/D");
	eTree->Branch("bbarY",&meTree.bbarY,"bbarY/D");

	eTree->Branch("eePt",&meTree.eePt,"eePt/D");
	eTree->Branch("eeEta",&meTree.eeEta,"eeEta/D");
	eTree->Branch("eePhi",&meTree.eePhi,"eePhi/D");
	eTree->Branch("eeRapidity",&meTree.eeRapidity,"eeRapidity/D");
	eTree->Branch("eeM",&meTree.eeM,"eeM/D");
	eTree->Branch("ePosParentGID",&meTree.ePosParentGID,"ePosParentGID/I");
	eTree->Branch("ePosPt",&meTree.ePosPt,"ePosPt/D");
	eTree->Branch("ePosEta",&meTree.ePosEta,"ePosEta/D");
	eTree->Branch("ePosPhi",&meTree.ePosPhi,"ePosPhi/D");
	eTree->Branch("eNegParentGID",&meTree.eNegParentGID,"eNegParentGID/I");
	eTree->Branch("eNegPt",&meTree.eNegPt,"eNegPt/D");
	eTree->Branch("eNegEta",&meTree.eNegEta,"eNegEta/D");
	eTree->Branch("eNegPhi",&meTree.eNegPhi,"eNegPhi/D");

	//======================== initial parameters ==========================

	//Initialize Pythia
	TPythia6 *myPythia6 = new TPythia6();

	//myPythia6->SetMSEL(1); //call PYGIVE('msel = 1')   ! pp min. bias. 
	//myPythia6->SetMSEL(4);//c trigger
	myPythia6->SetMSEL(5);//b trigger
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
	//myPythia6->SetCKIN(3,5);//Set the high pT trigger 4GeV/c    

	//tuned parameters:
	//myPythia6->SetMSTJ(11,3);//Peterson FF
	//myPythia6->SetPARJ(54,-1e-05);//mstj11 epsilon, make spectrum hard D=-0.05
	//myPythia6->SetPARJ(55,-1e-05);
	myPythia6->SetMSTP(33,1);//common kfactor stored in parp31
	myPythia6->SetMSTP(32,4);//Q^2 scale D:8
	myPythia6->SetMSTP(51,7);//CTEQ5L PDF
	myPythia6->SetPARP(31,3.5); //D:1.5
	myPythia6->SetPARP(91,1.5);//<kt>
	//myPythia6->SetPARP(93,10);//maximum pt broadening D:5
	myPythia6->SetPMAS(4,1,charmMass);//c-quark mass D=1.5
	myPythia6->SetPMAS(5,1,bottomMass);//b-quark mass  D=4.8
	//myPythia6->SetPARP(67,4);//mstp32*4 high pT tuned parameter

	//====================== particle decay mode ==========================

	//switch off non-electron decays
	for(Int_t i=673; i<=862; i++) myPythia6->SetMDME(i,1,0); //c

	for(Int_t i=869; i<=898; i++) myPythia6->SetMDME(i,1,0); //B0
	for(Int_t i=914; i<=943; i++) myPythia6->SetMDME(i,1,0); //B+
	for(Int_t i=959; i<=991; i++) myPythia6->SetMDME(i,1,0); //Bs

	//============================ initial run =============================

	myPythia6->SetMRPY(1,iseed);
	myPythia6->Initialize("CMS","p","p",54);

	myPythia6->Pystat(2);

	//============================ run events ==============================

	Int_t nevent = 0;
	for(Int_t i = 1;i<=Nevt;i++){
		myPythia6->GenerateEvent();

		if(debug){
			cout<<endl;
			cout<<"Woring on "<<i<<"-th event ..."<<endl;
			myPythia6->Pylist(1);
		}
		else if(i<10) myPythia6->Pylist(1);

		TObjArray *particles = myPythia6->GetListOfParticles();
		Int_t nParticles = particles->GetEntries();

		if(i%10000==0){
			cout<<"Woring on "<<i<<"-th event ..."<<endl;
		}

		Int_t nstring = 0;
		for(Int_t l=0;l<nParticles;l++){
			TMCParticle *mParticle = (TMCParticle*)particles->At(l);
			if(!mParticle) continue;

			Int_t stringpid = mParticle->GetKF();

			if(stringpid != 92) continue; //string

			Int_t string_parentID = mParticle->GetParent()-1;
			TMCParticle *mStringParent = (TMCParticle*)particles->At(string_parentID);
			if(!mStringParent) continue;

			Int_t bottompid = mStringParent->GetKF();
			if(TMath::Abs(bottompid) != 5) continue; //string comes from b-bbar

			nstring ++;
		}

		if(nstring != 2) continue;

		memset(&meTree,0,sizeof(meTree));

		Int_t nePos = 0;
		Int_t neNeg = 0;
		for(Int_t l=0;l<nParticles;l++){
			TMCParticle *mElectron = (TMCParticle*)particles->At(l);
			if(!mElectron) continue;

			Int_t bid  = mElectron->GetKF();
			Int_t b_ks = mElectron->GetKS();

			if(bid==5&&b_ks==21){
				Double_t bpx = mElectron->GetPx();
				Double_t bpy = mElectron->GetPy();
				Double_t bpz = mElectron->GetPz();
				Double_t bmass = mElectron->GetMass();
				TLorentzVector bFourMom(0,0,0,0);
				bFourMom.SetXYZM(bpx,bpy,bpz,bottomMass);
				meTree.bPt  = bFourMom.Pt();
				meTree.bEta = bFourMom.PseudoRapidity();
				meTree.bPhi = bFourMom.Phi();
				meTree.bY   = bFourMom.Rapidity();
			}
			if(bid==-5&&b_ks==21){
				Double_t bpx = mElectron->GetPx();
				Double_t bpy = mElectron->GetPy();
				Double_t bpz = mElectron->GetPz();
				Double_t bmass = mElectron->GetMass();
				TLorentzVector bbarFourMom(0,0,0,0);
				bbarFourMom.SetXYZM(bpx,bpy,bpz,bottomMass);
				meTree.bbarPt  = bbarFourMom.Pt();
				meTree.bbarEta = bbarFourMom.PseudoRapidity();
				meTree.bbarPhi = bbarFourMom.Phi();
				meTree.bbarY   = bbarFourMom.Rapidity();
			}

			Int_t epid = mElectron->GetKF();

			if(TMath::Abs(epid) != 11) continue; //e+: -11,  e-: 11

			Int_t e_parentID = mElectron->GetParent()-1;
			TMCParticle *meParent = (TMCParticle*)particles->At(e_parentID);
			if(!meParent) continue;

			Int_t e_parentpid = TMath::Abs(meParent->GetKF());
			if( e_parentpid != 511 &&
				e_parentpid != 521 &&
				e_parentpid != 531)  continue;

			Double_t px = mElectron->GetPx();
			Double_t py = mElectron->GetPy();
			Double_t pz = mElectron->GetPz();

			TVector3 Mom(px,py,pz);
			Double_t pt  = Mom.Perp();
			Double_t eta = Mom.PseudoRapidity();
			Double_t phi = Mom.Phi();

			if(epid == -11) {
				nePos++;
				meTree.ePosParentGID = e_parentpid;
				meTree.ePosPt        = pt;
				meTree.ePosEta       = eta;
				meTree.ePosPhi       = phi;
			}

			if(epid == 11) {
				neNeg++;
				meTree.eNegParentGID = e_parentpid;
				meTree.eNegPt        = pt;
				meTree.eNegEta       = eta;
				meTree.eNegPhi       = phi;
			}

		}//particles

		//if(nePos!=1 || neNeg!=1) continue;

		TLorentzVector ePosFourMom(0,0,0,0);
		TLorentzVector eNegFourMom(0,0,0,0);
		TLorentzVector eeFourMom(0,0,0,0);

		ePosFourMom.SetPtEtaPhiM(meTree.ePosPt,meTree.ePosEta,meTree.ePosPhi,eMass);
		eNegFourMom.SetPtEtaPhiM(meTree.eNegPt,meTree.eNegEta,meTree.eNegPhi,eMass);
		eeFourMom = ePosFourMom + eNegFourMom;

		meTree.eePt          = eeFourMom.Pt();
		meTree.eeEta         = eeFourMom.Eta();
		meTree.eePhi         = eeFourMom.Phi();
		meTree.eeRapidity    = eeFourMom.Rapidity();
		meTree.eeM           = eeFourMom.M();

		meTree.EventId       = nevent;

		eTree->Fill();

		nevent ++;
	}//event

	//============================= output ==============================

	char rootfilename[100];
	if(debug){
		sprintf(rootfilename,"test/pythiaevent_test.root");
	}
	else{
		sprintf(rootfilename,"output/pythiaevent%d%d%d.root",irun/100,(irun/10)%10,irun%10);
	}
	TFile* file = new TFile(rootfilename,"RECREATE");
	file->cd();
	eTree->Write();
	file->Close();

	myPythia6->Pystat(1);

	stopWatch->Stop();
	stopWatch->Print();
}



