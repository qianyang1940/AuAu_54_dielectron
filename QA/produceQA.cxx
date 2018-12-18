#include <cstdio> 
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "sys/types.h"
#include "dirent.h"
#include "math.h"
#include "string.h"

#ifndef __CINT__  
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "miniDst.h"
#include "StRefMultCorr.h"
#include "cuts.h"

using namespace std;
#endif


StRefMultCorr *refMultCorrUtil;

Int_t runIndex;
Int_t randomId;
map<Int_t,Int_t> mTotalRunId;

bool Init();
void bookHistograms();
bool passEvent(miniDst* event);
bool passTrack(miniDst* event, Int_t i);
void writeHistograms(char* outFile);

//define histograms
//inclusive QA
TH1D *hnEvts;
TH2F *hVyvsVx;
TH2F *hVPDVzvsTPCVz;
TH2F *hDeltaZvsTPCVz;
TH2F *hDeltaZvsVPDVz;
TH2F *hTPCVzvsRefMult;
TH2F *hVPDVzvsRefMult;
TH2F *hVrvsRefMult;
TH2F *hDeltaZvsRefMult;
TH2F *hZDCXvsRefMult;
TH2F *hBBCXvsRefMult;

TH2F *hVxvsRunIndex;
TH2F *hVyvsRunIndex;
TH2F *hTPCVzvsRunIndex;
TH2F *hVPDVzvsRunIndex;
TH2F *hDeltaZvsRunIndex;
TH2F *hRefMultvsRunIndex;
TH2F *hZDCXvsRunIndex;
TH2F *hBBCXvsRunIndex;

// track level QA
TH1F *hnHitsFit;
TH1F *hnHitsPoss;
TH1F *hnHitsDedx;
TH2F *hDcavsPt;
TH2F *hBetavsP;
TH2F *hBetavsEta;
TH2F *hBetavsPhi;
TH2F *hDedxvsP;
TH2F *hDedxvsEta;
TH2F *hDedxvsPhi;
TH2F *hnSigmaEvsP;
TH2F *hnSigmaEvsEta;
TH2F *hnSigmaEvsPhi;
TH2F *hMSquarevsP;
TH2F *hEtavsPhi;
TH2F *hEEtavsPhi;
TH2F *hEtavsPt;
TH2F *hEEtavsPt;
TH2F *hPhivsPt;
TH2F *hEPhivsPt;
//run by run QA
TH2F *hPtvsRunIndex;
TH2F *hEtavsRunIndex;
TH2F *hPhivsRunIndex;
TH2F *hDcavsRunIndex;
TH2F *hnHitsFitvsRunIndex;
TH2F *hnHitsPossvsRunIndex;
TH2F *hnHitsDedxvsRunIndex;
TH2F *hdEdxvsRunIndex;
TH2F *hBetavsRunIndex;
TH2F *hnSigmaEvsRunIndex;
TH2F *hEBetaDiffvsRunIndex;
TH2F *hKBetaDiffvsRunIndex;

int main(int argc, char** argv)
{
	if(argc!=1&&argc!=3) return -1;

	TString inFile="test.list";
	char outFile[1024];
	sprintf(outFile,"test/test");
	if(argc==3){
		inFile = argv[1];
		sprintf(outFile,"%s",argv[2]);
	}

	//+---------------------------------+
	//| open files and add to the chain |
	//+---------------------------------+
	TChain *chain = new TChain("miniDst");

	Int_t ifile=0;
	char filename[512];
	ifstream *inputStream = new ifstream;
	inputStream->open(inFile.Data());
	if (!(inputStream)) {
		printf("can not open list file\n");
		return 0;
	}
	for(;inputStream->good();){
		inputStream->getline(filename,512);
		if(inputStream->good()) {
			TFile *ftmp = new TFile(filename);
			if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
				cout<<"something wrong"<<endl;
			} else {
				cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
				chain->Add(filename);
				ifile++;
			}
			delete ftmp;
		}
	}
	delete inputStream;

	//intialization
	bookHistograms();
	if( !Init() ){
		cout<<"The initialization is failed !!!"<<endl;
		return 0;
	}
	//refMultCorrUtil = new StRefMultCorr("refmult");

	//+-------------+
	//| loop events |
	//+-------------+
	miniDst *event = new miniDst(chain);
	Int_t nEvts = chain->GetEntries();
	cout<<nEvts<<" events"<<endl;
	for(int i=0;i<nEvts;i++){

		//initialize the struct

		if(i%(nEvts/10)==0) cout << "begin " << i << "th entry...." << endl;
		event->GetEntry(i);

		Int_t runId = event->mRunId;
		map<Int_t,Int_t>::iterator iter = mTotalRunId.find(runId);
		if(iter != mTotalRunId.end())
			runIndex = iter->second;
		else{
			runIndex = -1;
			cout<<"Can not find the runNumber in the runNumber list"<<endl;
		}

		if(runIndex<0) continue;

		hnEvts->Fill(0);
		if(!passEvent(event)) continue; 
		hnEvts->Fill(1);

		Int_t nTrks = event->mNTrks;
		for(int j=0;j<nTrks;j++) passTrack(event,j);
	}

	writeHistograms(outFile);
	delete chain;

	cout<<"end of program"<<endl;
	return 0;
}
//________________________________________________________________
bool passEvent(miniDst* event)
{
	bool eventflag=kFALSE;
	for(Int_t i=0;i<event->mNTrigs;i++){
		if(event->mTrigId[i] == 580001) eventflag = kTRUE; 
		if(event->mTrigId[i] == 580011) eventflag = kTRUE; 
		if(event->mTrigId[i] == 580021) eventflag = kTRUE; 
	}
	if(!eventflag) return kFALSE;

	//Int_t runId = event->mRunId;
	Int_t zdcRate =event->mZDCRate;
	Int_t bbcRate =event->mBBCRate;
	Int_t refMult = event->mRefMult;
	Float_t vx = event->mVertexX;
	Float_t vy = event->mVertexY;
	Float_t vz = event->mVertexZ;
	Float_t vpdVz = event->mVpdVz;
	Float_t vzDiff = vz - vpdVz;
	Float_t vr = sqrt(vx*vx + vy*vy);

	//already have event cuts in pico production

	hVyvsVx->Fill(vx,vy);
	hVPDVzvsTPCVz->Fill(vz,vpdVz);
	hDeltaZvsTPCVz->Fill(vz,vzDiff);
	hDeltaZvsVPDVz->Fill(vpdVz,vzDiff);
	hTPCVzvsRefMult->Fill(refMult,vz);
	hVPDVzvsRefMult->Fill(refMult,vpdVz);
	hVrvsRefMult->Fill(refMult,vr);
	hDeltaZvsRefMult->Fill(refMult,vzDiff);
	hZDCXvsRefMult->Fill(refMult,zdcRate/1000.);
	hBBCXvsRefMult->Fill(refMult,bbcRate/1000.);

	hVxvsRunIndex->Fill(runIndex,vx); 
	hVyvsRunIndex->Fill(runIndex,vy); 
	hTPCVzvsRunIndex->Fill(runIndex,vz);
	hVPDVzvsRunIndex->Fill(runIndex,vpdVz);
	hDeltaZvsRunIndex->Fill(runIndex,vzDiff);
	hRefMultvsRunIndex->Fill(runIndex,refMult);
	hZDCXvsRunIndex->Fill(runIndex,zdcRate/1000.);
	hBBCXvsRunIndex->Fill(runIndex,bbcRate/1000.);

//	if(TMath::Abs(vz)<mVzCut){//make sure vz is in the range listed in the parameters file
//		refMultCorrUtil->init(runId);
//		refMultCorrUtil->initEvent(refMult,vz,zdcRate);
//		Double_t refMultCor = refMultCorrUtil->getRefMultCorr();//if you want to call the getRefMultCorr() with no argument, it must be called after initEvent()
//	}

	return kTRUE;
}
//______________________________________________________________
bool passTrack(miniDst* event, Int_t i)
{
	Char_t charge = event->mCharge[i];
	Char_t nHitsFit = event->mNHitsFit[i];
	Char_t nHitsDedx = event->mNHitsDedx[i];
	Char_t nHitsPoss = event->mNHitsPoss[i];
	Float_t dedx = event->mDedx[i];
	Float_t nSigmaE = event->mNSigmaE[i];
	Float_t dca = event->mDca[i];
	Float_t pt = event->mPt[i];
	Float_t eta = event->mEta[i];
	Float_t phi = event->mPhi[i];
	Float_t beta2TOF = event->mBeta2TOF[i];
	TVector3 mom;
	mom.SetPtEtaPhi(pt,eta,phi);
	Float_t p = mom.Mag();

	//already applied some track quality cuts in pico production

	hnHitsFit->Fill(charge*nHitsFit);
	hnHitsPoss->Fill(charge*nHitsPoss);
	hnHitsDedx->Fill(charge*nHitsDedx);
	hEtavsPhi->Fill(phi,eta);
	hEtavsPt->Fill(charge*pt,eta);
	hPhivsPt->Fill(charge*pt,phi);
	hDcavsPt->Fill(charge*pt, dca);
	if(beta2TOF>0.){
		hBetavsP->Fill(charge*p,1./beta2TOF);
		hBetavsEta->Fill(eta,1./beta2TOF);
		hBetavsPhi->Fill(phi,1./beta2TOF);
	}
	hDedxvsP->Fill(charge*p,dedx);
	hDedxvsEta->Fill(eta,dedx);
	hDedxvsPhi->Fill(phi,dedx);
	hnSigmaEvsP->Fill(charge*p,nSigmaE);
	hnSigmaEvsEta->Fill(eta,nSigmaE);
	hnSigmaEvsPhi->Fill(phi,nSigmaE);

	hPtvsRunIndex->Fill(runIndex,pt);
	hEtavsRunIndex->Fill(runIndex,eta);
	hPhivsRunIndex->Fill(runIndex,phi);
	hDcavsRunIndex->Fill(runIndex,dca);
	hnHitsFitvsRunIndex->Fill(runIndex,nHitsFit);
	hnHitsPossvsRunIndex->Fill(runIndex,nHitsPoss);
	hnHitsDedxvsRunIndex->Fill(runIndex,nHitsDedx);
	hdEdxvsRunIndex->Fill(runIndex,dedx);
	hnSigmaEvsRunIndex->Fill(runIndex,nSigmaE);
	if(beta2TOF>0.) hBetavsRunIndex->Fill(runIndex,1./beta2TOF);

	Float_t msquare = -999.;
	if(beta2TOF>0.) msquare = pow(p,2)*(1-pow(beta2TOF,2))/pow(beta2TOF,2);
	if(msquare>-990.) hMSquarevsP->Fill(p,msquare);

	Float_t expBeta2TOF = -999;

	Float_t mTpceNSigmaECutLow;
	if(p<1.){
		mTpceNSigmaECutLow = (mTpceNSigmaECut[0]+2)/(1.-mTpcePtCut[0])*(p-mTpcePtCut[0]) - 2;
	}else{
		mTpceNSigmaECutLow = mTpceNSigmaECut[0];
	}

	if(beta2TOF>0.
			&& TMath::Abs(1.-1./beta2TOF)<=mTpceBeta2TOFCut 
			&& nSigmaE>=mTpceNSigmaECutLow && nSigmaE<=mTpceNSigmaECut[1]){
		hEEtavsPhi->Fill(phi,eta);
		hEEtavsPt->Fill(charge*pt,eta);
		hEPhivsPt->Fill(charge*pt,phi);
		expBeta2TOF = p/sqrt(pow(Melectron,2)+pow(p,2));
		hEBetaDiffvsRunIndex->Fill(runIndex,1./beta2TOF-1./expBeta2TOF);
	}

	return kTRUE;
}
//____________________________________________________________
void bookHistograms()
{
	//inclusive QA
	hnEvts = new TH1D("hnEvts","hnEvts",5,0,5);

	hVyvsVx = new TH2F("hVyvsVx","hVyvsVx;Vx (cm); Vy (cm)",150,-1.5,1.5,150,-1.5,1.5);
	hVPDVzvsTPCVz = new TH2F("hVPDVzvsTPCVz","hVPDVzvsTPCVz; TPC Vz (cm); VPD Vz (cm)",400,-200,200,400,-200,200);
	hDeltaZvsTPCVz = new TH2F("hDeltaZvsTPCVz","hDeltaZvsTPCVz; TPC Vz (cm); Vz_{TPC} - Vz_{VPD} (cm)",400,-200,200,400,-20,20);
	hDeltaZvsVPDVz = new TH2F("hDeltaZvsVPDVz","hDeltaZvsVPDVz; VPD Vz (cm); Vz_{TPC} - Vz_{VPD} (cm)",400,-200,200,400,-20,20);
	hTPCVzvsRefMult = new TH2F("hTPCVzvsRefMult","hTPCVzvsRefMult; refMult; TPC Vz (cm)",200,0,1000,400,-200,200);
	hVPDVzvsRefMult = new TH2F("hVPDVzvsRefMult","hVPDVzvsRefMult; refMult; VPD Vz (cm)",200,0,1000,400,-200,200);
	hVrvsRefMult = new TH2F("hVrvsRefMult","hVrvsRefMult; refMult; Vr (cm)",200,0,1000,200,0,2);
	hDeltaZvsRefMult = new TH2F("hDeltaZvsRefMult","hDeltaZvsRefMult; refMult; Vz_{TPC} - Vz_{VPD} (cm)",200,0,1000,200,-10,10);
	hZDCXvsRefMult = new TH2F("hZDCXvsRefMult","hZDCXvsRefMult; refMult; zdcRate (kHz)",1000,0,1000,20,0,20);
	hBBCXvsRefMult = new TH2F("hBBCXvsRefMult","hBBCXvsRefMult; refMult; bbcRate (kHz)",1000,0,1000,2000,0,200);


	hVxvsRunIndex = new TH2F("hVxvsRunIndex","hVxvsRunIndex;runIndex; Vx (cm)",mTotalRun,0,mTotalRun,150,-1.5,1.5);
	hVyvsRunIndex = new TH2F("hVyvsRunIndex","hVyvsRunIndex;runIndex; Vy (cm)",mTotalRun,0,mTotalRun,150,-1.5,1.5);
	hTPCVzvsRunIndex = new TH2F("hTPCVzvsRunIndex","hTPCVzvsRunIndex;runIndex;TPC Vz (cm)",mTotalRun,0,mTotalRun,400,-200,200);
	hVPDVzvsRunIndex = new TH2F("hVPDVzvsRunIndex","hVPDVzvsRunIndex;runIndex;VPD Vz (cm)",mTotalRun,0,mTotalRun,400,-200,200);
	hDeltaZvsRunIndex = new TH2F("hDeltaZvsRunIndex","hDeltaZvsRunIndex;runIndex; Vz_{TPC} - Vz_{VPD} (cm)",mTotalRun,0,mTotalRun,400,-20,20);
	hRefMultvsRunIndex = new TH2F("hRefMultvsRunIndex","hRefMultvsRunIndex;runIndex; refMult",mTotalRun,0,mTotalRun,200,0,1000);
	hZDCXvsRunIndex = new TH2F("hZDCXvsRunIndex","hZDCXvsRunIndex;runIndex;zdcRate (KHz)",mTotalRun,0,mTotalRun,20,0,10);
	hBBCXvsRunIndex = new TH2F("hBBCXvsRunIndex","hBBCXvsRunIndex;runIndex;bbcRate (KHz)",mTotalRun,0,mTotalRun,2000,0,200);

	hnHitsFit = new TH1F("hnHitsFit","hnHitsFit;nHitsFit;Counts",100,-50,50);
	hnHitsPoss = new TH1F("hnHitsPoss","hnHitsPoss;nHitsPoss;Counts",100,-50,50);
	hnHitsDedx = new TH1F("hnHitsDedx","hnHitsDedx;nHitsDedx;Counts",100,-50,50);
	hDcavsPt = new TH2F("hDcavsPt","hDcavsPt;q*p_{T} (GeV/c);dca (cm)",200,-10,10,300,-1.e-6,3.-1.e-6);
	hBetavsP = new TH2F("hBetavsP","hBetavsP;q*p (GeV/c);1/#beta",200,-10,10,500,0,5);
	hBetavsEta = new TH2F("hBetavsEta","hBetavsEta;#eta;1/#beta",200,-1,1,500,0,5);
	hBetavsPhi = new TH2F("hBetavsPhi","hBetavsPhi;#phi;1/#beta",360,-TMath::Pi(),TMath::Pi(),500,0,5);
	hDedxvsP = new TH2F("hDedxvsP","hDedxvsP;q*p (GeV/c);dEdx (KeV/cm)",200,-10,10,300,-1.e-6,30-1.e-6);
	hDedxvsEta = new TH2F("hDedxvsEta","hDedxvsEta;#eta;dEdx (KeV/cm)",200,-1,1,300,-1.e-6,30-1.e-6);
	hDedxvsPhi = new TH2F("hDedxvsPhi","hDedxvsPhi;#phi;dEdx (KeV/cm)",360,-TMath::Pi(),TMath::Pi(),300,-1.e-6,30-1.e-6);
	hnSigmaEvsP = new TH2F("hnSigmaEvsP","hnSigmaEvsP;q*p (GeV/c);n#sigma_{e}",200,-10,10,600,-30-1.e-6,30-1.e-6);
	hnSigmaEvsEta = new TH2F("hnSigmaEvsEta","hnSigmaEvsEta;#eta;n#sigma_{e}",200,-1,1,600,-30-1.e-6,30-1.e-6);
	hnSigmaEvsPhi = new TH2F("hnSigmaEvsPhi","hnSigmaEvsPhi;#phi;n#sigma_{e}",360,-TMath::Pi(),TMath::Pi(),600,-30-1.e-6,30-1.e-6);
	hMSquarevsP = new TH2F("hMSquarevsP","hMSquarevsP; p (GeV/c);m^{2} ( (GeV/c^{2})^{2} )",100,0,10,1200,-0.2,1);
	hEtavsPhi = new TH2F("hEtavsPhi","hEtavsPhi;#phi;#eta",360,-TMath::Pi(),TMath::Pi(),200,-1,1);
	hEEtavsPhi = new TH2F("hEEtavsPhi","hEEtavsPhi;#phi;#eta",360,-TMath::Pi(),TMath::Pi(),200,-1,1);
	hEtavsPt = new TH2F("hEtavsPt","hEtavsPt; q*p_{T} (GeV/c); #eta",200,-10,10,200,-1,1);
	hEEtavsPt = new TH2F("hEEtavsPt","hEEtavsPt; q*p_{T} (GeV/c); #eta",200,-10,10,200,-1,1);
	hPhivsPt = new TH2F("hPhivsPt","hPhivsPt; q*p_{T} (GeV/c); #phi",2000,-10,10,1080,-TMath::Pi(),TMath::Pi());
	hEPhivsPt = new TH2F("hEPhivsPt","hEPhivsPt; q*p_{T} (GeV/c); #phi",2000,-10,10,1080,-TMath::Pi(),TMath::Pi());
	//run by run QA
	hPtvsRunIndex = new TH2F("hPtvsRunIndex","hPtvsRunIndex;runIndex;p_{T} (GeV/c)",mTotalRun,0,mTotalRun,100,0,10);
	hEtavsRunIndex = new TH2F("hEtavsRunIndex","hEtavsRunIndex;runIndex;#eta",mTotalRun,0,mTotalRun,300,-1.5,1.5);
	hPhivsRunIndex = new TH2F("hPhivsRunIndex","hPhivsRunIndex;runIndex;#phi",mTotalRun,0,mTotalRun,360,-TMath::Pi(),TMath::Pi());
	hDcavsRunIndex = new TH2F("hDcavsRunIndex","hDcavsRunIndex;runIndex;dca (cm)",mTotalRun,0,mTotalRun,300,-1.e-6,3-1.e-6);
	hnHitsFitvsRunIndex = new TH2F("hnHitsFitvsRunIndex","hnHitsFitvsRunIndex;runIndex;nHitsFit",mTotalRun,0, mTotalRun,50,0,50);
	hnHitsPossvsRunIndex = new TH2F("hnHitsPossvsRunIndex","hnHitsPossvsRunIndex;runIndex;nHitsPoss",mTotalRun,0, mTotalRun,50,0,50);
	hnHitsDedxvsRunIndex = new TH2F("hnHitsDedxvsRunIndex","hnHitsDedxvsRunIndex;runIndex;nHitsDedx",mTotalRun,0, mTotalRun,50,0,50);
	hdEdxvsRunIndex = new TH2F("hdEdxvsRunIndex","hdEdxvsRunIndex;runIndex;dE/dx (KeV/cm)",mTotalRun,0,mTotalRun,300,-1.e-6,30-1.e-6);
	hBetavsRunIndex = new TH2F("hBetavsRunIndex","hBetavsRunIndex;runIndex;1/#beta",mTotalRun,0,mTotalRun,500,0,5);
	hnSigmaEvsRunIndex = new TH2F("hnSigmaEvsRunIndex","hnSigmaEvsRunIndex;runIndex;n#sigma_{e}",mTotalRun,0,mTotalRun,600,-30-1.e-6,30-1.e-6);
	hEBetaDiffvsRunIndex = new TH2F("hEBetaDiffvsRunIndex","hEBetaDiffvsRunIndex;runIndex;1/#beta - 1/#beta_{exp}",mTotalRun,0,mTotalRun,100,-0.05,0.05);
	hKBetaDiffvsRunIndex = new TH2F("hKBetaDiffvsRunIndex","hKBetaDiffvsRunIndex;runIndex;1/#beta - 1/#beta_{exp}",mTotalRun,0,mTotalRun,100,-0.05,0.05);
}
//=======================================================================================
void writeHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.QAhisto.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	TFile *mFile = new TFile(buf,"recreate");
	mFile->cd();

	//inclusive QA
	hnEvts->Write();
	hVyvsVx->Write();
	hVPDVzvsTPCVz->Write();
	hDeltaZvsTPCVz->Write();
	hDeltaZvsVPDVz->Write();
	hTPCVzvsRefMult->Write();
	hVPDVzvsRefMult->Write();
	hVrvsRefMult->Write();
	hDeltaZvsRefMult->Write();
	hZDCXvsRefMult->Write();
	hBBCXvsRefMult->Write();

	hVxvsRunIndex->Write();
	hVyvsRunIndex->Write();
	hTPCVzvsRunIndex->Write();
	hVPDVzvsRunIndex->Write();
	hDeltaZvsRunIndex->Write();
	hRefMultvsRunIndex->Write();
	hZDCXvsRunIndex->Write();
	hBBCXvsRunIndex->Write();

	hnHitsFit->Write();
	hnHitsPoss->Write();
	hnHitsDedx->Write();
	hDcavsPt->Write();
	hBetavsP->Write();
	hBetavsEta->Write();
	hBetavsPhi->Write();
	hDedxvsP->Write();
	hDedxvsEta->Write();
	hDedxvsPhi->Write();
	hnSigmaEvsP->Write();
	hnSigmaEvsEta->Write();
	hnSigmaEvsPhi->Write();
	hMSquarevsP->Write();
	hEtavsPhi->Write();
	hEEtavsPhi->Write();
	hEtavsPt->Write();
	hEEtavsPt->Write();
	hPhivsPt->Write();
	hEPhivsPt->Write();

	//run by run QA
	hPtvsRunIndex->Write();
	hEtavsRunIndex->Write();
	hPhivsRunIndex->Write();
	hDcavsRunIndex->Write();
	hnHitsFitvsRunIndex->Write();
	hnHitsPossvsRunIndex->Write();
	hnHitsDedxvsRunIndex->Write();
	hdEdxvsRunIndex->Write();
	hBetavsRunIndex->Write();
	hnSigmaEvsRunIndex->Write();
	hEBetaDiffvsRunIndex->Write();
	hKBetaDiffvsRunIndex->Write();
}
//==============================================================================================
bool Init()
{
	ifstream indata;
	indata.open("/star/u/tc88qy/AuAu/run17/54GeV/QA/runList/output_all/GetRun/runList/mTotalRunList.dat");
	mTotalRunId.clear();
	if(indata.is_open()){
		cout<<"read in total run number list and recode run number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata>>oldId){
			mTotalRunId[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list !!!"<<endl;
		return kFALSE;
	}
	indata.close();
	for(map<Int_t,Int_t>::iterator iter=mTotalRunId.begin();iter!=mTotalRunId.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;
	cout<<endl;

	return kTRUE;
}
