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

#include "EVENT.h"
#include "StRefMultCorr.h"
#include "../../../cuts.h"

using namespace std;
#endif

StRefMultCorr *refMultCorrUtil;
Int_t mCentrality;

Bool_t studySysUnCertainty = kFALSE;
//Bool_t studySysUnCertainty = kTRUE;
const Double_t tmpCut4Per = 0.01;
const Double_t tmpCut4Cen = 0.0025;

map<Int_t,Int_t> mBadRunId;

const Double_t PI = TMath::Pi();

const Float_t nSigmaCut = 0.6;

const Int_t mMaxElectrons = 1000;
Int_t current_nEPlus;
Int_t current_nEMinus;
TLorentzVector current_ePlus[mMaxElectrons];
Int_t current_ePlusTrkId[mMaxElectrons];
Bool_t current_phePlusTag[mMaxElectrons];
TLorentzVector current_eMinus[mMaxElectrons];
Int_t current_eMinusTrkId[mMaxElectrons];
Bool_t current_pheMinusTag[mMaxElectrons];

bool Init();
void bookHistograms();
bool passEvent(EVENT* event);
bool passTrack(EVENT* event, Int_t i);
void photonE(EVENT* event);
void writeHistograms(char* outFile);

//***** constrain the bad dedx calibration geometry *****
TF1 *funPosHi;
TF1 *funPosLow;
TF1 *funNegHi;
TF1 *funNegLow;
Float_t par[4][4];
Float_t parErr[4][4];

//*************** book histograms ******************
TH3F *hDenPiPlusTofEff;
TH3F *hNumPiPlusTofEff;
TH3F *hDenPiMinusTofEff;
TH3F *hNumPiMinusTofEff;
TH3F *hDenKPlusTofEff;
TH3F *hNumKPlusTofEff;
TH3F *hDenKMinusTofEff;
TH3F *hNumKMinusTofEff;
TH3F *hDenPPlusTofEff;
TH3F *hNumPPlusTofEff;
TH3F *hDenPMinusTofEff;
TH3F *hNumPMinusTofEff;
TH2F *hULMvsPt;
TH2F *hLPosMvsPt;
TH2F *hLNegMvsPt;
TH2F *hPEPlusBetavsP;
TH2F *hPEMinusBetavsP;
TH3F *hDenPEPlusTofEff;
TH3F *hNumPEPlusTofEff;
TH3F *hDenPEMinusTofEff;
TH3F *hNumPEMinusTofEff;

TH3F *hDenPiPlusTofEffCen[16];
TH3F *hNumPiPlusTofEffCen[16];
TH3F *hDenPiMinusTofEffCen[16];
TH3F *hNumPiMinusTofEffCen[16];
TH3F *hDenKPlusTofEffCen[16];
TH3F *hNumKPlusTofEffCen[16];
TH3F *hDenKMinusTofEffCen[16];
TH3F *hNumKMinusTofEffCen[16];
TH3F *hDenPPlusTofEffCen[16];
TH3F *hNumPPlusTofEffCen[16];
TH3F *hDenPMinusTofEffCen[16];
TH3F *hNumPMinusTofEffCen[16];
TH3F *hULMvsPtCen;
TH3F *hLPosMvsPtCen;
TH3F *hLNegMvsPtCen;
TH3F *hPEPlusBetavsPCen;
TH3F *hPEMinusBetavsPCen;
TH3F *hDenPEPlusTofEffCen[16];
TH3F *hNumPEPlusTofEffCen[16];
TH3F *hDenPEMinusTofEffCen[16];
TH3F *hNumPEMinusTofEffCen[16];

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
	TChain *chain = new TChain("mEvent");

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

	refMultCorrUtil = new StRefMultCorr("refmult");

	//+-------------+
	//| loop events |
	//+-------------+
	EVENT *event = new EVENT(chain);
	Int_t nEvts = chain->GetEntries();
	cout<<nEvts<<" events"<<endl;
	for(int i=0;i<nEvts;i++){

		if(i%(nEvts/10)==0) cout << "begin " << i << "th entry...." << endl;
		event->GetEntry(i);

		Int_t runId = event->mRunId;
		map<Int_t,Int_t>::iterator iter = mBadRunId.find(runId);
		if(iter != mBadRunId.end()) continue;

		if(!passEvent(event)) continue; 

		current_nEPlus = 0;
		current_nEMinus = 0;
		Int_t nTrks = event->mNTrks;
		for(int j=0;j<nTrks;j++) passTrack(event,j);

		photonE(event);
	}

	writeHistograms(outFile);
	delete chain;

	cout<<"end of program"<<endl;
	return 0;
}
//________________________________________________________________
bool passEvent(EVENT* event)
{
	bool eventflag=kFALSE;
	for(Int_t i=0;i<event->mNTriggers;i++){
		if(event->mTriggerId[i] == 400005) eventflag=kTRUE; //vpd-zdce-tac-protected
		if(event->mTriggerId[i] == 400015) eventflag=kTRUE; //vpd-zdce-tac-protected
		if(event->mTriggerId[i] == 400025) eventflag=kTRUE; //vpd-zdce-tac-protected
		if(event->mTriggerId[i] == 400035) eventflag=kTRUE; //vpd-zdce-tac-protected
	}
	if(!eventflag) return kFALSE;

	Float_t vx = event->mVertexX;
	Float_t vy = event->mVertexY;
	Float_t vz = event->mVertexZ;
	Float_t vr = sqrt(vx*vx+vy*vy);
	Float_t vpdVz = event->mVpdVz;
	Float_t vzDiff = vz - vpdVz;

	if(TMath::Abs(vx)<1.e-5 && TMath::Abs(vy)<1.e-5 && TMath::Abs(vz)<1.e-5) return kFALSE;
	if(vr>=mVrCut) return kFALSE;
	if(TMath::Abs(vz)>=mVzCut) return kFALSE;//vz should also be in the range listed in the parameters file to do the refMult correction
	if(TMath::Abs(vzDiff)>=mVzDiffCut) return kFALSE;

	Int_t runId = event->mRunId;
	Int_t zdcRate = event->mZDCRate;
	Int_t refMult = event->mRefMult;
	refMultCorrUtil->init(runId);
	refMultCorrUtil->initEvent(refMult,vz,zdcRate);
	mCentrality = refMultCorrUtil->getCentralityBin16();

	if(mCentrality<0 || mCentrality>15) return kFALSE;

	return kTRUE;
}
//______________________________________________________________
bool passTrack(EVENT* event, Int_t i)
{
	Char_t charge = event->mCharge[i];
	Char_t nHitsFit = event->mNHitsFit[i];
	Char_t nHitsPoss = event->mNHitsPoss[i];
	Char_t nHitsDedx = event->mNHitsDedx[i];
	Float_t ratio = 1.0*nHitsFit/nHitsPoss;
	Float_t dedx = event->mDedx[i]/1000.;
	Float_t nSigmaE = event->mNSigmaE[i]/1000.;
	Float_t nSigmaPi = event->mNSigmaPi[i]/1000.;
	Float_t nSigmaK = event->mNSigmaK[i]/1000.;
	Float_t nSigmaP = event->mNSigmaP[i]/1000.;
	Float_t dca = event->mDca[i]/1000.;
	Float_t pt = event->mPt[i];
	Float_t eta = event->mEta[i];
	Float_t phi = event->mPhi[i];
	Char_t  matchFlag2TOF = event->mTOFMatchFlag[i];
	Float_t beta2TOF = event->mBeta2TOF[i];
	Float_t localY2TOF = event->mTOFLocalY[i];

	TVector3 mom;
	mom.SetPtEtaPhi(pt,eta,phi);
	Float_t p = mom.Mag();

	if(pt<mTpcePtCut[0] || pt>mTpcePtCut[1]) return kFALSE;
	if(nHitsFit<mTpceNHitsFitCut) return kFALSE;
	if(ratio<mTpceNHitsFitRatioCut) return kFALSE;
	if(nHitsDedx<mTpceNHitsDedxCut) return kFALSE;
	if(dca>mTpceDcaCut) return kFALSE;
	if(TMath::Abs(eta)>mTpceEtaCut) return kFALSE;

	if(charge*pt>0. && eta>0.
			&& phi<=funPosHi->Eval(charge*pt)
			&& phi>=funPosLow->Eval(charge*pt)
	  ) return kFALSE;

	if(charge*pt<0. && eta>0.
			&& phi<=funNegHi->Eval(charge*pt)
			&& phi>=funNegLow->Eval(charge*pt)
	  ) return kFALSE;

	if(phi<-PI+PI/12) phi += 2*PI;

	if(charge>0){
		if(TMath::Abs(nSigmaPi)<nSigmaCut){
			hDenPiPlusTofEff->Fill(pt,eta,phi);
			hDenPiPlusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			if(beta2TOF>0.){
				hNumPiPlusTofEff->Fill(pt,eta,phi);
				hNumPiPlusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			}
		}

		if(TMath::Abs(nSigmaK)<nSigmaCut){
			hDenKPlusTofEff->Fill(pt,eta,phi);
			hDenKPlusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			if(beta2TOF>0.){
				hNumKPlusTofEff->Fill(pt,eta,phi);
				hNumKPlusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			}
		}

		if(TMath::Abs(nSigmaP)<nSigmaCut){
			hDenPPlusTofEff->Fill(pt,eta,phi);
			hDenPPlusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			if(beta2TOF>0.){
				hNumPPlusTofEff->Fill(pt,eta,phi);
				hNumPPlusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			}
		}
	}else if(charge<0){
		if(TMath::Abs(nSigmaPi)<nSigmaCut){
			hDenPiMinusTofEff->Fill(pt,eta,phi);
			hDenPiMinusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			if(beta2TOF>0.){
				hNumPiMinusTofEff->Fill(pt,eta,phi);
				hNumPiMinusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			}
		}

		if(TMath::Abs(nSigmaK)<nSigmaCut){
			hDenKMinusTofEff->Fill(pt,eta,phi);
			hDenKMinusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			if(beta2TOF>0.){
				hNumKMinusTofEff->Fill(pt,eta,phi);
				hNumKMinusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			}
		}

		if(TMath::Abs(nSigmaP)<nSigmaCut){
			hDenPMinusTofEff->Fill(pt,eta,phi);
			hDenPMinusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			if(beta2TOF>0.){
				hNumPMinusTofEff->Fill(pt,eta,phi);
				hNumPMinusTofEffCen[mCentrality]->Fill(pt,eta,phi);
			}
		}
	}

	Float_t mTpceNSigmaECutLow;
	Float_t mTpceNSigmaECutHi;
	if(mCentrality<=3) { //60-80%
		mTpceNSigmaECutLow = mTpceNSigmaECut[0];  // old: -0.75
		mTpceNSigmaECutHi  = mTpceNSigmaECut[1];  // old: 1.25
	}
	else { //0-60%
		mTpceNSigmaECutLow = -1*nSigmaCut;
		mTpceNSigmaECutHi  = nSigmaCut;
	}
	if(nSigmaE<mTpceNSigmaECutLow+mNSigmaEShift || nSigmaE>mTpceNSigmaECutHi+mNSigmaEShift) return kFALSE;

	if(charge>0){
		current_ePlus[current_nEPlus].SetPtEtaPhiM(pt,eta,phi,Melectron);
		current_ePlusTrkId[current_nEPlus] = i;
		current_nEPlus++;
	}else if(charge<0){
		current_eMinus[current_nEMinus].SetPtEtaPhiM(pt,eta,phi,Melectron);
		current_eMinusTrkId[current_nEMinus] = i;
		current_nEMinus++;
	}

	return kTRUE;
}
//____________________________________________________________
void photonE(EVENT *event)
{
	for(Int_t i=0;i<current_nEPlus;i++) current_phePlusTag[i] = kFALSE;
	for(Int_t i=0;i<current_nEMinus;i++) current_pheMinusTag[i] = kFALSE;

	TLorentzVector pair(0,0,0,0);

	for(Int_t i=0;i<current_nEPlus;i++){
		for(Int_t j=0;j<current_nEMinus;j++){
			pair = current_ePlus[i] + current_eMinus[j];
			hULMvsPt->Fill(pair.Pt(),pair.M());
			hULMvsPtCen->Fill(mCentrality,pair.Pt(),pair.M());

			if(!studySysUnCertainty) {  // default selection criteria
				if(mCentrality<=3) { //60-80%
					if(pair.M()<mPheMassCutWTof){
						current_phePlusTag[i] = kTRUE;
						current_pheMinusTag[j] = kTRUE;
					}
				}
				else {  //0-60%
					if(pair.M()<mPheMassCutWoTof){
						current_phePlusTag[i] = kTRUE;
						current_pheMinusTag[j] = kTRUE;
					}
				}
			}
			else { // for systematic uncertainty study
				if(mCentrality<=3){ //60-80%
					if(pair.M()<tmpCut4Per){
						current_phePlusTag[i] = kTRUE;
						current_pheMinusTag[j] = kTRUE;
					}
				}
				else {  //0-60%
					if(pair.M()<tmpCut4Cen){
						current_phePlusTag[i] = kTRUE;
						current_pheMinusTag[j] = kTRUE;
					}
				}
			}
		}
	}

	for(Int_t i=0;i<current_nEPlus;i++){
		for(Int_t j=i+1;j<current_nEPlus;j++){
			pair = current_ePlus[i] + current_ePlus[j];
			hLPosMvsPt->Fill(pair.Pt(),pair.M());
			hLPosMvsPtCen->Fill(mCentrality,pair.Pt(),pair.M());

			if(!studySysUnCertainty) { // default selection criteria
				if(mCentrality<=3) { //60-80%
					if(pair.M()<mPheMassCutWTof){
						current_phePlusTag[i] = kFALSE;
						current_phePlusTag[j] = kFALSE;
					}
				}
				else { //0-60%
					if(pair.M()<mPheMassCutWoTof){
						current_phePlusTag[i] = kFALSE;
						current_phePlusTag[j] = kFALSE;
					}
				}
			}
			else { // for systematic uncertainty study
				if(mCentrality<=3) { //60-80%
					if(pair.M()<tmpCut4Per){
						current_phePlusTag[i] = kFALSE;
						current_phePlusTag[j] = kFALSE;
					}
				}
				else { //0-60%
					if(pair.M()<tmpCut4Cen){
						current_phePlusTag[i] = kFALSE;
						current_phePlusTag[j] = kFALSE;
					}
				}
			}
		}
	}

	for(Int_t i=0;i<current_nEMinus;i++){
		for(Int_t j=i+1;j<current_nEMinus;j++){
			pair = current_eMinus[i] + current_eMinus[j];
			hLNegMvsPt->Fill(pair.Pt(),pair.M());
			hLNegMvsPtCen->Fill(mCentrality,pair.Pt(),pair.M());

			if(!studySysUnCertainty) {  // default selection criteria
				if(mCentrality<=3) { //60-80%
					if(pair.M()<mPheMassCutWTof){
						current_pheMinusTag[i] = kFALSE;
						current_pheMinusTag[j] = kFALSE;
					}
				}
				else { //0-60%
					if(pair.M()<mPheMassCutWoTof){
						current_pheMinusTag[i] = kFALSE;
						current_pheMinusTag[j] = kFALSE;
					}
				}
			}
			else {  // for systematic uncertainty study
				if(mCentrality<=3) { //60-80%
					if(pair.M()<tmpCut4Per){
						current_pheMinusTag[i] = kFALSE;
						current_pheMinusTag[j] = kFALSE;
					}
				}
				else { //0-60%
					if(pair.M()<tmpCut4Cen){
						current_pheMinusTag[i] = kFALSE;
						current_pheMinusTag[j] = kFALSE;
					}
				}
			}
		}
	}

	for(Int_t i=0;i<current_nEPlus;i++){
		if(!current_phePlusTag[i]) continue;

		Float_t pt = current_ePlus[i].Pt();
		Float_t eta = current_ePlus[i].Eta();
		Float_t phi = current_ePlus[i].Phi();
		Float_t p = current_ePlus[i].P();
		Int_t trkId = current_ePlusTrkId[i];
		Float_t nSigmaE = event->mNSigmaE[trkId]/1000.;
		Float_t matchFlag2TOF = event->mTOFMatchFlag[trkId];
		Float_t beta2TOF = event->mBeta2TOF[trkId];
		Float_t localY2TOF = event->mTOFLocalY[trkId];

		if(phi<-PI+PI/12) phi += 2*PI;

		hDenPEPlusTofEff->Fill(pt,eta,phi);
		hDenPEPlusTofEffCen[mCentrality]->Fill(pt,eta,phi);
		if(beta2TOF>0.){
			hPEPlusBetavsP->Fill(p,1./beta2TOF);
			hPEPlusBetavsPCen->Fill(mCentrality,p,1./beta2TOF);
			hNumPEPlusTofEff->Fill(pt,eta,phi);
			hNumPEPlusTofEffCen[mCentrality]->Fill(pt,eta,phi);
		}
	}

	for(Int_t i=0;i<current_nEMinus;i++){
		if(!current_pheMinusTag[i]) continue;

		Float_t pt = current_eMinus[i].Pt();
		Float_t eta = current_eMinus[i].Eta();
		Float_t phi = current_eMinus[i].Phi();
		Float_t p = current_eMinus[i].P();
		Int_t trkId = current_eMinusTrkId[i];
		Float_t nSigmaE = event->mNSigmaE[trkId]/1000.;
		Float_t matchFlag2TOF = event->mTOFMatchFlag[trkId];
		Float_t beta2TOF = event->mBeta2TOF[trkId];
		Float_t localY2TOF = event->mTOFLocalY[trkId];

		if(phi<-PI+PI/12) phi += 2*PI;

		hDenPEMinusTofEff->Fill(pt,eta,phi);
		hDenPEMinusTofEffCen[mCentrality]->Fill(pt,eta,phi);
		if(beta2TOF>0.){
			hPEMinusBetavsP->Fill(p,1./beta2TOF);
			hPEMinusBetavsPCen->Fill(mCentrality,p,1./beta2TOF);
			hNumPEMinusTofEff->Fill(pt,eta,phi);
			hNumPEMinusTofEffCen[mCentrality]->Fill(pt,eta,phi);

		}
	}
}
//____________________________________________________________
void bookHistograms()
{
	hDenPiPlusTofEff = new TH3F("hDenPiPlusTofEff","hDenPiPlusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hNumPiPlusTofEff = new TH3F("hNumPiPlusTofEff","hNumPiPlusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hDenPiMinusTofEff = new TH3F("hDenPiMinusTofEff","hDenPiMinusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hNumPiMinusTofEff = new TH3F("hNumPiMinusTofEff","hNumPiMinusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hDenKPlusTofEff = new TH3F("hDenKPlusTofEff","hDenKPlusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hNumKPlusTofEff = new TH3F("hNumKPlusTofEff","hNumKPlusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hDenKMinusTofEff = new TH3F("hDenKMinusTofEff","hDenKMinusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hNumKMinusTofEff = new TH3F("hNumKMinusTofEff","hNumKMinusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hDenPPlusTofEff = new TH3F("hDenPPlusTofEff","hDenPPlusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hNumPPlusTofEff = new TH3F("hNumPPlusTofEff","hNumPPlusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hDenPMinusTofEff = new TH3F("hDenPMinusTofEff","hDenPMinusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hNumPMinusTofEff = new TH3F("hNumPMinusTofEff","hNumPMinusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hULMvsPt = new TH2F("hULMvsPt","hULMvsPt;p_{T} (GeV/c); M_{ee} (GeV/c^{2})",200,0,10,3000,0,0.3);
	hLPosMvsPt = new TH2F("hLPosMvsPt","hLPosMvsPt;p_{T} (GeV/c); M_{ee} (GeV/c^{2})",200,0,10,3000,0,0.3);
	hLNegMvsPt = new TH2F("hLNegMvsPt","hLNegMvsPt;p_{T} (GeV/c); M_{ee} (GeV/c^{2})",200,0,10,3000,0,0.3);
	hPEPlusBetavsP = new TH2F("hPEPlusBetavsP","hPEPlusBetavsP;p (GeV/c); 1/#beta",200,0,10,200,0.9,1.1);
	hPEMinusBetavsP = new TH2F("hPEMinusBetavsP","hPEMinusBetavsP;p (GeV/c); 1/#beta",200,0,10,200,0.9,1.1);
	hDenPEPlusTofEff = new TH3F("hDenPEPlusTofEff","hDenPEPlusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hNumPEPlusTofEff = new TH3F("hNumPEPlusTofEff","hNumPEPlusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hDenPEMinusTofEff = new TH3F("hDenPEMinusTofEff","hDenPEMinusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	hNumPEMinusTofEff = new TH3F("hNumPEMinusTofEff","hNumPEMinusTofEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);

	hULMvsPtCen = new TH3F("hULMvsPtCen","hULMvsPtCen;Centrality; p_{T} (GeV/c); M_{ee} (GeV/c^{2})",16,0,16,200,0,10,3000,0,0.3);
	hLPosMvsPtCen = new TH3F("hLPosMvsPtCen","hLPosMvsPtCen;Centrality; p_{T} (GeV/c); M_{ee} (GeV/c^{2})",16,0,16,200,0,10,3000,0,0.3);
	hLNegMvsPtCen = new TH3F("hLNegMvsPtCen","hLNegMvsPtCen;Centrality; p_{T} (GeV/c); M_{ee} (GeV/c^{2})",16,0,16,200,0,10,3000,0,0.3);
	hPEPlusBetavsPCen = new TH3F("hPEPlusBetavsPCen","hPEPlusBetavsPCen;Centrality; p (GeV/c); 1/#beta",16,0,16,200,0,10,200,0.9,1.1);
	hPEMinusBetavsPCen = new TH3F("hPEMinusBetavsPCen","hPEMinusBetavsPCen;Centrality; p (GeV/c); 1/#beta",16,0,16,200,0,10,200,0.9,1.1);

	for(Int_t i=0;i<16;i++){
		hDenPiPlusTofEffCen[i] = new TH3F(Form("hDenPiPlusTofEffCenBin%d",i),"hDenPiPlusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hNumPiPlusTofEffCen[i] = new TH3F(Form("hNumPiPlusTofEffCenBin%d",i),"hNumPiPlusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hDenPiMinusTofEffCen[i] = new TH3F(Form("hDenPiMinusTofEffCenBin%d",i),"hDenPiMinusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hNumPiMinusTofEffCen[i] = new TH3F(Form("hNumPiMinusTofEffCenBin%d",i),"hNumPiMinusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hDenKPlusTofEffCen[i] = new TH3F(Form("hDenKPlusTofEffCenBin%d",i),"hDenKPlusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hNumKPlusTofEffCen[i] = new TH3F(Form("hNumKPlusTofEffCenBin%d",i),"hNumKPlusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hDenKMinusTofEffCen[i] = new TH3F(Form("hDenKMinusTofEffCenBin%d",i),"hDenKMinusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hNumKMinusTofEffCen[i] = new TH3F(Form("hNumKMinusTofEffCenBin%d",i),"hNumKMinusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hDenPPlusTofEffCen[i] = new TH3F(Form("hDenPPlusTofEffCenBin%d",i),"hDenPPlusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hNumPPlusTofEffCen[i] = new TH3F(Form("hNumPPlusTofEffCenBin%d",i),"hNumPPlusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hDenPMinusTofEffCen[i] = new TH3F(Form("hDenPMinusTofEffCenBin%d",i),"hDenPMinusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hNumPMinusTofEffCen[i] = new TH3F(Form("hNumPMinusTofEffCenBin%d",i),"hNumPMinusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hDenPEPlusTofEffCen[i] = new TH3F(Form("hDenPEPlusTofEffCenBin%d",i),"hDenPEPlusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hNumPEPlusTofEffCen[i] = new TH3F(Form("hNumPEPlusTofEffCenBin%d",i),"hNumPEPlusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hDenPEMinusTofEffCen[i] = new TH3F(Form("hDenPEMinusTofEffCenBin%d",i),"hDenPEMinusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
		hNumPEMinusTofEffCen[i] = new TH3F(Form("hNumPEMinusTofEffCenBin%d",i),"hNumPEMinusTofEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI+PI/12,PI+PI/12);
	}
}
//____________________________________________________________
void writeHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.tofeff.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	TFile *mFile = new TFile(buf,"recreate");
	mFile->cd();

	for(Int_t i=0;i<16;i++){
		hDenPiPlusTofEffCen[i]->Write();
		hNumPiPlusTofEffCen[i]->Write();
		hDenPiMinusTofEffCen[i]->Write();
		hNumPiMinusTofEffCen[i]->Write();
		hDenKPlusTofEffCen[i]->Write();
		hNumKPlusTofEffCen[i]->Write();
		hDenKMinusTofEffCen[i]->Write();
		hNumKMinusTofEffCen[i]->Write();
		hDenPPlusTofEffCen[i]->Write();
		hNumPPlusTofEffCen[i]->Write();
		hDenPMinusTofEffCen[i]->Write();
		hNumPMinusTofEffCen[i]->Write();
		hDenPEPlusTofEffCen[i]->Write();
		hNumPEPlusTofEffCen[i]->Write();
		hDenPEMinusTofEffCen[i]->Write();
		hNumPEMinusTofEffCen[i]->Write();
	};
	hULMvsPtCen->Write();
	hLPosMvsPtCen->Write();
	hLNegMvsPtCen->Write();
	hPEPlusBetavsPCen->Write();
	hPEMinusBetavsPCen->Write();

	hDenPiPlusTofEff->Write();
	hNumPiPlusTofEff->Write();
	hDenPiMinusTofEff->Write();
	hNumPiMinusTofEff->Write();
	hDenKPlusTofEff->Write();
	hNumKPlusTofEff->Write();
	hDenKMinusTofEff->Write();
	hNumKMinusTofEff->Write();
	hDenPPlusTofEff->Write();
	hNumPPlusTofEff->Write();
	hDenPMinusTofEff->Write();
	hNumPMinusTofEff->Write();
	hULMvsPt->Write();
	hLPosMvsPt->Write();
	hLNegMvsPt->Write();
	hPEPlusBetavsP->Write();
	hPEMinusBetavsP->Write();
	hDenPEPlusTofEff->Write();
	hNumPEPlusTofEff->Write();
	hDenPEMinusTofEff->Write();
	hNumPEMinusTofEff->Write();
}
//____________________________________________________________
bool Init()
{
	cout<<endl;

	ifstream indata;

	indata.open("../../../produceRunList/runList/mBadRunList.dat");
	mBadRunId.clear();
	if(indata.is_open()){
		cout<<"read in the bad run number list ...";
		Int_t id;
		while(indata>>id){
			mBadRunId[id] = id;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the bad run number list !!!"<<endl;
		return kFALSE;
	}
	indata.close();
	for(map<Int_t,Int_t>::iterator iter=mBadRunId.begin();iter!=mBadRunId.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;
	cout<<endl;

	indata.open("/star/u/syang/run12/uu/minibias/badTPCSectorGeo/fitPars.dat");
	if(indata.is_open()){                                                      
		cout<<"read in bad TPC sector geometry parameters ...";

		Int_t idx = 0;
		Double_t tmp0,tmp1;
		while(indata>>tmp0>>tmp1){
			par[idx/4][idx%4] = tmp0;
			parErr[idx/4][idx%4] = tmp1;
			idx++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the bad TPC sector geometry parameters !!!"<<endl;
		return kFALSE;
	}
	indata.close();
	for(Int_t i=0;i<4;i++){
		for(Int_t j=0;j<4;j++){
			cout<<par[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
	funPosHi = new TF1("funPosHi","[0]*exp(-([1]/x)**[2])+[3]",0.1,60);
	funPosHi->SetParameters(par[0][0],par[0][1],par[0][2],par[0][3]);
	funPosLow = new TF1("funPosLow","[0]*exp(-([1]/x)**[2])+[3]",0.1,60);
	funPosLow->SetParameters(par[1][0],par[1][1],par[1][2],par[1][3]);
	funNegHi = new TF1("funNegHi","[0]*exp(-([1]/x)**[2])+[3]",-60,-0.1);
	funNegHi->SetParameters(par[2][0],par[2][1],par[2][2],par[2][3]);
	funNegLow = new TF1("funNegLow","[0]*exp(-([1]/x)**[2])+[3]",-60,-0.1);
	funNegLow->SetParameters(par[3][0],par[3][1],par[3][2],par[3][3]);

	cout<<"Initialization DONE !!!"<<endl;
	cout<<endl;

	return kTRUE;
}
