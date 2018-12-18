#include "/star/u/syang/Macro/headers.h"
#include "/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/MesonConstant.h"
#include "/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/FunUtil.h"
#include "/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/EffUtil.h"

#include "meTree.h"

const Double_t YCut = 1.;
const Double_t PtCut = 0.2;
const Double_t EtaCut = 1.;

Char_t  inFile[1024];
Char_t  outFile[1024];

Int_t   mCenIdx=0;
Bool_t  mUseScaleEff = kTRUE;

Double_t GetBR(Int_t id);

int main(int argc, char** argv)
{

	if(argc==1){
		sprintf(inFile,"test.list");
		sprintf(outFile,"test.root");
	}
	else if(argc==5){
		if (strcmp(argv[1],"080")==0) {
			mCenIdx = 0;
		}
		else if (strcmp(argv[1],"010")==0) {
			mCenIdx = 1;
		}
		else if (strcmp(argv[1],"1040")==0) {
			mCenIdx = 2;
		}
		else if (strcmp(argv[1],"4080")==0) {
			mCenIdx = 3;
		}
		else if (strcmp(argv[1],"4060")==0) {
			mCenIdx = 4;
		}
		else if (strcmp(argv[1],"6080")==0) {
			mCenIdx = 5;
		}
		else if (strcmp(argv[1],"6070")==0) {
			mCenIdx = 6;
		}
		else if (strcmp(argv[1],"7080")==0) {
			mCenIdx = 7;
		}
		else{
			cout<<"Centrality string is wrong !"<<endl;
			return -1;
		}

		system("mkdir -p output");
		system(Form("mkdir -p output/Cen%s",argv[1]));

		if(strcmp(argv[2],"Scale")== 0) mUseScaleEff = kTRUE;
		else if(strcmp(argv[2], "Individual")==0) mUseScaleEff = kFALSE; 
		else{
			cout<<"The input efficiency option is wrong! Should be \"Scale\" OR \"Individual\"!"<<endl;
			return -1;
		}

		sprintf(inFile,"%s",argv[3]);
		sprintf(outFile,"output/Cen%s/%s.histo.root",argv[1],argv[4]);
	}
	else{
		cout<<"The input parameters are wrong !"<<endl;
		return -1;
	}

	//----------------------------------
	// Open files and add to the chain
	//----------------------------------
	TChain *chain = new TChain("meTree");

	Int_t ifile=0;
	Char_t filename[1024];

	ifstream *inputStream = new ifstream;
	inputStream->open(inFile);
	if (!(inputStream)) {
		printf("can not open list file\n");
		return 0;
	}

	while(inputStream->good()){
		inputStream->getline(filename,512);

		if(!inputStream->good()) continue;

		TFile *ftmp = new TFile(filename);
		if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
			printf("file %s error in opening!!!\n",filename);
		} else {
			printf("read in file %s\n",filename);
		}

		if(ifile%1==0) cout<<"ifile="<<ifile<<" "<<filename<<endl;
		chain->Add(filename);

		ftmp->Delete();

		ifile++;
	}

	printf("%d files read in\n",ifile);
	delete inputStream;

	TH1D *hnEvts = new TH1D("hnEvts","hnEvts",5,0,5);
	hnEvts->GetXaxis()->SetBinLabel(1, "inclusiveEvts");
	hnEvts->GetXaxis()->SetBinLabel(2, "twoStringEvts");
	hnEvts->SetLabelSize(0.025,"X");
	hnEvts->SetLabelFont(62,"X");
	hnEvts->LabelsOption("d","X");

	TH2D *hMCFullYMvsPtcc2ee = new TH2D("hMCFullYMvsPtcc2ee","hMCFullYMvsPtcc2ee",1500,0.,15.,1200,0.,6.); //Without any cut
	TH2D *hMCAcc0MvsPtcc2ee = new TH2D("hMCAcc0MvsPtcc2ee","hMCAcc0MvsPtcc2ee",1500,0.,15.,1200,0.,6.); //|Yee|<=1.
	TH2D *hMCAcc1MvsPtcc2ee = new TH2D("hMCAcc1MvsPtcc2ee","hMCAcc1MvsPtcc2ee",1500,0.,15.,1200,0.,6.); //STAR Acceptance
	TH2D *hMCAcc2MvsPtcc2ee = new TH2D("hMCAcc2MvsPtcc2ee","hMCAcc2MvsPtcc2ee",1500,0.,15.,1200,0.,6.); //PHENIX Acceptance
	TH2D *hRCAcc1MvsPt3Dcc2ee = new TH2D("hRCAcc1MvsPt3Dcc2ee","hRCAcc1MvsPt3Dcc2ee",1500,0.,15.,1200,0.,6.); //3D-eff correction
	hMCFullYMvsPtcc2ee->Sumw2();
	hMCAcc0MvsPtcc2ee->Sumw2();
	hMCAcc1MvsPtcc2ee->Sumw2();
	hMCAcc2MvsPtcc2ee->Sumw2();
	hRCAcc1MvsPt3Dcc2ee->Sumw2();

	//****** Efficiency for single track ******
	InitializeEffFun(mCenIdx);

	//------------ additional resolution for RC tracks --------
	//Double_t mPtSmearPar[2] = {3.5e-3, 8.0e-03}; //untuned pt reso. - directly from run12 UU embedding fit
	Double_t mPtSmearPar[2] = {6.0e-3, 8.3e-03}; //tuned pt reso. - from published AuAu@200GeV PRC paper
	TF1 *funSmearPt = new TF1("funSmearPt","sqrt([0]*[0]*x*x+[1]*[1])",0.,10.);
	funSmearPt->SetParameters(mPtSmearPar);

	TF1 *momShape = new TF1("momShape",CrystalBall2,-1.,1.,7);
	//momShape->SetParameters(1., -1e-3, 0.01, 1.16, 1.84, 2.18, 1.92);//-3.9e-4
	momShape->SetParameters(1., -1e-3, 0.01, 1.29, 1.75, 2.92, 1.84);
	momShape->SetNpx(10000);

	meTree *t = new meTree(chain);

	Int_t nevents=(Int_t)chain->GetEntries();
	cout<<"== total entries : " <<nevents<<endl;
	for(int i=0; i<nevents; i++) {
		if(i%(nevents/10) == 0) cout<<"processing "<<i<<" events..."<<endl;

		t->GetEntry(i);

		hnEvts->Fill(0.5);

		Int_t ncString    = t->ncString;
		Int_t ncbarString = t->ncbarString;
		if(ncString!=1 || ncbarString!=1) continue; //this event must contain 2 strings form ccbar
		hnEvts->Fill(1.5);

		TLorentzVector smePos(0,0,0,0);
		TLorentzVector smeNeg(0,0,0,0);

		Int_t ietatpc, iphitpc, ipttpc;
		Int_t ietatof, iphitof, ipttof;
		Double_t ePosEff = 0;
		Double_t eNegEff = 0;

		Int_t nePos = t->nePos;
		Int_t neNeg = t->neNeg;
		for(Int_t ii=0; ii<nePos; ii++){
			Double_t ePosPt  = t->ePosPt[ii];
			Double_t ePosEta = t->ePosEta[ii];
			Double_t ePosPhi = t->ePosPhi[ii];
			Double_t ePosBR  = GetBR(t->ePosParentGID[ii]);

			Double_t ePosMomSig  = funSmearPt->Eval(ePosPt);
			Double_t ePosRcPt    = ePosPt*(1.+momShape->GetRandom()*ePosMomSig/0.01);
			smePos.SetPtEtaPhiM(ePosRcPt,ePosEta,ePosPhi,Masselectron);

			if(mUseScaleEff){
				tpcPtEtaPhi2Bin(0, ePosRcPt,ePosEta,ePosPhi, &ipttpc,&ietatpc,&iphitpc); //MB binning
				tofPtEtaPhi2Bin(0, ePosRcPt,ePosEta,ePosPhi, &ipttof,&ietatof,&iphitof); //MB binning
				ePosEff = hMBEff_Tpc_Pos[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*fTpcEffRatioToMB->Eval(ePosRcPt)*hMBEff_Tof_Pos[ietatof][iphitof]->GetBinContent(ipttof+1)*fMBEPiEff_Tof->Eval(ePosRcPt)*fTofEffRatioToMB->Eval(ePosRcPt)*funnSigEffP->Eval(smePos.P())*funBetaEffP->Eval(smePos.P())*funndEdxEffPt->Eval(ePosRcPt);
			}
			else{
				tpcPtEtaPhi2Bin(mCenIdx, ePosRcPt,ePosEta,ePosPhi, &ipttpc,&ietatpc,&iphitpc);
				tofPtEtaPhi2Bin(mCenIdx, ePosRcPt,ePosEta,ePosPhi, &ipttof,&ietatof,&iphitof);
				ePosEff = hEff_Tpc_Pos[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*hEff_Tof_Pos[ietatof][iphitof]->GetBinContent(ipttof+1)*fEPiEff_Tof->Eval(ePosRcPt)*funnSigEffP->Eval(smePos.P())*funBetaEffP->Eval(smePos.P())*funndEdxEffPt->Eval(ePosRcPt);
			}

			for(Int_t jj=0; jj<neNeg; jj++){
				Double_t eNegPt  = t->eNegPt[jj];
				Double_t eNegEta = t->eNegEta[jj];
				Double_t eNegPhi = t->eNegPhi[jj];
				Double_t eNegBR  = GetBR(t->eNegParentGID[jj]);

				Double_t eNegMomSig = funSmearPt->Eval(eNegPt);
				Double_t eNegRcPt   = eNegPt*(1.+momShape->GetRandom()*eNegMomSig/0.01);
				smeNeg.SetPtEtaPhiM(eNegRcPt,eNegEta,eNegPhi,Masselectron);

				if(mUseScaleEff){
					tpcPtEtaPhi2Bin(0, eNegRcPt,eNegEta,eNegPhi, &ipttpc,&ietatpc,&iphitpc); //MB binning
					tofPtEtaPhi2Bin(0, eNegRcPt,eNegEta,eNegPhi, &ipttof,&ietatof,&iphitof); //MB binning
					eNegEff = hMBEff_Tpc_Neg[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*fTpcEffRatioToMB->Eval(eNegRcPt)*hMBEff_Tof_Neg[ietatof][iphitof]->GetBinContent(ipttof+1)*fMBEPiEff_Tof->Eval(eNegRcPt)*fTofEffRatioToMB->Eval(eNegRcPt)*funnSigEffP->Eval(smeNeg.P())*funBetaEffP->Eval(smeNeg.P())*funndEdxEffPt->Eval(eNegRcPt);
				}
				else{
					tpcPtEtaPhi2Bin(mCenIdx, eNegRcPt,eNegEta,eNegPhi, &ipttpc,&ietatpc,&iphitpc);
					tofPtEtaPhi2Bin(mCenIdx, eNegRcPt,eNegEta,eNegPhi, &ipttof,&ietatof,&iphitof);
					eNegEff = hEff_Tpc_Neg[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*hEff_Tof_Neg[ietatof][iphitof]->GetBinContent(ipttof+1)*fEPiEff_Tof->Eval(eNegRcPt)*funnSigEffP->Eval(smeNeg.P())*funBetaEffP->Eval(smeNeg.P())*funndEdxEffPt->Eval(eNegRcPt);
				}

				TLorentzVector eePair(0,0,0,0);
				eePair = smePos + smeNeg;
				Double_t eePt = eePair.Pt();
				Double_t eeM  = eePair.M();

				hMCFullYMvsPtcc2ee->Fill(eePt, eeM, ePosBR*eNegBR);

				Int_t kPosPassFlag = PHENIXFilter(1,smePos.Px(),smePos.Py(),smePos.Pz(),0);
				Int_t kNegPassFlag = PHENIXFilter(-1,smeNeg.Px(),smeNeg.Py(),smeNeg.Pz(),0);
				if(kPosPassFlag>=0 && kNegPassFlag>=0) {
					hMCAcc2MvsPtcc2ee->Fill(eePt, eeM, ePosBR*eNegBR);
				}

				if(fabs(eePair.Rapidity())>YCut) continue;
				hMCAcc0MvsPtcc2ee->Fill(eePt, eeM, ePosBR*eNegBR);

				if(ePosRcPt<PtCut || eNegRcPt<PtCut ) continue;
				if(fabs(ePosEta)>EtaCut || fabs(eNegEta)>EtaCut) continue;
				hMCAcc1MvsPtcc2ee->Fill(eePt, eeM, ePosBR*eNegBR);
				hRCAcc1MvsPt3Dcc2ee->Fill(eePt, eeM, ePosBR*eNegBR*ePosEff*eNegEff);
			}
		}

	}

	TFile *fout = new TFile(outFile,"recreate");
	fout->cd();

	hnEvts->Write();
	hMCFullYMvsPtcc2ee->Write();
	hMCAcc0MvsPtcc2ee->Write();
	hMCAcc1MvsPtcc2ee->Write();
	hMCAcc2MvsPtcc2ee->Write();
	hRCAcc1MvsPt3Dcc2ee->Write();

	fout->Close();

	cout << "End!!" << endl;

	return 0;
}

Double_t GetBR(Int_t id){
	Double_t BR =0.;

	////PDG 2006
	//if(abs(id)==411)  BR = 0.172; // D+
	//if(abs(id)==421)  BR = 0.0671; // D0;
	//if(abs(id)==431)  BR = 0.08;  // Ds
	//if(abs(id)==4122) BR = 0.045; // Lc

	//PDG 2014
	if(abs(id)==411)  BR = 0.1607; // D+
	if(abs(id)==421)  BR = 0.0649; // D0;
	if(abs(id)==431)  BR = 0.065;  // Ds
	if(abs(id)==4122) BR = 0.045; // Lc

	return BR;
}
