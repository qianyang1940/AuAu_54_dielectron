#include "/star/u/syang/Macro/headers.h"
#include "/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/MesonConstant.h"
#include "/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/FunUtil.h"
#include "/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/EffUtil.h"

const Double_t YCut = 1.;
const Double_t PtCut = 0.2;
const Double_t EtaCut = 1.;

Char_t  inFile[1024];
Char_t  outFile[1024];

Int_t   mCenIdx=0;
Bool_t  mUseScaleEff = kTRUE;

int main(int argc, char** argv)
{

	if(argc==1){
		sprintf(inFile,"test.list");
		sprintf(outFile,"test.root");
	}
	else if(argc==4){

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
		sprintf(outFile,"output/Cen%s/dy2eeAll.root",argv[1]);
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

	TTree *t = chain;

	//TFile *f = new TFile("/star/u/huangbc/run9dilep/simulation/cocktail/pythia_hf/genTree/dy_fully/dy2ee_fully.root");
	//TTree *t = (TTree *)f->Get("meTree");

	Int_t   EventId;
	Double_t Z0Pt;
	Double_t Z0Eta;
	Double_t Z0Phi;
	Double_t Z0Y;
	Double_t eePt;
	Double_t eeEta;
	Double_t eePhi;
	Double_t eeRapidity;
	Double_t eeM;
	Double_t ePosPt;
	Double_t ePosEta;
	Double_t ePosPhi;
	Double_t eNegPt;
	Double_t eNegEta;
	Double_t eNegPhi;

	t->SetBranchAddress("EventId", &EventId);
	t->SetBranchAddress("Z0Pt", &Z0Pt);
	t->SetBranchAddress("Z0Eta", &Z0Eta);
	t->SetBranchAddress("Z0Phi", &Z0Phi);
	t->SetBranchAddress("Z0Y", &Z0Y);
	t->SetBranchAddress("eePt", &eePt);
	t->SetBranchAddress("eeEta", &eeEta);      
	t->SetBranchAddress("eePhi", &eePhi);      
	t->SetBranchAddress("eeRapidity", &eeRapidity);      
	t->SetBranchAddress("eeM", &eeM);      
	t->SetBranchAddress("ePosPt", &ePosPt);
	t->SetBranchAddress("ePosEta", &ePosEta);
	t->SetBranchAddress("ePosPhi", &ePosPhi);
	t->SetBranchAddress("eNegPt", &eNegPt);
	t->SetBranchAddress("eNegEta", &eNegEta);
	t->SetBranchAddress("eNegPhi", &eNegPhi);

	TH1D *hnParent = new TH1D("nParent","number of Parent",5,0,5);
	TH2D *hMCFullYMvsPtdy2ee = new TH2D("hMCFullYMvsPtdy2ee","hMCFullYMvsPtdy2ee",1500,0.,15.,1200,0.,6.);
	TH2D *hMCAcc0MvsPtdy2ee = new TH2D("hMCAcc0MvsPtdy2ee","hMCAcc0MvsPtdy2ee",1500,0.,15.,1200,0.,6.); //|Yee|<=1.
	TH2D *hMCAcc1MvsPtdy2ee = new TH2D("hMCAcc1MvsPtdy2ee","hMCAcc1MvsPtdy2ee",1500,0.,15.,1200,0.,6.); //STAR Acceptance
	TH2D *hMCAcc2MvsPtdy2ee = new TH2D("hMCAcc2MvsPtdy2ee","hMCAcc2MvsPtdy2ee",1500,0.,15.,1200,0.,6.); //PHENIX Acceptance
	TH2D *hRCAcc1MvsPt3Ddy2ee = new TH2D("hRCAcc1MvsPt3Ddy2ee","hRCAcc1MvsPt3Ddy2ee",1500,0.,15.,1200,0.,6.);

	hRCAcc1MvsPt3Ddy2ee->Sumw2();

	//================================ VPDMB setup =================================

	//****** Efficiency for single track ******
	InitializeEffFun(mCenIdx);

	//------------ additional resolution for RC tracks --------
	//Double_t mPtSmearPar[3] = {3.5e-3, 8.0e-03, 0.000511}; //untuned pt reso. - directly from run12 UU embedding fit
	Double_t mPtSmearPar[3] = {6.0e-3, 8.3e-03, 0.000511}; //tuned pt reso. - from published AuAu@200GeV PRC paper
	TF1 *funSmearPt = new TF1("funSmearPt","sqrt([0]*[0]*x*x+[1]*[1])",0.,10.);
	funSmearPt->SetParameters(mPtSmearPar);

	TF1 *momShape = new TF1("momShape",CrystalBall2,-1.,1.,7);
	//momShape->SetParameters(1., -1e-3, 0.01, 1.16, 1.84, 2.18, 1.92);//-3.9e-4
	momShape->SetParameters(1., -1e-3, 0.01, 1.29, 1.75, 2.92, 1.84);
	momShape->SetNpx(10000);

	Int_t nevents = t->GetEntries();
	cout<<"== total entries : " <<nevents<<endl;

	for(int i=0; i<nevents; i++) {

		if(i%(nevents/10) == 0) cout<<"processing "<<i<<" events..."<<endl;

		t->GetEntry(i);

		hnParent->Fill(0);

		TVector3 ePos(0,0,0);
		ePos.SetPtEtaPhi(ePosPt,ePosEta,ePosPhi);

		TVector3 eNeg(0,0,0);
		eNeg.SetPtEtaPhi(eNegPt,eNegEta,eNegPhi);

		Double_t emsig1 = funSmearPt->Eval(ePosPt);
		Double_t emsig2 = funSmearPt->Eval(eNegPt);

		double smrcpt1 = ePosPt*(1.+momShape->GetRandom()*emsig1/0.01);
		double smrcpt2 = eNegPt*(1.+momShape->GetRandom()*emsig2/0.01);

		TLorentzVector smePos(0,0,0,0);
		smePos.SetPtEtaPhiM(smrcpt1,ePosEta,ePosPhi,Masselectron);
		TLorentzVector smeNeg(0,0,0,0);
		smeNeg.SetPtEtaPhiM(smrcpt2,eNegEta,eNegPhi,Masselectron);

		TLorentzVector eePair(0,0,0,0);
		eePair = smePos + smeNeg;
		eePt = eePair.Pt();
		eeM  = eePair.M();

		float posEff = 0.;
		float negEff = 0.;

		int ietatpc,iphitpc,ipttpc;
		int ietatof,iphitof,ipttof;

		if(mUseScaleEff){
			tpcPtEtaPhi2Bin(0, smrcpt1,ePosEta,ePosPhi, &ipttpc,&ietatpc,&iphitpc); // MB binning 
			tofPtEtaPhi2Bin(0, smrcpt1,ePosEta,ePosPhi, &ipttof,&ietatof,&iphitof); // MB binning
			posEff = hMBEff_Tpc_Pos[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*fTpcEffRatioToMB->Eval(smrcpt1)*hMBEff_Tof_Pos[ietatof][iphitof]->GetBinContent(ipttof+1)*fMBEPiEff_Tof->Eval(smrcpt1)*fTofEffRatioToMB->Eval(smrcpt1)*funnSigEffP->Eval(smePos.P())*funBetaEffP->Eval(smePos.P())*funndEdxEffPt->Eval(smrcpt1);

			tofPtEtaPhi2Bin(0, smrcpt2,eNegEta,eNegPhi, &ipttof,&ietatof,&iphitof); // MB binning
			tpcPtEtaPhi2Bin(0, smrcpt2,eNegEta,eNegPhi, &ipttpc,&ietatpc,&iphitpc); // MB binning
			negEff = hMBEff_Tpc_Neg[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*fTpcEffRatioToMB->Eval(smrcpt2)*hMBEff_Tof_Neg[ietatof][iphitof]->GetBinContent(ipttof+1)*fMBEPiEff_Tof->Eval(smrcpt2)*fTofEffRatioToMB->Eval(smrcpt2)*funnSigEffP->Eval(smeNeg.P())*funBetaEffP->Eval(smeNeg.P())*funndEdxEffPt->Eval(smrcpt2);
		}
		else{
			tpcPtEtaPhi2Bin(mCenIdx, smrcpt1,ePosEta,ePosPhi, &ipttpc,&ietatpc,&iphitpc);
			tofPtEtaPhi2Bin(mCenIdx, smrcpt1,ePosEta,ePosPhi, &ipttof,&ietatof,&iphitof);
			posEff = hEff_Tpc_Pos[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*hEff_Tof_Pos[ietatof][iphitof]->GetBinContent(ipttof+1)*fEPiEff_Tof->Eval(smrcpt1)*funnSigEffP->Eval(smePos.P())*funBetaEffP->Eval(smePos.P())*funndEdxEffPt->Eval(smrcpt1);

			tofPtEtaPhi2Bin(mCenIdx, smrcpt2,eNegEta,eNegPhi, &ipttof,&ietatof,&iphitof);
			tpcPtEtaPhi2Bin(mCenIdx, smrcpt2,eNegEta,eNegPhi, &ipttpc,&ietatpc,&iphitpc);
			negEff = hEff_Tpc_Neg[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*hEff_Tof_Neg[ietatof][iphitof]->GetBinContent(ipttof+1)*fEPiEff_Tof->Eval(smrcpt2)*funnSigEffP->Eval(smeNeg.P())*funBetaEffP->Eval(smeNeg.P())*funndEdxEffPt->Eval(smrcpt2);
		}

		hMCFullYMvsPtdy2ee->Fill(eePt,eeM);

		int kPosPassFlag = PHENIXFilter(1,smePos.Px(),smePos.Py(),smePos.Pz(),0);
		int kNegPassFlag = PHENIXFilter(-1,smeNeg.Px(),smeNeg.Py(),smeNeg.Pz(),0);
		if(kPosPassFlag>=0&&kNegPassFlag>=0) {
			hMCAcc2MvsPtdy2ee->Fill(eePt,eeM);
		}

		if(fabs(eePair.Rapidity())>YCut) continue;
		hMCAcc0MvsPtdy2ee->Fill(eePt,eeM);

		if(smrcpt1<PtCut || smrcpt2<PtCut ) continue;
		if(fabs(ePosEta)>EtaCut || fabs(eNegEta)>EtaCut) continue;
		hMCAcc1MvsPtdy2ee->Fill(eePt,eeM);

		hRCAcc1MvsPt3Ddy2ee->Fill(eePt,eeM,posEff*negEff);
	}

	TFile *fout = new TFile(outFile,"recreate");	

	fout->cd();

	hnParent->Write();
	hMCFullYMvsPtdy2ee->Write();
	hMCAcc0MvsPtdy2ee->Write();
	hMCAcc1MvsPtdy2ee->Write();
	hMCAcc2MvsPtdy2ee->Write();
	hRCAcc1MvsPt3Ddy2ee->Write();

	fout->Close();

	cout << "End!!" << endl;
}
