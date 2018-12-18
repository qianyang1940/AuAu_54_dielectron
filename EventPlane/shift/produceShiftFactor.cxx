#include "/star/u/tc88qy/AuAu/run17/54GeV/EventPlane/func_macro/headers.h"
#include "/star/u/tc88qy/AuAu/run17/54GeV/EventPlane/func_macro/function.C"
#include "/star/u/tc88qy/AuAu/run17/54GeV/Analysis/cuts.h"

#include "miniDst.h"

TFile *fOutFile;

map<Int_t,Int_t> mTotalRunId;
map<Int_t,Int_t> mTotalDayId;
map<Int_t,Int_t> mBadRunId_001;
map<Int_t,Int_t> mBadRunId_021;

bool Init();
Double_t GetRefMultCorr(const Int_t RefMult, const Double_t z);
Double_t GetWeight(Double_t refMultCorr);
Int_t GetCentrality(Double_t refMultCorr);
void bookHistograms(char* outFile);
bool passEvent(miniDst const* const event);

TFile *fReCenter;
TProfile2D *etapluszplusQx;
TProfile2D *etapluszminusQx;
TProfile2D *etaminuszplusQx;
TProfile2D *etaminuszminusQx;
TProfile2D *etapluszplusQy;
TProfile2D *etapluszminusQy;
TProfile2D *etaminuszplusQy;
TProfile2D *etaminuszminusQy;

//define histograms
TH1D *hRawEventPlane;
TH1D *hReCenterEventPlane;
TProfile2D *shiftfactorcos[mArrayLength];
TProfile2D *shiftfactorsin[mArrayLength];

int main(int argc, char** argv)
{
	if(argc!=1&&argc!=3) return -1;

	TString inFile="test.list";
	char outFile[1024];
	sprintf(outFile,"test");
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
	if( !Init() ){
		cout<<"The initialization is failed !!!"<<endl;
		return 0;
	}
	bookHistograms(outFile);

	//+-------------+
	//| loop events |
	//+-------------+
	miniDst *event = new miniDst(chain);
	Int_t nEvts = chain->GetEntries();
	cout<<nEvts<<" events"<<endl;

	for(int i=0;i<nEvts;i++){

		if(i%(nEvts/10)==0) cout << "begin " << i << "th entry...." << endl;

		event->GetEntry(i);

		if(!passEvent(event)) continue;
	}

	if(fOutFile){
		fOutFile->cd();
		fOutFile->Write();
		fOutFile->Close();
	}

	delete chain;

	cout<<"end of program"<<endl;
	return 0;
}
//________________________________________________________________
bool passEvent(miniDst const* const event)
{
	Bool_t validTrig = kFALSE;
	Bool_t RefMVzCorFlag = kFALSE;
	Bool_t is001Trigger = kFALSE;
	Bool_t is021Trigger = kFALSE;
	for(Int_t i=0; i<event->mNTrigs; i++){
		if(event->mTrigId[i] == 580001) validTrig = kTRUE, is001Trigger = kTRUE, RefMVzCorFlag= kTRUE;
		//if(event->mTrigId[i] == 580011) validTrig = kTRUE;
		if(event->mTrigId[i] == 580021) validTrig = kTRUE, is021Trigger = kTRUE;
	}
	if(!validTrig){
		return kFALSE;
	}

	Int_t runId  = event->mRunId;
	if(runId<mMinRunId) return kFALSE;

    Int_t runIndex;
	map<Int_t, Int_t>::iterator iter = mTotalRunId.find(runId);
	if(iter != mTotalRunId.end()){
		runIndex = iter->second;
	}
	else{
		cout<<"Can not find the runNumber in the mTotalRunId list"<<endl;
		return kFALSE;
	}

	Int_t dayIndex;
	iter = mTotalDayId.find((runId/1000)%1000);
	if(iter != mTotalDayId.end()){
		dayIndex = iter->second;
	}
	else{
		cout<<"Can not find the dayNumber in the mTotalDayId list"<<endl;
		return kFALSE;
	}

	map<Int_t, Int_t>::iterator iter_001 = mBadRunId_001.find(runId);
	if(iter_001 != mBadRunId_001.end() && is001Trigger){
		//cout<<"bad run, continue"<<endl;
		return kFALSE;
	}

	map<Int_t, Int_t>::iterator iter_021 = mBadRunId_021.find(runId);
	if(iter_021 != mBadRunId_021.end() && is021Trigger){
		//cout<<"bad run, continue"<<endl;
		return kFALSE;
	}

	Int_t   mCentrality  = event->mCentrality;
	Double_t vx           = event->mVertexX;
	Double_t vy           = event->mVertexY;
	Double_t vz           = event->mVertexZ;
	Int_t  refMult       = event->mRefMult;
	Float_t vpdVz        = event->mVpdVz;
	Float_t vr           = sqrt(vx*vx + vy*vy);    
	Float_t vzDiff       = vz - vpdVz;

	if(vr>mVrCut)                     return kFALSE;
	if(TMath::Abs(vz)>mVzCut)         return kFALSE;
	if(TMath::Abs(vzDiff)>mVzDiffCut) return kFALSE;
	if(mCentrality<0)                 return kFALSE;

    Double_t RefMultCorr = refMult;
	if(RefMVzCorFlag)RefMultCorr = GetRefMultCorr(refMult, vz);	
    Double_t reweight = GetWeight(RefMultCorr);	
	mCentrality = GetCentrality(RefMultCorr);

	Float_t  mEtaPlusQx        = event->mEtaPlusQx;
	Float_t  mEtaPlusQy        = event->mEtaPlusQy;
	Float_t  mEtaPlusPtWeight  = event->mEtaPlusPtWeight;
	Short_t  mEtaPlusNTrks     = event->mEtaPlusNTrks;
	Float_t  mEtaMinusQx       = event->mEtaMinusQx;
	Float_t  mEtaMinusQy       = event->mEtaMinusQy;
	Float_t  mEtaMinusPtWeight = event->mEtaMinusPtWeight;
	Short_t  mEtaMinusNTrks    = event->mEtaMinusNTrks;

	Double_t mRawQx = mEtaPlusQx + mEtaMinusQx;
	Double_t mRawQy = mEtaPlusQy + mEtaMinusQy;

	TVector2 *mRawQ = new TVector2(mRawQx, mRawQy);
	if(mRawQ->Mod() > 0){
		Double_t mRawEventPlane = 0.5*mRawQ->Phi();
		if(mRawEventPlane<0.) mRawEventPlane += TMath::Pi();
		hRawEventPlane->Fill(mRawEventPlane);
	}

	// reCenter process
	Double_t mReCenterQx, mReCenterQy;
	if(vz>0){
		mReCenterQx = mRawQx - mEtaPlusNTrks*etapluszplusQx->GetBinContent(runIndex+1, mCentrality+1) - mEtaMinusNTrks*etaminuszplusQx->GetBinContent(runIndex+1, mCentrality+1);
		mReCenterQy = mRawQy - mEtaPlusNTrks*etapluszplusQy->GetBinContent(runIndex+1, mCentrality+1) - mEtaMinusNTrks*etaminuszplusQy->GetBinContent(runIndex+1, mCentrality+1);
	}
	else{
		mReCenterQx = mRawQx - mEtaPlusNTrks*etapluszminusQx->GetBinContent(runIndex+1, mCentrality+1) - mEtaMinusNTrks*etaminuszminusQx->GetBinContent(runIndex+1, mCentrality+1);
		mReCenterQy = mRawQy - mEtaPlusNTrks*etapluszminusQy->GetBinContent(runIndex+1, mCentrality+1) - mEtaMinusNTrks*etaminuszminusQy->GetBinContent(runIndex+1, mCentrality+1);
	}

    TVector2 *mReCenterQ = new TVector2(mReCenterQx, mReCenterQy);
	Double_t mReCenterEventPlane;
	if(mReCenterQ->Mod() > 0){
		mReCenterEventPlane = 0.5*mReCenterQ->Phi();
		if(mReCenterEventPlane<0.) mReCenterEventPlane += TMath::Pi();
		hReCenterEventPlane->Fill(mReCenterEventPlane);
	}

	for(Int_t j=0; j<mArrayLength; j++){
			shiftfactorcos[j]->Fill(dayIndex,mCentrality,cos(2*(j+1)*mReCenterEventPlane));
			shiftfactorsin[j]->Fill(dayIndex,mCentrality,sin(2*(j+1)*mReCenterEventPlane));
	}

	return kTRUE;
}
//____________________________________________________________
void bookHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.histo.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	fOutFile = new TFile(buf,"recreate");

	for(Int_t i=0;i<mArrayLength;i++){
		sprintf(buf,"shiftfactorcos_%d",i);
		shiftfactorcos[i] = new TProfile2D(buf,buf,mTotalDay,0,mTotalDay,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorsin_%d",i);
		shiftfactorsin[i] = new TProfile2D(buf,buf,mTotalDay,0,mTotalDay,mTotalCentrality,0,mTotalCentrality);
	}

	hRawEventPlane = new TH1D("hRawEventPlane","hRawEventPlane,Event Plane",360,0,TMath::Pi());
	hReCenterEventPlane = new TH1D("hReCenterEventPlane","hReCenterEventPlane,Event Plane",360,0,TMath::Pi());
}
//____________________________________________________________
bool Init()
{
	cout<<endl;

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

	indata.open("/star/u/tc88qy/AuAu/run17/54GeV/QA/runList/output_all/GetRun/runList/mTotalDayList.dat");
	mTotalDayId.clear();
	if(indata.is_open()){
		cout<<"read in day number list and recode day number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata>>oldId){
			mTotalDayId[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the day number list !!!"<<endl;
		return kFALSE;
	}
	indata.close();
	for(map<Int_t,Int_t>::iterator iter=mTotalDayId.begin();iter!=mTotalDayId.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;
	cout<<endl;


	//read in bad run for 580001 and 580021
	ifstream indata_001;
	indata_001.open("/star/u/tc88qy/AuAu/run17/54GeV/EventPlane/reCenter/Badrun/mBadRunList_001.dat");
	mBadRunId_001.clear();
	if(indata_001.is_open()){
		cout<<"read in total run number list and recode run number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata_001>>oldId){
			mBadRunId_001[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list !!!"<<endl;
		return kFALSE;
	}
	indata_001.close();
	
	//read in bad run for 580001 and 580021
	ifstream indata_021;
	indata_021.open("/star/u/tc88qy/AuAu/run17/54GeV/EventPlane/reCenter/Badrun/mBadRunList_021.dat");
	mBadRunId_021.clear();
	if(indata_021.is_open()){
		cout<<"read in total run number list and recode run number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata_021>>oldId){
			mBadRunId_021[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list !!!"<<endl;
		return kFALSE;
	}
	indata_021.close();

	cout<<"bad run for trigger 580001"<<endl;
	for(map<Int_t,Int_t>::iterator iter=mBadRunId_001.begin();iter!=mBadRunId_001.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;

	cout<<"bad run for trigger 580021"<<endl;
	for(map<Int_t,Int_t>::iterator iter=mBadRunId_021.begin();iter!=mBadRunId_021.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;


    fReCenter = TFile::Open("/star/u/tc88qy/AuAu/run17/54GeV/EventPlane/reCenter/output_all/reCenter.root");
	etapluszplusQx   = (TProfile2D *)fReCenter->Get("etapluszplusQx");
	etapluszminusQx  = (TProfile2D *)fReCenter->Get("etapluszminusQx");
	etaminuszplusQx  = (TProfile2D *)fReCenter->Get("etaminuszplusQx");
	etaminuszminusQx = (TProfile2D *)fReCenter->Get("etaminuszminusQx");
	etapluszplusQy   = (TProfile2D *)fReCenter->Get("etapluszplusQy");
	etapluszminusQy  = (TProfile2D *)fReCenter->Get("etapluszminusQy");
	etaminuszplusQy  = (TProfile2D *)fReCenter->Get("etaminuszplusQy");
	etaminuszminusQy = (TProfile2D *)fReCenter->Get("etaminuszminusQy");

	cout<<"Initialization DONE !!!"<<endl;
	cout<<endl;

	return kTRUE;
}
Double_t GetWeight(Double_t refMultCorr)
{
	Double_t Weight = 1.0;

	Double_t par0 =3.905;
	Double_t par1 = -204.4;
	Double_t par2 = 1.851;
	Double_t par3 = 24.32;
	Double_t par4 = -0.01746;
	Double_t A    = 0;
	Double_t par6 = 6405;
	Double_t par7 = 3.743e-05;

	// Additional z-vertex dependent correction
	//const Double_t A = ((1.27/1.21))/(30.0*30.0); // Don't ask...
	//const Double_t A = (0.05/0.21)/(30.0*30.0); // Don't ask...
	Double_t mVz = 30.0; // this has to be modified...

	if(
			refMultCorr != -(par3/par2) // avoid denominator = 0
			&& refMultCorr<=70         // normalization stop at refMultCorr=300
	  )
	{
		Weight = par0 + par1/(par2*refMultCorr + par3) + par4*(par2*refMultCorr + par3) + par6/pow(par2*refMultCorr+par3 ,2) + par7*pow(par2*refMultCorr+par3 ,2); // Parametrization of MC/data RefMult ratio
		//Weight = par0 + par1/(par2*refMultCorr + par3) + par4*(par2*refMultCorr + par3) ; // Parametrization of MC/data RefMult ratio
		Weight = Weight + (Weight-1.0)*(A*mVz*mVz); // z-dependent weight correction
	}

	return Weight;
}

Int_t GetCentrality(Double_t refMultCorr){
	if(refMultCorr > 361){
		return 9;
	}else if (refMultCorr > 299){
		return 8;
	}else if (refMultCorr > 205){
		return 7;
	}else if (refMultCorr > 138){
		return 6;
	}else if (refMultCorr > 89){
		return 5;
	}else if (refMultCorr > 54){
		return 4;
	}else if (refMultCorr > 31){
		return 3;
	}else if (refMultCorr > 16){
		return 2;
	}else if (refMultCorr > 7){
		return 1;
	}else{
		return 0;
	}
}

Double_t GetRefMultCorr(const Int_t RefMult, const Double_t z)
{
	const Double_t par_z[7] = {435.9, -0.02413, -0.003707, 0.0002204, 1.487e-5, -2.95e-7, -1.867e-8};
	const Double_t RefMult_ref = par_z[0];
	const Double_t  RefMult_z = par_z[0] + par_z[1]*z + par_z[2]*z*z + par_z[3]*z*z*z + par_z[4]*z*z*z*z + par_z[5]*z*z*z*z*z + par_z[6]*z*z*z*z*z*z; // Parametrization of mean RefMult vs. z_vertex position
	Double_t Hovno =1.0;
	if(RefMult_z > 0.0)
	{
		Hovno = (RefMult_ref + par_z[7])/RefMult_z;
	}
	Double_t RefMult_d = (Double_t)(RefMult) + gRandom->Rndm(); // random sampling over bin width -> avoid peak structures in corrected distribution
	return RefMult_d*Hovno;
}

