#include "headers.h"
#include "StMiniTreeMaker.h"


ClassImp(StMiniTreeMaker)

//_____________________________________________________________________________
StMiniTreeMaker::StMiniTreeMaker(const Char_t *name) : StMaker(name), mFillTree(0), mFillHisto(1), mPrintConfig(1), mPrintMemory(0), mPrintCpu(0), mStreamName(""), fOutFile(0), mOutFileName(""), mEvtTree(0), mDefaultVtx(1), mSelectVtxRank(1), mMaxVtxR(1.e4), mMaxVtxZ(1.e4), mMaxVzDiff(1.e4), mMinTrkPt(0.2), mMaxTrkEta(1.), mMinNHitsFit(15), mMinNHitsFitRatio(0.52), mMinNHitsDedx(10), mMaxDca(10.), mMaxnSigmaE(2.5), mMaxBeta2TOF(0.03)
{
	// default constructor

	// run15 st_mtd 
	mTriggerIDs.clear();
	mTriggerIDs.push_back(580001);     // AuAu@54  minbias
	mTriggerIDs.push_back(580011);     // AuAu@54  minbias 
	mTriggerIDs.push_back(580021);     // AuAu@54  minbias 
}
//_____________________________________________________________________________
StMiniTreeMaker::~StMiniTreeMaker()
{
	// default destructor
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::Init()
{
//	refMultCorr = CentralityMaker::instance()->getgRefMultCorr_VpdMBnoVtx()
//	refMultCorr = new StRefMultCorr("grefmult_VpdMBnoVtx");

	cout<<"StMiniTreeMaker::Init()"<<endl;
	if(!mOutFileName.Length()){
		LOG_ERROR << "StMiniTreeMaker:: no output file specified for tree and histograms." << endm;
		return kStERR;
	}
	fOutFile = new TFile(mOutFileName.Data(),"recreate");
	LOG_INFO << "StMiniTreeMaker:: create the output file to store the tree and histograms: " << mOutFileName.Data() << endm;

	if(mFillTree)    bookTree();
	if(mFillHisto)   bookHistos();

	return kStOK;
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::InitRun(const Int_t runnumber)
{
	LOG_INFO<<"Grab runnumber: "<<runnumber<<endm;
	return kStOK;
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::Finish()
{
	if(fOutFile){
		fOutFile->cd();
		fOutFile->Write();
		fOutFile->Close();
		LOG_INFO << "StMiniTreeMaker::Finish() -> write out tree in " << mOutFileName.Data() << endm;
	}
	if(mPrintConfig) printConfig();
	return kStOK;
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::Make()
{
	memset(&mEvtData, 0, sizeof(mEvtData)); //initial the mEvtData structure

	StTimer timer;
	if(mPrintMemory) StMemoryInfo::instance()->snapshot();
	if(mPrintCpu)    timer.start();

	mPicoDstMaker = (StPicoDstMaker *)GetMaker("picoDst");
	if(Debug()){
		LOG_INFO<<"PicoDstMaker pointer: "<<mPicoDstMaker<<endm;
	}

	if(mPicoDstMaker){
		if(Debug()) LOG_INFO<<"Use Pico file as input"<<endm;
		mPicoDst = mPicoDstMaker->picoDst();
		if(!mPicoDst){
			LOG_WARN<<"No PicoDst !"<<endm;
			return kStOK;
		}
	}
	else{
		LOG_WARN<<"No StPicoDstMaker !"<<endm;
		return kStOK;
	}


	// for event plane calculation
	vEtaPlusCosPart.clear();
	vEtaPlusSinPart.clear();
	vEtaPlusPtWeight.clear();
	vEtaMinusCosPart.clear();
	vEtaMinusSinPart.clear();
	vEtaMinusPtWeight.clear();

	if(!processPicoEvent())  return kStOK;

	if(mFillTree) mEvtTree->Fill();

	if(mPrintMemory){
		StMemoryInfo::instance()->snapshot();
		StMemoryInfo::instance()->print();
	}

	if(mPrintCpu){
		timer.stop();
		LOG_INFO << "CPU time for StMiniTreeMaker::Make(): " 
			<< timer.elapsedTime() << "sec " << endm;
	}

	return kStOK;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::processPicoEvent()
{
	if(mFillHisto) hEvent->Fill(0.5);

	StPicoEvent *picoEvent = mPicoDst->event();
	if(!picoEvent){
		LOG_WARN<<"No event level information !"<<endm;
		return kFALSE;
	}

	Bool_t validTrigger = kFALSE;
	Bool_t minbias      = kFALSE;

	Int_t nTrigs = 0;
	for(Int_t i=0;i<mTriggerIDs.size();i++){
		if(picoEvent->isTrigger(mTriggerIDs[i])){
			minbias = kTRUE;
			validTrigger = kTRUE;
			mEvtData.mTrigId[nTrigs] = mTriggerIDs[i];
			nTrigs++;
		}
	}

	mEvtData.mNTrigs = nTrigs;

	if(!validTrigger){
		return kFALSE;
	}

	if(mFillHisto){
		if(minbias)     hEvent->Fill(2.5);
	}

	mEvtData.mRunId          = picoEvent->runId();
	mEvtData.mEventId        = picoEvent->eventId();
	mEvtData.mRefMult        = picoEvent->refMult();
	mEvtData.mGRefMult       = picoEvent->grefMult();
	mEvtData.mBBCRate        = picoEvent->BBCx();
	mEvtData.mZDCRate        = picoEvent->ZDCx();
	mEvtData.mBField         = picoEvent->bField();
	mEvtData.mVpdVz          = picoEvent->vzVpd();

	StThreeVectorF vtxPos    = picoEvent->primaryVertex();
	mEvtData.mVertexX        = vtxPos.x();
	mEvtData.mVertexY        = vtxPos.y();
	mEvtData.mVertexZ        = vtxPos.z();
	if(Debug()){
		LOG_INFO<<"RunId: "<<mEvtData.mRunId<<endm;
		LOG_INFO<<"EventId: "<<mEvtData.mEventId<<endm;
		LOG_INFO<<"mZDCX: "<<mEvtData.mZDCRate<<endm;
		LOG_INFO<<"VPD Vz: "<<mEvtData.mVpdVz<<" \tTPC Vz: "<<mEvtData.mVertexZ<<endm;
	}

	if(mFillHisto){
		hVtxYvsVtxX->Fill(mEvtData.mVertexX, mEvtData.mVertexY);
		hVPDVzvsTPCVz->Fill(mEvtData.mVertexZ, mEvtData.mVpdVz);
		hVzDiff->Fill(mEvtData.mVertexZ - mEvtData.mVpdVz);
	}
	if(Debug()) LOG_INFO<<"gRefMult: "<<mEvtData.mGRefMult<<" \tgRefMultCorr: "<<mEvtData.mGRefMultCorr<<" \tmCentrality: "<<mEvtData.mCentrality<<endm;

	if(TMath::Abs(vtxPos.x())<1.e-5 && TMath::Abs(vtxPos.y())<1.e-5 && TMath::Abs(vtxPos.z())<1.e-5) return kFALSE;
	if(mFillHisto) hEvent->Fill(7.5);
	if(sqrt(vtxPos.x()*vtxPos.x()+vtxPos.y()*vtxPos.y())>=mMaxVtxR) return kFALSE;
	if(mFillHisto) hEvent->Fill(8.5);
	if(TMath::Abs(vtxPos.z())>=mMaxVtxZ) return kFALSE;
	if(mFillHisto) hEvent->Fill(9.5);
	if(TMath::Abs(mEvtData.mVertexZ - mEvtData.mVpdVz)>=mMaxVzDiff) return kFALSE;
	if(mFillHisto) hEvent->Fill(10.5);


	Int_t nNodes = mPicoDst->numberOfTracks();
	if(Debug()){
		LOG_INFO<<"# of global tracks in picoDst: "<<nNodes<<endm;
	}

	Short_t nTrks    = 0;
	for(Int_t i=0;i<nNodes;i++){
		StPicoTrack *pTrack = mPicoDst->track(i);
		if(!pTrack) continue;

		calQxQy(pTrack, vtxPos);
		if(!isValidTrack(pTrack, vtxPos)) continue;
		mEvtData.mTrkId[nTrks]            = i;  
		mEvtData.mTPCeTrkFlag[nTrks]      = kFALSE;

		mEvtData.mCharge[nTrks]           = pTrack->charge();

		StThreeVectorF pMom               = pTrack->pMom();
		StThreeVectorF gMom               = pTrack->gMom();
		StThreeVectorF origin             = pTrack->origin();

		mEvtData.mPt[nTrks]               = pMom.perp();
		mEvtData.mEta[nTrks]              = pMom.pseudoRapidity();
		mEvtData.mPhi[nTrks]              = pMom.phi();
		mEvtData.mgPt[nTrks]              = gMom.perp();
		mEvtData.mgEta[nTrks]             = gMom.pseudoRapidity();
		mEvtData.mgPhi[nTrks]             = gMom.phi();
		mEvtData.mgOriginX[nTrks]         = origin.x();
		mEvtData.mgOriginY[nTrks]         = origin.y();
		mEvtData.mgOriginZ[nTrks]         = origin.z();

		mEvtData.mNHitsFit[nTrks]         = pTrack->nHitsFit();
		mEvtData.mNHitsPoss[nTrks]        = pTrack->nHitsMax();
		mEvtData.mNHitsDedx[nTrks]        = pTrack->nHitsDedx();
		mEvtData.mDedx[nTrks]             = pTrack->dEdx(); 
		mEvtData.mNSigmaE[nTrks]          = pTrack->nSigmaElectron();
		mEvtData.mDca[nTrks]              = (pTrack->dcaPoint()-vtxPos).mag();

		if(mFillHisto){
			hdEdxvsP->Fill(pMom.mag(), mEvtData.mDedx[nTrks]);
			hnSigEvsP->Fill(pMom.mag(), mEvtData.mNSigmaE[nTrks]);
		}

		Int_t bTofPidTraitsIndex          = pTrack->bTofPidTraitsIndex();
		mEvtData.mTOFMatchFlag[nTrks]     = -1; 
		mEvtData.mTOFLocalY[nTrks]        = -999;
		mEvtData.mBeta2TOF[nTrks]         = -999;
		if(bTofPidTraitsIndex>=0){
			StPicoBTofPidTraits *btofPidTraits = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
			mEvtData.mTOFMatchFlag[nTrks] = btofPidTraits->btofMatchFlag(); 
			mEvtData.mTOFLocalY[nTrks]    = btofPidTraits->btofYLocal();
			mEvtData.mBeta2TOF[nTrks]     = btofPidTraits->btofBeta();

			if(mFillHisto) hBetavsP->Fill(pMom.mag(), 1./mEvtData.mBeta2TOF[nTrks]);
		}

		if(
				TMath::Abs(mEvtData.mNSigmaE[nTrks])<=mMaxnSigmaE
				&& mEvtData.mBeta2TOF[nTrks]>0.
				&& TMath::Abs(1.-1./mEvtData.mBeta2TOF[nTrks])<=mMaxBeta2TOF
		  )
			mEvtData.mTPCeTrkFlag[nTrks] = kTRUE;



		if(mEvtData.mTPCeTrkFlag[nTrks] ){
			nTrks++;
		}
	}

	//if(nTrks==0 ) return kFALSE;

	mEvtData.mNTrks       = nTrks;
	if(Debug()){
		LOG_INFO<<"# of primary tracks stored: "<<mEvtData.mNTrks<<endm;
	}

	fillEventPlane();

	return kTRUE;
}
//_____________________________________________________________________________
void StMiniTreeMaker::calQxQy(StPicoTrack *pTrack, StThreeVectorF vtxPos) const
{
	Float_t pt  = pTrack->pMom().perp();
	Float_t eta = pTrack->pMom().pseudoRapidity();
	Float_t phi = pTrack->pMom().phi();
	//Float_t dca = pTrack->dca();
	//Float_t dca = (pTrack->dca()-vtxPos).mag();
	Float_t dca = (pTrack->dcaPoint()-vtxPos).mag();

	if(pt<0.15 || pt>2.)                              return;
	if(TMath::Abs(eta)>1.)                            return;
	if(pTrack->nHitsFit()<15)                         return;
	if(pTrack->nHitsFit()*1./pTrack->nHitsMax()<0.52) return;
	if(dca>1.)                                        return;

	Double_t mCosPart = pt*TMath::Cos(2.*phi);
	Double_t mSinPart = pt*TMath::Sin(2.*phi);

	if(eta>0.){
		vEtaPlusCosPart.push_back(mCosPart);
		vEtaPlusSinPart.push_back(mSinPart);
		vEtaPlusPtWeight.push_back(pt);
	}else if(eta<0.){
		vEtaMinusCosPart.push_back(mCosPart);
		vEtaMinusSinPart.push_back(mSinPart);
		vEtaMinusPtWeight.push_back(pt);
	}

}
//_____________________________________________________________________________
void StMiniTreeMaker::fillEventPlane()
{
	mEvtData.mEtaPlusQx       = 0.;
	mEvtData.mEtaPlusQy       = 0.;
	mEvtData.mEtaPlusPtWeight = 0.;
	for(Int_t i=0;i<vEtaPlusCosPart.size();i++){
		mEvtData.mEtaPlusQx       += vEtaPlusCosPart[i];
		mEvtData.mEtaPlusQy       += vEtaPlusSinPart[i];
		mEvtData.mEtaPlusPtWeight += vEtaPlusPtWeight[i];
	}
	mEvtData.mEtaPlusNTrks   = vEtaPlusCosPart.size();

	mEvtData.mEtaMinusQx       = 0.;
	mEvtData.mEtaMinusQy       = 0.;
	mEvtData.mEtaMinusPtWeight = 0.;
	for(Int_t i=0;i<vEtaMinusCosPart.size();i++){
		mEvtData.mEtaMinusQx       += vEtaMinusCosPart[i];
		mEvtData.mEtaMinusQy       += vEtaMinusSinPart[i];
		mEvtData.mEtaMinusPtWeight += vEtaMinusPtWeight[i];
	}
	mEvtData.mEtaMinusNTrks  = vEtaMinusCosPart.size();

	if(mFillHisto){
		Double_t mRawQx = mEvtData.mEtaPlusQx+mEvtData.mEtaMinusQx;
		Double_t mRawQy = mEvtData.mEtaPlusQy+mEvtData.mEtaMinusQy;
		hRawQx->Fill(mRawQx);
		hRawQy->Fill(mRawQy);

		TVector2 *mRawQ = new TVector2(mRawQx, mRawQy);
		if(mRawQ->Mod()>0.){
			Double_t mRawEventPlane = 0.5*mRawQ->Phi();
			if(mRawEventPlane<0.) mRawEventPlane += TMath::Pi();
			hRawEventPlane->Fill(mRawEventPlane);
		}
		mRawQ->Delete();
	}
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::isValidTrack(StPicoTrack *pTrack, StThreeVectorF vtxPos) const
{
	Float_t pt  = pTrack->pMom().perp();
	Float_t eta = pTrack->pMom().pseudoRapidity();
	//Float_t dca = pTrack->dca();
	//Float_t dca = (pTrack->dca()-vtxPos).mag();
	Float_t dca = (pTrack->dcaPoint()-vtxPos).mag();

	if(pt<mMinTrkPt)                            return kFALSE;
	if(TMath::Abs(eta)>mMaxTrkEta)              return kFALSE;
	if(pTrack->nHitsFit()<mMinNHitsFit)         return kFALSE;
	if(pTrack->nHitsFit()*1./pTrack->nHitsMax()<mMinNHitsFitRatio) return kFALSE;
	if(pTrack->nHitsDedx()<mMinNHitsDedx)       return kFALSE;
	if(dca>mMaxDca)                             return kFALSE;

	return kTRUE;
}
//_____________________________________________________________________________
void StMiniTreeMaker::bookTree()
{
	LOG_INFO << "StMiniTreeMaker:: book the event tree to be filled." << endm;

	mEvtTree = new TTree("miniDst","miniDst");
	mEvtTree->SetAutoSave(100000000); // 100 MB

	// event information
	mEvtTree->Branch("mRunId", &mEvtData.mRunId, "mRunId/I");
	mEvtTree->Branch("mEventId", &mEvtData.mEventId, "mEventId/I");
	mEvtTree->Branch("mShouldHaveRejectEvent", &mEvtData.mShouldHaveRejectEvent, "mShouldHaveRejectEvent/B");
	mEvtTree->Branch("mNTrigs", &mEvtData.mNTrigs, "mNTrigs/I");
	mEvtTree->Branch("mTrigId", mEvtData.mTrigId, "mTrigId[mNTrigs]/I");
	mEvtTree->Branch("mRefMult", &mEvtData.mRefMult, "mRefMult/S");
	mEvtTree->Branch("mGRefMult", &mEvtData.mGRefMult, "mGRefMult/S");
	mEvtTree->Branch("mGRefMultCorr", &mEvtData.mGRefMultCorr, "mGRefMultCorr/F");
	mEvtTree->Branch("mEvtWeight", &mEvtData.mEvtWeight, "mEvtWeight/F");
	mEvtTree->Branch("mCentrality", &mEvtData.mCentrality, "mCentrality/S");
	mEvtTree->Branch("mBBCRate", &mEvtData.mBBCRate, "mBBCRate/F");
	mEvtTree->Branch("mZDCRate", &mEvtData.mZDCRate, "mZDCRate/F");
	mEvtTree->Branch("mBField", &mEvtData.mBField, "mBField/F");
	mEvtTree->Branch("mVpdVz", &mEvtData.mVpdVz, "mVpdVz/F");
	mEvtTree->Branch("mVertexX", &mEvtData.mVertexX, "mVertexX/F");
	mEvtTree->Branch("mVertexY", &mEvtData.mVertexY, "mVertexY/F");
	mEvtTree->Branch("mVertexZ", &mEvtData.mVertexZ, "mVertexZ/F");
	mEvtTree->Branch("mVertexRanking", &mEvtData.mVertexRanking, "mVertexRanking/F");
	mEvtTree->Branch("mEtaPlusQx", &mEvtData.mEtaPlusQx, "mEtaPlusQx/F");
	mEvtTree->Branch("mEtaPlusQy", &mEvtData.mEtaPlusQy, "mEtaPlusQy/F");
	mEvtTree->Branch("mEtaPlusPtWeight", &mEvtData.mEtaPlusPtWeight, "mEtaPlusPtWeight/F");
	mEvtTree->Branch("mEtaPlusNTrks", &mEvtData.mEtaPlusNTrks, "mEtaPlusNTrks/S");
	mEvtTree->Branch("mEtaMinusQx", &mEvtData.mEtaMinusQx, "mEtaMinusQx/F");
	mEvtTree->Branch("mEtaMinusQy", &mEvtData.mEtaMinusQy, "mEtaMinusQy/F");
	mEvtTree->Branch("mEtaMinusPtWeight", &mEvtData.mEtaMinusPtWeight, "mEtaMinusPtWeight/F");
	mEvtTree->Branch("mEtaMinusNTrks", &mEvtData.mEtaMinusNTrks, "mEtaMinusNTrks/S");

	//all tracks information
	mEvtTree->Branch("mNTrks", &mEvtData.mNTrks, "mNTrks/S");
	mEvtTree->Branch("mTrkId", mEvtData.mTrkId, "mTrkId[mNTrks]/S");
	mEvtTree->Branch("mTPCeTrkFlag", mEvtData.mTPCeTrkFlag, "mTPCeTrkFlag[mNTrks]/O");
	mEvtTree->Branch("mCharge", mEvtData.mCharge, "mCharge[mNTrks]/I");
	mEvtTree->Branch("mPt", mEvtData.mPt, "mPt[mNTrks]/F");
	mEvtTree->Branch("mEta", mEvtData.mEta, "mEta[mNTrks]/F");
	mEvtTree->Branch("mPhi", mEvtData.mPhi, "mPhi[mNTrks]/F");
	mEvtTree->Branch("mgPt", mEvtData.mgPt, "mgPt[mNTrks]/F");
	mEvtTree->Branch("mgEta", mEvtData.mgEta, "mgEta[mNTrks]/F");
	mEvtTree->Branch("mgPhi", mEvtData.mgPhi, "mgPhi[mNTrks]/F");
	mEvtTree->Branch("mgOriginX", mEvtData.mgOriginX, "mgOriginX[mNTrks]/F");
	mEvtTree->Branch("mgOriginY", mEvtData.mgOriginY, "mgOriginY[mNTrks]/F");
	mEvtTree->Branch("mgOriginZ", mEvtData.mgOriginZ, "mgOriginZ[mNTrks]/F");
	mEvtTree->Branch("mNHitsFit", mEvtData.mNHitsFit, "mNHitsFit[mNTrks]/I");
	mEvtTree->Branch("mNHitsPoss", mEvtData.mNHitsPoss, "mNHitsPoss[mNTrks]/I");
	mEvtTree->Branch("mNHitsDedx", mEvtData.mNHitsDedx, "mNHitsDedx[mNTrks]/I");
	mEvtTree->Branch("mDedx", mEvtData.mDedx, "mDedx[mNTrks]/F");
	mEvtTree->Branch("mNSigmaE", mEvtData.mNSigmaE, "mNSigmaE[mNTrks]/F");
	mEvtTree->Branch("mDca", mEvtData.mDca, "mDca[mNTrks]/F");
	mEvtTree->Branch("mTOFMatchFlag", mEvtData.mTOFMatchFlag, "mTOFMatchFlag[mNTrks]/I");
	mEvtTree->Branch("mTOFLocalY", mEvtData.mTOFLocalY, "mTOFLocalY[mNTrks]/F");
	mEvtTree->Branch("mBeta2TOF", mEvtData.mBeta2TOF, "mBeta2TOF[mNTrks]/F");

	//AddObj(mEvtTree,".hist");
}
//_____________________________________________________________________________
void StMiniTreeMaker::bookHistos()
{
	hEvent = new TH1D("hEvent","Event statistics",25,0,25);
	hEvent->GetXaxis()->SetBinLabel(1, "All events");
	hEvent->GetXaxis()->SetBinLabel(3, "minbias");
	hEvent->GetXaxis()->SetBinLabel(8, "None-Zero Vertex");
	hEvent->GetXaxis()->SetBinLabel(9,  Form("|V_{r}|<%1.2f cm",mMaxVtxR));
	hEvent->GetXaxis()->SetBinLabel(10, Form("|V_{z}|<%1.2f cm",mMaxVtxZ));
	hEvent->GetXaxis()->SetBinLabel(11, Form("|V_{z}Diff|<%1.2f cm",mMaxVzDiff));

	hVtxYvsVtxX = new TH2D("hVtxYvsVtxX","hVtxYvsVtxX; V_{x} (cm); V_{y} (cm)",120,-3,3,120,-3,3); 
	hVPDVzvsTPCVz = new TH2D("hVPDVzvsTPCVz","hVPDVzvsTPCVz; TPC V_{z} (cm); VPD V_{z} (cm)",200,-50,50,200,-50,50);
	hVzDiff = new TH1D("hVzDiff","hVzDiff; TPC V_{z} - VPD V_{z} (cm)",80,-20,20);
	hGRefMultvsGRefMultCorr = new TH2D("hGRefMultvsGRefMultCorr","hGRefMultvsGRefMultCorr; grefMultCorr; grefMult",1000,0,1000,1000,0,1000);
	hCentrality = new TH1D("hCentrality","hCentrality; mCentrality",16,0,16);

	hRawQx = new TH1D("hRawQx","hRawQx;Q_{x}",1000,-50,50);
	hRawQy = new TH1D("hRawQy","hRawQy;Q_{y}",1000,-50,50);
	hRawEventPlane = new TH1D("hRawEventPlane","hRawEventPlane,Event Plane",300,0,TMath::Pi());

	hdEdxvsP = new TH2D("hdEdxvsP","hdEdxvsP; p (GeV/c); dE/dx (KeV/cm)",300,0,15,400,0,20);
	hdNdxvsP = new TH2D("hdNdxvsP","hdNdxvsP; p (GeV/c); dN/dx",300,0,15,400,0,200);
	hnSigEvsP = new TH2D("hnSigEvsP","hnSigEvsP; p (GeV/c); n#sigma_{e}",300,0,15,700,-15,20);
	hBetavsP = new TH2D("hBetavsP","hBetavsP; p (GeV/c); 1/#beta",300,0,15,800,0,4);
}
//_____________________________________________________________________________
void StMiniTreeMaker::printConfig()
{
	const char *decision[2] = {"no","yes"};
	printf("=== Configuration for StMiniTreeMaker ===\n");
	printf("Fill the miniDst tree: %s\n",decision[mFillTree]);
	printf("Fill the QA histo: %s\n",decision[mFillHisto]);
	printf("Use default vertex: %s\n",decision[mDefaultVtx]);
	printf("Select positive vertex ranking: %s\n",decision[mSelectVtxRank]);
	printf("Maximum |Vr|: %1.2f\n",mMaxVtxR);
	printf("Maximum |Vz|: %1.2f\n",mMaxVtxZ);
	printf("Maximum |VzDiff|: %1.2f\n",mMaxVzDiff);
	printf("Minimum track pt: %1.2f\n",mMinTrkPt);
	printf("Maximum track |eta| : %1.2f\n",mMaxTrkEta);
	printf("Minimum number of fit hits: %d\n",mMinNHitsFit);
	printf("Minimum ratio of fit hits: %1.2f\n",mMinNHitsFitRatio);
	printf("Minimum number of dedx hits: %d\n",mMinNHitsDedx);
	printf("Maximum dca: %1.2f\n",mMaxDca);
	printf("Maximum |nSigmaE| for TPCe: %1.2f\n",mMaxnSigmaE);
	printf("Maximum |1-1/beta| for TPCe: %1.2f\n",mMaxBeta2TOF);
	printf("=======================================\n");
}
