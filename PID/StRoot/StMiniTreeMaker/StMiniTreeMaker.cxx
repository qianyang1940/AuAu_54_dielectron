#include "headers.h"
#include "StMiniTreeMaker.h"


ClassImp(StMiniTreeMaker)

//_____________________________________________________________________________
StMiniTreeMaker::StMiniTreeMaker(const Char_t *name) : StMaker(name),  mFillHisto(1), mPrintConfig(1), mPrintMemory(0), mPrintCpu(0), mStreamName(""), fOutFile(0), mOutFileName(""), mEvtTree(0), mDefaultVtx(1), mSelectVtxRank(1), mMaxVtxR(1.e4), mMaxVtxZ(1.e4), mMaxVzDiff(1.e4), mMinTrkPt(0.2), mMaxTrkEta(1.), mMinNHitsFit(15), mMinNHitsFitRatio(0.52), mMinNHitsDedx(10), mMaxDca(10.), mMaxnSigmaE(2.5),mMinnSigmaE(-0.75), mMaxBeta2TOF(0.03)
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

	if(!processPicoEvent())  return kStOK;

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
			nTrigs++;
		}
	}


	if(!validTrigger){
		return kFALSE;
	}

	if(mFillHisto){
		if(minbias)     hEvent->Fill(2.5);
	}

	StThreeVectorF vtxPos    = picoEvent->primaryVertex();
	double vpdvz            = picoEvent->vzVpd();

	if(mFillHisto){
		hVtxYvsVtxX->Fill(vtxPos.x(), vtxPos.y());
		hVPDVzvsTPCVz->Fill(vtxPos.z(), vpdvz);
		hVzDiff->Fill(vtxPos.z() - vpdvz);
	}


	if(TMath::Abs(vtxPos.x())<1.e-5 && TMath::Abs(vtxPos.y())<1.e-5 && TMath::Abs(vtxPos.z())<1.e-5) return kFALSE;
	if(mFillHisto) hEvent->Fill(7.5);
	if(sqrt(vtxPos.x()*vtxPos.x()+vtxPos.y()*vtxPos.y())>=mMaxVtxR) return kFALSE;
	if(mFillHisto) hEvent->Fill(8.5);
	if(TMath::Abs(vtxPos.z())>=mMaxVtxZ) return kFALSE;
	if(mFillHisto) hEvent->Fill(9.5);
	if( TMath::Abs(vtxPos.z() - vpdvz) >= mMaxVzDiff ) return kFALSE;
	if(mFillHisto) hEvent->Fill(10.5);


	Int_t nNodes = mPicoDst->numberOfTracks();
	if(Debug()){
		LOG_INFO<<"# of global tracks in picoDst: "<<nNodes<<endm;
	}

	Short_t nTrks    = 0;
	for(Int_t i=0;i<nNodes;i++){
		StPicoTrack *pTrack = mPicoDst->track(i);
		if(!pTrack) continue;

		if(!isValidTrack(pTrack, vtxPos)) continue;
		StThreeVectorF pMom               = pTrack->pMom();
		double p = pMom.mag();
		double pt1 = pMom.perp();
		double eta1 = pMom.pseudoRapidity();
		double phi1 = pMom.phi();
		double dedx = pTrack->dEdx();
		double nsigmaE1 = pTrack->nSigmaElectron();
		int charge1 = pTrack->charge();

		hdEdxvsP->Fill(p, dedx);
		hnSigEvsP->Fill(p, nsigmaE1);

		Int_t bTofPidTraitsIndex          = pTrack->bTofPidTraitsIndex();
		if(bTofPidTraitsIndex>=0){
			StPicoBTofPidTraits *btofPidTraits = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
			double beta  = btofPidTraits->btofBeta();
			hBetavsP->Fill(p, 1./beta);
			if(1./beta>0 && TMath::Abs(1.- 1./beta) < mMaxBeta2TOF){
				hnSigEvsPWTOF->Fill(p, nsigmaE1);
				if(!isElectron(pTrack))continue;

				for(int j =i+1; j< nNodes; j++){
					StPicoTrack *pNextTrack = mPicoDst->track(j);
					if(!pNextTrack)continue;
					if(!isValidTrack(pNextTrack, vtxPos))continue;
					if(!isElectron(pNextTrack))continue;
					Int_t bTofPidTraitsNextIndex          = pNextTrack->bTofPidTraitsIndex();

					StPicoBTofPidTraits *btofPidTraitsNext = mPicoDst->btofPidTraits(bTofPidTraitsNextIndex);
					double pt2 = pNextTrack->pMom().perp();
					double eta2 = pNextTrack->pMom().pseudoRapidity();
					double phi2 = pNextTrack->pMom().phi();
					double charge2 = pNextTrack->charge();

					TLorentzVector pair, ele1, ele2;
					ele1.SetPtEtaPhiM(pt1,eta1,phi1, Melectron);
					ele2.SetPtEtaPhiM(pt2,eta2,phi2, Melectron);
					pair = ele1 + ele2;

					double parimass = pair.M();
					double phiV = phiVangle(ele1, ele2, charge1, charge2 );
					int Tcharge  = charge2 + charge1;

					if(Tcharge == 0){	
						hULMvsphiV->Fill(parimass,phiV);

					}else if(Tcharge < 0){
						hLPosMvsphiV->Fill(parimass,phiV);
					}else{
						hLNegMvsphiV->Fill(parimass,phiV);
					}

				}//pair 

			}//ele1 tof 
		}//with tof

	}//



	return kTRUE;
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
	hVPDVzvsTPCVz = new TH2D("hVPDVzvsTPCVz","hVPDVzvsTPCVz; TPC V_{z} (cm); VPD V_{z} (cm)",800,-200,200,800,-200,200);
	hVzDiff = new TH1D("hVzDiff","hVzDiff; TPC V_{z} - VPD V_{z} (cm)",400,-20,20);


	hdEdxvsP = new TH2D("hdEdxvsP","hdEdxvsP; p (GeV/c); dE/dx (KeV/cm)",500,0,5,400,0,20);
	hnSigEvsP = new TH2D("hnSigEvsP","hnSigEvsP; p (GeV/c); n#sigma_{e}",3000,0,3,1200,-6,6);
	hnSigEvsPWTOF = new TH2D("hnSigEvsPWTOF","hnSigEvsP; p (GeV/c); n#sigma_{e}",3000,0,3,1200,-6,6);
	hBetavsP = new TH2D("hBetavsP","hBetavsP; p (GeV/c); 1/#beta",3000,0,3,800,0.5,1.3);
	hULMvsphiV = new TH2D("hULMvsphiV","UnLike-sign mass vs phiV;M_{primary} (GeV/c^{2});#phi_{V}",2000,0,0.2,180,0,TMath::Pi());
	hLPosMvsphiV = new TH2D("hLPosMvsphiV","Like-sign mass vs phiV;M_{primary} (GeV/c^{2});#phi_{V}",2000,0,0.2,180,0,TMath::Pi());
	hLNegMvsphiV = new TH2D("hLNegMvsphiV","Like-sign mass vs phiV;M_{primary} (GeV/c^{2});#phi_{V}",2000,0,0.2,180,0,TMath::Pi());
}
//_____________________________________________________________________________
void StMiniTreeMaker::printConfig()
{
	const char *decision[2] = {"no","yes"};
	printf("=== Configuration for StMiniTreeMaker ===\n");
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
	printf("Maximum min nSigmaE for TPCe: %1.2f\n",mMinnSigmaE);
	printf("Maximum max nSigmaE  for TPCe: %1.2f\n",mMaxnSigmaE);
	printf("Maximum |1-1/beta| for TPCe: %1.3f\n",mMaxBeta2TOF);
	printf("=======================================\n");
}
double StMiniTreeMaker::phiVangle(TLorentzVector e1, TLorentzVector e2, int q1, int q2)
{
	Double_t pt1 = e1.Pt();
	Double_t eta1 = e1.Eta();
	Double_t phi1 = e1.Phi();

	Double_t pt2 = e2.Pt();
	Double_t eta2 = e2.Eta();
	Double_t phi2 = e2.Phi();

	TRandom3 *myRandom = new TRandom3();;
	TVector3 e1Mom,e2Mom;
	if(q1>0&&q2<0){
		e2Mom.SetPtEtaPhi(pt1,eta1,phi1);//e+
		e1Mom.SetPtEtaPhi(pt2,eta2,phi2);//e-
	}else if(q1<0&&q2>0){
		e2Mom.SetPtEtaPhi(pt2,eta2,phi2);//e+
		e1Mom.SetPtEtaPhi(pt1,eta1,phi1);//e-
	}else if(q1==q2&&TMath::Abs(q1)==1){
		Double_t ran = myRandom->Uniform(-1,1);
		if(ran>0){
			e2Mom.SetPtEtaPhi(pt1,eta1,phi1);
			e1Mom.SetPtEtaPhi(pt2,eta2,phi2);
		}
		else{
			e2Mom.SetPtEtaPhi(pt2,eta2,phi2);
			e1Mom.SetPtEtaPhi(pt1,eta1,phi1);
		}
	}else return -1;

	TVector3 pu=e1Mom+e2Mom;
	TVector3 pv=e1Mom.Cross(e2Mom);
	TVector3 pw=pu.Cross(pv);
	TVector3 pnz(0.,0.,-1);
	TVector3 pwc=pu.Cross(pnz);

	Double_t angleV = pw.Angle(pwc);

	return angleV;

}
bool StMiniTreeMaker::isElectron(StPicoTrack *pTrack){

	Int_t bTofPidTraitsIndex          = pTrack->bTofPidTraitsIndex();
	if(bTofPidTraitsIndex<0)return kFALSE;
	StPicoBTofPidTraits *btofPidTraits = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
	double beta  = btofPidTraits->btofBeta();
	if(beta<=0 || TMath::Abs(1.- 1./beta) > mMaxBeta2TOF)return kFALSE;
	StThreeVectorF pMom               = pTrack->pMom();
	double p = pMom.mag();
	Float_t mTpceNSigmaECutLow;
	if(p<1.){
		mTpceNSigmaECutLow = (mMinnSigmaE+2)/(1.- mMinnSigmaE)*(p - mMinnSigmaE) - 2;
	}else{
		mTpceNSigmaECutLow = mMinnSigmaE;
	}
	double nsigmaE = pTrack->nSigmaElectron();
	if(nsigmaE < mTpceNSigmaECutLow || nsigmaE >= mMaxnSigmaE )return kFALSE;
	return kTRUE;
}

