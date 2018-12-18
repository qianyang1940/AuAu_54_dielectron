const Int_t mMax = 2000;
struct StEvtData
{
	// event information
	Int_t    mRunId;
	Int_t    mEventId;
	Char_t   mShouldHaveRejectEvent;
	Char_t   mNTrigs;
	Int_t    mTrigId[64];
	Short_t  mRefMult;
	Short_t  mGRefMult;
	Float_t  mGRefMultCorr;
	Float_t  mEvtWeight;
	Short_t  mCentrality;
	Float_t  mBBCRate;
	Float_t  mZDCRate;
	Float_t  mBField;
	Float_t  mVpdVz;
	Float_t  mVertexX;
    Float_t  mVertexY;
    Float_t  mVertexZ;
	Float_t  mVertexRanking;
	Float_t  mEtaPlusQx;
	Float_t  mEtaPlusQy;
	Float_t  mEtaPlusPtWeight;
	Short_t  mEtaPlusNTrks;
	Float_t  mEtaMinusQx;
	Float_t  mEtaMinusQy;
	Float_t  mEtaMinusPtWeight;
	Short_t  mEtaMinusNTrks;

	//track information
	Short_t  mNTrks;
	Short_t  mTrkId[mMax];
	Bool_t   mTPCeTrkFlag[mMax];
	Short_t  mBEMCTraitsIndex[mMax];
	Short_t  mMTDTraitsIndex[mMax];
	Int_t   mCharge[mMax];
	Float_t  mPt[mMax];
	Float_t  mEta[mMax];
	Float_t  mPhi[mMax];
	Float_t  mgPt[mMax];
	Float_t  mgEta[mMax];
	Float_t  mgPhi[mMax];
	Float_t  mgOriginX[mMax];
	Float_t  mgOriginY[mMax];
	Float_t  mgOriginZ[mMax];
	Int_t   mNHitsFit[mMax];
	Int_t   mNHitsPoss[mMax];
	Int_t   mNHitsDedx[mMax];
	Float_t  mDedx[mMax];
	Float_t  mDndx[mMax];
	Float_t  mDndxError[mMax];
	Float_t  mNSigmaE[mMax];
	Float_t  mDca[mMax];
	Float_t  mChi2[mMax];
	Float_t  mChi2Prob[mMax];
	Bool_t   mIsHFTTrk[mMax];
	Bool_t   mHasHFT4Layers[mMax];
	Int_t   mTOFMatchFlag[mMax];
	Float_t  mTOFLocalY[mMax];
	Float_t  mBeta2TOF[mMax];

};
