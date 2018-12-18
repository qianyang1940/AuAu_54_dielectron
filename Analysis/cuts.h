const Float_t Mmuon = 0.105658367;
const Float_t Mpion = 0.13957018;
const Float_t Mkaon = 0.493677;
const Float_t Mproton = 0.93827231;
const Float_t Melectron = 0.00051099907;
const Float_t c_light = 29.9792458; //cm/ns
const Float_t v_signal = 56.; //ps/cm,MTD signal velocity 

const Int_t   mTotalRun = 580;
const Int_t   mTotalDay = 20;
const Int_t   mTotalCentrality = 9; //default value is 16;
const Int_t   mArrayLength = 20;

const Int_t    mMinRunId = 18153052;

//event cuts
const Float_t mVrCut = 2.;
const Float_t mVzCut = 60.;
const Float_t mVzDiffCut = 3.;
const Int_t   mTrigId[3] ={580001, 580011, 580021}; 

//cuts for miniTree production and eventPlane calculation
const Int_t   mNHitsFitCut = 20;
const Int_t   mNHitsDedxCut = 15;
const Float_t mNHitsFitRatioCut[2] = {0.52, 1.05};
const Float_t mPtCut[2] = {0.2, 2.};
const Float_t mEtaCut = 1.;
const Float_t mDcaCut  = 2.;
const Float_t mBeta2TOFCut = 0.05;
const Float_t mNSigmaECut = 3.;

//electron cuts, using TPC+TOF to do the EID
const Int_t   mTpceNHitsFitCut = 20;
const Int_t   mTpceNHitsDedxCut = 15;
const Float_t mTpceNHitsFitRatioCut = 0.52;
const Float_t mTpcePtCut[2] = {0.2, 30.};
const Float_t mTpceEtaCut = 1.;
const Float_t mTpceDcaCut = 1.;
const Float_t mTpceBeta2TOFCut = 0.025;
const Float_t mTpceNSigmaECut[2] = {-0.75, 2.0}; //assume the mean of electron nSigmaE is 0
const Float_t mNSigmaEShift = 0.; //the mean of electron nsigmaE shift to -0.34

//electron cuts, using TPC+BEMC to do the EID
const Int_t   mEmceNHitsFitCut = 20;
const Int_t   mEmceNHitsDedxCut = 15;
const Float_t mEmceNHitsFitRatioCut = 0.52;
const Float_t mEmcePtCut[2] = {2.5, 30.0};
const Float_t mEmceEtaCut = 1.;
const Float_t mEmceDcaCut = 1.;
const Float_t mEmceBeta2TOFCut = 0.025;
const Float_t mEmceNSigmaECut[2] = {-0.75, 2.}; //assume the mean of electron nsigmaE is 0
const Float_t mEmcePECut[2] = {0.3, 1.5};
const Int_t   mEmceAdcCut = 50;
const Int_t   mEmceNEtaCut = 1;
const Int_t   mEmceNPhiCut = 1;
const Float_t mEmcePhiDistCut = 0.04;
const Float_t mEmceZDistCut = 12.;

//muon cuts, MTD constants and MuID
const Int_t   mBackLegs = 30;
const Int_t   mTrays = 5;
const Int_t   mCells = 12;
const Int_t   mChannels = 24;
const Float_t mStripLength = 87.0; //cm
const Float_t mStripWidth  = 3.8 ; //cm
const Float_t mStripGap    = 0.6 ; //cm
const Int_t   mMuNHitsFitCut = 20;
const Int_t   mMuNHitsDedxCut = 15;
const Float_t mMuNHitsFitRatioCut = 0.52;
const Float_t mMuPtCut[2] = {2., 30.};
const Float_t mMuEtaCut = 0.65;
const Float_t mMuDcaCut = 2.0;
const Float_t mMuNSigmaPiCut[2] = {-1., 3.};
const Float_t mMuDeltaZCut = 20.;
const Float_t mMuTofDiffCut[2] = {-0.35, 0.1};
const Float_t mMuBetaCut[2] = {0.97, 1.02};

//pair cuts
const Float_t mPairYCut = 1.;
const Float_t mPheMassCutWoTof = 0.005;//photonic electron rejection
const Float_t mPheMassCutWTof = 0.015;//photonic electron rejection
