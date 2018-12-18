const Float_t Mmuon = 0.105658367;
const Float_t Mpion = 0.13957018;
const Float_t Mkaon = 0.493677;
const Float_t Mproton = 0.93827231;
const Float_t Melectron = 0.00051099907;
const Float_t c_light = 29.9792458; //cm/ns
const Float_t v_signal = 56.; //ps/cm,MTD signal velocity 

const Int_t   mTotalRun = 580;
const Int_t   mTotalDay = 20;
const Int_t   mTotalCentrality = 16; //default value is 16;
const Int_t   mArrayLength = 20;

//event cuts
const Float_t mVrCut = 2.;
const Float_t mVzCut = 30.;
const Float_t mVzDiffCut = 3.;

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
const Float_t mNSigmaEShift = -0.34; //the mean of electron nsigmaE shift to -0.34

//pair cuts
const Float_t mPairYCut = 1.;
const Float_t mPheMassCutWoTof = 0.005;//photonic electron rejection
const Float_t mPheMassCutWTof = 0.015;//photonic electron rejection
