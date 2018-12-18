#include "TMath.h"

const Double_t Mmuon = 0.105658367;
const Double_t Mpion = 0.13957018;
const Double_t Mkaon = 0.493677;
const Double_t Mproton = 0.93827231;
const Double_t Melectron = 0.00051099907;

const Int_t    mTotalCentrality = 16;
const Int_t    mTotalRun = 1654;
const Int_t    mTotalDay = 93;
const Int_t    mArrayLength = 20;

const Int_t    mMinRunId = 15107008;

//event cuts
const Double_t mVzCut = 6.;
const Double_t mVzDiffCut = 3.;
const Double_t mVrCut = 2.;

//electron cuts, using TPC+TOF to do the EID
const Double_t mTpcePtCut[2] = {0.2, 30.};
const Double_t mTpceEtaCut = 1.;
const Int_t    mTpceNHitsFitCut = 20;
const Double_t mTpceNHitsFitRatioCut = 0.52;
const Int_t    mTpceNHitsDedxCut = 15;
const Double_t mTpceDcaCutWoHFT = 1.5;
const Double_t mMinTpceDcaCutWHFT = 0.004;
const Double_t mMaxTpceDcaCutWHFT = 0.1;
const Double_t mTpceBeta2TOFCut = 0.03;
const Double_t mTpceTOFLocalYCut = 1.8;
const Double_t mTpceNSigmaECut[2] = {-1., 2.}; //assume the mean of electron nSigmaE is 0
const Double_t mNSigmaEShift = -0.3; //the mean of electron nsigmaE shift to -0.34

//gamma cut
const Double_t mGammaMassCut = 0.005;
const Double_t mGammaAngleCut = TMath::Pi()/10;
const Double_t mGammaDcaDaughtersCut = 0.5;
const Double_t mMinGammaDecayLengthCut = 2;
const Double_t mMaxGammaDecayLengthCut = 20.;
const Double_t mGammaDcaCut = 0.5;

//muon cuts, MTD constants and MuID
const Int_t    mMuNHitsFitCut = 20;
const Double_t mMuNHitsFitRatioCut = 0.52;
const Int_t    mMuNHitsDedxCut = 15;
const Double_t mLowMuPtCut = 1.2; 
const Double_t mMuEtaCut = 0.8;
const Double_t mMinMuDcaCutWHFT = 0.002;
const Double_t mMaxMuDcaCutWHFT = 0.1;
const Double_t mMuDcaCutWoHFT = 1.0;
const Double_t mMuNSigmaPiCut[2] = {-1., 3.};
const Double_t mMuNSigma[2] = {2.5, 3}; //pT<3 use 2.; pT>=3 use 2.5
const Double_t mMuDeltaTCut = 0.75;

const Double_t mEMuDcaCut = 0.004;

const Double_t mPairYCut = 1.;
