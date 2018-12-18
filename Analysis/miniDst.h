//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 10 23:38:47 2018 by ROOT version 5.34/30
// from TTree miniDst/miniDst
// found on file: /star/u/tc88qy/AuAu/run17/54GeV/Sample/picoDst/test/test.root
//////////////////////////////////////////////////////////

#ifndef miniDst_h
#define miniDst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class miniDst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           mRunId;
   Int_t           mEventId;
   Char_t          mShouldHaveRejectEvent;
   Int_t           mNTrigs;
   Int_t           mTrigId[3];   //[mNTrigs]
   Short_t         mRefMult;
   Short_t         mGRefMult;
   Float_t         mGRefMultCorr;
   Float_t         mEvtWeight;
   Short_t         mCentrality;
   Float_t         mBBCRate;
   Float_t         mZDCRate;
   Float_t         mBField;
   Float_t         mVpdVz;
   Float_t         mVertexX;
   Float_t         mVertexY;
   Float_t         mVertexZ;
   Float_t         mVertexRanking;
   Float_t         mEtaPlusQx;
   Float_t         mEtaPlusQy;
   Float_t         mEtaPlusPtWeight;
   Short_t         mEtaPlusNTrks;
   Float_t         mEtaMinusQx;
   Float_t         mEtaMinusQy;
   Float_t         mEtaMinusPtWeight;
   Short_t         mEtaMinusNTrks;
   Short_t         mNTrks;
   Short_t         mTrkId[100];   //[mNTrks]
   Bool_t          mTPCeTrkFlag[100];   //[mNTrks]
   Int_t           mCharge[100];   //[mNTrks]
   Float_t         mPt[100];   //[mNTrks]
   Float_t         mEta[100];   //[mNTrks]
   Float_t         mPhi[100];   //[mNTrks]
   Float_t         mgPt[100];   //[mNTrks]
   Float_t         mgEta[100];   //[mNTrks]
   Float_t         mgPhi[100];   //[mNTrks]
   Float_t         mgOriginX[100];   //[mNTrks]
   Float_t         mgOriginY[100];   //[mNTrks]
   Float_t         mgOriginZ[100];   //[mNTrks]
   Int_t           mNHitsFit[100];   //[mNTrks]
   Int_t           mNHitsPoss[100];   //[mNTrks]
   Int_t           mNHitsDedx[100];   //[mNTrks]
   Float_t         mDedx[100];   //[mNTrks]
   Float_t         mNSigmaE[100];   //[mNTrks]
   Float_t         mDca[100];   //[mNTrks]
   Int_t           mTOFMatchFlag[100];   //[mNTrks]
   Float_t         mTOFLocalY[100];   //[mNTrks]
   Float_t         mBeta2TOF[100];   //[mNTrks]

   // List of branches
   TBranch        *b_mRunId;   //!
   TBranch        *b_mEventId;   //!
   TBranch        *b_mShouldHaveRejectEvent;   //!
   TBranch        *b_mNTrigs;   //!
   TBranch        *b_mTrigId;   //!
   TBranch        *b_mRefMult;   //!
   TBranch        *b_mGRefMult;   //!
   TBranch        *b_mGRefMultCorr;   //!
   TBranch        *b_mEvtWeight;   //!
   TBranch        *b_mCentrality;   //!
   TBranch        *b_mBBCRate;   //!
   TBranch        *b_mZDCRate;   //!
   TBranch        *b_mBField;   //!
   TBranch        *b_mVpdVz;   //!
   TBranch        *b_mVertexX;   //!
   TBranch        *b_mVertexY;   //!
   TBranch        *b_mVertexZ;   //!
   TBranch        *b_mVertexRanking;   //!
   TBranch        *b_mEtaPlusQx;   //!
   TBranch        *b_mEtaPlusQy;   //!
   TBranch        *b_mEtaPlusPtWeight;   //!
   TBranch        *b_mEtaPlusNTrks;   //!
   TBranch        *b_mEtaMinusQx;   //!
   TBranch        *b_mEtaMinusQy;   //!
   TBranch        *b_mEtaMinusPtWeight;   //!
   TBranch        *b_mEtaMinusNTrks;   //!
   TBranch        *b_mNTrks;   //!
   TBranch        *b_mTrkId;   //!
   TBranch        *b_mTPCeTrkFlag;   //!
   TBranch        *b_mCharge;   //!
   TBranch        *b_mPt;   //!
   TBranch        *b_mEta;   //!
   TBranch        *b_mPhi;   //!
   TBranch        *b_mgPt;   //!
   TBranch        *b_mgEta;   //!
   TBranch        *b_mgPhi;   //!
   TBranch        *b_mgOriginX;   //!
   TBranch        *b_mgOriginY;   //!
   TBranch        *b_mgOriginZ;   //!
   TBranch        *b_mNHitsFit;   //!
   TBranch        *b_mNHitsPoss;   //!
   TBranch        *b_mNHitsDedx;   //!
   TBranch        *b_mDedx;   //!
   TBranch        *b_mNSigmaE;   //!
   TBranch        *b_mDca;   //!
   TBranch        *b_mTOFMatchFlag;   //!
   TBranch        *b_mTOFLocalY;   //!
   TBranch        *b_mBeta2TOF;   //!

   miniDst(TTree *tree=0);
   virtual ~miniDst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef miniDst_cxx
miniDst::miniDst(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/star/u/tc88qy/AuAu/run17/54GeV/Sample/picoDst/test/test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/star/u/tc88qy/AuAu/run17/54GeV/Sample/picoDst/test/test.root");
      }
      f->GetObject("miniDst",tree);

   }
   Init(tree);
}

miniDst::~miniDst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t miniDst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t miniDst::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void miniDst::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mRunId", &mRunId, &b_mRunId);
   fChain->SetBranchAddress("mEventId", &mEventId, &b_mEventId);
   fChain->SetBranchAddress("mShouldHaveRejectEvent", &mShouldHaveRejectEvent, &b_mShouldHaveRejectEvent);
   fChain->SetBranchAddress("mNTrigs", &mNTrigs, &b_mNTrigs);
   fChain->SetBranchAddress("mTrigId", mTrigId, &b_mTrigId);
   fChain->SetBranchAddress("mRefMult", &mRefMult, &b_mRefMult);
   fChain->SetBranchAddress("mGRefMult", &mGRefMult, &b_mGRefMult);
   fChain->SetBranchAddress("mGRefMultCorr", &mGRefMultCorr, &b_mGRefMultCorr);
   fChain->SetBranchAddress("mEvtWeight", &mEvtWeight, &b_mEvtWeight);
   fChain->SetBranchAddress("mCentrality", &mCentrality, &b_mCentrality);
   fChain->SetBranchAddress("mBBCRate", &mBBCRate, &b_mBBCRate);
   fChain->SetBranchAddress("mZDCRate", &mZDCRate, &b_mZDCRate);
   fChain->SetBranchAddress("mBField", &mBField, &b_mBField);
   fChain->SetBranchAddress("mVpdVz", &mVpdVz, &b_mVpdVz);
   fChain->SetBranchAddress("mVertexX", &mVertexX, &b_mVertexX);
   fChain->SetBranchAddress("mVertexY", &mVertexY, &b_mVertexY);
   fChain->SetBranchAddress("mVertexZ", &mVertexZ, &b_mVertexZ);
   fChain->SetBranchAddress("mVertexRanking", &mVertexRanking, &b_mVertexRanking);
   fChain->SetBranchAddress("mEtaPlusQx", &mEtaPlusQx, &b_mEtaPlusQx);
   fChain->SetBranchAddress("mEtaPlusQy", &mEtaPlusQy, &b_mEtaPlusQy);
   fChain->SetBranchAddress("mEtaPlusPtWeight", &mEtaPlusPtWeight, &b_mEtaPlusPtWeight);
   fChain->SetBranchAddress("mEtaPlusNTrks", &mEtaPlusNTrks, &b_mEtaPlusNTrks);
   fChain->SetBranchAddress("mEtaMinusQx", &mEtaMinusQx, &b_mEtaMinusQx);
   fChain->SetBranchAddress("mEtaMinusQy", &mEtaMinusQy, &b_mEtaMinusQy);
   fChain->SetBranchAddress("mEtaMinusPtWeight", &mEtaMinusPtWeight, &b_mEtaMinusPtWeight);
   fChain->SetBranchAddress("mEtaMinusNTrks", &mEtaMinusNTrks, &b_mEtaMinusNTrks);
   fChain->SetBranchAddress("mNTrks", &mNTrks, &b_mNTrks);
   fChain->SetBranchAddress("mTrkId", mTrkId, &b_mTrkId);
   fChain->SetBranchAddress("mTPCeTrkFlag", mTPCeTrkFlag, &b_mTPCeTrkFlag);
   fChain->SetBranchAddress("mCharge", mCharge, &b_mCharge);
   fChain->SetBranchAddress("mPt", mPt, &b_mPt);
   fChain->SetBranchAddress("mEta", mEta, &b_mEta);
   fChain->SetBranchAddress("mPhi", mPhi, &b_mPhi);
   fChain->SetBranchAddress("mgPt", mgPt, &b_mgPt);
   fChain->SetBranchAddress("mgEta", mgEta, &b_mgEta);
   fChain->SetBranchAddress("mgPhi", mgPhi, &b_mgPhi);
   fChain->SetBranchAddress("mgOriginX", mgOriginX, &b_mgOriginX);
   fChain->SetBranchAddress("mgOriginY", mgOriginY, &b_mgOriginY);
   fChain->SetBranchAddress("mgOriginZ", mgOriginZ, &b_mgOriginZ);
   fChain->SetBranchAddress("mNHitsFit", mNHitsFit, &b_mNHitsFit);
   fChain->SetBranchAddress("mNHitsPoss", mNHitsPoss, &b_mNHitsPoss);
   fChain->SetBranchAddress("mNHitsDedx", mNHitsDedx, &b_mNHitsDedx);
   fChain->SetBranchAddress("mDedx", mDedx, &b_mDedx);
   fChain->SetBranchAddress("mNSigmaE", mNSigmaE, &b_mNSigmaE);
   fChain->SetBranchAddress("mDca", mDca, &b_mDca);
   fChain->SetBranchAddress("mTOFMatchFlag", mTOFMatchFlag, &b_mTOFMatchFlag);
   fChain->SetBranchAddress("mTOFLocalY", mTOFLocalY, &b_mTOFLocalY);
   fChain->SetBranchAddress("mBeta2TOF", mBeta2TOF, &b_mBeta2TOF);
   Notify();
}

Bool_t miniDst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void miniDst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t miniDst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef miniDst_cxx
