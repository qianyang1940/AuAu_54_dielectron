//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 27 15:06:36 2017 by ROOT version 5.34/30
// from TTree meTree/eTree
// found on file: /star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/heavyflavor/cc2ee/genTree/output/pythiaevent36.root
//////////////////////////////////////////////////////////

#ifndef meTree_h
#define meTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class meTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ncString;
   Int_t           ncbarString;
   Int_t           nc;
   Bool_t          isStringC[100];   //[nc]
   Double_t        cPt[100];   //[nc]
   Double_t        cEta[100];   //[nc]
   Double_t        cPhi[100];   //[nc]
   Double_t        cY[100];   //[nc]
   Int_t           ncbar;
   Bool_t          isStringCbar[100];   //[ncbar]
   Double_t        cbarPt[100];   //[ncbar]
   Double_t        cbarEta[100];   //[ncbar]
   Double_t        cbarPhi[100];   //[ncbar]
   Double_t        cbarY[100];   //[ncbar]
   Int_t           nePos;
   Int_t           ePosParentGID[100];   //[nePos]
   Double_t        ePosPt[100];   //[nePos]
   Double_t        ePosEta[100];   //[nePos]
   Double_t        ePosPhi[100];   //[nePos]
   Int_t           neNeg;
   Int_t           eNegParentGID[100];   //[neNeg]
   Double_t        eNegPt[100];   //[neNeg]
   Double_t        eNegEta[100];   //[neNeg]
   Double_t        eNegPhi[100];   //[neNeg]

   // List of branches
   TBranch        *b_ncString;   //!
   TBranch        *b_ncbarString;   //!
   TBranch        *b_nc;   //!
   TBranch        *b_isStringC;   //!
   TBranch        *b_cPt;   //!
   TBranch        *b_cEta;   //!
   TBranch        *b_cPhi;   //!
   TBranch        *b_cY;   //!
   TBranch        *b_ncbar;   //!
   TBranch        *b_isStringCbar;   //!
   TBranch        *b_cbarPt;   //!
   TBranch        *b_cbarEta;   //!
   TBranch        *b_cbarPhi;   //!
   TBranch        *b_cbarY;   //!
   TBranch        *b_nePos;   //!
   TBranch        *b_ePosParentGID;   //!
   TBranch        *b_ePosPt;   //!
   TBranch        *b_ePosEta;   //!
   TBranch        *b_ePosPhi;   //!
   TBranch        *b_neNeg;   //!
   TBranch        *b_eNegParentGID;   //!
   TBranch        *b_eNegPt;   //!
   TBranch        *b_eNegEta;   //!
   TBranch        *b_eNegPhi;   //!

   meTree(TTree *tree=0);
   virtual ~meTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef meTree_cxx
meTree::meTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/heavyflavor/cc2ee/genTree/output/pythiaevent36.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/heavyflavor/cc2ee/genTree/output/pythiaevent36.root");
      }
      f->GetObject("meTree",tree);

   }
   Init(tree);
}

meTree::~meTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t meTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t meTree::LoadTree(Long64_t entry)
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

void meTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("ncString", &ncString, &b_ncString);
   fChain->SetBranchAddress("ncbarString", &ncbarString, &b_ncbarString);
   fChain->SetBranchAddress("nc", &nc, &b_nc);
   fChain->SetBranchAddress("isStringC", isStringC, &b_isStringC);
   fChain->SetBranchAddress("cPt", cPt, &b_cPt);
   fChain->SetBranchAddress("cEta", cEta, &b_cEta);
   fChain->SetBranchAddress("cPhi", cPhi, &b_cPhi);
   fChain->SetBranchAddress("cY", cY, &b_cY);
   fChain->SetBranchAddress("ncbar", &ncbar, &b_ncbar);
   fChain->SetBranchAddress("isStringCbar", isStringCbar, &b_isStringCbar);
   fChain->SetBranchAddress("cbarPt", cbarPt, &b_cbarPt);
   fChain->SetBranchAddress("cbarEta", cbarEta, &b_cbarEta);
   fChain->SetBranchAddress("cbarPhi", cbarPhi, &b_cbarPhi);
   fChain->SetBranchAddress("cbarY", cbarY, &b_cbarY);
   fChain->SetBranchAddress("nePos", &nePos, &b_nePos);
   fChain->SetBranchAddress("ePosParentGID", ePosParentGID, &b_ePosParentGID);
   fChain->SetBranchAddress("ePosPt", ePosPt, &b_ePosPt);
   fChain->SetBranchAddress("ePosEta", ePosEta, &b_ePosEta);
   fChain->SetBranchAddress("ePosPhi", ePosPhi, &b_ePosPhi);
   fChain->SetBranchAddress("neNeg", &neNeg, &b_neNeg);
   fChain->SetBranchAddress("eNegParentGID", eNegParentGID, &b_eNegParentGID);
   fChain->SetBranchAddress("eNegPt", eNegPt, &b_eNegPt);
   fChain->SetBranchAddress("eNegEta", eNegEta, &b_eNegEta);
   fChain->SetBranchAddress("eNegPhi", eNegPhi, &b_eNegPhi);
   Notify();
}

Bool_t meTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void meTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t meTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef meTree_cxx
