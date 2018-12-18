TTree *eTree;

struct trkTree
{
	Int_t    EventId;

	Double_t bPt;
	Double_t bEta;
	Double_t bPhi;
	Double_t bY;

	Double_t bbarPt;
	Double_t bbarEta;
	Double_t bbarPhi;
	Double_t bbarY;

	Double_t eePt;
	Double_t eeEta;
	Double_t eePhi;
	Double_t eeRapidity;
	Double_t eeM;
	Int_t    ePosParentGID;
	Double_t ePosPt;
	Double_t ePosEta;
	Double_t ePosPhi;
	Int_t    eNegParentGID;
	Double_t eNegPt;
	Double_t eNegEta;
	Double_t eNegPhi;

}; trkTree meTree;
