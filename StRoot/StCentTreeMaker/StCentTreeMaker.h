#ifndef STAR_StCentTreeMaker
#define STAR_StCentTreeMaker
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TString.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"

#include "StRoot/StCFMult/StCFMult.h"
#include "StRoot/TpcShiftTool/TpcShiftTool.h"
#include "StRoot/TriggerTool/TriggerTool.h"
#include "StRoot/MeanDcaTool/MeanDcaTool.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoDstMaker;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TTree;
class TH2D;
class TriggerTool;
class MeanDcaTool;
#include "TVector3.h"

class StCentTreeMaker : public StMaker
{
public:
	StCentTreeMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName = ".root");
	virtual ~StCentTreeMaker();

	virtual Int_t Init();
	virtual Int_t Make();
	virtual void Clear(Option_t *opt = "");
	virtual Int_t Finish();

private:
	StPicoDstMaker *mPicoDstMaker;
	StPicoDst *mPicoDst;
	StPicoEvent *mPicoEvent;
	StPicoTrack *picoTrack;

	TpcShiftTool* mtShift;
	StCFMult* mtMult;
	TriggerTool* mtTrg;
	MeanDcaTool* mtDca;

	TString mOutputName;
	TFile *mFileOut;
	TTree *stree;

	ClassDef(StCentTreeMaker, 1)
};

#endif
