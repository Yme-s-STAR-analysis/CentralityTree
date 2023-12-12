#include "StCentTreeMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StThreeVectorF.hh"
#include "TFile.h"
#include "Stiostream.h"
#include <TMath.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include <fstream>
#include <vector>
#include <algorithm>
#include "TProfile.h"
#include "TTree.h"

ClassImp(StCentTreeMaker)

struct Event {

	Int_t RunId;
	Int_t TriggerId;

	Double_t Vz;
	Double_t ZDCx;

	Int_t refMult;
	Int_t refMult3;
	Int_t refMult3X;
	Int_t nTofMatch;
	Int_t nTofMatchZ;
	Int_t nTofBeta;
	Int_t tofMult;
};
Event event;

StCentTreeMaker::StCentTreeMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName) : StMaker(name) {
	mOutputName = outName;
	mPicoDstMaker = picoMaker;
	mPicoDst = 0;
}

StCentTreeMaker::~StCentTreeMaker() {}

// ------
Int_t StCentTreeMaker::Init() {
	mFileOut = new TFile(mOutputName, "recreate");
	stree = new TTree("centTree", "the reduced tree for centality defination");
	stree->Branch("TriggerId", &event.TriggerId, "TriggerId/I");
	stree->Branch("RunId", &event.RunId, "RunId/I");
	stree->Branch("Vz", &event.Vz, "Vz/D");
	stree->Branch("ZDCx", &event.ZDCx, "ZDCx/D");
	stree->Branch("refMult", &event.refMult, "refMult/I");
	stree->Branch("refMult3", &event.refMult3, "refMult3/I");
	stree->Branch("refMult3X", &event.refMult3X, "refMult3X/I");
	stree->Branch("nTofMatch", &event.nTofMatch, "nTofMatch/I");
	stree->Branch("nTofMatchZ", &event.nTofMatchZ, "nTofMatchZ/I");
	stree->Branch("nTofBeta", &event.nTofBeta, "nTofBeta/I");
	stree->Branch("tofMult", &event.tofMult, "tofMult/I");

	trg = new TriggerTool();
	// Mean DCA tool needs to be set up here
	// the parameters come from MakeMeanDcaCut
	// here the examples are from 14.6 GeV
	dcatool = new MeanDcaTool();
	dcatool->SetUpperCurveParZ(-0.0531713, 1.71842, 0.465328);
	dcatool->SetLowerCurveParZ(0.0532244, -1.72135, 0.464709);
	dcatool->SetUpperCurveParXY(0.045491, 2.14648, 0.558145);
	dcatool->SetLowerCurveParXY(-0.102939, -2.14641, 0.54303);

	return kStOK;
}

// ==== init shift util
void StCentTreeMaker::SetShiftTool(const char* fname, const char* h1name, const char* h2name) {
	shift = new TpcShiftTool();
	shift->Init(fname, h1name, h2name);
	cnter = new StCFMult();
	cnter->ImportShiftTool(shift);
}

//---------------------------------------------------------
Int_t StCentTreeMaker::Finish() {
	mFileOut->cd();
	stree->Write();
	mFileOut->Close();
	return kStOK;
}

void StCentTreeMaker::Clear(Option_t *opt) {}

//---------------------------------------------------------------
Int_t StCentTreeMaker::Make() {
	if (!mPicoDstMaker) {
		LOG_WARN << " No PicoDstMaker! Skip! " << endm;
		return kStWarn;
	}

	mPicoDst = mPicoDstMaker->picoDst();
	if (!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStWarn;
	}

	// Load event
	mPicoEvent = (StPicoEvent *)mPicoDst->event();
	if (!mPicoEvent) {
		cerr << "Error opening picoDst Event, skip!" << endl;
		return kStOk;
	}


	TVector3 pVtx = mPicoEvent->primaryVertex();
	Double_t vx = pVtx.X();
	Double_t vy = pVtx.Y();
	Double_t vz = pVtx.Z();

	if (TMath::Abs(vx) < 1.e-5 && TMath::Abs(vy) < 1.e-5 && TMath::Abs(vz) < 1.e-5) {
		return kStOk;
	}
	if (sqrt(vx * vx + vy * vy) >= 2.0) {
		return kStOk;
	}
	if (fabs(vz) > 70.) {
		return kStOk;
	}

	// others
	auto trgid = trg->GetTriggerID(mPicoEvent);
	if (trgid < 0) {
		return kStOK;
	}
	event.RunId = mPicoEvent->runId();
	event.Vz = vz;
	event.ZDCx = mPicoEvent->ZDCx();
	event.TriggerId = trgid;

	// do mean dca cut
	if (!dcatool->Make(mPicoDst)) {
		return kStOK;
	}

	if (dcatool->IsBadMeanDcaZEvent(mPicoDst)) {
		return kStOK;
	}
	if (dcatool->IsBadMeanDcaXYEvent(mPicoDst)) {
		return kStOK;
	}

	// get multiplicities
	cnter->make(mPicoDst);

	event.refMult = cnter->mRefMult;
	event.refMult3 = cnter->mRefMult3;
	event.refMult3X = cnter->mRefMult3X;
	event.nTofMatch = cnter->mNTofMatch;
	event.nTofMatchZ = cnter->mNTofMatchZ;
	event.nTofBeta = cnter->mNTofBeta;
	event.tofMult = cnter->mTofMult;

	stree->Fill();

	return kStOK;
}
