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
	Int_t refMult3S;
	Int_t refMult3E;
	Int_t nTofMatch;
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
	stree->Branch("refMult3S", &event.refMult3S, "refMult3S/I");
	stree->Branch("refMult3E", &event.refMult3E, "refMult3E/I");
	stree->Branch("nTofMatch", &event.nTofMatch, "nTofMatch/I");
	stree->Branch("nTofBeta", &event.nTofBeta, "nTofBeta/I");
	stree->Branch("tofMult", &event.tofMult, "tofMult/I");

	mtTrg = new TriggerTool();

	mtDca = new MeanDcaTool();
	mtDca->ReadParams();

	mtShift = new TpcShiftTool();
	mtShift->Init();

	mtMult = new StCFMult();
	mtMult->ImportShiftTool(mtShift);
	// in centrality tree, RefMult3E just uses 1.0
	// we can scan it later and find the proper efficiency value
	mtMult->SetEfficiencyRefMult3E(1.0); 

	mtVtx = new VtxShiftTool();

	return kStOK;
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

	if (fabs(vx) < 1.e-5 && 
		fabs(vy) < 1.e-5 &&
		fabs(vz) < 1.e-5) {
		return kStOK;
	}

	if (fabs(vz) > 50) { return kStOK; }
	auto vr = mtVtx->GetShiftedVr(vx, vy);
	if (vr > 1) { return kStOK; }

	// others
	auto trgid = mtTrg->GetTriggerID(mPicoEvent);
	if (trgid < 0) { return kStOK; }
	event.RunId = mPicoEvent->runId();
	event.Vz = vz;
	event.ZDCx = mPicoEvent->ZDCx();
	event.TriggerId = trgid;

	// do mean dca cut
	if (!mtDca->Make(mPicoDst)) { return kStOK; }

	if (mtDca->IsBadMeanDcaZEvent(mPicoDst)) { return kStOK; }
	if (mtDca->IsBadMeanDcaXYEvent(mPicoDst)) { return kStOK; }

	// get multiplicities
	mtMult->make(mPicoDst);

	event.refMult = mtMult->mRefMult;
	event.refMult3 = mtMult->mRefMult3;
	event.refMult3X = mtMult->mRefMult3X;
	event.refMult3S = mtMult->mRefMult3S;
	event.refMult3E = mtMult->mRefMult3E;
	event.nTofMatch = mtMult->mNTofMatch;
	event.nTofBeta = mtMult->mNTofBeta;
	event.tofMult = mtMult->mTofMult;

	stree->Fill();

	return kStOK;
}
