#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;
class StCentTreeMaker; 

StChain *chain;
void readPicoDst(const Char_t *inputFile = "file.list", TString JobIdName ="CentTree") {
    Int_t nEvents = 15000000;
    
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    
    gSystem->Load("StUtilities");    
    gSystem->Load("StPicoEvent");
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StCentTreeMaker");
    gSystem->Load("TpcShiftTool");
    gSystem->Load("MeanDcaTool");
    gSystem->Load("StCFMult");
    gSystem->Load("TriggerTool");
    
    chain = new StChain();
    
    TString Name = JobIdName ;  Name.Append(".root") ;
    StPicoDstMaker *picoMaker = new StPicoDstMaker(2, inputFile, "picoDst");
    StCentTreeMaker *anaMaker = new StCentTreeMaker("ana",picoMaker, Name);
    anaMaker->SetShiftTool(
        "/star/u/yghuang/Work/DataAnalysis/BES2/19p6/yqa/7ShiftFile/U.shift.root", 
        "ProtonShift1D", "ProtonShift2D"
    );

    chain->Init();
    cout << "chain->Init();" << endl;
    int total = picoMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;
    if (nEvents > total) nEvents = total;
    for (Int_t i = 0; i < nEvents; i++) {
        if (i % 10000 == 0)
            cout << "Working on eventNumber " << i << endl;
        
        chain->Clear();
        int iret = chain->Make(i);
        
        if (iret) {
            cout << "Bad return code!" << iret << endl;
            break;
        }
        total++;
    }
    
    chain->Finish();
    cout << ">>>>>>>> " << endl;
    cout << "All done! Total number of events  " << nEvents << endl;
    
    delete chain;
    
}
