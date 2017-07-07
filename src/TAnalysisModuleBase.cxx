#include "TAnalysisModuleBase.hxx"

ClassImp(CP::TAnalysisModuleBase);

CP::TAnalysisModuleBase::TAnalysisModuleBase():
    TNamed("Undefined: Problem in derived class!", 
           "Undefined: Problem in derived class!"),
    fOutputTree(NULL),
    fIsEnabled(kFALSE),
    fIsUsedForPreselection(kFALSE),
    fBufferSize(32000),
    fSplitLevel(99) { }

CP::TAnalysisModuleBase::~TAnalysisModuleBase() {}

Bool_t CP::TAnalysisModuleBase::Configure(const std::string &option) {
    CaptError("No configurable options defined for this class");
    return kFALSE;
}


std::string CP::TAnalysisModuleBase::GetDirectoryName() const {
    switch (GetTreeType()) {
    case kHeader: return "HeaderDir";
    case kTruth: return "TruthDir";
    case kRecon: return "ReconDir";
    default: throw CP::EUndefinedTreeType();
    }
}


void CP::TAnalysisModuleBase::Initialize(TTree *tree) {
    fOutputTree = tree;
    
    fOutputTree->SetAutoFlush(3000000);
    
    InitializeModule();
    fOutputTree->Branch("RunId", &fRunId, "RunId/I", fBufferSize);
    fOutputTree->Branch("EventId", &fEventId, "EventId/I", fBufferSize);
    fOutputTree->Branch("Preselected", &fPreselected, "Preselected/B", 
                        fBufferSize);
    fOutputTree->Branch("SubrunId", &fSubrunId, "SubrunId/I", fBufferSize);
    InitializeBranches();
}


bool CP::TAnalysisModuleBase::Process(CP::TEvent& event) {
    fRunId = event.GetRunId();
    fEventId = event.GetEventId();
    fSubrunId = event.GetContext().GetSubRun();
    if (!FillTree(event)) {
        throw EAnalysisFailure();
    }
    fPreselected = SaveFullEvent(event);
    fOutputTree->Fill();
    return fPreselected;
}

bool CP::TAnalysisModuleBase::SaveFullEvent(CP::TEvent& event){
    // Do not save the main event tree entry by default, just the analysis one
    return false; 
}

Bool_t CP::TAnalysisModuleBase::ProcessFirstEvent(CP::TEvent&) {
    return true;
}


void CP::TAnalysisModuleBase::PrintMessage() {
    // Override this function if more details need to be printed out
    CaptLog("Module: " << GetName() 
             << (IsEnabled() ? " (Enabled)  " : " (Disabled) ")
             << (IsUsedForPreselection() ? " (Preselecting)  "
                 : " (Not Preselecting) "));
}
