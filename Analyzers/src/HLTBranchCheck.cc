#include "HLTBranchCheck.h"

HLTBranchCheck::HLTBranchCheck()
    : targetTrigger("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350"),
      passTargetTrigger(false),
      hasTargetTrigger(false) {}

HLTBranchCheck::~HLTBranchCheck() {}

void HLTBranchCheck::initializeAnalyzer() {
    if (!fChain) {
        cerr << "[HLTBranchCheck::initializeAnalyzer] fChain is null" << endl;
        exit(EXIT_FAILURE);
    }

    hasTargetTrigger = (fChain->GetBranch(targetTrigger) != nullptr);
    cout << "[HLTBranchCheck::initializeAnalyzer] targetTrigger = " << targetTrigger << endl;
    cout << "[HLTBranchCheck::initializeAnalyzer] hasTargetTrigger = " << hasTargetTrigger << endl;

    if (!hasTargetTrigger) {
        cerr << "[HLTBranchCheck::initializeAnalyzer] Branch not found: " << targetTrigger << endl;
        exit(EXIT_FAILURE);
    }

    // This minimal check only needs the target HLT bit.  EXONanoAOD files also
    // contain CMSSW/EDM metadata branches without dictionaries in the plain LCG
    // ROOT environment; leaving all branches active can crash at file boundaries.
    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus(targetTrigger, 1);
    fChain->SetBranchAddress(targetTrigger, &passTargetTrigger);
}

void HLTBranchCheck::executeEvent() {
    const int decision = passTargetTrigger ? 1 : 0;
    FillHist("HLT/Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350", decision, 1.0, 2, -0.5, 1.5);
}
