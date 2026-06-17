#ifndef HLTBranchCheck_h
#define HLTBranchCheck_h

#include "AnalyzerCore.h"

class HLTBranchCheck : public AnalyzerCore {
public:
    HLTBranchCheck();
    ~HLTBranchCheck();

    void initializeAnalyzer() override;
    void executeEvent() override;

private:
    TString targetTrigger;
    Bool_t passTargetTrigger;
    bool hasTargetTrigger;
};

#endif
