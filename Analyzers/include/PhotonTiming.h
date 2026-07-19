#ifndef PhotonTiming_h
#define PhotonTiming_h

#include "AnalyzerCore.h"

#include <vector>

class PhotonTiming : public AnalyzerCore {
public:
    PhotonTiming();
    ~PhotonTiming();

    void initializeAnalyzer() override;
    void executeEvent() override;

private:
    static constexpr int kMaxPho = 1000;
    static constexpr int kMaxPhoECALIdx = 1000000;
    static constexpr int kMaxEcalRechit = 1000000;

    void RequireBranch(const TString &branchName, void *address);
    void OptionalBranch(const TString &branchName, void *address);

    Int_t npho;
    Int_t nphoECALIdx;
    Int_t necalRechit;

    std::vector<Float_t> pho_pt;
    std::vector<Float_t> pho_eta;
    std::vector<Float_t> pho_phi;
    std::vector<Float_t> pho_seedRechitTime;
    std::vector<Float_t> pho_seedRechitEnergy;
    std::vector<Int_t> pho_seedRechitIdx;
    std::vector<Int_t> pho_ecalRechitIdxStart;
    std::vector<Int_t> pho_nEcalRechits;
    std::vector<UChar_t> pho_isOOT;
    std::vector<UChar_t> pho_passHLTEG20HEUnseeded;
    std::vector<UChar_t> pho_passHLTEG30HEUnseeded;

    std::vector<Int_t> phoECALIdx_EcalRechitIndex;

    std::vector<Float_t> ecalRechit_time;
    std::vector<Float_t> ecalRechit_energy;
    std::vector<UChar_t> ecalRechit_flagOOT;
    std::vector<UChar_t> ecalRechit_gainSwitch1;
    std::vector<UChar_t> ecalRechit_gainSwitch6;

    bool hasHLT20FilterBranch;
    bool hasHLT30FilterBranch;
};

#endif
