#ifndef LLPMuonSystemAcceptance_h
#define LLPMuonSystemAcceptance_h

#include "AnalyzerCore.h"

#include <string>

class LLPMuonSystemAcceptance : public AnalyzerCore {
public:
    LLPMuonSystemAcceptance();
    ~LLPMuonSystemAcceptance();

    void initializeAnalyzer();
    void executeEvent();
    void WriteHist() override;

private:
    static constexpr int LLP_PDGID = 1000022;
    static constexpr int DECAY_PHOTON_PDGID = 22;
    static constexpr float MAX_LLP_ABS_ETA = 2.4;
    static constexpr float MIN_ABS_Z_CM = 400.0;
    static constexpr float MAX_ABS_Z_CM = 1100.0;
    static constexpr float MAX_R_CM = 695.5;
    static constexpr float MAX_CSC_MATCH_DR = 0.3;
    static constexpr int TARGET_RUN = 1;
    static constexpr int TARGET_LUMI = 21;
    static constexpr Long_t TARGET_EVENT = 20525;

    Int_t ncscMDSCluster = 0;
    RVec<Float_t> cscMDSCluster_eta;
    RVec<Float_t> cscMDSCluster_phi;
    Bool_t HLT_CscCluster50_Photon20Unseeded = false;

    TTree *matchedEventTree = nullptr;
    Int_t matchedRun = 0;
    Int_t matchedLumi = 0;
    Long64_t matchedEvent = 0;
    Int_t matchedPhotonIndex = -1;
    Int_t matchedClusterIndex = -1;
    Bool_t matchedPhotonMatchedToCSC = false;
    Bool_t matchedHLT_CscCluster50_Photon20Unseeded = false;
    Float_t matchedDeltaR = -1.;
    Float_t matchedPhotonPt = -1.;
    Float_t matchedPhotonEta = -999.;
    Float_t matchedPhotonPhi = -999.;
    Float_t matchedPhotonVx = -999.;
    Float_t matchedPhotonVy = -999.;
    Float_t matchedPhotonVz = -999.;
    Float_t matchedPhotonDecayR = -1.;
    Float_t matchedPhotonDecayLxyz = -1.;
    Bool_t matchedPhotonPassMuonSystemGeometry = false;
    std::string matchedInputFile;

    TTree *targetPhotonTree = nullptr;
    Bool_t targetEventPrinted = false;
    Int_t targetRun = 0;
    Int_t targetLumi = 0;
    Long64_t targetEvent = 0;
    Int_t targetPhotonIndex = -1;
    Int_t targetPhotonMotherIndex = -1;
    Float_t targetPhotonPt = -1.;
    Float_t targetPhotonEta = -999.;
    Float_t targetPhotonPhi = -999.;
    Float_t targetPhotonVx = -999.;
    Float_t targetPhotonVy = -999.;
    Float_t targetPhotonVz = -999.;
    std::string targetPhotonInputFile;

    bool PassMuonSystemGeometry(int gen_index) const;
    int MotherIndex(int gen_index) const;
    float DeltaR(float eta1, float phi1, float eta2, float phi2) const;
};

#endif
