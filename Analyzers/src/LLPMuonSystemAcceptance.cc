#include "LLPMuonSystemAcceptance.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <vector>

LLPMuonSystemAcceptance::LLPMuonSystemAcceptance() {}
LLPMuonSystemAcceptance::~LLPMuonSystemAcceptance() {}

void LLPMuonSystemAcceptance::initializeAnalyzer() {
    cout << "[LLPMuonSystemAcceptance::initializeAnalyzer] LLP PDG ID = " << LLP_PDGID << endl;
    cout << "[LLPMuonSystemAcceptance::initializeAnalyzer] Decay daughter PDG ID = " << DECAY_PHOTON_PDGID << endl;
    cout << "[LLPMuonSystemAcceptance::initializeAnalyzer] Muon-system geometry: |eta_LLP| < "
         << MAX_LLP_ABS_ETA << ", " << MIN_ABS_Z_CM << " < |z_decay| < "
         << MAX_ABS_Z_CM << " cm, r_decay < " << MAX_R_CM << " cm" << endl;

    if (!fChain->GetBranch("GenPart_vx") || !fChain->GetBranch("GenPart_vy") || !fChain->GetBranch("GenPart_vz")) {
        cerr << "[LLPMuonSystemAcceptance::initializeAnalyzer] Missing GenPart_vx/vy/vz branches. "
             << "This analyzer requires EXONanoAOD generator vertex branches." << endl;
        exit(EXIT_FAILURE);
    }
    if (!fChain->GetBranch("ncscMDSCluster") || !fChain->GetBranch("cscMDSCluster_eta") ||
        !fChain->GetBranch("cscMDSCluster_phi")) {
        cerr << "[LLPMuonSystemAcceptance::initializeAnalyzer] Missing cscMDSCluster branches. "
             << "This analyzer requires EXONanoAOD CSC MDS cluster branches." << endl;
        exit(EXIT_FAILURE);
    }

    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("nGenPart", 1);
    fChain->SetBranchStatus("GenPart_eta", 1);
    fChain->SetBranchStatus("GenPart_phi", 1);
    fChain->SetBranchStatus("GenPart_pt", 1);
    fChain->SetBranchStatus("GenPart_pdgId", 1);
    fChain->SetBranchStatus("GenPart_genPartIdxMother", 1);
    fChain->SetBranchStatus("GenPart_vx", 1);
    fChain->SetBranchStatus("GenPart_vy", 1);
    fChain->SetBranchStatus("GenPart_vz", 1);
    fChain->SetBranchStatus("ncscMDSCluster", 1);
    fChain->SetBranchStatus("cscMDSCluster_eta", 1);
    fChain->SetBranchStatus("cscMDSCluster_phi", 1);
    fChain->SetBranchStatus("run", 1);
    fChain->SetBranchStatus("luminosityBlock", 1);
    fChain->SetBranchStatus("event", 1);

    const int max_csc_clusters = std::max(1, static_cast<int>(std::ceil(fChain->GetMaximum("ncscMDSCluster"))));
    cscMDSCluster_eta.resize(max_csc_clusters);
    cscMDSCluster_phi.resize(max_csc_clusters);
    fChain->SetBranchAddress("ncscMDSCluster", &ncscMDSCluster);
    fChain->SetBranchAddress("cscMDSCluster_eta", cscMDSCluster_eta.data());
    fChain->SetBranchAddress("cscMDSCluster_phi", cscMDSCluster_phi.data());

    if (fChain->GetBranch("HLT_CscCluster50_Photon20Unseeded")) {
        fChain->SetBranchStatus("HLT_CscCluster50_Photon20Unseeded", 1);
        fChain->SetBranchAddress("HLT_CscCluster50_Photon20Unseeded", &HLT_CscCluster50_Photon20Unseeded);
    } else {
        cerr << "[LLPMuonSystemAcceptance::initializeAnalyzer] Missing HLT_CscCluster50_Photon20Unseeded branch. "
             << "HLT-split CSC matching histograms will be filled as HLT-fail only." << endl;
    }

    GetOutfile()->cd();
    matchedEventTree = new TTree("MatchedCSCPhotonEvents", "Events with LLP photon matched to CSC MDS cluster");
    matchedEventTree->Branch("run", &matchedRun, "run/I");
    matchedEventTree->Branch("lumi", &matchedLumi, "lumi/I");
    matchedEventTree->Branch("event", &matchedEvent, "event/L");
    matchedEventTree->Branch("photonIndex", &matchedPhotonIndex, "photonIndex/I");
    matchedEventTree->Branch("cscMDSClusterIndex", &matchedClusterIndex, "cscMDSClusterIndex/I");
    matchedEventTree->Branch("deltaR", &matchedDeltaR, "deltaR/F");
    matchedEventTree->Branch("photonMatchedToCSC", &matchedPhotonMatchedToCSC, "photonMatchedToCSC/O");
    matchedEventTree->Branch("photonPt", &matchedPhotonPt, "photonPt/F");
    matchedEventTree->Branch("photonEta", &matchedPhotonEta, "photonEta/F");
    matchedEventTree->Branch("photonPhi", &matchedPhotonPhi, "photonPhi/F");
    matchedEventTree->Branch("photonVx", &matchedPhotonVx, "photonVx/F");
    matchedEventTree->Branch("photonVy", &matchedPhotonVy, "photonVy/F");
    matchedEventTree->Branch("photonVz", &matchedPhotonVz, "photonVz/F");
    matchedEventTree->Branch("photonDecayR", &matchedPhotonDecayR, "photonDecayR/F");
    matchedEventTree->Branch("photonDecayLxyz", &matchedPhotonDecayLxyz, "photonDecayLxyz/F");
    matchedEventTree->Branch("photonPassMuonSystemGeometry",
                             &matchedPhotonPassMuonSystemGeometry,
                             "photonPassMuonSystemGeometry/O");
    matchedEventTree->Branch("inputFile", &matchedInputFile);
    matchedEventTree->Branch("HLT_CscCluster50_Photon20Unseeded",
                             &matchedHLT_CscCluster50_Photon20Unseeded,
                             "HLT_CscCluster50_Photon20Unseeded/O");

    targetPhotonTree = new TTree("TargetEventGenPhotons",
                                 "Gen photons from neutralinos in the requested event");
    targetPhotonTree->Branch("run", &targetRun, "run/I");
    targetPhotonTree->Branch("lumi", &targetLumi, "lumi/I");
    targetPhotonTree->Branch("event", &targetEvent, "event/L");
    targetPhotonTree->Branch("photonIndex", &targetPhotonIndex, "photonIndex/I");
    targetPhotonTree->Branch("motherIndex", &targetPhotonMotherIndex, "motherIndex/I");
    targetPhotonTree->Branch("pt", &targetPhotonPt, "pt/F");
    targetPhotonTree->Branch("eta", &targetPhotonEta, "eta/F");
    targetPhotonTree->Branch("phi", &targetPhotonPhi, "phi/F");
    targetPhotonTree->Branch("vx", &targetPhotonVx, "vx/F");
    targetPhotonTree->Branch("vy", &targetPhotonVy, "vy/F");
    targetPhotonTree->Branch("vz", &targetPhotonVz, "vz/F");
    targetPhotonTree->Branch("inputFile", &targetPhotonInputFile);
}

void LLPMuonSystemAcceptance::executeEvent() {
    if (IsDATA) return;

    FillHist("Event_Counts", 0., 1., 4, 0., 4.);

    std::set<int> counted_llp_indices;
    int n_llp_decays = 0;
    int n_llp_in_muon_system = 0;
    bool has_csc_photon_match = false;
    struct DecayedPhotonInfo {
        int photon_index;
        int cluster_index;
        float delta_r;
        bool matched_to_csc;
    };
    std::vector<DecayedPhotonInfo> decayed_photons;
    const bool is_target_event = (RunNumber == TARGET_RUN &&
                                  LumiBlock == TARGET_LUMI &&
                                  EventNumber == TARGET_EVENT);
    int target_event_n_photons = 0;

    for (int i = 0; i < nGenPart; ++i) {
        if (std::abs(GenPart_pdgId[i]) != DECAY_PHOTON_PDGID) continue;

        const int mother_idx = MotherIndex(i);
        if (mother_idx < 0 || mother_idx >= nGenPart) continue;
        if (std::abs(GenPart_pdgId[mother_idx]) != LLP_PDGID) continue;

        FillHist("GenPhotonFromNeutralino_Pt", GenPart_pt[i], 1., 100, 0., 2000.);
        FillHist("GenPhotonFromNeutralino_Eta", GenPart_eta[i], 1., 100, -5., 5.);

        float min_csc_dr = std::numeric_limits<float>::max();
        int matched_cluster_index = -1;
        for (int j = 0; j < ncscMDSCluster; ++j) {
            const float dr = DeltaR(GenPart_eta[i], GenPart_phi[i], cscMDSCluster_eta[j], cscMDSCluster_phi[j]);
            if (dr < min_csc_dr) {
                min_csc_dr = dr;
                matched_cluster_index = j;
            }
        }

        const bool photon_matched_to_csc = min_csc_dr < MAX_CSC_MATCH_DR;
        decayed_photons.push_back({i, matched_cluster_index, min_csc_dr, photon_matched_to_csc});

        if (photon_matched_to_csc) {
            has_csc_photon_match = true;
            FillHist("CSCPhotonMatch_DeltaR", min_csc_dr, 1., 60, 0., MAX_CSC_MATCH_DR);
            if (HLT_CscCluster50_Photon20Unseeded) {
                FillHist("CSCPhotonMatch_DeltaR_HLT_CscCluster50_Photon20Unseeded", min_csc_dr, 1., 60, 0., MAX_CSC_MATCH_DR);
            } else {
                FillHist("CSCPhotonMatch_DeltaR_NoHLT_CscCluster50_Photon20Unseeded", min_csc_dr, 1., 60, 0., MAX_CSC_MATCH_DR);
            }
        }

        if (is_target_event && targetPhotonTree) {
            ++target_event_n_photons;
            targetRun = RunNumber;
            targetLumi = LumiBlock;
            targetEvent = static_cast<Long64_t>(EventNumber);
            targetPhotonIndex = i;
            targetPhotonMotherIndex = mother_idx;
            targetPhotonPt = GenPart_pt[i];
            targetPhotonEta = GenPart_eta[i];
            targetPhotonPhi = GenPart_phi[i];
            targetPhotonVx = GenPart_vx[i];
            targetPhotonVy = GenPart_vy[i];
            targetPhotonVz = GenPart_vz[i];
            TFile *current_file = fChain->GetCurrentFile();
            targetPhotonInputFile = current_file ? current_file->GetName() : "";
            targetPhotonTree->Fill();

            if (!targetEventPrinted) {
                cout << "[LLPMuonSystemAcceptance] Gen photon from neutralino in target event "
                     << "run:lumi:event = " << TARGET_RUN << ":" << TARGET_LUMI << ":" << TARGET_EVENT
                     << ", photonIndex = " << targetPhotonIndex
                     << ", motherIndex = " << targetPhotonMotherIndex
                     << ", pt = " << targetPhotonPt
                     << ", eta = " << targetPhotonEta
                     << ", phi = " << targetPhotonPhi
                     << ", decay vertex (vx, vy, vz) cm = ("
                     << targetPhotonVx << ", " << targetPhotonVy << ", " << targetPhotonVz << ")"
                     << ", inputFile = " << targetPhotonInputFile << endl;
            }
        }

        if (counted_llp_indices.count(mother_idx)) continue;
        counted_llp_indices.insert(mother_idx);

        const float vx = GenPart_vx[i];
        const float vy = GenPart_vy[i];
        const float vz = GenPart_vz[i];
        const float r = std::hypot(vx, vy);
        const float lxyz = std::sqrt(vx * vx + vy * vy + vz * vz);
        const bool pass_ms = PassMuonSystemGeometry(i);

        ++n_llp_decays;
        if (pass_ms) ++n_llp_in_muon_system;

        FillHist("LLP_DecayVertex_R", r, 1., 100, 0., 1500.);
        FillHist("LLP_DecayVertex_absZ", std::abs(vz), 1., 100, 0., 2000.);
        FillHist("LLP_DecayVertex_Lxyz", lxyz, 1., 100, 0., 2500.);
        FillHist("LLP_Eta", GenPart_eta[mother_idx], 1., 100, -5., 5.);
        FillHist("LLP_Pt", GenPart_pt[mother_idx], 1., 100, 0., 2000.);
        FillHist("LLP_DecayVertex_R_vs_absZ", std::abs(vz), r, 1., 100, 0., 2000., 100, 0., 1500.);

        if (pass_ms) {
            FillHist("LLPInMS_DecayVertex_R", r, 1., 100, 0., 1500.);
            FillHist("LLPInMS_DecayVertex_absZ", std::abs(vz), 1., 100, 0., 2000.);
            FillHist("LLPInMS_Eta", GenPart_eta[mother_idx], 1., 100, -5., 5.);
        }

    }

    const bool has_llp_decay = n_llp_decays > 0;
    const bool has_llp_in_muon_system = n_llp_in_muon_system > 0;

    if (has_llp_decay) FillHist("Event_Counts", 1., 1., 4, 0., 4.);
    if (has_llp_in_muon_system) FillHist("Event_Counts", 2., 1., 4, 0., 4.);
    if (n_llp_in_muon_system >= 2) FillHist("Event_Counts", 3., 1., 4, 0., 4.);

    FillHist("Event_NLLPDecay", n_llp_decays, 1., 5, 0., 5.);
    FillHist("Event_NLLPInMuonSystem", n_llp_in_muon_system, 1., 5, 0., 5.);

    if (is_target_event && !targetEventPrinted) {
        if (target_event_n_photons == 0) {
            cout << "[LLPMuonSystemAcceptance] Target event "
                 << TARGET_RUN << ":" << TARGET_LUMI << ":" << TARGET_EVENT
                 << " has no gen photons with neutralino mother." << endl;
        }
        targetEventPrinted = true;
    }

    if (has_csc_photon_match && matchedEventTree) {
        TFile *current_file = fChain->GetCurrentFile();
        matchedInputFile = current_file ? current_file->GetName() : "";
        matchedRun = RunNumber;
        matchedLumi = LumiBlock;
        matchedEvent = static_cast<Long64_t>(EventNumber);
        matchedHLT_CscCluster50_Photon20Unseeded = HLT_CscCluster50_Photon20Unseeded;

        for (const auto &photon : decayed_photons) {
            matchedPhotonIndex = photon.photon_index;
            matchedClusterIndex = photon.cluster_index;
            matchedDeltaR = photon.cluster_index >= 0 ? photon.delta_r : -1.;
            matchedPhotonMatchedToCSC = photon.matched_to_csc;
            matchedPhotonPt = GenPart_pt[photon.photon_index];
            matchedPhotonEta = GenPart_eta[photon.photon_index];
            matchedPhotonPhi = GenPart_phi[photon.photon_index];
            matchedPhotonVx = GenPart_vx[photon.photon_index];
            matchedPhotonVy = GenPart_vy[photon.photon_index];
            matchedPhotonVz = GenPart_vz[photon.photon_index];
            matchedPhotonDecayR = std::hypot(matchedPhotonVx, matchedPhotonVy);
            matchedPhotonDecayLxyz = std::sqrt(matchedPhotonVx * matchedPhotonVx +
                                               matchedPhotonVy * matchedPhotonVy +
                                               matchedPhotonVz * matchedPhotonVz);
            matchedPhotonPassMuonSystemGeometry = PassMuonSystemGeometry(photon.photon_index);
            matchedEventTree->Fill();
        }
    }
}

bool LLPMuonSystemAcceptance::PassMuonSystemGeometry(int gen_index) const {
    const int mother_idx = MotherIndex(gen_index);
    if (mother_idx < 0 || mother_idx >= nGenPart) return false;

    const float vx = GenPart_vx[gen_index];
    const float vy = GenPart_vy[gen_index];
    const float vz = GenPart_vz[gen_index];
    const float r = std::hypot(vx, vy);
    const float abs_z = std::abs(vz);

    return std::abs(GenPart_eta[mother_idx]) < MAX_LLP_ABS_ETA
        && abs_z < MAX_ABS_Z_CM
        && abs_z > MIN_ABS_Z_CM
        && r < MAX_R_CM;
}

int LLPMuonSystemAcceptance::MotherIndex(int gen_index) const {
    if (Run == 3) return GenPart_genPartIdxMother[gen_index];
    return GenPart_genPartIdxMother_RunII[gen_index];
}

float LLPMuonSystemAcceptance::DeltaR(float eta1, float phi1, float eta2, float phi2) const {
    float dphi = std::remainder(phi1 - phi2, 2.0f * static_cast<float>(M_PI));
    const float deta = eta1 - eta2;
    return std::sqrt(deta * deta + dphi * dphi);
}

void LLPMuonSystemAcceptance::WriteHist() {
    if (matchedEventTree) {
        GetOutfile()->cd();
        matchedEventTree->Write();
    }
    if (targetPhotonTree) {
        GetOutfile()->cd();
        targetPhotonTree->Write();
    }
    AnalyzerCore::WriteHist();
}
