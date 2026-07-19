#include "PhotonTiming.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

PhotonTiming::PhotonTiming()
    : npho(0),
      nphoECALIdx(0),
      necalRechit(0),
      hasHLT20FilterBranch(false),
      hasHLT30FilterBranch(false) {
    pho_pt.resize(kMaxPho);
    pho_eta.resize(kMaxPho);
    pho_phi.resize(kMaxPho);
    pho_seedRechitTime.resize(kMaxPho);
    pho_seedRechitEnergy.resize(kMaxPho);
    pho_seedRechitIdx.resize(kMaxPho);
    pho_ecalRechitIdxStart.resize(kMaxPho);
    pho_nEcalRechits.resize(kMaxPho);
    pho_isOOT.resize(kMaxPho);
    pho_passHLTEG20HEUnseeded.resize(kMaxPho);
    pho_passHLTEG30HEUnseeded.resize(kMaxPho);

    phoECALIdx_EcalRechitIndex.resize(kMaxPhoECALIdx);

    ecalRechit_time.resize(kMaxEcalRechit);
    ecalRechit_energy.resize(kMaxEcalRechit);
    ecalRechit_flagOOT.resize(kMaxEcalRechit);
    ecalRechit_gainSwitch1.resize(kMaxEcalRechit);
    ecalRechit_gainSwitch6.resize(kMaxEcalRechit);
}

PhotonTiming::~PhotonTiming() {}

void PhotonTiming::RequireBranch(const TString &branchName, void *address) {
    if (!fChain || !fChain->GetBranch(branchName)) {
        cerr << "[PhotonTiming::initializeAnalyzer] Required branch not found: "
             << branchName << endl;
        exit(EXIT_FAILURE);
    }

    fChain->SetBranchStatus(branchName, 1);
    fChain->SetBranchAddress(branchName, address);
}

void PhotonTiming::OptionalBranch(const TString &branchName, void *address) {
    if (!fChain || !fChain->GetBranch(branchName)) {
        cout << "[PhotonTiming::initializeAnalyzer] Optional branch not found: "
             << branchName << endl;
        return;
    }

    fChain->SetBranchStatus(branchName, 1);
    fChain->SetBranchAddress(branchName, address);
}

void PhotonTiming::initializeAnalyzer() {
    if (!fChain) {
        cerr << "[PhotonTiming::initializeAnalyzer] fChain is null" << endl;
        exit(EXIT_FAILURE);
    }

    RequireBranch("npho", &npho);
    RequireBranch("pho_pt", pho_pt.data());
    RequireBranch("pho_eta", pho_eta.data());
    RequireBranch("pho_phi", pho_phi.data());
    RequireBranch("pho_seedRechitTime", pho_seedRechitTime.data());
    RequireBranch("pho_seedRechitEnergy", pho_seedRechitEnergy.data());
    RequireBranch("pho_seedRechitIdx", pho_seedRechitIdx.data());
    RequireBranch("pho_ecalRechitIdxStart", pho_ecalRechitIdxStart.data());
    RequireBranch("pho_nEcalRechits", pho_nEcalRechits.data());
    RequireBranch("pho_isOOT", pho_isOOT.data());

    RequireBranch("nphoECALIdx", &nphoECALIdx);
    RequireBranch("phoECALIdx_EcalRechitIndex", phoECALIdx_EcalRechitIndex.data());

    RequireBranch("necalRechit", &necalRechit);
    RequireBranch("ecalRechit_time", ecalRechit_time.data());
    RequireBranch("ecalRechit_energy", ecalRechit_energy.data());
    RequireBranch("ecalRechit_flagOOT", ecalRechit_flagOOT.data());
    RequireBranch("ecalRechit_gainSwitch1", ecalRechit_gainSwitch1.data());
    RequireBranch("ecalRechit_gainSwitch6", ecalRechit_gainSwitch6.data());

    hasHLT20FilterBranch = fChain->GetBranch("pho_passHLTEG20HEUnseeded");
    hasHLT30FilterBranch = fChain->GetBranch("pho_passHLTEG30HEUnseeded");
    OptionalBranch("pho_passHLTEG20HEUnseeded", pho_passHLTEG20HEUnseeded.data());
    OptionalBranch("pho_passHLTEG30HEUnseeded", pho_passHLTEG30HEUnseeded.data());

    cout << "[PhotonTiming::initializeAnalyzer] Reading photon timing from pho_seedRechitTime"
         << " and ecalRechit_time through phoECALIdx_EcalRechitIndex" << endl;
}

void PhotonTiming::executeEvent() {
    if (npho > kMaxPho || nphoECALIdx > kMaxPhoECALIdx || necalRechit > kMaxEcalRechit) {
        cerr << "[PhotonTiming::executeEvent] Branch size exceeds local buffer: "
             << "npho=" << npho << "/" << kMaxPho
             << ", nphoECALIdx=" << nphoECALIdx << "/" << kMaxPhoECALIdx
             << ", necalRechit=" << necalRechit << "/" << kMaxEcalRechit << endl;
        exit(EXIT_FAILURE);
    }

    const int safeNPho = npho;
    const int safeNPhoECALIdx = nphoECALIdx;
    const int safeNEcalRechit = necalRechit;

    FillHist("Event/npho", safeNPho, 1., 20, -0.5, 19.5);
    FillHist("Event/necalRechit", safeNEcalRechit, 1., 500, 0., 50000.);

    for (int ipho = 0; ipho < safeNPho; ++ipho) {
        const bool isOOT = pho_isOOT[ipho];
        const bool passHLT20 = hasHLT20FilterBranch && pho_passHLTEG20HEUnseeded[ipho];
        const bool passHLT30 = hasHLT30FilterBranch && pho_passHLTEG30HEUnseeded[ipho];
        const float seedTime = pho_seedRechitTime[ipho];
        const float weight = 1.;

        FillHist("Photon/pt", pho_pt[ipho], weight, 200, 0., 1000.);
        FillHist("Photon/eta", pho_eta[ipho], weight, 120, -6., 6.);
        FillHist("Photon/nEcalRechits", pho_nEcalRechits[ipho], weight, 200, -0.5, 199.5);

        if (seedTime > -900.) {
            FillHist("PhotonTiming/seedRechitTime_all", seedTime, weight, 200, -50., 50.);
            if (isOOT) FillHist("PhotonTiming/seedRechitTime_OOTPhoton", seedTime, weight, 200, -50., 50.);
            else FillHist("PhotonTiming/seedRechitTime_inTimePhoton", seedTime, weight, 200, -50., 50.);
            if (passHLT20) FillHist("PhotonTiming/seedRechitTime_passHLTEG20HEUnseeded", seedTime, weight, 200, -50., 50.);
            if (passHLT30) FillHist("PhotonTiming/seedRechitTime_passHLTEG30HEUnseeded", seedTime, weight, 200, -50., 50.);
        }

        const int start = pho_ecalRechitIdxStart[ipho];
        const int nRechits = pho_nEcalRechits[ipho];
        if (start < 0 || nRechits <= 0) continue;

        double energySum = 0.;
        double energyTimeSum = 0.;
        int nValidRechits = 0;

        for (int iidx = start; iidx < start + nRechits; ++iidx) {
            if (iidx < 0 || iidx >= safeNPhoECALIdx) continue;

            const int rechitIdx = phoECALIdx_EcalRechitIndex[iidx];
            if (rechitIdx < 0 || rechitIdx >= safeNEcalRechit) continue;

            const float rechitTime = ecalRechit_time[rechitIdx];
            const float rechitEnergy = ecalRechit_energy[rechitIdx];
            if (rechitTime < -900.) continue;

            FillHist("PhotonTiming/associatedRechitTime_all", rechitTime, weight, 200, -50., 50.);
            if (isOOT) FillHist("PhotonTiming/associatedRechitTime_OOTPhoton", rechitTime, weight, 200, -50., 50.);
            else FillHist("PhotonTiming/associatedRechitTime_inTimePhoton", rechitTime, weight, 200, -50., 50.);
            if (ecalRechit_flagOOT[rechitIdx]) FillHist("PhotonTiming/associatedRechitTime_flagOOT", rechitTime, weight, 200, -50., 50.);
            if (ecalRechit_gainSwitch1[rechitIdx]) FillHist("PhotonTiming/associatedRechitTime_gainSwitch1", rechitTime, weight, 200, -50., 50.);
            if (ecalRechit_gainSwitch6[rechitIdx]) FillHist("PhotonTiming/associatedRechitTime_gainSwitch6", rechitTime, weight, 200, -50., 50.);

            if (rechitEnergy > 0.) {
                energySum += rechitEnergy;
                energyTimeSum += rechitEnergy * rechitTime;
            }
            nValidRechits++;
        }

        FillHist("PhotonTiming/nValidAssociatedRechits", nValidRechits, weight, 200, -0.5, 199.5);
        if (energySum > 0.) {
            FillHist("PhotonTiming/energyWeightedAssociatedRechitTime", energyTimeSum / energySum, weight, 200, -50., 50.);
        }
    }
}
