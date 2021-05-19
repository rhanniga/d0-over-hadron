#include "AliAnalysisTaskD0HadronRatio.h"

//BOTH OF THESE ARE WEIRD, BUT APPARENTLY NECESSARRY
class AliAnalysisTaskD0HadronRatio;
ClassImp(AliAnalysisTaskD0HadronRatio);


AliAnalysisTaskD0HadronRatio::AliAnalysisTaskD0HadronRatio() :
    AliAnalysisTaskSE(),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHD0{0},
    fDphiHD0Rotated{0},
    fDphiHD0RotatedKaon{0},
    fDphiHD0RotatedPi{0},
    fDphiHD0Flipped{0},
    fDphiHH{0},
    fDphiHD0LS{0},
    fDphiHD0Mixed{0},
    fDphiHHMixed{0},
    fDphiHD0LSMixed{0},
    fCorPoolMgr{0}
    // fPid{0},
    // fSignalAnalysis{0}
{
    MULT_LOW = 0;
    MULT_HIGH = 20;
    CENT_ESTIMATOR = "V0A";
    DAUGHTER_TRK_BIT = AliAODTrack::kTrkGlobalNoDCA; // = 16
    ASSOC_TRK_BIT = 1024; // global tracks with tight pt dependent dca cut (selecting primaries)
    TRIG_TRK_BIT = AliAODTrack::kIsHybridGCG; // = 2^20
}

AliAnalysisTaskD0HadronRatio::AliAnalysisTaskD0HadronRatio(const char *name) :
    AliAnalysisTaskSE(name),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHD0{0},
    fDphiHD0Rotated{0},
    fDphiHD0RotatedKaon{0},
    fDphiHD0RotatedPi{0},
    fDphiHD0Flipped{0},
    fDphiHH{0},
    fDphiHD0LS{0},
    fDphiHD0Mixed{0},
    fDphiHHMixed{0},
    fDphiHD0LSMixed{0},
    fCorPoolMgr{0}
    // fPid{0},
    // fSignalAnalysis{0}
{

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    MULT_LOW = 0;
    MULT_HIGH = 20;
    CENT_ESTIMATOR = "V0A";
    DAUGHTER_TRK_BIT = AliAODTrack::kTrkGlobalNoDCA; // = 16
    ASSOC_TRK_BIT = 1024; // Global tracks with tight, pt dependent DCA
    TRIG_TRK_BIT = AliAODTrack::kIsHybridGCG; // = 2^20

}

AliAnalysisTaskD0HadronRatio::~AliAnalysisTaskD0HadronRatio()
{

    if(fOutputList) delete fOutputList;

}

void AliAnalysisTaskD0HadronRatio::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(true);

    //Generating the mixed event pools:
    int poolSize = 500;
    int trackDepth = 1000;

    int numMultBins = 1;
    double multBins[2] = {MULT_LOW, MULT_HIGH};

    int numzVtxBins = 10;
    double zVtxBins[11] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};

    fCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);

    fTriggersAndD0sPerEvent_All = new TH2D("fTriggersAndD0sPerEvent_All", "Triggers and D0s per event (all p_{T})", 10, 0, 10, 10, 0, 10);
    fOutputList->Add(fTriggersAndD0sPerEvent_All);

    fTriggersAndD0sPerEvent_2_4 = new TH2D("fTriggersAndD0sPerEvent_2_4", "Triggers and D0s per event (2-4 p_{T})", 10, 0, 10, 10, 0, 10);
    fOutputList->Add(fTriggersAndD0sPerEvent_2_4);


    //Distribution axes are: Pt, Phi, Eta, zVtx
    int dist_bins[4] = {100, 16, 20, 10};
    double dist_mins[4] = {0, 0, -1, -10};
    double dist_maxes[4] = {15, 6.28, 1, 10};

    fLooseDist = new THnSparseF("fLooseDist", "All Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fLooseDist);

    fTriggerDist = new THnSparseF("fTriggerDist", "Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fTriggerDist);

    fAssociatedHDist = new THnSparseF("fAssociatedHDist", "Associated Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fAssociatedHDist);

    fD0Dist = new THnSparseF("fD0Dist", "D0 Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fD0Dist);


    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Inv Mass, Zvtx
    int hd_cor_bins[6] = {8, 10, 16, 20, 400, 10};
    double hd_cor_mins[6] = {4.0, 1, -1.0*TMath::Pi()/2.0, -2.0, 0.4, -10};
    double hd_cor_maxes[6] = {12.0, 6, 3.0*TMath::Pi()/2.0, 2.0, 2.4, 10};

    fDphiHD0 = new THnSparseF("fDphiHD0", "Hadron-D0 Correlation Histogram", 6, hd_cor_bins, hd_cor_mins, hd_cor_maxes);
    fOutputList->Add(fDphiHD0);

    fDphiHD0V0 = new THnSparseF("fDphiHD0V0", "Hadron-D0 (using V0) Correlation Histogram", 6, hd_cor_bins, hd_cor_mins, hd_cor_maxes);
    fOutputList->Add(fDphiHD0V0);

    fDphiHD0Rotated = new THnSparseF("fDphiHD0Rotated", "Hadron-D0 (rotated) Correlation Histogram", 6, hd_cor_bins, hd_cor_mins, hd_cor_maxes);
    fOutputList->Add(fDphiHD0Rotated);

    fDphiHD0RotatedKaon = new THnSparseF("fDphiHD0RotatedKaon", "Hadron-D0 (rotated K) Correlation Histogram", 6, hd_cor_bins, hd_cor_mins, hd_cor_maxes);
    fOutputList->Add(fDphiHD0RotatedKaon);

    fDphiHD0RotatedPi = new THnSparseF("fDphiHD0RotatedPi", "Hadron-D0 (pion rotated) Correlation Histogram", 6, hd_cor_bins, hd_cor_mins, hd_cor_maxes);
    fOutputList->Add(fDphiHD0RotatedPi);

    fDphiHD0Flipped = new THnSparseF("fDphiHD0Flipped", "Hadron-D0 (flipped) Correlation Histogram", 6, hd_cor_bins, hd_cor_mins, hd_cor_maxes);
    fOutputList->Add(fDphiHD0Flipped);

    fDphiHD0LS = new THnSparseF("fDphiHD0LS", "Hadron-D0 LS Correlation Histogram", 6, hd_cor_bins, hd_cor_mins, hd_cor_maxes);
    fOutputList->Add(fDphiHD0LS);

    fDphiHD0Mixed = new THnSparseF("fDphiHD0Mixed", "Mixed Hadron-D0 Correlation Histogram", 6, hd_cor_bins, hd_cor_mins, hd_cor_maxes);
    fOutputList->Add(fDphiHD0Mixed);

    fDphiHD0LSMixed = new THnSparseF("fDphiHD0LSMixed", "Mixed Hadron-D0 LS Correlation Histogram", 6, hd_cor_bins, hd_cor_mins, hd_cor_maxes);
    fOutputList->Add(fDphiHD0LSMixed);


    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Zvtx
    int hh_cor_bins[5] = {20, 20, 16, 20, 10};
    double hh_cor_mins[5] = {2, 2, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double hh_cor_maxes[5] = {12, 12, 3.0*TMath::Pi()/2.0, 2.0, 10};

    fDphiHH = new THnSparseF("fDphiHH", "Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fOutputList->Add(fDphiHH);

    fDphiTriggerTrigger = new THnSparseF("fDphiTriggerTrigger", "Trigger-Trigger Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fOutputList->Add(fDphiTriggerTrigger);

    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fOutputList->Add(fDphiHHMixed);

    fDphiTriggerTriggerMixed = new THnSparseF("fDphiTriggerTriggerMixed", "MixedTrigger-Trigger Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fOutputList->Add(fDphiTriggerTriggerMixed);


    //axes are pion 0: dca pion 1: dca Kaon
    //              2: pt pion  3: pt Kaon
    //              4: pt D0    5: mass D0
    int dca_bins[6] = {100, 100, 20, 20, 20, 400};
    double dca_mins[6] = {-2.4, -2.4, 0, 0, 0, 0.4};
    double dca_maxes[6] = {2.4, 2.4, 10, 10, 10, 2.4};

    fD0DaughterDCA = new THnSparseF("fD0DaughterDCA", "#D0^{0} daughter DCA dist", 6, dca_bins, dca_mins, dca_maxes);
    fOutputList->Add(fD0DaughterDCA);

    PostData(1, fOutputList);

}


void AliAnalysisTaskD0HadronRatio::FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist)
{

    double dist_points[4]; //Pt, Phi, Eta, zVtx
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i];
        dist_points[0] = particle->Pt();
        dist_points[1] = particle->Phi();
        dist_points[2] = particle->Eta();
        dist_points[3] = zVtx;
        fDist->Fill(dist_points);
    }

}

void AliAnalysisTaskD0HadronRatio::FillSingleParticleDist(std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> particle_list, double zVtx, THnSparse* fDist)
{

    double dist_points[4]; //Pt, Phi, Eta, zVtx
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].particle;
        dist_points[0] = particle.Pt();
        dist_points[1] = particle.Phi();
        dist_points[2] = particle.Eta();
        dist_points[3] = zVtx;
        fDist->Fill(dist_points);
    }

}

AliAnalysisTaskD0HadronRatio::AliMotherContainer AliAnalysisTaskD0HadronRatio::DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{

    AliAnalysisTaskD0HadronRatio::AliMotherContainer mom;
    mom.particle.SetPx(track1->Px() + track2->Px());
    mom.particle.SetPy(track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}

AliAnalysisTaskD0HadronRatio::AliMotherContainer AliAnalysisTaskD0HadronRatio::RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle)
{

    AliAnalysisTaskD0HadronRatio::AliMotherContainer mom;
    // Rotating track1
    TVector3 track1Vector(track1->Px(), track1->Py(), track1->Pz());
    track1Vector.RotateZ(angle);
    mom.particle.SetPx(track1Vector(0) + track2->Px());
    mom.particle.SetPy(track1Vector(1) + track2->Py());
    mom.particle.SetPz(track1Vector(2) + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}

AliAnalysisTaskD0HadronRatio::AliMotherContainer AliAnalysisTaskD0HadronRatio::FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{

    AliAnalysisTaskD0HadronRatio::AliMotherContainer mom;
    // Flipping track1
    mom.particle.SetPx(-track1->Px() + track2->Px());
    mom.particle.SetPy(-track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}

void AliAnalysisTaskD0HadronRatio::MakeSameHD0Correlations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list, THnSparse* fDphi, double zVtx, bool eff)
{

    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)D0_list.size(); i++) {
            auto D0 = D0_list[i];

            //Make sure trigger isn't one of the daughters of D0
            if((trigger->GetID() == D0.daughter1ID) || (trigger->GetID() == D0.daughter2ID)) continue;

            dphi_point[1] = D0.particle.Pt();
            dphi_point[2] = trigger->Phi() - D0.particle.Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - D0.particle.Eta();
            dphi_point[4] = D0.particle.M();
            dphi_point[5] = zVtx;
            if(eff) {
                // double triggerScale = 1.0/ftrigEff->Eval(trigger->Pt());
                // double D0Scale = 1.0/fD0Eff->Eval(D0.particle.Pt());
                // double totalScale = triggerScale*D0Scale;
                // fDphi->Fill(dphi_point, totalScale);
                fDphi->Fill(dphi_point);
            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }

}

void AliAnalysisTaskD0HadronRatio::MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
{

    double dphi_point[5];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)associated_h_list.size(); i++) {
            auto associate = associated_h_list[i];

            dphi_point[1] = associate->Pt();
            dphi_point[2] = trigger->Phi() - associate->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - associate->Eta();
            dphi_point[4] = zVtx;
            if(eff) {
                // double triggerScale = 1.0/ftrigEff->Eval(trigger->Pt());
                // double assocatedScale = 1.0/fassociatedEff->Eval(associate->Pt());
                // double totalScale = triggerScale*associatedScale;
                fDphi->Fill(dphi_point);
            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }

}

void AliAnalysisTaskD0HadronRatio::MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx, bool eff)
{

    double dphi_point[5];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = j+1; i < (int)trigger_list.size(); i++) {
            auto associate = trigger_list[i];

            dphi_point[1] = associate->Pt();
            dphi_point[2] = trigger->Phi() - associate->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - associate->Eta();
            dphi_point[4] = zVtx;
            if(eff) {
                // double triggerScale = 1.0/ftrigEff->Eval(trigger->Pt());
                // double assocatedScale = 1.0/ftrigEff->Eval(associate->Pt());
                // double totalScale = triggerScale*associatedScale;
                fDphi->Fill(dphi_point);
            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }

}

void AliAnalysisTaskD0HadronRatio::MakeMixedHD0Correlations(AliEventPool* fPool, std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list , THnSparse* fDphi, double zVtx, bool eff)
{

    double dphi_point[6];
    int numEvents = fPool->GetCurrentNEvents();
    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)D0_list.size(); j++) {
                auto D0 = D0_list[j];

                dphi_point[1] = D0.particle.Pt();
                dphi_point[2] = trigger->Phi() - D0.particle.Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - D0.particle.Eta();
                dphi_point[4] = D0.particle.M();
                dphi_point[5] = zVtx;
                if(eff) {
                    // double triggerScale = 1.0/ftrigEff->Eval(trigger->Pt());
                    // double D0Scale = 1.0/fD0Eff->Eval(D0.particle.Pt());
                    // double totalScale = triggerScale*D0Scale;
                    fDphi->Fill(dphi_point);
                }
                else{
                    fDphi->Fill(dphi_point);
                }
            }
        }
    }
}

void AliAnalysisTaskD0HadronRatio::MakeMixedHHCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
{

    double dphi_point[5];

    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)associated_h_list.size(); j++) {
                auto associate = associated_h_list[j];

                dphi_point[1] = associate->Pt();
                dphi_point[2] = trigger->Phi() - associate->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - associate->Eta();
                dphi_point[4] = zVtx;
                if(eff) {
                    // double triggerScale = 1.0/ftrigEff->Eval(trigger->Pt());
                    // double assocatedScale = 1.0/fassociatedEff->Eval(associate->Pt());
                    // double totalScale = triggerScale*associatedScale;
                    fDphi->Fill(dphi_point);
                }
                else{
                    fDphi->Fill(dphi_point);
                }
            }
        }
    }

}

bool AliAnalysisTaskD0HadronRatio::PassDaughterCuts(AliAODTrack *track){

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && (track->TestFilterMask(DAUGHTER_TRK_BIT));

    // pass = pass && (track->IsOn(AliAODTrack::kTPCrefit));

    // pass = pass && (track->GetTPCCrossedRows() > 70);

    // float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    // pass = pass && (ratio > 0.8);

    return pass;
}

bool AliAnalysisTaskD0HadronRatio::PassAssociatedCuts(AliAODTrack *track){

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && track->TestFilterMask(ASSOC_TRK_BIT);

    return pass;
}

Bool_t AliAnalysisTaskD0HadronRatio::PassTriggerCuts(AliAODTrack *track){

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && track->TestBit(TRIG_TRK_BIT);

    return pass;
}

void AliAnalysisTaskD0HadronRatio::UserExec(Option_t*)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        std::cout << "THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!" << std::endl;
        std::abort();
    }


    fpidResponse = fInputHandler->GetPIDResponse();


    //Event cuts
    TString cent_estimator = CENT_ESTIMATOR;
    double multPercentile = 0;

    fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(fMultSelection) multPercentile = fMultSelection->GetMultiplicityPercentile(cent_estimator.Data());
    else return;

    if(multPercentile < MULT_LOW || multPercentile > MULT_HIGH) return;

    AliVVertex *prim = fAOD->GetPrimaryVertex();
    int NcontV = prim->GetNContributors();
    if(NcontV < 3) return;

    double primZ = prim->GetZ();
    if(primZ < -10 || primZ > 10) return;


    int numTracks = fAOD->GetNumberOfTracks();

    std::vector<AliAODTrack*> unlikelyKaon_list;
    std::vector<AliAODTrack*> unlikelyPion_list;
    std::vector<AliAODTrack*> kPlus_list;
    std::vector<AliAODTrack*> kMinus_list;
    std::vector<AliAODTrack*> piPlus_list;
    std::vector<AliAODTrack*> piMinus_list;
    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> associated_h_list;
    std::vector<AliAODTrack*> all_hadron_list;
    std::vector<AliAODTrack*> k_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;

        //List for comparison with cuts/filter bits
        all_hadron_list.push_back(track);

        //Filter for trigger particles
        if(PassTriggerCuts(track)) {
            trigger_list.push_back(track);
            AliCFParticle *triggerPart = new AliCFParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0);
            fMixedTrackObjArray->Add(triggerPart);
        }

        if(PassAssociatedCuts(track)) {
            associated_h_list.push_back(track);
        }

        if(PassDaughterCuts(track)) {
            double TPCNSigmaPion = 1000;
            double TOFNSigmaPion = 1000;

            TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
            TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);

            if(TOFNSigmaPion != 1000 && track->Charge() != 1) unlikelyPion_list.push_back(track);

            if(TMath::Abs(TPCNSigmaPion) <= 3 && (TMath::Abs(TOFNSigmaPion) <= 3 || TOFNSigmaPion == 1000)) {

                if(track->Charge() == 1){
                    piPlus_list.push_back(track);
                }
                else {
                    piMinus_list.push_back(track);
                }
            }

            double TPCNSigmaKaon = 1000;
            double TOFNSigmaKaon = 1000;


            TPCNSigmaKaon = fpidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
            TOFNSigmaKaon = fpidResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

            if(TOFNSigmaKaon != 1000 && track->Charge() == 1) unlikelyKaon_list.push_back(track);

            if(TMath::Abs(TPCNSigmaKaon) <= 3 && (TMath::Abs(TOFNSigmaKaon) <= 3 || TOFNSigmaKaon == 1000)) {

                if(track->Charge() == 1){
                    kPlus_list.push_back(track);
                }
                else {
                    kMinus_list.push_back(track);
                }
            }
        }
    }

    //Making list of possible D0s (have to do +/- for kaon or pi):

    std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list;
    std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list_v0;
    std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list_signal_region;
    std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list_signal_region_2_4;
    std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list_RotatedPion;
    std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list_RotatedKaon;
    std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list_RotatedPi; //pi radians, not the particle
    std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list_Flipped;
    std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list_LS;

    // for(int i = 0; i < (int)unlikelyPion_list.size(); i++) {
    //     for(int j = 0; j < (int) unlikelyKaon_list.size(); j++) {

    //         double TPCNSigmaPion = 1000;
    //         double TOFNSigmaPion = 1000;
    //         double TPCNSigmaKaon = 1000;
    //         double TOFNSigmaKaon = 1000;

    //         TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(unlikelyPion_list[i], AliPID::kPion);
    //         TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(unlikelyPion_list[i], AliPID::kPion);
    //         TPCNSigmaKaon = fpidResponse->NumberOfSigmasTPC(unlikelyKaon_list[j], AliPID::kKaon);
    //         TOFNSigmaKaon = fpidResponse->NumberOfSigmasTOF(unlikelyKaon_list[j], AliPID::kKaon);

    //         AliMotherContainer mother = DaughtersToMother(unlikelyPion_list[i], unlikelyKaon_list[j], 0.1396, 0.493677);
    //         double mass = mother.particle.M();
    //         double pt = mother.particle.Pt();

    //         double signal_array[6] = {TPCNSigmaPion, TOFNSigmaPion, TPCNSigmaKaon, TOFNSigmaKaon, mass, pt};
    //         fSignalAnalysis->Fill(signal_array);
    //     }
    // }


    for(int i = 0; i < (int)piMinus_list.size(); i++) {
        for(int j = 0; j < (int) kPlus_list.size(); j++) {
            auto pion = piMinus_list[i];
            auto kaon = kPlus_list[j];

            AliMotherContainer D0 = DaughtersToMother(pion, kaon, 0.1396, 0.493677);
            
            double pion_dz[2];
            double pion_covar[3];

            double kaon_dz[2];
            double kaon_covar[3];

            bool is_pionDCA = piMinus_list[i]->PropagateToDCA(prim, fAOD->GetMagneticField(), 20., pion_dz, pion_covar);
            bool is_kaonDCA = kPlus_list[j]->PropagateToDCA(prim, fAOD->GetMagneticField(), 20., kaon_dz, kaon_covar);

            if(is_pionDCA && is_kaonDCA) {
                double fillArray[6] = {pion_dz[0], kaon_dz[0], pion->Pt(), kaon->Pt(), D0.particle.Pt(), D0.particle.M()};
                fD0DaughterDCA->Fill(fillArray);
            }

            AliMotherContainer D0_RotatedPi = RotatedDaughtersToMother(pion, kaon, 0.1396, 0.493677, TMath::Pi());
            AliMotherContainer D0_Flipped = FlippedDaughtersToMother(pion, kaon, 0.1396, 0.493677);
            D0_list.push_back(D0);
            D0_list_RotatedPi.push_back(D0_RotatedPi);
            D0_list_Flipped.push_back(D0_Flipped);

            AliMotherContainer D0_Rotated;
            AliMotherContainer D0_RotatedKaon;
            for(int k = 1; k < 12; k++) {
                D0_Rotated = RotatedDaughtersToMother(pion, kaon, 0.1396, 0.493677, (2*TMath::Pi()*k)/12);
                D0_RotatedKaon = RotatedDaughtersToMother(kaon, pion, 0.493677, 0.1396, (2*TMath::Pi()*k)/12);
                D0_list_RotatedPion.push_back(D0_Rotated);
                D0_list_RotatedKaon.push_back(D0_RotatedKaon);

            }
        }
    }

    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = 0; j < (int) kMinus_list.size(); j++) {
            auto pion = piPlus_list[i];
            auto kaon = kMinus_list[j];

            AliMotherContainer D0 = DaughtersToMother(pion, kaon, 0.1396, 0.493677);
            AliMotherContainer D0_RotatedPi = RotatedDaughtersToMother(pion, kaon, 0.1396, 0.493677, TMath::Pi());
            AliMotherContainer D0_Flipped = FlippedDaughtersToMother(pion, kaon, 0.1396, 0.493677);
            D0_list.push_back(D0);
            D0_list_RotatedPi.push_back(D0_RotatedPi);
            D0_list_Flipped.push_back(D0_Flipped);

            AliMotherContainer D0_Rotated;
            AliMotherContainer D0_RotatedKaon;
            for(int k = 1; k < 12; k++) {
                D0_Rotated = RotatedDaughtersToMother(pion, kaon, 0.1396, 0.493677, (2*TMath::Pi()*k)/12);
                D0_RotatedKaon = RotatedDaughtersToMother(kaon, pion, 0.493677, 0.1396, (2*TMath::Pi()*k)/12);
                D0_list_RotatedPion.push_back(D0_Rotated);
                D0_list_RotatedKaon.push_back(D0_RotatedKaon);

            }

        }
    }

    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = 0; j < (int) kPlus_list.size(); j++) {
            if(piPlus_list[i]->GetID() == kPlus_list[j]->GetID()) continue;
            AliMotherContainer D0 = DaughtersToMother(piPlus_list[i], kPlus_list[j], 0.1396, 0.493677);
            D0_list_LS.push_back(D0);
        }
    }

    for(int i = 0; i < (int)piMinus_list.size(); i++) {
        for(int j = 0; j < (int) kMinus_list.size(); j++) {
            if(piMinus_list[i]->GetID() == kMinus_list[j]->GetID()) continue;
            AliMotherContainer D0 = DaughtersToMother(piMinus_list[i], kMinus_list[j], 0.1396, 0.493677);
            D0_list_LS.push_back(D0);
        }
    }


    for(int i = 0; i < (int)D0_list.size(); i++) {
        if(D0_list[i].particle.M() < 1.95 && D0_list[i].particle.M() > 1.8) {
            D0_list_signal_region.push_back(D0_list[i]);
            if(D0_list[i].particle.Pt() < 4 && D0_list[i].particle.Pt() > 2) {
                D0_list_signal_region_2_4.push_back(D0_list[i]);
            }
        }
    }


    // V0 SECTION

    int numV0s = fAOD->GetNumberOfV0s();
    for(int i = 0; i < numV0s; i++) {
        AliAODv0 *v0 = fAOD->GetV0(i);

        AliAODTrack* posTrack = (AliAODTrack*) v0->GetDaughter(0);
        AliAODTrack* negTrack = (AliAODTrack*) v0->GetDaughter(1);

        // Occasionally returns null, not quite sure why...
        if(!posTrack || !negTrack) continue;
        if(!(PassDaughterCuts(posTrack) && PassDaughterCuts(negTrack))) continue;


        double TPCNSigmaKaon = 1000;
        double TOFNSigmaKaon = 1000;
        double TPCNSigmaPion = 1000;
        double TOFNSigmaPion = 1000;

        TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(posTrack, AliPID::kPion);
        TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(posTrack, AliPID::kPion);
        TPCNSigmaKaon = fpidResponse->NumberOfSigmasTPC(negTrack, AliPID::kKaon);
        TOFNSigmaKaon = fpidResponse->NumberOfSigmasTOF(negTrack, AliPID::kKaon);

        bool isPosTrackPion = TMath::Abs(TPCNSigmaPion) <= 3 && (TMath::Abs(TOFNSigmaPion) <= 3 || TOFNSigmaPion == 1000);
        bool isNegTrackKaon = TMath::Abs(TPCNSigmaKaon) <= 3 && (TMath::Abs(TOFNSigmaKaon) <= 3 || TOFNSigmaKaon == 1000);

        if(isPosTrackPion && isNegTrackKaon) {
            auto D0 = DaughtersToMother(posTrack, negTrack, 0.1396, 0.493677);
            D0_list_v0.push_back(D0);
        }
    }

    // Filling all of our single particle distribution histograms:
    FillSingleParticleDist(trigger_list, primZ, fTriggerDist);
    FillSingleParticleDist(associated_h_list, primZ, fAssociatedHDist);
    FillSingleParticleDist(all_hadron_list, primZ, fLooseDist);
    FillSingleParticleDist(D0_list_signal_region, primZ, fD0Dist);

    // Filling all of our correlation histograms
    MakeSameHD0Correlations(trigger_list, D0_list, fDphiHD0, primZ);
    MakeSameHD0Correlations(trigger_list, D0_list_v0, fDphiHD0V0, primZ);
    MakeSameHD0Correlations(trigger_list, D0_list_RotatedPi, fDphiHD0RotatedPi, primZ);
    MakeSameHD0Correlations(trigger_list, D0_list_Flipped, fDphiHD0Flipped, primZ);
    MakeSameHD0Correlations(trigger_list, D0_list_RotatedPion, fDphiHD0Rotated, primZ);
    MakeSameHD0Correlations(trigger_list, D0_list_RotatedKaon, fDphiHD0RotatedKaon, primZ);
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHH, primZ);
    MakeSameTriggerTriggerCorrelations(trigger_list, fDphiTriggerTrigger, primZ);
    MakeSameHD0Correlations(trigger_list, D0_list_LS, fDphiHD0LS, primZ);


    fTriggersAndD0sPerEvent_All->Fill(trigger_list.size(), D0_list_signal_region.size());
    fTriggersAndD0sPerEvent_2_4->Fill(trigger_list.size(), D0_list_signal_region_2_4.size());


    if(D0_list.size() > 0 && associated_h_list.size() > 0) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }


        else {
            if(fCorPool->IsReady()) {
                MakeMixedHD0Correlations(fCorPool, D0_list, fDphiHD0Mixed, primZ);
                MakeMixedHD0Correlations(fCorPool, D0_list_LS, fDphiHD0LSMixed, primZ);
                MakeMixedHHCorrelations(fCorPool, associated_h_list, fDphiHHMixed, primZ);
                MakeMixedHHCorrelations(fCorPool, trigger_list, fDphiTriggerTriggerMixed, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }

    PostData(1, fOutputList);

}

void AliAnalysisTaskD0HadronRatio::Terminate(Option_t *option)
{
}
