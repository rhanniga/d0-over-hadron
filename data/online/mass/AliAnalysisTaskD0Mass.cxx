#include "AliAnalysisTaskD0Mass.h"

//BOTH OF THESE ARE WEIRD, BUT APPARENTLY NECESSARRY
class AliAnalysisTaskD0Mass;
ClassImp(AliAnalysisTaskD0Mass);


AliAnalysisTaskD0Mass::AliAnalysisTaskD0Mass() :
    AliAnalysisTaskSE(),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fD0Dist{0}, 
    fHFD0Dist{0},
    fTriggeredD0Dist{0}
{
    MULT_LOW = 0;
    MULT_HIGH = 20;
    CENT_ESTIMATOR = "V0A";
    DAUGHTER_TRK_BIT = AliAODTrack::kTrkGlobalNoDCA; // = 16, not used
}

AliAnalysisTaskD0Mass::AliAnalysisTaskD0Mass(const char *name) :
    AliAnalysisTaskSE(name),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fD0Dist{0}, 
    fHFD0Dist{0},
    fTriggeredD0Dist{0}
{

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    MULT_LOW = 0;
    MULT_HIGH = 20;
    CENT_ESTIMATOR = "V0A";
    DAUGHTER_TRK_BIT = AliAODTrack::kTrkGlobalNoDCA; // = 16
}

AliAnalysisTaskD0Mass::~AliAnalysisTaskD0Mass()
{

    if(fOutputList) delete fOutputList;

}

void AliAnalysisTaskD0Mass::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(true);



    //Distribution axes are: Pt, Phi, Eta, InvMass
    int dist_bins[4] = {60, 16, 20, 400};
    double dist_mins[4] = {0, -3.14, -2, 0.4};
    double dist_maxes[4] = {15, 3.14, 2, 2.4};

    fD0Dist = new THnSparseF("fD0Dist", "D0 Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fD0Dist);

    fHFD0Dist = new THnSparseF("fHFD0Dist", "D0 Distribution from HF vertex thing", 4, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fHFD0Dist);

    //Distribution axes are: Pt, Phi, Eta, InvMass
    int triggered_dist_bins[5] = {60, 16, 20, 400, 50};
    double triggered_dist_mins[5] = {0, -3.14, -2, 0.4, 0};
    double triggered_dist_maxes[5] = {15, 3.14, 2, 2.4, 10};

    fTriggeredD0Dist = new THnSparseF("fTriggeredD0Dist", "D0 Distribution with max trigger pt axis", 5, triggered_dist_bins, triggered_dist_mins, triggered_dist_maxes);
    fOutputList->Add(fTriggeredD0Dist);

    fTriggeredHFD0Dist = new THnSparseF("fTriggeredHFD0Dist", "D0 Distribution from HF vertex thing with max trigger pt axis", 5, triggered_dist_bins, triggered_dist_mins, triggered_dist_maxes);
    fOutputList->Add(fTriggeredHFD0Dist);

    PostData(1, fOutputList);

}

void AliAnalysisTaskD0Mass::FillSingleParticleDist(std::vector<AliAnalysisTaskD0Mass::AliMotherContainer> particle_list, THnSparse* fDist)
{

    double dist_points[4]; //Pt, Phi, Eta, Mass
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].particle;
        dist_points[0] = particle.Pt();
        dist_points[1] = particle.Phi();
        dist_points[2] = particle.Eta();
        dist_points[3] = particle.M();
        fDist->Fill(dist_points);
    }

}

void AliAnalysisTaskD0Mass::FillHFD0Dist(TClonesArray* hf_d0s, THnSparse* fDist)
{

    double dist_points[4]; //Pt, Phi, Eta, Mass
    int numD0s = hf_d0s->GetEntriesFast();

    for(int iD0 = 0; iD0 < numD0s; iD0++) {

        AliAODRecoDecayHF2Prong* d0 = (AliAODRecoDecayHF2Prong *)hf_d0s->UncheckedAt(iD0);

        dist_points[0] = d0->Pt();
        dist_points[1] = d0->Phi();
        dist_points[2] = d0->Eta();

        double daughter0_TPCNSigmaPion = 1000;
        double daughter0_TOFNSigmaPion = 1000;
        double daughter1_TPCNSigmaPion = 1000;
        double daughter1_TOFNSigmaPion = 1000;

        daughter0_TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC((AliAODTrack*)d0->GetDaughter(0), AliPID::kPion);
        daughter0_TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF((AliAODTrack*)d0->GetDaughter(0), AliPID::kPion);
        daughter1_TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC((AliAODTrack*)d0->GetDaughter(1), AliPID::kPion);
        daughter1_TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF((AliAODTrack*)d0->GetDaughter(1), AliPID::kPion);

        double daughter0_TPCNSigmaKaon = 1000;
        double daughter0_TOFNSigmaKaon = 1000;
        double daughter1_TPCNSigmaKaon = 1000;
        double daughter1_TOFNSigmaKaon = 1000;

        daughter0_TPCNSigmaKaon = fpidResponse->NumberOfSigmasTPC((AliAODTrack*)d0->GetDaughter(0), AliPID::kKaon);
        daughter0_TOFNSigmaKaon = fpidResponse->NumberOfSigmasTOF((AliAODTrack*)d0->GetDaughter(0), AliPID::kKaon);
        daughter1_TPCNSigmaKaon = fpidResponse->NumberOfSigmasTPC((AliAODTrack*)d0->GetDaughter(1), AliPID::kKaon);
        daughter1_TOFNSigmaKaon = fpidResponse->NumberOfSigmasTOF((AliAODTrack*)d0->GetDaughter(1), AliPID::kKaon);
        
        bool posPion = TMath::Abs(daughter0_TPCNSigmaPion) <= 3 && (TMath::Abs(daughter0_TOFNSigmaPion) <= 3 || daughter0_TOFNSigmaPion == 1000)

        if(TMath::Abs(TPCNSigmaPion) <= 3 && (TMath::Abs(TOFNSigmaPion) <= 3 || TOFNSigmaPion == 1000)) {

            if(track->Charge() == 1){
                piPlus_list.push_back(track);
            }
            else {
                piMinus_list.push_back(track);
            }
        }

        if(TMath::Abs(TPCNSigmaKaon) <= 2 && (TMath::Abs(TOFNSigmaKaon) <= 2 || TOFNSigmaKaon == 1000)) {

            if(track->Charge() == 1){
                kPlus_list.push_back(track);
            }
            else {
                kMinus_list.push_back(track);
            }
        }

        dist_points[3] = d0->InvMassD0();
        fDist->Fill(dist_points);
    }

}

void AliAnalysisTaskD0Mass::FillTriggeredSingleParticleDist(std::vector<AliAnalysisTaskD0Mass::AliMotherContainer> particle_list, THnSparse* fDist, double maxTriggerPt)
{

    double dist_points[5]; //Pt, Phi, Eta, Mass, max trigger pt in event
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].particle;
        dist_points[0] = particle.Pt();
        dist_points[1] = particle.Phi();
        dist_points[2] = particle.Eta();
        dist_points[3] = particle.M();
        dist_points[4] = maxTriggerPt;
        fDist->Fill(dist_points);
    }

}

AliAnalysisTaskD0Mass::AliMotherContainer AliAnalysisTaskD0Mass::DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{

    AliAnalysisTaskD0Mass::AliMotherContainer mom;
    mom.particle.SetPx(track1->Px() + track2->Px());
    mom.particle.SetPy(track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}

bool AliAnalysisTaskD0Mass::PassDaughterCuts(AliAODTrack *track){

    bool pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && (track->IsOn(AliAODTrack::kTPCrefit));

    pass = pass && (track->GetTPCCrossedRows() > 70);

    float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    pass = pass && (ratio > 0.8);

    return pass;
}

bool AliAnalysisTaskD0Mass::PassTriggerCuts(AliAODTrack *track){

    bool pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && track->TestBit(AliAODTrack::kIsHybridGCG);

    return pass;
}

void AliAnalysisTaskD0Mass::UserExec(Option_t*)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        std::cout << "THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!" << std::endl;
        std::abort();
    }

    // std::cout << fAOD->GetList()->FindObject("D0toKpi") << std::endl;

    AliAnalysisManager* currentMgr = AliAnalysisManager::GetAnalysisManager();
    if(!currentMgr) {
        AliFatal("NO MANAGER CONNECTED, EXITING");
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

    std::vector<AliAODTrack*> kPlus_list;
    std::vector<AliAODTrack*> kMinus_list;
    std::vector<AliAODTrack*> piPlus_list;
    std::vector<AliAODTrack*> piMinus_list;
    std::vector<AliAODTrack*> trigger_list;

    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;


        if(PassTriggerCuts(track)) trigger_list.push_back(track);

        if(PassDaughterCuts(track)) {
            double TPCNSigmaPion = 1000;
            double TOFNSigmaPion = 1000;

            TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
            TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);

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
    std::vector<AliAnalysisTaskD0Mass::AliMotherContainer> D0_list;



    for(int i = 0; i < (int)piMinus_list.size(); i++) {
        for(int j = 0; j < (int) kPlus_list.size(); j++) {
            auto pion = piMinus_list[i];
            auto kaon = kPlus_list[j];

            AliMotherContainer D0 = DaughtersToMother(pion, kaon, 0.1396, 0.493677);
            D0_list.push_back(D0);
        }
    }

    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = 0; j < (int) kMinus_list.size(); j++) {
            auto pion = piPlus_list[i];
            auto kaon = kMinus_list[j];

            AliMotherContainer D0 = DaughtersToMother(pion, kaon, 0.1396, 0.493677);
            D0_list.push_back(D0);
        }
    }

    // double maxTriggerPt = FindMaxTriggerPt(D0_list, trigger_list);
    double maxTriggerPt = 3;

    FillSingleParticleDist(D0_list, fD0Dist);
    FillTriggeredSingleParticleDist(D0_list, fD0Dist, maxTriggerPt);
    
    // THIS IS THE HF VERTEXING SECTION!!!!!!!!!!!!!!!!!!!

    TClonesArray* hf_D0_array = 0x0;
    TString arrayName = "D0toKpi";
    hf_D0_array = (TClonesArray*)fAOD->GetList()->FindObject(arrayName.Data());
    if(!hf_D0_array) {
        std::cout << "NO D0 ARRAY FOUND" << std::endl;
        return;
    }

    // int nD0s_hf = hf_D0_array->GetEntriesFast();

    // for(int iD0 = 0; iD0 < nD0s_hf; iD0++) {
    //     AliAODRecoDecayHF2Prong* D0 = (AliAODRecoDecayHF2Prong *)hf_D0_array->UncheckedAt(iD0);
    //     AliAODTrack* daughterFirst = (AliAODTrack*)D0->GetDaughter(0);
    //     AliAODTrack* daughterSecond = (AliAODTrack*)D0->GetDaughter(1);
    //     std::cout << daughterFirst->GetMostProbablePID() << ", " << daughterFirst->Charge() << " is D1 (pid, charge)" << std::endl;
    //     std::cout << daughterSecond->GetMostProbablePID() << ", " << daughterSecond->Charge() << " is D2 (pid, charge)" << std::endl;
    // }


    //Making list of possible D0s (from hf finder):
    std::vector<AliAnalysisTaskD0Mass::AliMotherContainer> hf_D0_list;


    //  maxTriggerPt = FindMaxTriggerPt(hf_D0_list, trigger_list);

    // FillSingleParticleDist(hf_D0_list, fHFD0Dist);
    FillHFD0Dist(hf_D0_array, fHFD0Dist);
    FillTriggeredSingleParticleDist(hf_D0_list, fTriggeredHFD0Dist, maxTriggerPt);

    PostData(1, fOutputList);

}

void AliAnalysisTaskD0Mass::Terminate(Option_t *option)
{
}