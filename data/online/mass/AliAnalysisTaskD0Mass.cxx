#include "AliAnalysisTaskD0Mass.h"

//BOTH OF THESE ARE WEIRD, BUT APPARENTLY NECESSARRY
class AliAnalysisTaskD0Mass;
ClassImp(AliAnalysisTaskD0Mass);


AliAnalysisTaskD0Mass::AliAnalysisTaskD0Mass() :
    AliAnalysisTaskSE(),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fD0Dist{0}
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
    fD0Dist{0}
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

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && (track->IsOn(AliAODTrack::kTPCrefit));

    pass = pass && (track->GetTPCCrossedRows() > 70);

    float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    pass = pass && (ratio > 0.8);

    return pass;
}

void AliAnalysisTaskD0Mass::UserExec(Option_t*)
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

    std::vector<AliAODTrack*> kPlus_list;
    std::vector<AliAODTrack*> kMinus_list;
    std::vector<AliAODTrack*> piPlus_list;
    std::vector<AliAODTrack*> piMinus_list;

    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;

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

            if(TMath::Abs(TPCNSigmaKaon) <= 2 && (TMath::Abs(TOFNSigmaKaon) <= 2 || TOFNSigmaKaon == 1000)) {

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

    FillSingleParticleDist(D0_list, fD0Dist);
    PostData(1, fOutputList);

}

void AliAnalysisTaskD0Mass::Terminate(Option_t *option)
{
}