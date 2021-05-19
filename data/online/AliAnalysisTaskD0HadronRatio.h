#ifndef AliAnalysisTaskD0HadronRatio_H
#define AliAnalysisTaskD0HadronRatio_H

//All relevant AliPhysics includes (this list will continue to grow):
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1F.h"
#include "TH3D.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliPID.h"
#include "TChain.h"
#include "TVector.h"
#include "AliEventPoolManager.h"
#include "AliCFParticle.h"

//These includes probably aren't necessary but I feel like it's proper etiquette
#include <vector>
#include <iostream>

class AliAnalysisTaskD0HadronRatio : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskD0HadronRatio();
  AliAnalysisTaskD0HadronRatio(const char *name);
  virtual ~AliAnalysisTaskD0HadronRatio();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);

  struct AliMotherContainer {
    TLorentzVector particle;
    int daughter1ID;
    int daughter2ID;
  };

 private:

  float MULT_LOW;
  float MULT_HIGH;
  float DAUGHTER_TRK_BIT;
  float ASSOC_TRK_BIT;
  float TRIG_TRK_BIT;

  TString CENT_ESTIMATOR;
  AliAODEvent* fAOD; //!>! input event
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager

  TH2D* fTriggersAndD0sPerEvent_All; //!>! triggers and all D0s per event
  TH2D* fTriggersAndD0sPerEvent_2_4; //!>! triggers and 2-4 GeV D0s per event

  THnSparseF* fLooseDist;  //!>! single particle all hadron dist (no cuts at all)
  THnSparseF* fTriggerDist;  //!>! single particle trigger dist
  THnSparseF* fAssociatedHDist;  //!>! single particle associated hadron dist
  THnSparseF* fD0Dist;  //!>! single particle D0 dist

  THnSparseF* fDphiHD0;  //!>! hadron-D0 correlation hist
  THnSparseF* fDphiHD0V0;  //!>! hadron-D0 correlation hist (using v0 finder for D0)
  THnSparseF* fDphiHD0Rotated;  //!>! hadron-D0 correlation hist with rotated pion
  THnSparseF* fDphiHD0RotatedPi;  //!>! hadron-D0 correlation hist with daughter rotated by pi
  THnSparseF* fDphiHD0RotatedKaon;  //!>! hadron-D0 correlation hist with rotated Kaon
  THnSparseF* fDphiHD0Flipped;  //!>! hadron-D0 correlation hist with flipped pion
  THnSparseF* fDphiHH;   //!>! hadron-hadron correlation hist
  THnSparseF* fDphiTriggerTrigger;   //!>! trigger-trigger correlation hist
  THnSparseF* fDphiHD0LS; //!>! hadron-Kaon+pion like sign correlation hist
  THnSparseF* fDphiHD0Mixed; //!>! hadron-D0 mixed correlation hist
  THnSparseF* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparseF* fDphiHD0LSMixed; //!>! hadron-Kaon+pion like sign mixed correlation hist
  THnSparseF* fDphiTriggerTriggerMixed;   //!>! mixed trigger-trigger correlation hist

  THnSparseF* fD0DaughterDCA;

  // THnSparseF* fPid; //!>! histogram to visualize pid cuts
  // THnSparseF* fSignalAnalysis; //!>! histogram to analyze signal with nsigma cuts

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection

  //hand written functions:

  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  AliMotherContainer RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle);
  AliMotherContainer FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist);
  void FillSingleParticleDist(std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> particle_list, double zVtx, THnSparse* fDist);
  void MakeSameHD0Correlations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHD0Correlations(AliEventPool *fPool, std::vector<AliAnalysisTaskD0HadronRatio::AliMotherContainer> D0_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, bool eff=true);
  bool PassDaughterCuts(AliAODTrack *track);
  bool PassTriggerCuts(AliAODTrack *track);
  bool PassAssociatedCuts(AliAODTrack *track);

  ClassDef(AliAnalysisTaskD0HadronRatio, 0);

};
#endif