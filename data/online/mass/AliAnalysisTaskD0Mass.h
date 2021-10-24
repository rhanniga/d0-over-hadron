#ifndef AliAnalysisTaskD0Mass_H
#define AliAnalysisTaskD0Mass_H

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
#include "AliAODHandler.h"
#include "AliPID.h"
#include "TChain.h"
#include "TVector.h"
#include "AliEventPoolManager.h"
#include "AliCFParticle.h"
#include "AliAnalysisManager.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliPhysicsSelectionTask.h"
#include "TSystem.h"
#include "AliAnalysisAlien.h"

//These includes probably aren't necessary but I feel like it's proper etiquette
#include <vector>
#include <iostream>

class AliAnalysisTaskD0Mass : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskD0Mass();
  AliAnalysisTaskD0Mass(const char *name);
  virtual ~AliAnalysisTaskD0Mass();
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

  THnSparseF* fD0Dist;  //!>! single particle D0 dist
  THnSparseF* fD0V0Dist;  //!>! single particle D0 dist using V0 finder
  THnSparseF* fTriggeredD0Dist;  //!>! single particle D0 dist with trigger pt axis
  THnSparseF* fHFD0Dist;  //!>! single particle D0 dist (from HF)
  THnSparseF* fTriggeredHFD0Dist;  //!>! single particle D0 dist (from HF) with trigger pt axis

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection

  //hand written functions:
  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliMotherContainer> particle_list, THnSparse* fDist);
  void FillHFD0Dist(TClonesArray* hf_d0s, THnSparse* fDist);
  void FillTriggeredSingleParticleDist(std::vector<AliMotherContainer> particle_list, THnSparse* fDist, double maxTriggerPt);
  double FindMaxTriggerPt(std::vector<AliMotherContainer> particle_list, std::vector<AliAODTrack*> trigger_list);
  bool PassDaughterCuts(AliAODTrack *track);
  bool PassTriggerCuts(AliAODTrack *track);

  ClassDef(AliAnalysisTaskD0Mass, 0);

};
#endif