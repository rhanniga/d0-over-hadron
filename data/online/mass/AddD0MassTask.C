#include <iostream>

AliAnalysisTaskD0Mass* AddD0MassTask(TString name = "D0Mass") {

  AliAnalysisManager *manage = AliAnalysisManager::GetAnalysisManager();

  if (!manage) return 0x0;

  if(!manage->GetInputEventHandler()) return 0x0;



  TString file_name = AliAnalysisManager::GetCommonFileName();

  AliAnalysisTaskD0Mass* task = new AliAnalysisTaskD0Mass(name.Data());

  if(!task) return 0x0;

  manage->AddTask(task);

  manage->ConnectInput(task, 0, manage->GetCommonInputContainer());
  manage->ConnectOutput(task, 1, manage->CreateContainer("h-D0", TList::Class(), AliAnalysisManager::kOutputContainer, file_name.Data()));

  return task;

}
