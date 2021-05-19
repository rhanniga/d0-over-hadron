#include <iostream>

AliAnalysisTaskD0HadronRatio* AddD0HadronRatioTask(TString name = "D0HadronRatio") {

  AliAnalysisManager *manage = AliAnalysisManager::GetAnalysisManager();

  if (!manage) return 0x0;

  if(!manage->GetInputEventHandler()) return 0x0;



  TString file_name = AliAnalysisManager::GetCommonFileName();

  AliAnalysisTaskD0HadronRatio* task = new AliAnalysisTaskD0HadronRatio(name.Data());

  if(!task) return 0x0;

  manage->AddTask(task);

  manage->ConnectInput(task, 0, manage->GetCommonInputContainer());
  manage->ConnectOutput(task, 1, manage->CreateContainer("h-D0", TList::Class(), AliAnalysisManager::kOutputContainer, file_name.Data()));

  return task;

}
