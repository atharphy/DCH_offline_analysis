#pragma once

#include <string>
#include <vector>

const int NlepMax = 4;

const bool APPLY_OFFICIAL_MET_CORRECTION = true;
const bool APPLY_ZPT_REWEIGHTING = true;
const bool APPLY_DY_RECOIL_CORRECTION = true;
const bool APPLY_WJ_RECOIL_CORRECTION = false;
const bool USE_CUSTOM_RECOIL_CORRECTIONS = true;

const std::vector<std::string> finalStates = {
    "ee","em","et","mm","mt","tt","eee","eem","eet","eme","emm","emt","ete","etm","ett","mme","mmm","mmt","mte","mtm","mtt","tte","ttm","ttt",
    "eeee","eeem","eeet","eemm","eemt","eett","emem","emet","emmm","emmt","emtt","etet","etmm","etmt","ettt","mmmm","mmmt","mmtt","mtmt","mttt","tttt"
};
