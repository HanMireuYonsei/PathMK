#ifndef __CALCACTION__
#define __CALCACTION__

#include <tuple>
#include <vector>
#include <TTree.h>

std::tuple<double, double> CalcAction(const std::vector<double>& taus, const std::vector<double>& zets, const std::vector<double>& thes);
std::tuple<double, double> CalcAction(TTree* tree);

#endif