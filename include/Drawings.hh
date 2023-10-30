#ifndef __DRAWINGS__
#define __DRAWINGS__

#include <vector>
#include <TGraph.h>
#include <TString.h>

void DrawInit();

void GraphSetting(TGraph& graph);

void DrawPath(std::vector<double> zets, std::vector<double> thes, TString title);
void DrawTaZe(std::vector<double> taus, std::vector<double> zets, TString title);
void DrawTaTh(std::vector<double> taus, std::vector<double> thes, TString title);

void DrawAction(std::vector<double> indices, std::vector<double> actions_ori, std::vector<double> actions_mod);

#endif