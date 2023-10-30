#include <vector>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TMath.h>
#include "Drawings.hh"
#include "PathVariables.hh"

void DrawInit()
{
    gSystem->Unlink("Action.png");
    gSystem->Unlink("Path.gif");
    gSystem->Unlink("Zeta.gif");
    gSystem->Unlink("Theta.gif");
}

void GraphSetting(TGraph& graph)
{
    graph.SetMarkerStyle(8);
    graph.SetMarkerSize(0.8);
    graph.SetMarkerColor(kBlack);
    graph.SetLineWidth(2);
    graph.SetLineColor(kRed);
}

void DrawPath(std::vector<double> zets, std::vector<double> thes, TString title)
{
    static TCanvas canv("Path", "Path", 1000, 1000);
    canv.cd();

    const double maxwindow = 1.1 * TMath::Max(zet_ini, zet_fin);
    static TH2D axis("axis_path", "axis_path", 10, -maxwindow, maxwindow, 10, -maxwindow, maxwindow);
    axis.SetTitle(title);
    axis.SetStats(false);
    axis.Draw();

    std::vector<double> xs, ys;
    for (int i = 0; i < nbin + 1; i++)
    {
        xs.push_back(zets[i] * TMath::Cos(thes[i]));
        ys.push_back(zets[i] * TMath::Sin(thes[i]));
    }

    TGraph graph(nbin + 1, &xs[0], &ys[0]);
    GraphSetting(graph);
    graph.Draw("lp");
    
    canv.Modified();
    canv.Update();
    canv.Print("Path.gif+");
}

void DrawTaZe(std::vector<double> taus, std::vector<double> zets, TString title)
{
    static TCanvas canv("Zeta", "Zeta", 1500, 1000);
    canv.cd();

    TGraph graph(nbin + 1, &taus[0], &zets[0]);
    graph.SetTitle(title);
    GraphSetting(graph);
    graph.Draw("alp");

    canv.Modified();
    canv.Update();
    canv.Print("Zeta.gif+");
}

void DrawTaTh(std::vector<double> taus, std::vector<double> thes, TString title)
{
    static TCanvas canv("Theta", "Theta", 1500, 1000);
    canv.cd();

    TGraph graph(nbin + 1, &taus[0], &thes[0]);
    graph.SetTitle(title);
    GraphSetting(graph);
    graph.Draw("alp");

    canv.Modified();
    canv.Update();
    canv.Print("Theta.gif+");
}

void DrawAction(std::vector<double> indices, std::vector<double> actions_ori, std::vector<double> actions_mod)
{
    static TCanvas canv("Action", "Action", 1500, 1000);
    canv.cd();

    TLegend legend(0.7, 0.85, 0.95, 0.95);

    TGraph graph_ori(indices.size(), &indices[0], &actions_ori[0]);
    TGraph graph_mod(indices.size(), &indices[0], &actions_mod[0]);

    graph_mod.SetTitle("Action;try;Action");

    graph_ori.SetLineWidth(1);
    graph_ori.SetLineColor(kRed);
    graph_mod.SetLineWidth(1);
    graph_mod.SetLineColor(kBlue);

    graph_mod.Draw("al");
    graph_ori.Draw("l");

    legend.AddEntry(&graph_mod, "Current Action Value");
    legend.AddEntry(&graph_ori, "Least Action Value");
    legend.Draw();

    canv.Modified();
    canv.Update();
    canv.Print("Action.png");
}