#include <tuple>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TString.h>
#include "CalcAction.hh"
#include "Drawings.hh"
#include "PathVariables.hh"

int main(const int argc, const char *argv[])
{
    gRandom->SetSeed(100);

    TF1 linear_tau("linear_tau", "[0] * x + [1]", tau_ini, tau_fin);
    TF1 linear_zet("linear_zet", "[0] * x + [1]", zet_ini, zet_fin);
    TF1 linear_the("linear_the", "[0] * x + [1]", the_ini, the_fin);
    linear_tau.SetParameter(0, (tau_fin - tau_ini) / nbin);
    linear_tau.SetParameter(1, tau_ini);
    linear_zet.SetParameter(0, (zet_fin - zet_ini) / nbin);
    linear_zet.SetParameter(1, zet_ini);
    linear_the.SetParameter(0, (the_fin - the_ini) / nbin);
    linear_the.SetParameter(1, the_ini);

    std::vector<double> taus, zets_ori, thes_ori, zets_mod, thes_mod, actions_ori, actions_mod, indices;

    for (int i = 0; i < nbin + 1; i++)
    {
        taus.push_back(linear_tau.Eval(i));
        zets_ori.push_back(linear_zet.Eval(i));
        thes_ori.push_back(linear_the.Eval(i));
        zets_mod.push_back(linear_zet.Eval(i));
        thes_mod.push_back(linear_the.Eval(i));
    }

    std::tuple<double, double> act_tup_ori = CalcAction(taus, zets_ori, thes_ori);
    std::tuple<double, double> act_tup_mod = CalcAction(taus, zets_mod, thes_mod);
    double action_ori = std::get<0>(act_tup_ori) + std::get<1>(act_tup_ori);
    double action_mod = std::get<0>(act_tup_mod) + std::get<1>(act_tup_mod);
    actions_ori.push_back(action_ori);
    actions_mod.push_back(action_mod);
    indices.push_back(0);

    DrawInit();

    DrawPath(zets_mod, thes_mod, TString::Format("try 0: action = %f + %f = %f;x;y"      , std::get<0>(act_tup_mod), std::get<1>(act_tup_mod), action_mod));
    DrawTaZe(taus    , zets_mod, TString::Format("try 0: action = %f + %f = %f;tau;zeta" , std::get<0>(act_tup_mod), std::get<1>(act_tup_mod), action_mod));
    DrawTaTh(taus    , thes_mod, TString::Format("try 0: action = %f + %f = %f;tau;theta", std::get<0>(act_tup_mod), std::get<1>(act_tup_mod), action_mod));

    std::cout << "==========MK Start==========" << std::endl;
    int ntried = 0;
    int nretried = 0;
    while (nretried < limretry)
    {
        ntried++;
        nretried++;

        if (action_ori <= action_mod)
        {
            for (int i = 0; i < nbin + 1; i++)
            {
                zets_mod[i] = zets_ori[i] + gRandom->Gaus(0., sigma * TMath::Power(fowdec, nretried) * double(i) * double(nbin - i) * 4. / double(nbin) / double(nbin));
                thes_mod[i] = thes_ori[i] + gRandom->Gaus(0., sigma * TMath::Power(fowdec, nretried) * double(i) * double(nbin - i) * 4. / double(nbin) / double(nbin));
            }
        }
        else
        {
            for (int i = 0; i < nbin + 1; i++)
            {
                double zet_buf, the_buf;

                zet_buf = zets_mod[i] + (zets_mod[i] - zets_ori[i]);
                the_buf = thes_mod[i] + (thes_mod[i] - thes_ori[i]);

                zets_ori[i] = zets_mod[i];
                thes_ori[i] = thes_mod[i];
                
                zets_mod[i] = zet_buf;
                thes_mod[i] = the_buf;
            }
        }

        act_tup_ori = CalcAction(taus, zets_ori, thes_ori);
        act_tup_mod = CalcAction(taus, zets_mod, thes_mod);
        action_ori = std::get<0>(act_tup_ori) + std::get<1>(act_tup_ori);
        action_mod = std::get<0>(act_tup_mod) + std::get<1>(act_tup_mod);
        actions_ori.push_back(action_ori);
        actions_mod.push_back(action_mod);
        indices.push_back(ntried);

        if (action_ori > action_mod)
        {
            std::cout << "try " << nretried << " / " << ntried << std::endl;

            DrawPath(zets_mod, thes_mod, TString::Format("try %d: action = %f + %f = %f;x;y"      , ntried, std::get<0>(act_tup_mod), std::get<1>(act_tup_mod), action_mod));
            DrawTaZe(taus    , zets_mod, TString::Format("try %d: action = %f + %f = %f;tau;zeta" , ntried, std::get<0>(act_tup_mod), std::get<1>(act_tup_mod), action_mod));
            DrawTaTh(taus    , thes_mod, TString::Format("try %d: action = %f + %f = %f;tau;theta", ntried, std::get<0>(act_tup_mod), std::get<1>(act_tup_mod), action_mod));

            nretried = 0;
        }
    }

    DrawAction(indices, actions_ori, actions_mod);
    
    TTree tree("Path", "Path");
    double tau, zet, the;
    tree.Branch("tau", &tau);
    tree.Branch("zeta", &zet);
    tree.Branch("theta", &the);

    for (int i = 0; i < nbin + 1; i++)
    {
        tau = taus[i];
        zet = zets_ori[i];
        the = thes_ori[i];

        tree.Fill();
    }

    TFile file("Path.root", "recreate");
    tree.Write();
    
    return 0;
}