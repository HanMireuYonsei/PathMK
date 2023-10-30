#include <tuple>
#include <vector>
#include <TTree.h>
#include <TMath.h>
#include "CalcAction.hh"

std::tuple<double, double> CalcAction(const std::vector<double>& taus, const std::vector<double>& zets, const std::vector<double>& thes)
{
    double dtau, dzet, dthe;
    double kin = 0.;
    double pot = 0.;

    for (size_t i = 1; i < taus.size(); i++)
    {
        dtau = taus[i] - taus[i - 1];
        dzet = zets[i] - zets[i - 1];
        dthe = thes[i] - thes[i - 1];

        kin += 0.5 * (TMath::Power(dzet / dtau, 2) + TMath::Power(zets[i - 1] * dthe / dtau, 2)) * dtau;
        pot += dtau / zets[i - 1];
    }
    
    return std::tuple<double, double>(kin, pot);
}

std::tuple<double, double> CalcAction(TTree* tree)
{
    double tau_old, zet_old, the_old, tau, zet, the, dtau, dzet, dthe;
    double kin = 0.;
    double pot = 0.;

    tree->SetBranchAddress("tau"  , &tau);
    tree->SetBranchAddress("zeta" , &zet);
    tree->SetBranchAddress("theta", &the);
    tree->GetEntry(0);

    for (Long64_t i = 1; i < tree->GetEntries(); i++)
    {
        tau_old = tau;
        zet_old = zet;
        the_old = the;
        tree->GetEntry(i);
        dtau = tau - tau_old;
        dzet = zet - zet_old;
        dthe = the - the_old;

        kin += 0.5 * (TMath::Power(dzet / dtau, 2) + TMath::Power(zet_old * dthe / dtau, 2)) * dtau;
        pot += dtau / zet_old;
    }

    return std::tuple<double, double>(kin, pot);
}
