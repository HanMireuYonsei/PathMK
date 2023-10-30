#ifndef __PATHVARIABLES__
#define __PATHVARIABLES__

#include <TMath.h>

constexpr int nbin = 15;
constexpr int limretry = 1000;
constexpr double sigma = 0.01;
constexpr double fowdec = 0.99;

constexpr double tau_ini = 0.;
constexpr double tau_fin = 1.0188 * TMath::Pi();
constexpr double tau_half = (tau_fin - tau_ini) / 2.;
constexpr double zet_ini = 0.9;
constexpr double zet_fin = 1.125;
constexpr double the_ini = 0.;
constexpr double the_fin = TMath::Pi();

#endif