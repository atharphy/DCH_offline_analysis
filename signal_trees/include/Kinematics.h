#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "TLorentzVector.h"
#include "TMath.h"
#include <vector>
#include <set>
#include <cmath>
#include <string>
#include <algorithm>

using std::string;
using std::vector;

vector<double> SortPt();
string applyHEMveto(string cat);
void applyTauES(string cat);
double ST(string cat);
int remaining_idx(vector<std::pair<int,int>>& pair, string cat_name);
double dPhi(double phi1, double phi2);
double deltaPhi(const TLorentzVector& v1, const TLorentzVector& v2);
double getDR(double eta1, double phi1, double eta2, double phi2);
double deltaR(const TLorentzVector& v1, const TLorentzVector& v2);
double calculateMT(const TLorentzVector& V, const TLorentzVector& MET);
double calculateMTtot(const TLorentzVector& l1, const TLorentzVector& l2);

struct Lepton {
    double pt, eta, phi, mass;
    int q;
    double d0, dZ, iso;
    Lepton(double pt, double eta, double phi, double mass, int q, double d0, double dZ, double iso)
        : pt(pt), eta(eta), phi(phi), mass(mass), q(q), d0(d0), dZ(dZ), iso(iso) {}
};

bool isDuplicate(const Lepton& l1, const Lepton& l2);
string pairFunc(int m, int n, string cat, double Zwindow);
vector<int> ZCandMaker(string cat, double Zwindow);
vector<int> ZVetoMaker(string cat, double Zwindow);
TLorentzVector* ZCandMaker_pair(string cat, TLorentzVector l1, TLorentzVector l2,
                                TLorentzVector l3, TLorentzVector l4, double Zwindow);
TLorentzVector* ZVetoMaker_pair(string cat, TLorentzVector l1, TLorentzVector l2,
                                TLorentzVector l3, TLorentzVector l4, double Zwindow);

#endif
