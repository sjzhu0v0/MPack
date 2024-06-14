#ifndef MPHYSICS_H
#define MPHYSICS_H

#include "MSystem.h"

const double mass_proton = 0.938272;
const double mass_kaon = 0.493677;
const double mass_pion = 0.139570;
const double mass_muon = 0.105658;
const double mass_electron = 0.000511;
const double mass_d0 = 1.86483;
const double mass_lambda_c = 2.2846;
const double mass_jpsi = 3.0969;

namespace TrentoKinematics {
TLorentzVector p4BeamHadron;
TLorentzVector p4BeamElectron;
double mass_nucleon = mass_proton;

double Pt_hadron(TLorentzVector p4ScatteredLepton, TLorentzVector p4Hadron) {
  TLorentzVector q4 = p4BeamElectron - p4ScatteredLepton;
  double gamma = mass_nucleon * q4.Mag() / (p4Hadron.Dot(q4));
  double q4_p4Hadron = q4.Dot(p4Hadron);
  double p4Hadron_p4BeamHadron = p4Hadron.Dot(p4BeamHadron);
  double q4_p4BeamHadron = q4.Dot(p4BeamHadron);

  double part1 = p4Hadron.Mag2();
  double part2 = -2. * q4_p4Hadron * p4Hadron_p4BeamHadron / (q4_p4BeamHadron) /
                 (1 + gamma * gamma);
  double part3 = -(gamma * gamma / (1 + gamma * gamma)) *
                 (q4_p4Hadron * q4_p4Hadron / q4.Mag2() +
                  p4Hadron_p4BeamHadron * p4Hadron_p4BeamHadron / mass_nucleon /
                      mass_nucleon);
  double pt_hadron = sqrt(-1 * (part1 + part2 + part3));
  return pt_hadron;
}
} // namespace TrentoKinematics


#endif
