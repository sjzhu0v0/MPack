#ifndef __MFormula_Flow__
#define __MFormula_Flow__

#include "MMath.h"

namespace MFormula {
namespace Flow {
double (*c2_2_1)(vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[0].Rho();
};
double (*v2_2_1)(vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return sqrt(vec_val[0].Rho());
};
double (*c2_4_1)(vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[1].Rho() - 2 * vec_val[0].Rho2();
};
double (*v2_4_1)(vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return sqrt(sqrt(-1. * (vec_val[1].Rho() - 2 * vec_val[0].Rho2())));
};
double (*c2_6_1)(vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[2].Rho() - 9 * vec_val[1].Rho() * vec_val[0].Rho() +
         12 * vec_val[0].Rho2() * vec_val[0].Rho();
};
double (*v2_6_1)(vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  double c = vec_val[2].Rho() - 9 * vec_val[1].Rho() * vec_val[0].Rho() +
             12 * vec_val[0].Rho2() * vec_val[0].Rho();
  return pow(c / 4., 1. / 6.);
};

double (*c2_4_3)(vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[2].Rho() - 2 * vec_val[1].Rho() * vec_val[0].Rho();
};
double (*v2_4_3)(vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return sqrt(sqrt(-1. *
              (vec_val[2].Rho() - 2 * vec_val[1].Rho() * vec_val[0].Rho())));
};

double (*Cp_mf_2) (vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[0].Rho();
};

double (*Cf_mp_2) (vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[0].Rho();
};

double (*Cp_mf_4) (vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[2].Rho() - 2 * vec_val[0].Rho()*vec_val[1].Rho();
};

double (*Cf_mp_4) (vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[2].Rho() - 2 * vec_val[0].Rho()*vec_val[1].Rho();
};

double (*Cmp_ff_mf_4) (vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[2].Rho() - 2 * vec_val[0].Rho()*vec_val[1].Rho();
};

double (*Cm_pf_mf_4) (vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[4].Rho() - vec_val[3].Rho()*vec_val[0].Rho() - vec_val[1].Rho()*vec_val[2].Rho();
};

double (*C_m_ff_mp_4) (vector<TComplex> vec_val) = [](vector<TComplex> vec_val) {
  return vec_val[2].Rho() - 2 * vec_val[0].Rho()*vec_val[1].Rho();
};
} // namespace Flow
}; // namespace MFormula

#endif