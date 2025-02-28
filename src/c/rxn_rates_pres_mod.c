#include <math.h>
#include "header.h"
#include "rates.h"

void get_rxn_pres_mod (const double T, const double pres, const double * __restrict__ C, double * __restrict__ pres_mod) {
  // third body variable declaration
  double thd;

  // pressure dependence variable declarations
  double k0;
  double kinf;
  double Pr;

  // troe variable declarations
  double logFcent;
  double A;
  double B;

  double logT = log(T);
  double m = pres / (8.31446210e+03 * T);

  // reaction 8;
  pres_mod[0] = m + 1.5 * C[1] + 11.0 * C[5];

  // reaction 9;
  pres_mod[1] = m + 1.5 * C[1] + 11.0 * C[5];

  // reaction 10;
  pres_mod[2] = m + 1.5 * C[1] + 11.0 * C[5] - 0.17000000000000004 * C[8];

  // reaction 11;
  pres_mod[3] = m + 1.5 * C[1] + 11.0 * C[5] - 0.17000000000000004 * C[8];

  // reaction 12;
  pres_mod[4] = m + 1.5 * C[1] + 11.0 * C[5] - 0.25 * C[8];

  // reaction 13;
  pres_mod[5] = m + 1.5 * C[1] + 11.0 * C[5] - 0.25 * C[8];

  // reaction 14;
  pres_mod[6] = m - 0.27 * C[1] + 11.0 * C[5] - 0.62 * C[8];

  // reaction 15;
  pres_mod[7] = m - 0.27 * C[1] + 11.0 * C[5] - 0.62 * C[8];

  // reaction 16;
  thd = m + 0.30000000000000004 * C[1] + 13.0 * C[5] - 0.32999999999999996 * C[8];
  k0 = exp(2.4273457771286456e+01 - 0.411 * logT - (-5.6108981481796650e+02 / T));
  kinf = exp(2.1111923826738195e+01 + 0.6 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.00000000e-01 * exp(-T / 1.00000000e-30) + 5.00000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+100 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[8] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 29;
  thd = m + 1.5 * C[1] + 11.0 * C[5] - 0.36 * C[8];
  k0 = exp(3.2420178138029655e+01 - (2.2896490201091907e+04 / T));
  kinf = exp(3.3318335397877441e+01 - (2.4370923526129252e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.00000000e-01 * exp(-T / 1.00000000e-30) + 5.00000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+100 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[9] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

} // end get_rxn_pres_mod

