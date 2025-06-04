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

  // reaction 6;
  pres_mod[0] = m + 1.5 * C[0] + 11.0 * C[6] - 1.0 * C[1];

  // reaction 7;
  pres_mod[1] = m + 1.5 * C[0] + 11.0 * C[6];

  // reaction 8;
  pres_mod[2] = m + 1.5 * C[0] + 11.0 * C[6] - 1.0 * C[1];

  // reaction 9;
  pres_mod[3] = m + 2.0 * C[0] - 1.0 * C[6] + 1.0 * C[25];

  // reaction 11;
  thd = m + 1.0 * C[0] + 13.0 * C[6] - 0.21999999999999997 * C[1];
  k0 = exp(3.4087162630776540e+01 - 1.72 * logT - (2.6408962763808859e+02 / T));
  kinf = exp(2.2270828345662423e+01 + 0.44 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.00000000e-01 * exp(-T / 1.00000000e-30) + 5.00000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[4] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 24;
  thd = m + 6.5 * C[6] + 0.5 * C[25] + 0.19999999999999996 * C[1] + 6.7 * C[7] + 2.7 * C[0];
  k0 = exp(4.9270577684749114e+01 - 2.3 * logT - (2.4531450567319323e+04 / T));
  kinf = exp(2.8324168296488494e+01 + 0.9 * logT - (2.4531450567319323e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.70000000e-01 * exp(-T / 1.00000000e-30) + 4.30000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[5] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 30;
  pres_mod[6] = m;

  // reaction 55;
  pres_mod[7] = m;

  // reaction 82;
  pres_mod[8] = m + 9.0 * C[6];

  // reaction 83;
  pres_mod[9] = m;

  // reaction 93;
  thd = m + 0.6000000000000001 * C[25];
  k0 = exp(1.9296149481306266e+01 + 0.206 * logT - (-7.7999032553170230e+02 / T));
  kinf = exp(2.8036486224036711e+01 - 0.41 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(1.80000000e-01 * exp(-T / 1.00000000e-30) + 8.20000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[10] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 101;
  thd = m;
  k0 = exp(4.2998340473490288e+01 - 2.87 * logT - (7.7999032553170230e+02 / T));
  kinf = exp(2.7893385380396040e+01 - 0.75 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.50000000e-01 * exp(-T / 1.00000000e+03) + 7.50000000e-01 * exp(-T / 1.00000000e+05) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[11] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 102;
  pres_mod[12] = m;

  // reaction 108;
  thd = m;
  k0 = exp(4.0365366298828434e+01 - 2.5 * logT);
  kinf = exp(2.5423746202738826e+01 - 0.3 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.50000000e-01 * exp(-T / 1.00000000e-30) + 7.50000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[13] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 111;
  thd = m;
  k0 = exp(4.4826845844638555e+01 - 3.0 * logT);
  kinf = 30000000000.0;
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.00000000e-01 * exp(-T / 1.00000000e-30) + 4.00000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[14] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 118;
  thd = m + 0.7 * C[25] + 0.3999999999999999 * C[1] + 11.0 * C[6];
  k0 = exp(2.7218531392883420e+01 - (2.8935124979401859e+04 / T));
  kinf = exp(2.5318385687081001e+01 - (2.9166605979237072e+04 / T));
  Pr = k0 * thd / kinf;
  pres_mod[15] =  Pr / (1.0 + Pr);

  // reaction 131;
  thd = m;
  k0 = exp(6.4942386233079020e+01 - 5.49 * logT - (9.9989727537515637e+02 / T));
  kinf = exp(2.7051202620675607e+01 - 0.414 * logT - (3.3212491280704739e+01 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.90000000e-01 * exp(-T / 1.00000000e-30) + 3.10000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[16] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 169;
  pres_mod[17] = m;

  // reaction 170;
  pres_mod[18] = m;

} // end get_rxn_pres_mod

