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
  pres_mod[0] = m + 1.5 * C[2] + 11.0 * C[8] - 1.0 * C[1] - 1.0 * C[0] - 1.0 * C[3];

  // reaction 9;
  pres_mod[1] = m + 1.5 * C[2] + 11.0 * C[8] - 1.0 * C[1] - 1.0 * C[0];

  // reaction 12;
  pres_mod[2] = m + 1.5 * C[2] + 11.0 * C[8] - 0.25 * C[1] - 0.25 * C[0] - 1.0 * C[3];

  // reaction 13;
  pres_mod[3] = m + 2.0 * C[2] - 1.0 * C[8] + 0.10000000000000009 * C[0] + 1.0 * C[32] - 1.0 * C[3];

  // reaction 15;
  thd = m + 1.0 * C[2] + 13.0 * C[8] - 0.21999999999999997 * C[3] - 0.32999999999999996 * C[1] - 0.19999999999999996 * C[0];
  k0 = exp(3.4087162630776540e+01 - 1.72 * logT - (2.6408962763808859e+02 / T));
  kinf = exp(2.2270828345662423e+01 + 0.44 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.00000000e-01 * exp(-T / 1.00000000e-30) + 5.00000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[4] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 28;
  thd = m + 6.5 * C[8] + 0.5 * C[32] + 0.19999999999999996 * C[3] - 0.35 * C[0] + 6.7 * C[9] + 2.7 * C[2];
  k0 = exp(4.9270577684749114e+01 - 2.3 * logT - (2.4531450567319323e+04 / T));
  kinf = exp(2.8324168296488494e+01 + 0.9 * logT - (2.4531450567319323e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.70000000e-01 * exp(-T / 1.00000000e-30) + 4.30000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[5] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 34;
  pres_mod[6] = m;

  // reaction 59;
  pres_mod[7] = m;

  // reaction 86;
  thd = m;
  k0 = exp(7.9974292115367788e+01 - 5.96 * logT - (3.3606512199989469e+04 / T));
  kinf = exp(4.6388174096502127e+01 - 1.31 * logT - (3.2246309716175150e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.90000000e-01 * exp(-T / 1.00000000e-30) + 3.10000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[8] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 99;
  pres_mod[9] = m + 9.0 * C[8];

  // reaction 100;
  pres_mod[10] = m + 9.0 * C[8];

  // reaction 101;
  pres_mod[11] = m;

  // reaction 111;
  pres_mod[12] = m + 9.0 * C[8];

  // reaction 122;
  thd = m + 0.6000000000000001 * C[32];
  k0 = exp(1.9296149481306266e+01 + 0.206 * logT - (-7.7999032553170230e+02 / T));
  kinf = exp(2.8036486224036711e+01 - 0.41 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(1.80000000e-01 * exp(-T / 1.00000000e-30) + 8.20000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[13] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 130;
  thd = m;
  k0 = exp(4.2998340473490288e+01 - 2.87 * logT - (7.7999032553170230e+02 / T));
  kinf = exp(2.7893385380396040e+01 - 0.75 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.50000000e-01 * exp(-T / 1.00000000e+03) + 7.50000000e-01 * exp(-T / 1.00000000e+05) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[14] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 131;
  pres_mod[15] = m;

  // reaction 139;
  thd = m;
  k0 = exp(4.0365366298828434e+01 - 2.5 * logT);
  kinf = exp(2.5423746202738826e+01 - 0.3 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.50000000e-01 * exp(-T / 1.00000000e-30) + 7.50000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[16] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 148;
  thd = m;
  k0 = exp(3.5670178506401783e+01 - (1.5851416293063627e+04 / T));
  kinf = exp(3.3152482033790797e+01 - (1.6253991944950956e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(-1.49000000e-01 * exp(-T / 1.00000000e-30) + 1.14900000e+00 * exp(-T / 3.12500000e+03) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[17] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 152;
  thd = m;
  k0 = exp(3.3152482033790797e+01 - 1.5 * logT);
  kinf = exp(2.1976028805441779e+01 + 0.24 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.90000000e-01 * exp(-T / 1.00000000e-30) + 7.10000000e-01 * exp(-T / 1.70000000e+03) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[18] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 158;
  thd = m;
  k0 = exp(4.4826845844638555e+01 - 3.0 * logT);
  kinf = 30000000000.0;
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.00000000e-01 * exp(-T / 1.00000000e-30) + 4.00000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[19] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 163;
  thd = m + 0.7 * C[32] + 0.3999999999999999 * C[3] + 11.0 * C[8];
  k0 = exp(2.7218531392883420e+01 - (2.8935124979401859e+04 / T));
  kinf = exp(2.5318385687081001e+01 - (2.9166605979237072e+04 / T));
  Pr = k0 * thd / kinf;
  pres_mod[20] =  Pr / (1.0 + Pr);

  // reaction 178;
  thd = m;
  k0 = exp(6.4942386233079020e+01 - 5.49 * logT - (9.9989727537515637e+02 / T));
  kinf = exp(2.7051202620675607e+01 - 0.414 * logT - (3.3212491280704739e+01 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.90000000e-01 * exp(-T / 1.00000000e-30) + 3.10000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[21] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 217;
  pres_mod[22] = m;

  // reaction 218;
  pres_mod[23] = m;

  // reaction 219;
  pres_mod[24] = m;

} // end get_rxn_pres_mod

