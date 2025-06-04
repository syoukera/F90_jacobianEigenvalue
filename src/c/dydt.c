#include "header.h"
#include "chem_utils.h"
#include "rates.h"

#if defined(CONP)

void dydt (const double t, const double pres, const double * __restrict__ y, double * __restrict__ dy) {

  // species molar concentrations
  double conc[26];
  double y_N;
  double mw_avg;
  double rho;
  eval_conc (y[0], pres, &y[1], &y_N, &mw_avg, &rho, conc);

  // local arrays holding reaction rates
  double fwd_rates[171];
  double rev_rates[168];
  eval_rxn_rates (y[0], pres, conc, fwd_rates, rev_rates);

  // get pressure modifications to reaction rates
  double pres_mod[19];
  get_rxn_pres_mod (y[0], pres, conc, pres_mod);

  // evaluate species molar net production rates
  double dy_N;
  eval_spec_rates (fwd_rates, rev_rates, pres_mod, &dy[1], &dy_N);
  // local array holding constant pressure specific heat
  double cp[26];
  eval_cp (y[0], cp);

  // constant pressure mass-average specific heat
  double cp_avg = (cp[0] * y[1]) + (cp[1] * y[2]) + (cp[2] * y[3]) + (cp[3] * y[4])
              + (cp[4] * y[5]) + (cp[5] * y[6]) + (cp[6] * y[7]) + (cp[7] * y[8])
              + (cp[8] * y[9]) + (cp[9] * y[10]) + (cp[10] * y[11]) + (cp[11] * y[12])
              + (cp[12] * y[13]) + (cp[13] * y[14]) + (cp[14] * y[15]) + (cp[15] * y[16])
              + (cp[16] * y[17]) + (cp[17] * y[18]) + (cp[18] * y[19]) + (cp[19] * y[20])
              + (cp[20] * y[21]) + (cp[21] * y[22]) + (cp[22] * y[23]) + (cp[23] * y[24])
              + (cp[24] * y[25]) + (cp[25] * y_N);

  // local array for species enthalpies
  double h[26];
  eval_h(y[0], h);
  // rate of change of temperature
  dy[0] = (-1.0 / (rho * cp_avg)) * ((dy[1] * h[0] * 2.0158800000000001e+00)
        + (dy[2] * h[1] * 3.1998799999999999e+01) + (dy[3] * h[2] * 1.0079400000000001e+00)
        + (dy[4] * h[3] * 1.5999400000000000e+01) + (dy[5] * h[4] * 1.7007339999999999e+01)
        + (dy[6] * h[5] * 3.3006740000000001e+01) + (dy[7] * h[6] * 1.8015280000000001e+01)
        + (dy[8] * h[7] * 3.4014679999999998e+01) + (dy[9] * h[8] * 1.4006740000000001e+01)
        + (dy[10] * h[9] * 1.7030560000000001e+01) + (dy[11] * h[10] * 1.6022620000000000e+01)
        + (dy[12] * h[11] * 1.5014680000000000e+01) + (dy[13] * h[12] * 2.9021420000000003e+01)
        + (dy[14] * h[13] * 3.0006140000000002e+01) + (dy[15] * h[14] * 4.4012880000000003e+01)
        + (dy[16] * h[15] * 3.1014080000000000e+01) + (dy[17] * h[16] * 3.1014080000000000e+01)
        + (dy[18] * h[17] * 3.2022019999999998e+01) + (dy[19] * h[18] * 4.6005539999999996e+01)
        + (dy[20] * h[19] * 4.7013480000000001e+01) + (dy[21] * h[20] * 6.3012879999999996e+01)
        + (dy[22] * h[21] * 3.0029360000000000e+01) + (dy[23] * h[22] * 3.0029360000000000e+01)
        + (dy[24] * h[23] * 3.2045240000000000e+01) + (dy[25] * h[24] * 3.1037300000000002e+01)
        + (dy_N * h[25] * 2.8013480000000001e+01));

  // calculate rate of change of species mass fractions
  dy[1] *= (2.0158800000000001e+00 / rho);
  dy[2] *= (3.1998799999999999e+01 / rho);
  dy[3] *= (1.0079400000000001e+00 / rho);
  dy[4] *= (1.5999400000000000e+01 / rho);
  dy[5] *= (1.7007339999999999e+01 / rho);
  dy[6] *= (3.3006740000000001e+01 / rho);
  dy[7] *= (1.8015280000000001e+01 / rho);
  dy[8] *= (3.4014679999999998e+01 / rho);
  dy[9] *= (1.4006740000000001e+01 / rho);
  dy[10] *= (1.7030560000000001e+01 / rho);
  dy[11] *= (1.6022620000000000e+01 / rho);
  dy[12] *= (1.5014680000000000e+01 / rho);
  dy[13] *= (2.9021420000000003e+01 / rho);
  dy[14] *= (3.0006140000000002e+01 / rho);
  dy[15] *= (4.4012880000000003e+01 / rho);
  dy[16] *= (3.1014080000000000e+01 / rho);
  dy[17] *= (3.1014080000000000e+01 / rho);
  dy[18] *= (3.2022019999999998e+01 / rho);
  dy[19] *= (4.6005539999999996e+01 / rho);
  dy[20] *= (4.7013480000000001e+01 / rho);
  dy[21] *= (6.3012879999999996e+01 / rho);
  dy[22] *= (3.0029360000000000e+01 / rho);
  dy[23] *= (3.0029360000000000e+01 / rho);
  dy[24] *= (3.2045240000000000e+01 / rho);
  dy[25] *= (3.1037300000000002e+01 / rho);

} // end dydt

#elif defined(CONV)

void dydt (const double t, const double rho, const double * __restrict__ y, double * __restrict__ dy) {

  // species molar concentrations
  double conc[26];
  double y_N;
  double mw_avg;
  double pres;
  eval_conc_rho (y[0]rho, &y[1], &y_N, &mw_avg, &pres, conc);

  // local arrays holding reaction rates
  double fwd_rates[171];
  double rev_rates[168];
  eval_rxn_rates (y[0], pres, conc, fwd_rates, rev_rates);

  // get pressure modifications to reaction rates
  double pres_mod[19];
  get_rxn_pres_mod (y[0], pres, conc, pres_mod);

  // evaluate species molar net production rates
  double dy_N;  eval_spec_rates (fwd_rates, rev_rates, pres_mod, &dy[1], &dy_N);

  double cv[26];
  eval_cv(y[0], cv);

  // constant volume mass-average specific heat
  double cv_avg = (cv[0] * y[1]) + (cv[1] * y[2]) + (cv[2] * y[3]) + (cv[3] * y[4])
              + (cv[4] * y[5]) + (cv[5] * y[6]) + (cv[6] * y[7]) + (cv[7] * y[8])
              + (cv[8] * y[9]) + (cv[9] * y[10]) + (cv[10] * y[11]) + (cv[11] * y[12])
              + (cv[12] * y[13]) + (cv[13] * y[14]) + (cv[14] * y[15]) + (cv[15] * y[16])
              + (cv[16] * y[17]) + (cv[17] * y[18]) + (cv[18] * y[19]) + (cv[19] * y[20])
              + (cv[20] * y[21]) + (cv[21] * y[22]) + (cv[22] * y[23]) + (cv[23] * y[24])
              + (cv[24] * y[25])(cv[25] * y_N);

  // local array for species internal energies
  double u[26];
  eval_u (y[0], u);

  // rate of change of temperature
  dy[0] = (-1.0 / (rho * cv_avg)) * ((dy[1] * u[0] * 2.0158800000000001e+00)
        + (dy[2] * u[1] * 3.1998799999999999e+01) + (dy[3] * u[2] * 1.0079400000000001e+00)
        + (dy[4] * u[3] * 1.5999400000000000e+01) + (dy[5] * u[4] * 1.7007339999999999e+01)
        + (dy[6] * u[5] * 3.3006740000000001e+01) + (dy[7] * u[6] * 1.8015280000000001e+01)
        + (dy[8] * u[7] * 3.4014679999999998e+01) + (dy[9] * u[8] * 1.4006740000000001e+01)
        + (dy[10] * u[9] * 1.7030560000000001e+01) + (dy[11] * u[10] * 1.6022620000000000e+01)
        + (dy[12] * u[11] * 1.5014680000000000e+01) + (dy[13] * u[12] * 2.9021420000000003e+01)
        + (dy[14] * u[13] * 3.0006140000000002e+01) + (dy[15] * u[14] * 4.4012880000000003e+01)
        + (dy[16] * u[15] * 3.1014080000000000e+01) + (dy[17] * u[16] * 3.1014080000000000e+01)
        + (dy[18] * u[17] * 3.2022019999999998e+01) + (dy[19] * u[18] * 4.6005539999999996e+01)
        + (dy[20] * u[19] * 4.7013480000000001e+01) + (dy[21] * u[20] * 6.3012879999999996e+01)
        + (dy[22] * u[21] * 3.0029360000000000e+01) + (dy[23] * u[22] * 3.0029360000000000e+01)
        + (dy[24] * u[23] * 3.2045240000000000e+01) + (dy[25] * u[24] * 3.1037300000000002e+01)
        + (dy_N * u[25] * 2.8013480000000001e+01));

  // calculate rate of change of species mass fractions
  dy[1] *= (2.0158800000000001e+00 / rho);
  dy[2] *= (3.1998799999999999e+01 / rho);
  dy[3] *= (1.0079400000000001e+00 / rho);
  dy[4] *= (1.5999400000000000e+01 / rho);
  dy[5] *= (1.7007339999999999e+01 / rho);
  dy[6] *= (3.3006740000000001e+01 / rho);
  dy[7] *= (1.8015280000000001e+01 / rho);
  dy[8] *= (3.4014679999999998e+01 / rho);
  dy[9] *= (1.4006740000000001e+01 / rho);
  dy[10] *= (1.7030560000000001e+01 / rho);
  dy[11] *= (1.6022620000000000e+01 / rho);
  dy[12] *= (1.5014680000000000e+01 / rho);
  dy[13] *= (2.9021420000000003e+01 / rho);
  dy[14] *= (3.0006140000000002e+01 / rho);
  dy[15] *= (4.4012880000000003e+01 / rho);
  dy[16] *= (3.1014080000000000e+01 / rho);
  dy[17] *= (3.1014080000000000e+01 / rho);
  dy[18] *= (3.2022019999999998e+01 / rho);
  dy[19] *= (4.6005539999999996e+01 / rho);
  dy[20] *= (4.7013480000000001e+01 / rho);
  dy[21] *= (6.3012879999999996e+01 / rho);
  dy[22] *= (3.0029360000000000e+01 / rho);
  dy[23] *= (3.0029360000000000e+01 / rho);
  dy[24] *= (3.2045240000000000e+01 / rho);
  dy[25] *= (3.1037300000000002e+01 / rho);

} // end dydt

#endif
