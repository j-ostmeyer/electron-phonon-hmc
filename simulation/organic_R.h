#ifndef ORGANIC_R_INTERFACE
#define ORGANIC_R_INTERFACE

static R_INLINE double complex toC99(const Rcomplex *x);
static R_INLINE void SET_C99_COMPLEX(Rcomplex *x, R_xlen_t i, double complex value);

SEXP simulate_hmc(SEXP L1, SEXP L2, SEXP Nt, SEXP hop, SEXP dHop, SEXP mu, SEXP omega, SEXP steps, SEXP traj_length, SEXP therm, SEXP meas, SEXP meas_freq, SEXP obs, SEXP corr);

#endif
