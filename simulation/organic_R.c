#include <R.h>
#include <Rinternals.h>
#include <complex.h>

#include "organic_aux.h"
#include "organic_hmc.h"

static R_INLINE double complex toC99(const Rcomplex *x){
#if __GNUC__
	double complex ans = (double complex) 0; /* -Wall */
	__real__ ans = x->r;
	__imag__ ans = x->i;
	return ans;
#else
	return x->r + x->i * I;
#endif
}

static R_INLINE void SET_C99_COMPLEX(Rcomplex *x, R_xlen_t i, double complex value){
	Rcomplex *ans = x+i;
	ans->r = creal(value);
	ans->i = cimag(value);
}

SEXP simulate_hmc(SEXP L1, SEXP L2, SEXP Nt, SEXP hop, SEXP dHop, SEXP mu, SEXP omega, SEXP steps, SEXP traj_length, SEXP therm, SEXP meas, SEXP meas_freq, SEXP obs, SEXP corr){
	const unsigned l1 = asInteger(L1), l2 = asInteger(L2);
	const unsigned nn = 6, nn2 = nn/2;
	const unsigned res_dim = asInteger(meas)/asInteger(meas_freq) * (nn2*nn2*asInteger(Nt) + NUM_RES);

	SEXP out = PROTECT(allocVector(REALSXP, res_dim)); // need space for results
	double *res = calloc(res_dim, sizeof(double));

	unsigned length = run_hmc(res, l1, l2, asInteger(Nt), REAL(hop), REAL(dHop), asReal(mu), asReal(omega), NULL, asInteger(steps), asReal(traj_length), asInteger(therm), asInteger(meas), asInteger(meas_freq), CHAR(asChar(obs)), CHAR(asChar(corr)));

	if(length != res_dim)
		printf("ERROR: Result vector expected of length %d, is %d!\n", res_dim, length);

	memcpy(REAL(out), res, length*sizeof(double));

	free(res);
	UNPROTECT(1);

	return out;
}
