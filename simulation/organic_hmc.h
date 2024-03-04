#ifndef ORGANIC_FFTW3
#define ORGANIC_FFTW3
#include <fftw3.h>
#endif

#ifndef ORGANIC_HMC
#define ORGANIC_HMC

double check_force(double *x, double *p, double *p_dot, double *M_hat, double *w, int *perm, double *t, double *kappa, double mu, double omega, unsigned flavors, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, organic_flags *mode);

double energy_momenta(double *p, double complex *pc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, organic_flags *mode);
double hamilton(double *x, double *p, double logDetM, double energyP, double omega, unsigned flavors, unsigned ns, unsigned nn, unsigned nt);

void force(double *x, double *p_dot, double *M_hat, double *QM1P, double *w, int *perm, double *t, double *kappa, double mu, double omega, unsigned flavors, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, double logDetM, organic_flags *mode);
void update_x(double *x, double *p, double *p_dot, double complex *xc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, double h, double half, organic_flags *mode);

void leap_frog(double *x, double *p, double *p_dot, double complex *xc, double *M_hat, double *QM1P, double *w, int *perm, double *t, double *kappa, double mu, double omega, unsigned flavors, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, unsigned steps, double traj_length, const fftw_plan *fft, organic_flags *mode);

short trajectory(double *x, double complex *xc, int *perm, double *t, double *kappa, double mu, double omega, unsigned flavors, double *contr_mat, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, unsigned steps, double traj_length, double *results, double *greens, int iteration, unsigned meas_freq, const fftw_plan *fft, double *old, FILE *res_out, FILE *corr_out, organic_flags *mode);

unsigned run_hmc(double *greens, unsigned l1, unsigned l2, unsigned nt, double *t, double *kappa, double mu, double omega, unsigned flavors, double *contr_mat, unsigned steps, double traj_length, unsigned therm, unsigned meas, unsigned meas_freq, const char *res_name, const char *corr_name, organic_flags *mode);

#endif
