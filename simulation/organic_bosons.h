#ifndef ORGANIC_FFTW3
#define ORGANIC_FFTW3
#include <fftw3.h>
#endif

#ifndef ORGANIC_BOSONS
#define ORGANIC_BOSONS

void hot_start(double *x, double omega, unsigned ns, unsigned nn, unsigned nt);

double higher_Matsubaras_norm(double *x, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt);

void sample_free_bosons(double *x, double complex *xc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft);
void sample_fourier_momenta(double *p, double complex *pc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, organic_flags *mode);

double energy_fourier_momenta(double complex *pc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft);

void evolve_bosons(double complex *xc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, double h);
void evolve_bosons_approx_step(double complex *xc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, double h, organic_flags *mode);
void boson_force(double *x, double *p_dot, double omega, unsigned ns, unsigned nn, unsigned nt, int overwrite);

#endif
