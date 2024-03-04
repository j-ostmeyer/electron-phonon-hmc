#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>

#include "mt19937-64.h"
#include "organic_random.h"

#include "organic_flags.h"
#include "organic_aux.h"
#include "organic_dense_matrix.h"

#include "organic_bosons.h"

void hot_start(double *x, double omega, unsigned ns, unsigned nn, unsigned nt){
	// Similar to random_vector, but constant over time
	// and with std dev governed by phonon number
	const unsigned loc_dim = ns*nn/2;
	const double sd = sqrt(0.5 / tanh(0.5 * nt*omega));

	//printf("sd = %g\n", sd);

	for(unsigned k = 0; k < loc_dim; k += 2){
		const double u1 = sd * sqrt(-2 * log(genrand64_real2()));
		const double u2 = M_2PI * genrand64_real2();
		for(unsigned i = 0; i < nt; i++){
			
			x[k + i*loc_dim] = u1 * cos(u2);
			x[(k+1)%loc_dim + i*loc_dim] = u1 * sin(u2);
			
		}
	}
}

double higher_Matsubaras_norm(double *x, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt)
{
	const unsigned nn2 = nn/2;
	double xd2 = 0.0;
	
	for(unsigned i=0; i<ns; i++)
		for(unsigned k=0; k<nn2; k++)
		{
			double xm = 0.0;
			for(unsigned tau=0; tau<nt; tau++)
				xm += x[tau*ns*nn2 + i*nn2 + k];
			xm /= (double)nt;
			for(unsigned tau=0; tau<nt; tau++)
				xd2 += pow(x[tau*ns*nn2 + i*nn2 + k] - xm, 2);
		};
	
	return xd2;
}

void sample_free_bosons(double *x, double complex *xc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft){
	const unsigned loc_dim = ns*nn/2, compl_dim = nt/2+1, dim = nt*loc_dim;
	const double matsubara = M_PI / nt;

	random_vector(x, dim);
	fftw_execute(fft[0]);

#pragma omp parallel for
	for(unsigned tau = 0; tau < compl_dim; tau++){
		const double freq = sqrt(omega + pow(2*sin(matsubara * tau), 2) / omega);
		const double ampl = 1./freq / nt; // factor nt because fft accumulates it from there-back trafo

		for(unsigned i = 0; i < loc_dim; i++){
			const unsigned pos = tau*loc_dim + i;
			xc[pos] *= ampl;
		}
	}

	fftw_execute(fft[2]);
}

void sample_fourier_momenta(double *p, double complex *pc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, organic_flags *mode){
	const unsigned loc_dim = ns*nn/2, compl_dim = nt/2+1, dim = nt*loc_dim;
	const double matsubara = M_PI / nt;

	random_vector(p, dim);

	if(mode->no_fourier_acc) return;
	if(mode->classical_fourier_acc)
		omega = sqrt(pow(mode->cfa_mass, 2) + pow(omega, 2));

	fftw_execute(fft[1]);

#pragma omp parallel for
	for(unsigned tau = 0; tau < compl_dim; tau++){
		const double freq = sqrt(omega + pow(2*sin(matsubara * tau), 2) / omega);
		const double ampl = freq / nt; // factor nt because fft accumulates it from there-back trafo

		for(unsigned i = 0; i < loc_dim; i++){
			const unsigned pos = tau*loc_dim + i;
			pc[pos] *= ampl;
		}
	}

	fftw_execute(fft[3]);
}

double energy_fourier_momenta(double complex *pc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft){
	const unsigned loc_dim = ns*nn/2, compl_dim = nt/2+1, centerT = (nt+1)/2;
	const double matsubara = M_PI / nt;
	double en = 0;

	fftw_execute(fft[1]);

#pragma omp parallel for reduction(+:en)
	for(unsigned tau = 0; tau < compl_dim; tau++){
		const double freq2 = omega + pow(2*sin(matsubara * tau), 2) / omega;
		double weight = .5 / nt / freq2; // factor nt because fft accumulates it from there-back trafo

		// elements t=1...Nt/2-1 occur twice (as complex conj pairs), but are stored only once
		if(tau > 0 && tau < centerT) weight *= 2;

		const unsigned shift = tau*loc_dim;
		for(unsigned i = 0; i < loc_dim; i++){
			en += norm(pc[shift + i]) * weight;
		}
	}

	return en;
}

void evolve_bosons(double complex *xc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, double h){
	const unsigned loc_dim = ns*nn/2, compl_dim = nt/2+1, dim = compl_dim*loc_dim;
	const double matsubara = M_PI / nt;
	const double c = cos(h) / nt, s = sin(h) / nt; // factor nt because fft accumulates it from there-back trafo
	double complex *pc = xc + dim;

	fftw_execute(fft[0]);
	fftw_execute(fft[1]);

#pragma omp parallel for
	for(unsigned tau = 0; tau < compl_dim; tau++){
		const double freq2 = omega + pow(2*sin(matsubara * tau), 2) / omega;
		const double sw = s * freq2, s1w = s / freq2;

		for(unsigned i = 0; i < loc_dim; i++){
			const unsigned pos = tau*loc_dim + i;
			const double complex x0 = xc[pos], p0 = pc[pos];
			xc[pos] = c * x0 + s1w * p0;
			pc[pos] = c * p0 - sw  * x0;
		}
	}

	fftw_execute(fft[2]);
	fftw_execute(fft[3]);
}

void evolve_bosons_approx_step(double complex *xc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, double h, organic_flags *mode){
	const unsigned loc_dim = ns*nn/2, compl_dim = nt/2+1, dim = compl_dim*loc_dim;
	const double matsubara = M_PI / nt;
	const double norm = 1. / nt; // factor nt because fft accumulates it from there-back trafo
	double complex *pc = xc + dim;

	if(mode->classical_fourier_acc)
		omega = sqrt(pow(mode->cfa_mass, 2) + pow(omega, 2));

	fftw_execute(fft[0]);
	fftw_execute(fft[1]);

#pragma omp parallel for
	for(unsigned tau = 0; tau < compl_dim; tau++){
		const double freq2 = omega + pow(2*sin(matsubara * tau), 2) / omega;
		const double inv_mass = h / freq2;

		for(unsigned i = 0; i < loc_dim; i++){
			const unsigned pos = tau*loc_dim + i;
			xc[pos] += inv_mass * pc[pos];
			xc[pos] *= norm;
		}
	}

	fftw_execute(fft[2]);
}

void boson_force(double *x, double *p_dot, double omega, unsigned ns, unsigned nn, unsigned nt, int overwrite){
	const unsigned loc_dim = ns*nn/2, dim = nt*loc_dim;

	if(overwrite)
		for(unsigned i = 0; i < dim; i++) p_dot[i] = 0;

	for(unsigned tau = 0; tau < nt; tau++){
		const unsigned shift = tau*loc_dim,
			  shiftP = ((tau+1) % nt) * loc_dim, shiftM = ((tau+nt-1) % nt) * loc_dim;
		for(unsigned i = 0; i < loc_dim; i++){
			const unsigned pos = shift + i;
			p_dot[pos] -= omega*x[pos];
			p_dot[pos] += (x[shiftP + i] - 2*x[pos] + x[shiftM + i]) / omega;
		}
	}
}
