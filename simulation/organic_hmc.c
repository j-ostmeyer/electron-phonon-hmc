#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <time.h>
#include <omp.h>

#include "mt19937-64.h"
#include "organic_random.h"

#include "organic_flags.h"
#include "organic_aux.h"
#include "organic_lattice.h"
#include "organic_dense_matrix.h"
#include "organic_bosons.h"

#include "organic_hmc.h"

double check_force(double *x, double *p, double *p_dot, double *M_hat, double *w, int *perm, double *t, double *kappa, double mu, double omega, unsigned flavors, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, organic_flags *mode){
	const double delta = pow(DBL_EPSILON, 1./3);
	const unsigned dim = nt*ns*nn/2;
	double det, e_plus, e_minus, force;
	double error = 0, norm = 0;
	//double norm_p = 0, norm_f = 0;

	for(unsigned i = 0; i < dim; i++){
		const double save_x = x[i];

		x[i] = save_x + delta;
		construct_M_hat(x, M_hat, w, t, kappa, mu, nnt, ns, nn, nt, mode);
		det = logDetM_and_inv_M_hat(M_hat, w, perm, ns, flavors, mode);
		e_plus = -det*flavors; //hamilton(x, p, det, omega, flavors, ns, nn, nt);

		x[i] = save_x - delta;
		construct_M_hat(x, M_hat, w, t, kappa, mu, nnt, ns, nn, nt, mode);
		det = logDetM_and_inv_M_hat(M_hat, w, perm, ns, flavors, mode);
		e_minus = -det*flavors; //hamilton(x, p, det, omega, flavors, ns, nn, nt);

		force = -(e_plus - e_minus) / (2*delta);
		
		error += pow(p_dot[i] - force, 2);
		norm += pow(fabs(p_dot[i]) + fabs(force), 2);

		//norm_p += p_dot[i]*p_dot[i];
		//norm_f += force*force;

		x[i] = save_x;
	}

	//printf("|p_dot| = %g,\t|force| = %g\n", sqrt(norm_p), sqrt(norm_f));
	//printf("abs err = %g,\tnorm = %g\n", sqrt(error), sqrt(norm));
	return sqrt(error/norm); // should not be much larger than 1e-10
}

double energy_momenta(double *p, double complex *pc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, organic_flags *mode){
	if(mode->no_fourier_acc){
		const unsigned loc_dim = ns*nn/2, dim = nt*loc_dim;
		return .5 * dim * average_sq(p, dim);
	}else{
		if(mode->classical_fourier_acc)
			omega = sqrt(pow(mode->cfa_mass, 2) + pow(omega, 2));
		return energy_fourier_momenta(pc, omega, ns, nn, nt, fft);
	}
}

double hamilton(double *x, double *p, double logDetM, double energyP, double omega, unsigned flavors, unsigned ns, unsigned nn, unsigned nt){
	const unsigned nn2 = nn/2;
	double energyPot = 0, energyKin = 0;

#pragma omp parallel for reduction(+:energyPot,energyKin)
	for(unsigned tau = 0; tau < nt; tau++){
		for(unsigned i = 0; i < ns; i++){
			const unsigned shift = (tau*ns + i) * nn2;
			const unsigned shiftT = ((tau+1)%nt * ns + i) * nn2;
			for(unsigned k = 0; k < nn2; k++){
				energyPot += pow(x[shift + k], 2);
				energyKin += pow(x[shiftT + k] - x[shift + k], 2);
			}
		}
	}

	return -logDetM*flavors + energyP + (energyKin/omega + energyPot*omega)/2;
}

void force(double *x, double *p_dot, double *M_hat, double *QM1P, double *w, int *perm, double *t, double *kappa, double mu, double omega, unsigned flavors, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, double logDetM, organic_flags *mode){
	//printf("calculating force with %u flavors\n", flavors);
	const unsigned nn2 = nn/2, mat_dim = ns*ns;

	//double force_ferm = 0;

	if(mode->constant_force){
		// Construct representative derivative dM/dx for 0th time slice
		construct_QM1P(x, M_hat, QM1P, w, t, kappa, mu, nnt, ns, nn, nt, 0, mode);
	}else{
		// Construct the derivatives dM/dx for every time slice
		construct_all_QM1P(x, M_hat, QM1P, w, perm, t, kappa, mu, nnt, ns, nn, nt, mode);
	}

	for(unsigned tau = 0; tau < nt; tau++){
		if(!mode->constant_force){
			refine_QM1P(x, QM1P, w, t, kappa, mu, nnt, ns, nn, tau);
			//refine_QM1P_approx(x, QM1P, w, t, kappa, mu, nnt, ns, nn, tau); // only useful for high Taylor orders >= 6
		}

//#pragma omp parallel for
		for(unsigned i = 0; i < ns; i++){
			const unsigned shift = (tau*ns + i)*nn2;
			const unsigned shiftN = i*nn, shiftS = i*ns;

			for(unsigned k = 0; k < nn2; k++){
				const unsigned pos = shift + k;
				const unsigned neighbour = nnt[shiftN + 2*k];
				double projection = QM1P[shiftS + neighbour] + QM1P[neighbour*ns + i];

				if(mode->schur_form_inv && fabs(projection) > 2) // crude regularisation of the force to avoid inf
					projection = projection > 0? 2: -2;

				p_dot[pos] = -kappa[k] * projection * flavors;

				if(mode->single_particle) p_dot[pos] *= exp(-logDetM);

				//force_ferm += pow(p_dot[pos], 2);
				//force_ferm = force_ferm > fabs(projection)? force_ferm: fabs(projection);
			}
		}

		if(!mode->constant_force) QM1P += mat_dim;
	}

	if(mode->no_fourier_acc) boson_force(x, p_dot, omega, ns, nn, nt, 0);

	//printf("|F_ferm| = %g\n", force_ferm);
}

void update_x(double *x, double *p, double *p_dot, double complex *xc, double omega, unsigned ns, unsigned nn, unsigned nt, const fftw_plan *fft, double h, double half, organic_flags *mode){
	const unsigned loc_dim = ns*nn/2, dim = nt*loc_dim;

	if(mode->no_fourier_acc){
		h *= fabs(half);
		for(unsigned i = 0; i < dim; i++) x[i] += h*p[i];
	}else if(mode->classical_fourier_acc){
		h /= mode->cfa_steps;
		if(half != .5){
			boson_force(x, p_dot, omega, ns, nn, nt, 1);
			for(unsigned i = 0; i < dim; i++) p[i] += .5*h * p_dot[i];
			for(unsigned k = 0; k < mode->cfa_steps/2; k++){
				evolve_bosons_approx_step(xc, omega, ns, nn, nt, fft, h, mode);
				boson_force(x, p_dot, omega, ns, nn, nt, 1);
				for(unsigned i = 0; i < dim; i++) p[i] += h * p_dot[i];
			}
			evolve_bosons_approx_step(xc, omega, ns, nn, nt, fft, .5*h, mode);
		}
		if(half != -.5){
			evolve_bosons_approx_step(xc, omega, ns, nn, nt, fft, .5*h, mode);
			for(unsigned k = 0; k < mode->cfa_steps/2; k++){
				boson_force(x, p_dot, omega, ns, nn, nt, 1);
				for(unsigned i = 0; i < dim; i++) p[i] += h * p_dot[i];
				evolve_bosons_approx_step(xc, omega, ns, nn, nt, fft, h, mode);
			}
			boson_force(x, p_dot, omega, ns, nn, nt, 1);
			for(unsigned i = 0; i < dim; i++) p[i] += .5*h * p_dot[i];
		}
	}else evolve_bosons(xc, omega, ns, nn, nt, fft, fabs(half)*h);
}

void leap_frog(double *x, double *p, double *p_dot, double complex *xc, double *M_hat, double *QM1P, double *w, int *perm, double *t, double *kappa, double mu, double omega, unsigned flavors, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, unsigned steps, double traj_length, const fftw_plan *fft, organic_flags *mode){
	double h = steps? traj_length/steps: 2*traj_length;
	double half = .5;
	const unsigned loc_dim = ns*nn/2, dim = nt*loc_dim;

	update_x(x, p, p_dot, xc, omega, ns, nn, nt, fft, h, half, mode);
	half = 1;

	for(unsigned s = 1; s <= steps; s++){
		//print_config(x, ns, nn, nt);
		construct_M_hat(x, M_hat, w, t, kappa, mu, nnt, ns, nn, nt, mode);
		double logDetM = logDetM_and_inv_M_hat(M_hat, w, perm, ns, flavors, mode);

		force(x, p_dot, M_hat, QM1P, w, perm, t, kappa, mu, omega, flavors, nnt, ns, nn, nt, logDetM, mode);
		//double error = check_force(x, p, p_dot, M_hat, w, perm, t, kappa, mu, omega, flavors, nnt, ns, nn, nt, mode);
		//printf("force err = %g\n", error); fflush(stdout);

#pragma omp parallel for
		for(unsigned i = 0; i < dim; i++) p[i] += h*p_dot[i];

		if(s == steps) half = -.5; // only half step in the end
		update_x(x, p, p_dot, xc, omega, ns, nn, nt, fft, h, half, mode);
	}
}

short trajectory(double *x, double complex *xc, int *perm, double *t, double *kappa, double mu, double omega, unsigned flavors, double *contr_mat, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, unsigned steps, double traj_length, double *results, double *greens, int iteration, unsigned meas_freq, const fftw_plan *fft, double *old, FILE *res_out, FILE *corr_out, organic_flags *mode){
	const unsigned nn2 = nn/2, dim = nt*ns*nn2, mat_dim = ns*ns;
	// total dimension of x >= 4*dim + 8*mat_dim  + mat_dim*nt
	double *p = x + dim; // dim
	double *p_dot = p + dim; // dim
	double *x_old = p_dot + dim; // dim
	double *M_hat = x_old + dim; // mat_dim
	double *M_hat_old = M_hat + mat_dim; // mat_dim
	double *QM1P = M_hat_old + mat_dim; // mat_dim*nt
	double *R_0t = QM1P + mat_dim*nt; // mat_dim
	double *w  = R_0t + mat_dim; // 5*mat_dim

	const unsigned loc_dim = ns*nn2, compl_dim = loc_dim * (nt/2 + 1);
	double complex *pc = xc + compl_dim;

	double logDetM = 0, energy, energyP;
	double boltzmann = 1;
	short acc = 1;

	double logDetM_old = old[0], energy_old;

	if(!mode->zero_filling){
#pragma omp parallel for
		for(unsigned i = 0; i < dim; i++) x_old[i] = x[i];
#pragma omp parallel for
		for(unsigned i = 0; i < mat_dim; i++) M_hat_old[i] = M_hat[i];

		sample_fourier_momenta(p, pc, omega, ns, nn, nt, fft, mode);
		energyP = energy_momenta(p, pc, omega, ns, nn, nt, fft, mode);
		energy_old = hamilton(x, p, logDetM_old, energyP, omega, flavors, ns, nn, nt);

		leap_frog(x, p, p_dot, xc, M_hat, QM1P, w, perm, t, kappa, mu, omega, flavors, nnt, ns, nn, nt, steps, traj_length, fft, mode);

		construct_M_hat(x, M_hat, w, t, kappa, mu, nnt, ns, nn, nt, mode);
		logDetM = logDetM_and_inv_M_hat(M_hat, w, perm, ns, flavors, mode);
		energyP = energy_momenta(p, pc, omega, ns, nn, nt, fft, mode);
		energy = hamilton(x, p, logDetM, energyP, omega, flavors, ns, nn, nt);

		//printf("logDetM: %g -> %g\n", logDetM_old, logDetM);
		//printf("energies: %g -> %g ,  dE = %g\n", energy_old, energy, energy-energy_old);

		boltzmann = exp(energy_old-energy);
		if(!(boltzmann > genrand64_real2())){ // cumbersome expression to avoid problems with NaN
			// reject (nothing to be done if accepted)
#pragma omp parallel for
			for(unsigned i = 0; i < dim; i++) x[i] = x_old[i];
#pragma omp parallel for
			for(unsigned i = 0; i < mat_dim; i++) M_hat[i] = M_hat_old[i];
			logDetM = logDetM_old;
			energy = energy_old;
			acc = 0;
		}

		old[0] = logDetM;
	}else{
		//This is the zero filling part
		sample_free_bosons(x, xc, omega, ns, nn, nt, fft);

		construct_M_hat(x, M_hat, w, t, kappa, mu, nnt, ns, nn, nt, mode);
		logDetM = 0;
		energy = hamilton(x, NULL, 0, 0, omega, flavors, ns, nn, nt);
		//End of zero filling part
	}

	if(iteration % meas_freq == 0){
		//print_config(x, ns, nn, nt);

		if(greens){
			cc_correlator(x, M_hat, QM1P, R_0t, w, t, kappa, mu, nnt, ns, nn, nt, greens, logDetM, mode);
			//phonon_correlator(x, ns, nn, nt, greens);
			fprint_corr(corr_out, greens, contr_mat, nn, nt);
		}

//Measure the deviation of x from its time average


		if(results){
			// collection of different results, current dimension = NUM_RES = 8+3
			results[0] = average(x, dim); // effective rescaling of hopping
			results[1] = average_sq(x, dim) - 0.5; // phonon number expectation
			results[2] = energy; // HMC energy
			results[3] = -logDetM; // fermionic energy
			results[4] = higher_Matsubaras_norm(x, nnt, ns, nn, nt);//omega * dim * (results[1]+.5); // phonon energy = E_kin + E_pot, where E_kin = E_pot through the Virial theorem
			if(!mode->zero_filling && !mode->single_particle)
				results[5] = 1 - trace(M_hat, ns) / ns; // average electron occupation
			else
				results[5] = trace(M_hat, ns) / ns; // average electron occupation to 1st order
			results[6] = acc; // acceptance
			results[7] = boltzmann; // exp(-dH), should average to 1
			average_links(x, kappa, ns, nn/2, nt, results+8); // <t_a \lambda_a x_a>, i.e. additive rescaling of the hopping

			if(greens) results[11] = contract_greens(greens, contr_mat, nn, nt/2) * nt*nt/M_PI; // low-frequency conductivity in units of e^2/hbar

			//fprint_config("configs.txt", x, ns, nn, nt, iteration? "a":"w");
			//fwrite_config("configs.bin", x, dim, iteration? "a":"w");
			fprint_results(res_out, results, 1, NUM_RES);
		}

		//print_config(p, ns, nn, nt);
	}

	return acc;
}

unsigned run_hmc(double *greens, unsigned l1, unsigned l2, unsigned nt, double *t, double *kappa, double mu, double omega, unsigned flavors, double *contr_mat, unsigned steps, double traj_length, unsigned therm, unsigned meas, unsigned meas_freq, const char *res_name, const char *corr_name, organic_flags *mode){

	//printf("Hi from run_hmc with %u flavors!!!\n", flavors); fflush(stdout);
	unsigned ns, nn;

	unsigned *nnt = construct_triangular(l1, l2, &ns, &nn);

	const unsigned nn2 = nn/2, dim = nt*ns*nn2, mat_dim = ns*ns, greens_dim = nt*nn2*nn2;
	const unsigned loc_dim = ns*nn2, compl_dim = loc_dim * (nt/2 + 1);

	double *results = greens + greens_dim*meas/meas_freq;
	double *x = malloc((4*dim + (8+nt) * mat_dim) * sizeof(double));
	fftw_complex *xc = (fftw_complex*) fftw_malloc(2*compl_dim * sizeof(fftw_complex));
	int *perm = malloc(max(ns, 2*nt) * sizeof(int));

	FILE *res_out = res_name? fopen(res_name, "a"): NULL;
	FILE *corr_out = corr_name? fopen(corr_name, "a"): NULL;

	init_genrand64(time(NULL));
	set_taylor_pair_order();

	int fft_dim[1] = {nt};
	fftw_plan fft[4]; // x/p forward/backward in every combination
	fft[0] = fftw_plan_many_dft_r2c(1, fft_dim, loc_dim, x, NULL, loc_dim, 1, xc, NULL, loc_dim, 1, FFTW_MEASURE);
	fft[1] = fftw_plan_many_dft_r2c(1, fft_dim, loc_dim, x+dim, NULL, loc_dim, 1, xc+compl_dim, NULL, loc_dim, 1, FFTW_MEASURE);
	fft[2] = fftw_plan_many_dft_c2r(1, fft_dim, loc_dim, xc, NULL, loc_dim, 1, x, NULL, loc_dim, 1, FFTW_MEASURE);
	fft[3] = fftw_plan_many_dft_c2r(1, fft_dim, loc_dim, xc+compl_dim, NULL, loc_dim, 1, x+dim, NULL, loc_dim, 1, FFTW_MEASURE);

	//hot_start(x, omega, ns, nn, nt);
	sample_free_bosons(x, xc, omega, ns, nn, nt, fft);
	
	//printf("This is after hot_start!!!\n"); fflush(stdout);

	double *M_hat = x + 4*dim;
	double *w = x + 4*dim + (nt+3)*mat_dim;

	construct_M_hat(x, M_hat, w, t, kappa, mu, nnt, ns, nn, nt, mode);
	double logDetM = 0;
	if(!mode->zero_filling)
		logDetM = logDetM_and_inv_M_hat(M_hat, w, perm, ns, flavors, mode);

	for(unsigned i = 0; i < therm; i++)
		trajectory(x, xc, perm, t, kappa, mu, omega, flavors, contr_mat, nnt, ns, nn, nt, steps, traj_length, NULL, NULL, -1, meas_freq, fft, &logDetM, NULL, NULL, mode);

	for(unsigned i = 0; i < meas; i++)
		trajectory(x, xc, perm, t, kappa, mu, omega, flavors, contr_mat, nnt, ns, nn, nt, steps, traj_length, results + i/meas_freq*NUM_RES, greens + i/meas_freq*greens_dim, i, meas_freq, fft, &logDetM, res_out, corr_out, mode);

	free(nnt);
	free(x);
	fftw_free(xc);
	free(perm);
	if(res_out) fclose(res_out);
	if(corr_out) fclose(corr_out);

	for(unsigned i = 0; i < 4; i++)	fftw_destroy_plan(fft[i]);

	return meas/meas_freq * (greens_dim + NUM_RES); // dimension of results array for sanity check
}
