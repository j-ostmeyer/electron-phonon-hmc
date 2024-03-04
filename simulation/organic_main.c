#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>

#include "mycpu.h"

#include "organic_flags.h"
#include "organic_aux.h"
#include "organic_lattice.h"
#include "organic_dense_matrix.h"

#include "organic_hmc.h"
#include "organic_trans_loc.h"

#define ROOM_TEMPERATURE 0.025

int main(int argc, char **argv){
	const char *corr_name;

	switch(argc){
		case 3:
			corr_name = NULL;
			break;
		case 4:
			corr_name = argv[3];
			break;
		default:
			printf("Error: One input and one or two output files needed!\n");
			return 0;
	}

#pragma omp parallel
	{
		//printf("Process %d of %d on CPU %d\n", omp_get_thread_num(), omp_get_max_threads(), findmycpu_());
	}

	const unsigned nn = 6, nn2 = nn/2;

	FILE *in = fopen(argv[1], "r");

	unsigned l1, l2, nt;
	if(fscanf(in, "%u %u\n%u\n", &l1, &l2, &nt) != 3){
		printf("Error: Need L1 L2, Nt!\n");
		return 0;
	}

	double beta;
	if(fscanf(in, "%lg\n", &beta) != 1){
		printf("Error: Need beta!\n");
		return 0;
	}

	double t[3];
	double kappa[3];
	if(fscanf(in, "%lg %lg %lg\n%lg %lg %lg\n", t, t+1, t+2, kappa, kappa+1, kappa+2) != nn){
		printf("Error: Need 3 x t, 3 x kappa!\n");
		return 0;
	}

	double mu0, mu, omega;
	if(fscanf(in, "%lg\n%lg\n", &mu0, &omega) != 2){
		printf("Error: Need mu0, omega!\n");
		return 0;
	}
	mu = mu0;

	unsigned steps;
	double traj_length;
	if(fscanf(in, "%u %lg\n", &steps, &traj_length) != 2){
		printf("Error: Need steps, traj-length!\n");
		return 0;
	}

	unsigned therm, meas, meas_freq;
	if(fscanf(in, "%u\n%u %u\n", &therm, &meas, &meas_freq) != 3){
		printf("Error: Need therm, meas, meas-freq!\n");
		return 0;
	}

	double unit_vec_a[2];
	double unit_vec_b[2];
	double *contr_mat = NULL;
	if(fscanf(in, "%lg %lg\n%lg %lg\n", unit_vec_a, unit_vec_a+1, unit_vec_b, unit_vec_b+1) == 4){
		contr_mat = construct_contraction_mat(unit_vec_a, unit_vec_b, nn2);
	}

	unsigned flavors;
	if(fscanf(in, "%u\n", &flavors) != 1){
		flavors = 1;
	}

	organic_flags mode[1];
	set_flags(in, mode);

	fclose(in);

	const double x_avg = mode->static_disorder? 1: sqrt(.5 / tanh(.5 * omega / ROOM_TEMPERATURE));

	if(!mode->static_disorder) beta /= nt;
	else beta *= -1;

	mu *= beta;
	omega *= beta;
	if(mode->classical_fourier_acc) mode->cfa_mass *= beta;

	for(unsigned k = 0; k < nn2; k++){
		t[k] *= beta;
		kappa[k] *= t[k] / x_avg;
		mu -= 2*fabs(t[k]);
	}

	if(traj_length == 0){
		traj_length = HALF_M_PI; // de-correlate all phonon modes exactly with EFA
	}

	//printf("Params:\n%u %u\n%u\n%g\n%g %g %g\n%g %g %g\n%g\n%g\n%u %g\n%u\n%u %u\n",\
	//		l1, l2, nt, beta,\
	//		t[0], t[1], t[2], kappa[0], kappa[1], kappa[2],\
	//		mu, omega, steps, traj_length, therm, meas, meas_freq);

	unsigned greens_dim, res_dim, length;
	if(!mode->static_disorder){
		greens_dim = nt*nn2*nn2;
		res_dim = meas/meas_freq * (greens_dim + NUM_RES);
	}else{
		greens_dim = nt; // here nt is the number of bins in the conductivity
		res_dim = meas + greens_dim; // all samples independent, measure every time
	}

	double *greens = calloc(res_dim, sizeof(double));

	if(!mode->static_disorder)
		length = run_hmc(greens, l1, l2, nt, t, kappa, mu, omega, flavors, contr_mat, steps, traj_length, therm, meas, meas_freq, argv[2], corr_name, mode);
	else
		length = run_trans_loc(greens, l1, l2, nt, t, kappa, mu, contr_mat, meas, argv[2], corr_name);

	if(length != res_dim)
		printf("ERROR: Result vector expected of length %d, is %d!\n", res_dim, length);

	if(contr_mat) free(contr_mat);
	free(greens);

	return 0;
}
