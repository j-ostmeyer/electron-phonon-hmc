#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>

#include "organic_flags.h"
#include "organic_aux.h"

double sigma_of_w(double w, double t, double omega, double width){
	// w: integration variable; omega: external frequency
	const double effW = -(w*width + omega);

	t *= effW;

	return exp(t) / (1 + w*w) * (1 + exp(effW - 2*t)) / (1 + exp(effW));
}

double calc_integral_sigma(double t, double omega, double width){
	if(t > .5) t = 1-t;
	if(t == 0) return M_PI;

	const double h = 0.03;
	const unsigned max_iter = 200;

	double sum = HALF_M_PI * sigma_of_w(0, t, omega, width);

	for(unsigned i = 1; i < max_iter; i++){
		double x = i*h;
		double sh_x = HALF_M_PI * sinh(x);
		double ch_x = HALF_M_PI * cosh(x);
		double ch_sh_x = cosh(sh_x);
		double w = sinh(sh_x);

		double sigma_val = sigma_of_w(w, t, omega, width) + sigma_of_w(w, t, -omega, width);
		double deriv = ch_x * ch_sh_x;

		double int_val = sigma_val * deriv;

		if(int_val >= DBL_EPSILON * sum) sum += int_val;
		else{
			//printf("It took %d iterations.\n", i);
			break;
		}
	}

	return h * sum;
}

double matsubara_integral(double t, double omega, double width, unsigned nt){
	const double complex effW = width + I*omega;
	double sum = creal(ccos(effW*(t-.5)) / ccos(.5*effW));

	for(unsigned i = 0; i < nt; i++){
		const double freq = (2*i+1) * M_PI;
		const double complex diff = I*freq - omega;
		sum += creal(4*width * sin(freq*t) / (width*width + diff*diff));
	}

	return sum;
}

double eval_GK_kernel(double t, double omega, double width, unsigned nt){
	if(nt){
		return matsubara_integral(t, omega, width, nt);
	}else{
		if(width == 0)
			return sigma_of_w(0, t, omega, 1);
		else
			return calc_integral_sigma(t, omega, width) / M_PI;
	}
}

double lorentz(double x, double width){
	return width / (M_PI*(width*width + x*x));
	//return (atan((x+h*.5)/width) - atan((x-h*.5)/width)) / (h*M_PI);
}

double gaussian(double x, double width){
	x /= width;
	return exp(-.5*x*x) / (sqrt(M_2PI) * width);
}

double kernel(double x, double x0, double par){
	return lorentz(x-x0, par) + lorentz(x+x0, par);
	//return gaussian(x-x0, par) + gaussian(x+x0, par);
}

double smeared_sigma_w(double omega, double *sigma, double width, double h, unsigned n){
	if(width == 0){
		const unsigned bin = (unsigned) (omega/h + .5);
		if(bin == 0) return sigma[bin]; // get rid of integral weight 1/2
		else return M_PI * .5*tanhc(.5*omega) * sigma[bin];
	}

	double sum = 0;
	
	for(unsigned i = 0; i < n; i++){
		sum += sigma[i] * kernel(omega, i*h, width);
	}

	return M_PI * .5*tanhc(.5*omega) * h * sum; // factor pi for compatibility
}

double smeared_sigma_t(double t, double *sigma, double width, double h, unsigned n, unsigned nt){
	double sum = 0;
	
	for(unsigned i = 0; i < n; i++){
		sum += sigma[i] * eval_GK_kernel(t, i*h, width, nt);
	}

	return h * sum;
}

int main(int argc, char **argv){
	if(argc != 4 && argc != 5){
		printf("Error! 3 or 4 files needed: greens data, list of smearing widths, frequency output, opt. time output.\n");
		return 0;
	}

	unsigned n, n_smear, nt;
	double h, beta;
	double *sigma, *omega;

	sigma = fread_greens(argv[1], &h, &n);
	omega = fscan_array(argv[2], &n_smear, &beta, &nt);

	//printf("read %d values for sigma and %d for omega\n", n, n_smear);

	for(unsigned k = 0; k < n_smear; k++) omega[k] *= beta;

	FILE *out_w = fopen(argv[3], "w");

	for(unsigned i = 0; i < n; i++){
		const double freq = i * h;

		fprintf(out_w, "%g\t", freq/beta);
		for(unsigned k = 0; k < n_smear; k++)
			fprintf(out_w, "%g\t", smeared_sigma_w(freq, sigma, omega[k], h, n));
		fprintf(out_w, "\n");
	}

	fclose(out_w);

	if(argc > 4){
		FILE *out_t = fopen(argv[4], "w");

		for(unsigned i = 0; i <= nt/2; i++){
			const double t = (double) i / nt;

			fprintf(out_t, "%g\t", t);
			for(unsigned k = 0; k < n_smear; k++)
				fprintf(out_t, "%g\t", smeared_sigma_t(t, sigma, omega[k], h, n, 0));
			fprintf(out_t, "\n");
		}

		fclose(out_t);
	}

	free(sigma);
	free(omega);

	return 0;
}
