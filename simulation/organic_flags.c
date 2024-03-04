#include <stdio.h>
#include <string.h>

#include "organic_flags.h"

void set_flags(FILE *in, organic_flags *mode){
	// set default values (all off)
	mode->schur_form_inv = 0;
	mode->zero_filling = 0;
	mode->static_disorder = 0;
	mode->single_particle = 0;
	mode->constant_force = 0;
	mode->no_fourier_acc = 0;
	mode->classical_fourier_acc = 0;

	char flag[200];

	while(fscanf(in, "%s\n", flag) != EOF){
		//printf("%s\n", flag);

		if(!strcmp(flag, "SCHUR_FORM_INV")){
			mode->schur_form_inv = 1;
		}else if(!strcmp(flag, "ZERO_FILLING")){
			mode->zero_filling = 1;
		}else if(!strcmp(flag, "STATIC_DISORDER")){
			mode->static_disorder = 1;
		}else if(!strcmp(flag, "SINGLE_PARTICLE")){
			mode->single_particle = 1;
		}else if(!strcmp(flag, "CONSTANT_FORCE")){
			mode->constant_force = 1;
		}else if(!strcmp(flag, "NO_FOURIER_ACC")){
			mode->no_fourier_acc = 1;
		}else if(!strcmp(flag, "CLASSICAL_FOURIER_ACC")){
			mode->classical_fourier_acc = 1;
			if(fscanf(in, "%u %lg\n", &mode->cfa_steps, &mode->cfa_mass) != 2){
				mode->cfa_steps = 1;
				mode->cfa_mass = 0;
			}
		}
	}

	//printf("sfi = %d\n", mode->schur_form_inv);
	//printf("zf = %d\n", mode->zero_filling);
	//printf("sd = %d\n", mode->static_disorder);
	//printf("sp = %d\n", mode->single_particle);
	//printf("cf = %d\n", mode->constant_force);
	//printf("nfa = %d\n", mode->no_fourier_acc);
	//printf("cfa = %d\n", mode->classical_fourier_acc);
}
