CC = gcc
MYFILES=lattice flags aux random dense_matrix bosons hmc trans_loc
RNG=mt19937-64 mycpu
PARAM_FILES=taylor-coeff_pairs.h

CC += -std=gnu11 -march=native
MYFLAGS=-Wall -pedantic -O3
MYFLAGS += -Wno-unknown-pragmas
CC += -Wno-unused-result
MYFLAGS += -Wno-comment

MYLIBS=-lfftw3 -lm

MYFLAGS += -DTAYLOR_ORDER_4
#MYFLAGS += -DTAYLOR_JACKSON_6_4
#MYFLAGS += -DCHEBYSHEV_ORDER_6


findany = $(strip $(foreach W,$1,$(findstring $W,$2)))

ifneq (,$(call findany,login viz,$(HOSTNAME)))
#  #We are on DIaL
  CC = icc -std=c11
  #MYFLAGS += -qopenmp
  MYFLAGS += -no-multibyte-chars -mkl
  #MYFLAGS += -fopenmp
  MYFLAGS += -xCORE-AVX512 #-xHost #-ipo -fp-model fast=2
  #MYFLAGS += -I/opt/gridware/depots/e2b91392/el7/pkg/apps/R/4.1.0/gcc-5.5.0+lapack-3.5.0+blas-3.6.0/lib64/R/include
  #MYFLAGS += -I/opt/gridware/depots/e2b91392/el7/pkg/apps/R/3.6.3/gcc-5.5.0+lapack-3.5.0+blas-3.6.0/lib64/R/include
  MYFLAGS += -I$$FFTW3_DOUBLEINCLUDE
  MYFLAGS += -DMKL_LAPACKE
  MYLIBS += -lmkl_core -liomp5 -lmkl_intel_thread -lmkl_rt -lmkl_intel_lp64 #-lpthread
  MYLIBS += -L$$FFTW3_DOUBLELIB
  #MYFLAGS += -I$$LAPACKINCLUDE
  #MYLIBS += -L$$LAPACKLIB
  #MYLIBS += -L$$OPENBLASLIB -lopenblas
else ifeq ($(USER),pavel)
  MYFLAGS += -I/opt/OpenBLAS/include -L/opt/OpenBLAS/lib
  #MYFLAGS += -g -Og -fsanitize=address -static-libasan
  MYLIBS += -lopenblas #-nostartfiles
else
  MYFLAGS += -I/usr/share/R/include -flto -fPIC
  #MYFLAGS += -g -Og -fsanitize=address -static-libasan
  MYLIBS += -llapacke -llapack -lblas
endif

all: organic_main organic_conv

# Take file identifiers and change empty prefix by organic_ % empty suffix by .o
# Compile all dependencies $^ to target $@
organic_main: $(MYFILES:%=organic_%.o) $(RNG:%=%.o) organic_main.o
	$(CC) $(MYFLAGS) -o $@ $^ $(MYLIBS)

organic_conv: $(MYFILES:%=organic_%.o) $(RNG:%=%.o) organic_conv.o
	$(CC) $(MYFLAGS) -o $@ $^ $(MYLIBS)

organic_hmc.so: $(MYFILES:%=organic_%.o) $(RNG:%=%.o) organic_R.o
	$(CC) $(MYFLAGS) -shared -o $@ $^ $(MYLIBS)

min_test: minimal_test.o
	$(CC) $(MYFLAGS) -o $@ $^ $(MYLIBS)

# For every object .o check for .c and .h file with same name
# Compile first dependency $< (the .c file)
%.o: %.c %.h Makefile $(PARAM_FILES)
	$(CC) $(MYFLAGS) -c $<

print-%  : ; @echo $* = $($*)

clean:
	rm -f *.o *.so exact organic_hmc organic_main organic_conv

distclean:
	rm -f *.o *.so *.exe *~ tags exact organic_hmc organic_main organic_conv
