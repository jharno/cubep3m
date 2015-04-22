fort=mpif90
xflags=-fpp -O3 -shared-intel -openmp -mt_mpi -mcmodel=medium -traceback
oflags=-c $(xflags)
dflags=-O0 -C -B -check all -fp-stack-check -CB -fpe0

MOD_FILES := $(wildcard *.f90)
OBJ_FILES := $(addprefix ,$(notdir $(MOD_FILES:.f90=.o)))

PENCIL_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include
PENCIL_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib

all: Modules

Modules: $(OBJ_FILES)
	@echo "Finished making modules: $(value MOD_FILES)"

include Dependencies

%.o: %.f90 
	$(fort) $(oflags) -o $@ $< -I$(PENCIL_INC) -L$(PENCIL_LIB) -lp3dfft -lfftw3f

clean:
	rm *.mod *.o *.x *~
