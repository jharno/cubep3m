fort=mpif90
xflags=-fpp -O3 -shared-intel -openmp -mt_mpi -mcmodel=medium -traceback
oflags=-c $(xflags)
dflags=-O0 -C -B -check all -fp-stack-check -CB -fpe0

F90_FILES := $(wildcard *.f90)
OBJ_FILES := $(addprefix ,$(notdir $(F90_FILES:.f90=.o)))
EXC_FILES := $(addprefix ,$(notdir $(OBJ_FILES:.o=.x)))

MOD_FILES := $(wildcard ../Modules/*.f90)
LNK_FILES := $(wildcard ../Modules/*.o)

PENCIL_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include
PENCIL_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib

CUR_DIR=`basename "$(CURID)"`
MOD_DIR=../Modules

all: Modules $(EXC_FILES)
	@echo "Finished making modules: $(value MOD_FILES)"

Modules:
	cd $(MOD_DIR); make; cd $(CUR_DIR)

%.x: %.o 
	@echo "Files to be linked: $(value LNK_FILES)"
	$(fort) $(xflags)  -I$(SCINET_FFTW_INC) -I$(PENCIL_INC) -o $@ $< $(value LNK_FILES) -L$(SCINET_FFTW_LIB) -L$(PENCIL_LIB) -lp3dfft -lfftw3f

%.o: %.f90 
	$(fort) $(oflags) -I$(MOD_DIR) -I$(SCINET_FFTW_INC) -I$(PENCIL_INC) -o $@ $< -L$(PENCIL_LIB) -L$(SCINET_FFTW_LIB) -lp3dfft -lfftw3f

clean:
	rm *.mod *.o *.x *~
