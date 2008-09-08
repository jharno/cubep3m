SHELL = /bin/sh

FC=mpif77

## Mckenzie 
#FFTWLIB=/home/merz/lib/fftw-2.1.5_intel8/lib
#FFTWINC=/home/merz/lib/fftw-2.1.5_intel8/include
## Lobsters
FFTWLIB=/scratch/merz/fftw/lib
FFTWINC=/scratch/merz/fftw/include

## Mckenzie
#LDLIBS= -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
## Lobsters
LDLIBS= -lrfftw_mpi -lrfftw -lfftw_mpi -lfftw -lm

## Mckenzie
#FFLAGS=-fpp -g -w -O3 -axN -DDEBUG #-DDEBUG_LOW #-openmp 
## Lobsters
FFLAGS=-fpp -g -w -O3 -tpp1 -DDEBUG -DIA64 #-DDEBUG_LOW #-openmp 
# OPTIONAL --- -DFIXED_DIS :: if I've fixed the displacement routine to do np=nc
#              -DDEBUG_LOW :: debugging fftw problems
#              -DDEBUG :: general debugging 

OBJS= checknoise.o denstats.o displacement.o dist_init.o dmpowerspectrum.o fftw3d.o functions.o kernel.o mpi_initialization.o noisemap.o  phi_buffer.o powerspectrum.o transferfnc.o zeldovich.o 

run: /tmp/dist_init 

/tmp/dist_init: $(OBJS)
	$(FC) $(FFLAGS) $^ -o $@ -L$(FFTWLIB) $(LDLIBS) 

dist_init.o: dist_init.f90
	$(FC) $(FFLAGS) -c $<  

checknoise.o: checknoise.f90
	$(FC) $(FFLAGS) -c $<

denstats.o: denstats.f90
	$(FC) $(FFLAGS) -c $<

displacement.o: displacement.f90
	$(FC) $(FFLAGS) -c $<

dmpowerspectrum.o: dmpowerspectrum.f90
	$(FC) $(FFLAGS) -c $<

fftw3d.o: fftw3d.f90
	$(FC) $(FFLAGS) -I$(FFTWINC) -c $<

functions.o: functions.f90
	$(FC) $(FFLAGS) -c $<

kernel.o: kernel.f90
	$(FC) $(FFLAGS) -c $<

mpi_initialization.o: mpi_initialization.f90
	$(FC) $(FFLAGS) -c $<

noisemap.o: noisemap.f90
	$(FC) $(FFLAGS) -c $<

phi_buffer.o: phi_buffer.f90
	$(FC) $(FFLAGS) -c $<

powerspectrum.o: powerspectrum.f90
	$(FC) $(FFLAGS) -c $<

transferfnc.o: transferfnc.f90
	$(FC) $(FFLAGS) -c $<

zeldovich.o: zeldovich.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.pc dist_init 

tar:
	tar -cvjf dist_init.tar.bz2 *.par *.fh *.f90 Makefile tf.wmap 
