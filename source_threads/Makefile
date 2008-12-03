SHELL = /bin/sh

FC=mpif90

#TACC Lonestar	
FFTWLIB=$(TACC_FFTW2_LIB) #/home/merz/lib/fftw-2.1.5_intel8/lib
FFTWINC=$(TACC_FFTW2_INC) #/home/merz/lib/fftw-2.1.5_intel8/include

#Mckenzie	
LDLIBS= -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

#Mckenzie
FFLAGS=-fpp -g -w -O3 -axN -DDIAG #-DNGP -DMHD #-check bounds -DDEBUG #-DDEBUG_LOW #-DDEBUG_VEL -DDEBUG_CCIC -DDEBUG_RHOC -openmp -DDEBUG_LOW

#OPTIONAL FLAGS :: -DMHD -DDEBUG_CCIC -DDEBUG_CRHO -DDEBUG_RHOC -DPOINTSRC ##-DBFTEST

OBJS= mpi_tvd_mhd.o gas_density_buffer.o fine_gas_mass.o checkpoint.o coarse_cic_mass.o coarse_cic_mass_buffer.o coarse_force.o coarse_force_buffer.o coarse_mass.o coarse_max_dt.o coarse_mesh.o coarse_power.o coarse_velocity.o cubepm.o delete_particles.o fftw3ds.o fine_cic_mass.o fine_cic_mass_buffer.o fine_mesh.o fine_ngp_mass.o fine_ngp_mass_halos.o fine_velocity.o halofind.o fftw2.o kernel_initialization.o link_list.o mpi_initialization.o particle_initialization.o particle_mesh.o particle_pass.o projection.o report_pair.o set_pair.o indexedsort.o timers.o timestep.o update_position.o variable_initialization.o 

run: cubep3m

cubep3m: $(OBJS)
	$(FC) $(FFLAGS) $^ -o $@ -L$(FFTWLIB) $(LDLIBS) 

mpi_tvd_mhd.o: mpi_tvd_mhd.f90
	$(FC) $(FFLAGS) -c $<

gas_density_buffer.o: gas_density_buffer.f90 mpi_tvd_mhd.o
	$(FC) $(FFLAGS) -c $<

fine_gas_mass.o: fine_gas_mass.f90 mpi_tvd_mhd.o
	$(FC) $(FFLAGS) -c $<

cubepm.o: cubepm.f90 mpi_tvd_mhd.o
	$(FC) $(FFLAGS) -c $<  

checkpoint.o: checkpoint.f90 mpi_tvd_mhd.o
	$(FC) $(FFLAGS) -c $<

coarse_cic_mass.o: coarse_cic_mass.f90
	$(FC) $(FFLAGS) -c $<

coarse_cic_mass_buffer.o: coarse_cic_mass_buffer.f90
	$(FC) $(FFLAGS) -c $<

coarse_force.o: coarse_force.f90
	$(FC) $(FFLAGS) -c $<

coarse_force_buffer.o: coarse_force_buffer.f90
	$(FC) $(FFLAGS) -c $<

coarse_mass.o: coarse_mass.f90 mpi_tvd_mhd.o
	$(FC) $(FFLAGS) -c $<

coarse_max_dt.o: coarse_max_dt.f90
	$(FC) $(FFLAGS) -c $<

coarse_mesh.o: coarse_mesh.f90
	$(FC) $(FFLAGS) -c $<

coarse_power.o: coarse_power.f90
	$(FC) $(FFLAGS) -c $<

coarse_velocity.o: coarse_velocity.f90 mpi_tvd_mhd.o
	$(FC) $(FFLAGS) -c $<

delete_particles.o: delete_particles.f90
	$(FC) $(FFLAGS) -c $<

fftw3ds.o: fftw3ds.f90
	$(FC) $(FFLAGS) -I$(FFTWINC) -c $<

fine_cic_mass.o: fine_cic_mass.f90
	$(FC) $(FFLAGS) -c $<

fine_cic_mass_buffer.o: fine_cic_mass_buffer.f90
	$(FC) $(FFLAGS) -c $<

fine_mesh.o: fine_mesh.f90 mpi_tvd_mhd.o
	$(FC) $(FFLAGS) -c $<

fine_ngp_mass.o: fine_ngp_mass.f90
	$(FC) $(FFLAGS) -c $<

fine_ngp_mass_halos.o: fine_ngp_mass_halos.f90
	$(FC) $(FFLAGS) -c $<

fine_velocity.o: fine_velocity.f90 mpi_tvd_mhd.o
	$(FC) $(FFLAGS) -c $<

halofind.o: halofind.f90
	$(FC) $(FFLAGS) -c $<

fftw2.o: fftw2.f90
	$(FC) $(FFLAGS) -I$(FFTWINC) -c $<

kernel_initialization.o: kernel_initialization.f90
	$(FC) $(FFLAGS) -c $<

link_list.o: link_list.f90
	$(FC) $(FFLAGS) -c $<

mpi_initialization.o: mpi_initialization.f90
	$(FC) $(FFLAGS) -c $<

particle_initialization.o: particle_initialization.f90
	$(FC) $(FFLAGS) -c $<

particle_mesh.o: particle_mesh.f90
	$(FC) $(FFLAGS) -c $<

particle_pass.o: particle_pass.f90
	$(FC) $(FFLAGS) -c $<

projection.o: projection.f90
	$(FC) $(FFLAGS) -c $<

report_pair.o: report_pair.f90
	$(FC) $(FFLAGS) -c $<

set_pair.o: set_pair.f90
	$(FC) $(FFLAGS) -c $<

indexedsort.o: indexedsort.f90
	$(FC) $(FFLAGS) -c $<

timers.o: timers.f90
	$(FC) $(FFLAGS) -c $<

timestep.o: timestep.f90 mpi_tvd_mhd.o
	$(FC) $(FFLAGS) -c $<

update_position.o: update_position.f90
	$(FC) $(FFLAGS) -c $<

variable_initialization.o: variable_initialization.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.pc cubep3m mpi_tvd_mhd.mod

tar:
	tar -cvjf cubepm.tar.bz2 *.c *.par *.fh *.f90 Makefile kernels
