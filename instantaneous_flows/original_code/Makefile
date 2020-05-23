FFTW =  -lfftw3 -L/apps/Modules/developers/libraries/fftw/3.3.5
F90 = gfortran -I. -ffree-line-length-0
dipole_bound :  dipole_bound.F90 subs.F90
	$(F90) -c subs.F90
	$(F90) -o dipole_bound dipole_bound.F90 subs.o $(FFTW) 
