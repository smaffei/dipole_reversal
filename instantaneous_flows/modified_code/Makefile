F90 = gfortran -I. -ffree-line-length-0
FFTW = -lfftw3 -L /Users/sam/Documents/University/phd/fftw-3.3.8/fftw_installed/lib
dipole_bound :  dipole_bound_mod.F90 subs.F90
	$(F90) -c subs.F90
	$(F90) -o dipole_bound dipole_bound_mod.F90 subs.o $(FFTW)
