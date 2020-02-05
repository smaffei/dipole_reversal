F90 = ifort # gfortran -ffree-line-length-0
timestepping_induction :  timestepping_induction.F90 /nfs/see-fs-01_users/earsmaf/SWIGS/OptimalFlow/subs.F90
	$(F90) -c   /nfs/see-fs-01_users/earsmaf/SWIGS/OptimalFlow/subs.F90
	$(F90) -o timestepping_induction timestepping_induction.F90 subs.o  -lfftw3 
    
