PROGRAM TIMESTEPPING_INDUCTION
USE SUBS
IMPLICIT NONE

! Reads in data for both B and u and calculates d g10/dt.

! The file for B is the usual Gauss coefficients, in nT
! The file for u are the Gauss coefficients describing the flow on the unit sphere in yr^(-1). They must be multiplied by the CMB radius to give dimensions of m/yr (or km/yr).

TYPE (HARMONIC_STRUCTURE), ALLOCATABLE, DIMENSION(:) :: HARMONICS

REAL( KIND = 8) ::  B_LOCATION(3), THETA, PHI, PHI_OBS
INTEGER :: L, M, SINCOS, POLTOR, LMAX_B_OBS, LMAX_U, IOS, LMAX_FILE, NP, i, HARMONIC, LMAX_SV, NTHETA_GRID, NPHI_GRID, I_PHI, I_THETA, NUM_LEGENDRE_SV, INDEX_PLM, NUMBER_REALISATIONS
CHARACTER(300) :: FILENAME_B, FILENAME_U, JUNK, FILENAME_SV, FORM_B, FORM_SV, STRNB, STRNSV
CHARACTER(1) :: STRING_HARM_TYPE2

INTEGER :: J, LMAX_U_FILE, INDEX, REVERSE
REAL( KIND = 8), ALLOCATABLE, DIMENSION(:) :: TOR_U, POL_U, POL_B_DOT, LEGENDRE_GAUW, COSTHETA_GRID
REAL( KIND = EXTRA_LONG_REAL), ALLOCATABLE, DIMENSION(:,:)  ::  TEMP_THETA_TRANSFORM_DALF, TEMP_THETA_TRANSFORM_ALF
REAL( KIND = 8), ALLOCATABLE :: GRAD_H_B_R(:,:,:), B_R(:,:), DIV_H_U_H(:,:), U_H(:,:,:),ALF_LOCATION(:,:), DALF_LOCATION(:,:), &
ALF(:,:), DALF(:,:), NONLIN_BR_CMB(:,:), SPEC(:), LEGENDRE_INV(:,:), ONE_DIV_SIN_THETA(:), FLOW_L_SPEC(:), SPEC_L(:), GAUSS(:), SPEC_U_TOR(:), SPEC_U_POL(:)
REAL( KIND = 8) :: COSTHETA_LOC, SINTHETA_LOC, PHI_DEP, DERV_PHI_DEP, B_DOT_LOCATION(3), F_LOC, MAX_INTENSITY_RATE, NORM_FACTOR, KE, TOTAL_RMS, SCALE_FACTOR, VALUE(1:4), DIPOLE_TILT_RATE, NORM_G
REAL( KIND = EXTRA_LONG_REAL) :: LOCATIONS(1)

REAL( KIND = 8), ALLOCATABLE :: GAUSS_N(:), SPEC_P(:)
REAL( KIND = 8) :: alpha, M1, M2, t0, tf, dt, eta, RHS
INTEGER :: N_IT

PRINT*, 'ENTER MAX DEGREE FOR U'
READ*, LMAX_U

PRINT*, 'ENTER MAX DEGREE FOR B'
READ*, LMAX_B_OBS

PRINT*, 'ENTER FILENAME FOR U'
READ*, FILENAME_U

PRINT*, 'ENTER FILENAME FOR B'
READ*, FILENAME_B

! this actually seems to be doing nothing
!PRINT*, 'ENTER FILENAME FOR OUTPUT SV'
!READ*, FILENAME_SV

PRINT*, 'MAGNETIC DIFFUSITVITY (in m^2/s)' 
READ*, eta

PRINT*, 'DIFFUSION TIMESTEPPING METHOD PARAMETER? (1/2 = Crank-Nicholson, 0 = Euler)' 
READ*, alpha

PRINT*, 'INITIAL YEAR' 
READ*, t0

PRINT*, 'FINAL YEAR' 
READ*, tf

PRINT*, 'TIMESTEP (in yr)' 
READ*, dt

PRINT*, 'CHANGE THE SIGN OF THE FLOW? (0 = NO, 1 = YES)'
READ*, REVERSE


PRINT*, ''


! Read in Gauss coefficients for magnetic field.
PRINT*, 'Read in Gauss coefficients for initial magnetic field'

 OPEN(11, FILE = FILENAME_B, STATUS = 'OLD', FORM = 'FORMATTED', & 
                             IOSTAT = IOS, ACTION = 'READ')
    IF( IOS .NE. 0) THEN
    PRINT*, 'ERROR IN OPENING FILE ', FILENAME_B
    STOP
    ENDIF
   

     NP = LMAX_B_OBS*(LMAX_B_OBS +2)
     ALLOCATE( GAUSS(1:NP ), GAUSS_N(1:NP ))
     GAUSS(:) = 0.0_8

     READ(11, *) JUNK
     J = 1      
    DO
      READ(11,*, END = 1001) L, M,  VALUE(1), VALUE(2)
      IF( M .NE. 0) THEN
	GAUSS(J) = VALUE(1)
	GAUSS(J+1) = VALUE(2)
	J = J + 2
      ELSE
	GAUSS(J) = VALUE(1)
	J = J + 1
      ENDIF
      IF( J > NP ) EXIT
    ENDDO


1001 CONTINUE


! Read in coefficients for flow.
PRINT*, 'Read in coefficients for flow'

  OPEN(11, FILE = FILENAME_U, STATUS = 'OLD', FORM = 'FORMATTED', &
                             IOSTAT = IOS, ACTION = 'READ')
    IF( IOS .NE. 0) THEN
    PRINT*, 'ERROR IN OPENING FILE ', FILENAME_U
    STOP
    ENDIF

NP = LMAX_U * (LMAX_U +2)
ALLOCATE( SPEC_U_TOR(1:NP), SPEC_U_POL(1:NP)   )
SPEC_U_TOR(:) = 0.0_LONG_REAL
SPEC_U_POL(:) = 0.0_LONG_REAL

READ(11, *) JUNK

     J = 1
     DO
      READ(11,*, END = 1002) L, M,  VALUE(1), VALUE(2), VALUE(3), VALUE(4)
      IF( M .NE. 0) THEN
       SPEC_U_TOR(J) = (-1)**REVERSE * VALUE(3)
       SPEC_U_TOR(J+1) = (-1)**REVERSE * VALUE(4)
       SPEC_U_POL(J) = (-1)**REVERSE * VALUE(1)
       SPEC_U_POL(J+1) = (-1)**REVERSE * VALUE(2)
       J = J + 2
      ELSE
       SPEC_U_TOR(J) = (-1)**REVERSE * VALUE(3)
       SPEC_U_POL(J) = (-1)**REVERSE * VALUE(1)
       J = J + 1
       ENDIF
IF( J > NP ) EXIT
       ENDDO



1002 CONTINUE

     CLOSE(11)


! Multiply the spectral coefficients so that u is given in km/yr
! ******************
! THINK ABOUT THIS
! ******************
SPEC_U_TOR(:) = SPEC_U_TOR(:) * CMB_RADIUS / 1000.0_LONG_REAL
SPEC_U_POL(:) = SPEC_U_POL(:) * CMB_RADIUS / 1000.0_LONG_REAL

! Transform eta so that everything is in Km, yr units
eta = eta * (3600*24*365) / (1000.0_LONG_REAL)**2

PRINT*, 'Preparing Spehrical Harmonics'
! only one set of harmonic coefficients are required - the SV coefficients which have the highest truncation. Those for B[obs] and u are simply just the first section of the coefficients.      
      LMAX_SV = LMAX_U + LMAX_B_OBS
      ALLOCATE(HARMONICS(1: LMAX_SV * (LMAX_SV + 2) ) )
      HARMONIC = 1       
      DO L = 1, LMAX_SV
      DO M = 0, L
      DO SINCOS = COSINE_HARMONIC, SINE_HARMONIC  !cos is 1, sin is 2.
      IF( M .eq. 0 .AND. SINCOS .eq. 2) CYCLE
      HARMONICS(HARMONIC)%M = M
      HARMONICS(HARMONIC)%SINCOS = SINCOS
      HARMONICS(HARMONIC)%L = L
      HARMONIC= HARMONIC+1
      ENDDO
      ENDDO
      ENDDO

      OPEN(21, FILE = 'HARMONICS.DAT', FORM = 'FORMATTED', STATUS = 'REPLACE')
      DO I=1, LMAX_SV * (LMAX_SV + 2)
      IF(HARMONICS(I)%SINCOS .eq. COSINE_HARMONIC) WRITE(STRING_HARM_TYPE2,'(A)') 'C'
      IF(HARMONICS(I)%SINCOS .eq. SINE_HARMONIC) WRITE(STRING_HARM_TYPE2,'(A)') 'S'
   
      WRITE(21,'(i3,3x, i3, 3x, A1)')  & 
      HARMONICS(I)%L, HARMONICS(I)%M, STRING_HARM_TYPE2
      ENDDO
      CLOSE(21)
      
     
   NUM_LEGENDRE_SV = (LMAX_SV+1) * (LMAX_SV+2) / 2
   NTHETA_GRID =  LMAX_SV + 1  !The maximum degree of the colatitude-integrand is Lmax_B_OBS + Lmax_u + Lmax_sv = 2 * Lmax_sv. So number of pts is 1/2 ( 2 * Lmax_sv ) + 1.
   NPHI_GRID = 3 * LMAX_SV + 1 !(using the usual 3/2 rule).

   ALLOCATE( FFT_TRANSFORM_ARRAY(1:NTHETA_GRID,0:NPHI_GRID-1) )

! FFTs in PHI, REAL to Half complex
      CALL DFFTW_PLAN_MANY_R2R(PLAN_FFT_R2HC,       &  ! PLAN  
      1,                                            &  ! Rank of transforms
      NPHI_GRID,                                    &  ! size of each transform
      NTHETA_GRID,                                  &  ! how many
      FFT_TRANSFORM_ARRAY,                          &  ! Input array
      SIZE(FFT_TRANSFORM_ARRAY),                    &  ! subset size of transform (whole set)
      NTHETA_GRID,                                  &  ! Stride for transform
      1,                                            &  ! Gap between transforms in array
      FFT_TRANSFORM_ARRAY,                          &  ! Output array
      SIZE(FFT_TRANSFORM_ARRAY),                    &  ! subset size of output array
      NTHETA_GRID,                                  &  ! output stride
      1,                                            &  ! Gap between outputs for transforms.
      FFTW_R2HC, TRANSFORM_OPTIMISATION)                        ! Transform type; flags
 


    
  ALLOCATE( COSTHETA_GRID(NTHETA_GRID), LEGENDRE_GAUW(NTHETA_GRID), ALF(NTHETA_GRID, NUM_LEGENDRE_SV), DALF(NTHETA_GRID, NUM_LEGENDRE_SV), &
            TEMP_THETA_TRANSFORM_DALF(NTHETA_GRID, NUM_LEGENDRE_SV), TEMP_THETA_TRANSFORM_ALF(NTHETA_GRID, NUM_LEGENDRE_SV), LEGENDRE_INV(NTHETA_GRID, NUM_LEGENDRE_SV) , &
             ONE_DIV_SIN_THETA(1:NTHETA_GRID)  )
             
! Evaluate B_r and its horizontal derivative on (theta, phi) grid
   CALL GAUWTS ( NTHETA_GRID, COSTHETA_GRID, LEGENDRE_GAUW )
      
! Get the Associated Legendre Functions (ALF) at the theta grid points, with their theta derivatives (DALF), of the ALF's up to degree LMAX_PRECOMP
   CALL GET_LEGENDRE_FUNCTIONS( REAL(COSTHETA_GRID, KIND = EXTRA_LONG_REAL), LMAX_SV, LMAX_SV, TEMP_THETA_TRANSFORM_ALF, TEMP_THETA_TRANSFORM_DALF) 
   ALF = REAL(  TEMP_THETA_TRANSFORM_ALF, KIND = 8)
   DALF = REAL( TEMP_THETA_TRANSFORM_DALF, KIND = 8)
   ONE_DIV_SIN_THETA(:) = 1.0_LONG_REAL / SQRT( 1.0_LONG_REAL - COSTHETA_GRID(:)**2 )

   DO I = 1, NUM_LEGENDRE_SV
   LEGENDRE_INV(1:NTHETA_GRID,I) = ALF(1:NTHETA_GRID,I) * LEGENDRE_GAUW(1:NTHETA_GRID) 
   ENDDO
   
  ALLOCATE( NONLIN_BR_CMB(NTHETA_GRID, 0:NPHI_GRID-1), SPEC(1:LMAX_SV * (LMAX_SV + 2)), SPEC_P(1:LMAX_SV * (LMAX_SV + 2)) )
  

    ALLOCATE( B_R(1:NTHETA_GRID, 0:NPHI_GRID-1), GRAD_H_B_R(1:NTHETA_GRID,0:NPHI_GRID-1,2) )

      ALLOCATE( U_H(1:NTHETA_GRID, 0:NPHI_GRID-1,2), DIV_H_U_H(1:NTHETA_GRID, 0:NPHI_GRID-1) )

ALLOCATE( SPEC_L(1: LMAX_SV))
 ALLOCATE( FLOW_L_SPEC(1:LMAX_U) )

SPEC_L(:) = 0.0_8
NONLIN_BR_CMB(:,:) = 0.0_LONG_REAL
GRAD_H_B_R(:,:,:) = 0.0_LONG_REAL
U_H(:,:,:) = 0.0_LONG_REAL
DIV_H_U_H(:,:) = 0.0_LONG_REAL

OPEN(12, FILE = 'MF_COEFFS.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')
OPEN(13, FILE = 'SV_COEFFS.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')
! SV coefficients do not include diffusion

! write initial Gauss coefficients and SV
PRINT*, 'Write initial Gauss and SV coefficients'

WRITE(12,'(F10.3)',advance="no")  t0
WRITE(STRNB, '(I5)') LMAX_B_OBS * (LMAX_B_OBS + 2)
FORM_B = '(' // TRIM(STRNB) // 'ES23.15)'
WRITE(12,FORM_B)  (GAUSS(I), I=1,LMAX_B_OBS * (LMAX_B_OBS + 2))

! Initial SV coefficients (make it a subroutine?)
! Calculate Br on the CMB, on (theta,phi) grid
CALL EVALUATE_B_R_GRID( NTHETA_GRID, NPHI_GRID, GAUSS(:), HARMONICS, LMAX_B_OBS, LMAX_SV, B_R, GRAD_H_B_R, ALF, DALF, ONE_DIV_SIN_THETA )
! Evaluate u_h and its horizontal divergence on (theta, phi) grid
CALL EVALUATE_U_H( NTHETA_GRID, NPHI_GRID, SPEC_U_TOR(:), SPEC_U_POL(:),  U_H, DIV_H_U_H, ALF, DALF, HARMONICS, LMAX_U, LMAX_SV, ONE_DIV_SIN_THETA  )
! Evaluate nonlinear term in Br_dot (at CMB) on the physical grid
DO I_THETA = 1, NTHETA_GRID
DO I_PHI = 0, NPHI_GRID-1
  NONLIN_BR_CMB(I_THETA, I_PHI) = - U_H(I_THETA, I_PHI,1) * GRAD_H_B_R(I_THETA, I_PHI,1) &
				  - U_H(I_THETA, I_PHI,2) * GRAD_H_B_R(I_THETA, I_PHI,2) &
			          - B_R(I_THETA, I_PHI) * DIV_H_U_H(I_THETA, I_PHI)
ENDDO
ENDDO
! NONLIN_BR_CMB is produced from a flow in units of km/yr. Need to multiply by 1000 to get a magnetic field in units of nT/yr.
! *********************
! Do I want this?
! *********************
NONLIN_BR_CMB(:,:) = NONLIN_BR_CMB(:,:) * 1000.0_LONG_REAL
! Calculate SH coefficients of the nonlinear term
CALL REAL_2_SPEC( NONLIN_BR_CMB, SPEC, LMAX_SV, LMAX_SV, HARMONICS, NPHI_GRID, NTHETA_GRID , LEGENDRE_INV )
! Transform to Gauss-type coefficients:
DO J=1, LMAX_SV * (LMAX_SV + 2)
  SPEC(J) = SPEC(J) * (CMB_RADIUS / ES_RADIUS )**(HARMONICS(J)%L+2) / REAL( HARMONICS(J)%L + 1, KIND = LONG_REAL)
ENDDO
! *********************
! SHOULD ADD DIFFUSION
! *********************
WRITE(13,'(F10.3)',advance="no")  t0
WRITE(STRNSV, '(I5)') LMAX_SV * (LMAX_SV + 2)
FORM_SV = '(' // TRIM(STRNSV) // 'ES23.15)'
WRITE(13,FORM_SV)  (SPEC(I), I=1,LMAX_SV * (LMAX_SV + 2))


PRINT*, 'Calculating integration coefficients'
! number of iterations
N_IT = 1
! Time-independent matrices
DO J=1, LMAX_B_OBS * (LMAX_B_OBS + 2)
    ! M1 = 1/dt
    M1 = 1/dt
    ! M2 = eta*l*(l+1)/c^2
    M2 =  eta * (HARMONICS(J)%L + 1) * HARMONICS(J)%L / CMB_RADIUS**(2)
ENDDO


! ***********************************
! 	big loop should start here
! ***********************************
PRINT*, ''
PRINT*, 'BEGIN OF TIME INTEGRATION'
PRINT*, ''
DO

  IF ( t0 + dt*N_IT .GT. tf) EXIT
  
  PRINT*, '*********************************'
  PRINT*, 'TIME: ', t0 + dt*(N_IT)
  PRINT*, 'time integration'
  DO J=1, LMAX_B_OBS * (LMAX_B_OBS + 2)
    IF ( N_IT .EQ. 1 ) THEN
      ! first timestep (AB1 = Euler's method)
      ! (M1 + alpha M2) GAUSS(t+dt) = (M - (1- alpha) M2) GAUSS(t) + SPEC(t)
      RHS = (M1 -(1-alpha)*M2) * GAUSS(J) + SPEC(J)
    ELSE
      ! Further timesteps (AB2)
      ! (M1 + alpha M2) GAUSS(t+dt) = (M - (1- alpha) M2) GAUSS(t) + 3/2 SPEC(t) - 1/2 SPEC(t-dt)
      RHS = (M1 -(1-alpha)*M2) * GAUSS(J) + (3.0/2.0)*SPEC(J) - (1.0/2.0)*SPEC_P(J)
    ENDIF
    GAUSS_N(J) = RHS / (M1 + alpha*M2)
  ENDDO

  
  PRINT*, 'Write updated Gauss coefficients'
  WRITE(12,'(F10.3)',advance="no")  t0 + N_IT*dt
  WRITE(12,FORM_B)  (GAUSS_N(I), I=1,LMAX_B_OBS * (LMAX_B_OBS + 2))

  PRINT*, 'update for the next timestep'
  SPEC_P = SPEC
  ! Update SV coefficients (make it a subroutine?)
  ! Calculate Br on the CMB, on (theta,phi) grid
  CALL EVALUATE_B_R_GRID( NTHETA_GRID, NPHI_GRID, GAUSS_N(:), HARMONICS, LMAX_B_OBS, LMAX_SV, B_R, GRAD_H_B_R, ALF, DALF, ONE_DIV_SIN_THETA )
  ! Evaluation of u_h is not needed anymore, flow is taken to be constant for now.
  ! Evaluate nonlinear term in Br_dot (at CMB) on the physical grid
  DO I_THETA = 1, NTHETA_GRID
  DO I_PHI = 0, NPHI_GRID-1
    NONLIN_BR_CMB(I_THETA, I_PHI) = - U_H(I_THETA, I_PHI,1) * GRAD_H_B_R(I_THETA, I_PHI,1) &
				    - U_H(I_THETA, I_PHI,2) * GRAD_H_B_R(I_THETA, I_PHI,2) &
				    - B_R(I_THETA, I_PHI) * DIV_H_U_H(I_THETA, I_PHI)
  ENDDO
  ENDDO
  ! NONLIN_BR_CMB is produced from a flow in units of km/yr. Need to multiply by 1000 to get a magnetic field in units of nT/yr.
  ! *********************
  ! Do I want this?
  ! *********************
  NONLIN_BR_CMB(:,:) = NONLIN_BR_CMB(:,:) * 1000.0_LONG_REAL
  ! Calculate SH coefficients of the nonlinear term
  CALL REAL_2_SPEC( NONLIN_BR_CMB, SPEC, LMAX_SV, LMAX_SV, HARMONICS, NPHI_GRID, NTHETA_GRID , LEGENDRE_INV )
  ! Transform to Gauss-type coefficients:
  DO J=1, LMAX_SV * (LMAX_SV + 2)
    SPEC(J) = SPEC(J) * (CMB_RADIUS / ES_RADIUS )**(HARMONICS(J)%L+2) / REAL( HARMONICS(J)%L + 1, KIND = LONG_REAL)
  ENDDO
  ! *********************
  ! SHOULD ADD DIFFUSION
  ! *********************
  PRINT*, 'Write updated SV coefficients'
  WRITE(13,'(F10.3)',advance="no")  t0 + N_IT*dt
  WRITE(13,FORM_SV)  (SPEC(I), I=1,LMAX_SV * (LMAX_SV + 2))

  
  GAUSS = GAUSS_N
  N_IT = N_IT +1

ENDDO



! OUTPUTS and WRAP UP

CLOSE(12)
CLOSE(13)


! WRITE FINAL MF COEFFICIENTS SEPARATELY
OPEN(14, FILE = 'MF_COEFFS_FINAL.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')

  WRITE(14,'(A,I3)') 'FINAL MF WITH L_MAX = ', LMAX_B_OBS
  J = 1      
  DO
    IF( HARMONICS(J)%M .NE. 0) THEN
      WRITE(14,'(I3,x,I3,2(x,ES15.5,x))') HARMONICS(J)%L, HARMONICS(J)%M, GAUSS_N(J), GAUSS_N(J+1)
      J = J + 2
    ELSE
      WRITE(14,'(I3,x,I3,2(x,ES15.5,x))') HARMONICS(J)%L, HARMONICS(J)%M, GAUSS_N(J), 0
      J = J + 1
    ENDIF
    IF( J > LMAX_B_OBS*(LMAX_B_OBS +2) ) EXIT
  ENDDO
    
CLOSE(14)

! *****************************
! COMMENT OUT THIS LAST PART
! *****************************
GOTO 101 

! Calculate spectrum of flow
!       OPEN(23, FILE = 'FLOW_L_SPEC.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')
      
       FLOW_L_SPEC(:) = 0.0_LONG_REAL
       DO I=1, LMAX_U * (LMAX_U + 2)
       NORM_FACTOR = HARMONICS(I)%L * (HARMONICS(I)%L + 1) * 4.0_LONG_REAL * Pi / REAL( 2 * HARMONICS(I)%L + 1, KIND = LONG_REAL )
       FLOW_L_SPEC(HARMONICS(I)%L) = FLOW_L_SPEC(HARMONICS(I)%L) + (SPEC_U_TOR(I)**2 + SPEC_U_POL(I)**2) * NORM_FACTOR
       ENDDO
       
!       DO I = 1, LMAX_U
!       WRITE(23,*) I, FLOW_L_SPEC(I)/(4.0_LONG_REAL * Pi)
!       ENDDO
!       CLOSE(23)

       WRITE(6,'(A,x,x,F8.3)') 'RMS FLOW IS (KM/YR)', SQRT( SUM(FLOW_L_SPEC) / (4.0_8 * Pi) ) 
       WRITE(6,'(A,x,x,ES10.3)') 'RMS FLOW IS (M/S)', SQRT( SUM(FLOW_L_SPEC) / (4.0_8 * Pi) ) * 1000 / (3600 * 24 * 365)




 PRINT*, 'd/dt g^1_0 =', SPEC(1)
! the harmonics have already been transported to the Earth's surface 
 NORM_G = (GAUSS(1)**2 + GAUSS(2)**2 + GAUSS(3)**2) * SQRT( GAUSS(2)**2 + GAUSS(3)**2 )
 DIPOLE_TILT_RATE =  (SPEC(2) * GAUSS(1) * GAUSS(2) + SPEC(3) * GAUSS(1) * GAUSS(3) &
    - SPEC(1) * (GAUSS(2)**2 +GAUSS(3)**2) )  / NORM_G
 PRINT*, 'd/dt (dipole tilt) [rad/yr] =', DIPOLE_TILT_RATE
 PRINT*, 'd/dt (dipole tilt) [deg/yr] =', DIPOLE_TILT_RATE*180/Pi
 
 
! Calculate Lowes-spectrum on CMB

DO I=1, SIZE(HARMONICS)
SPEC_L(HARMONICS(I)%L) = SPEC_L(HARMONICS(I)%L)  + (ES_RADIUS/CMB_RADIUS)**(2*HARMONICS(I)%L + 4)  * REAL( HARMONICS(I)%L + 1, KIND = LONG_REAL) * SPEC(I)**2
ENDDO

PRINT*, 'SV COEFF'
DO I =1, LMAX_SV * (LMAX_SV + 2)
PRINT*, HARMONICS(I)%L, HARMONICS(I)%M, SPEC(I)

ENDDO
PRINT*, 'LOWES SPEC ON CORE SURFACE OF SV'
DO I =1, LMAX_SV
PRINT*, I, SPEC_L(I)
ENDDO

PRINT*, ''

PRINT*, 'LOWES SPEC ON EARTH''S SURFACE OF SV'
DO I =1, LMAX_SV
PRINT*, I, SPEC_L(I) / (ES_RADIUS/CMB_RADIUS)**(2*I + 4)
ENDDO

101 CONTINUE


       CALL DFFTW_DESTROY_PLAN(PLAN_FFT_R2HC)
       DEALLOCATE(FFT_TRANSFORM_ARRAY)    

       
STOP
     
END PROGRAM TIMESTEPPING_INDUCTION