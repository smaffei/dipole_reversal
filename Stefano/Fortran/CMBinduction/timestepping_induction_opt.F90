PROGRAM TIMESTEPPING_INDUCTION_OPT
USE SUBS
IMPLICIT NONE

! Timesteps the induction equation and optises the velocity field at every timestep.
!
! The file for B is the usual Gauss coefficients, in nT
! The file for u are the Gauss coefficients describing the flow on the unit sphere in yr^(-1). They must be multiplied by the CMB radius to give dimensions of m/yr (or km/yr).

TYPE (HARMONIC_STRUCTURE), ALLOCATABLE, DIMENSION(:) :: HARMONICS

REAL( KIND = EXTRA_LONG_REAL) ::  B_LOCATION(3), THETA, PHI
INTEGER :: L, M, SINCOS, POLTOR, LMAX_B_OBS, LMAX_U, IOS, LMAX_FILE, NP, i, HARMONIC, LMAX_SV, NTHETA_GRID, NPHI_GRID, I_PHI, I_THETA, NUM_LEGENDRE_SV, INDEX_PLM, NUMBER_REALISATIONS
CHARACTER(300) :: FILENAME_B, FILENAME_U, JUNK, FILENAME_SV, FORM_B, FORM_SV, FORM_U, STRNB, STRNSV, STRNU
CHARACTER(1) :: STRING_HARM_TYPE2
INTEGER :: MIN_LOC(2)
INTEGER :: J, LMAX_U_FILE, INDEX, REVERSE, FLAG_U_INIT
REAL( KIND = EXTRA_LONG_REAL), ALLOCATABLE, DIMENSION(:) :: TOR_U, POL_U, POL_B_DOT, LEGENDRE_GAUW, COSTHETA_GRID
REAL( KIND = EXTRA_LONG_REAL), ALLOCATABLE, DIMENSION(:,:)  ::  TEMP_THETA_TRANSFORM_DALF, TEMP_THETA_TRANSFORM_ALF
REAL( KIND = EXTRA_LONG_REAL), ALLOCATABLE :: GRAD_H_B_R(:,:,:), B_R(:,:), DIV_H_U_H(:,:), U_H(:,:,:),ALF_LOCATION(:,:), DALF_LOCATION(:,:), &
ALF(:,:), DALF(:,:), NONLIN_BR_CMB(:,:), SPEC(:), LEGENDRE_INV(:,:), ONE_DIV_SIN_THETA(:), FLOW_L_SPEC(:), SPEC_L(:), B_R_SURF(:,:), B_T_SURF(:,:), B_P_SURF(:,:), F_SURF(:,:)
REAL( KIND = EXTRA_LONG_REAL) :: COSTHETA_LOC, SINTHETA_LOC, PHI_DEP, DERV_PHI_DEP, B_DOT_LOCATION(3), F_LOC, MAX_INTENSITY_RATE, NORM_FACTOR, KE, TOTAL_RMS, VALUE(1:4), DIPOLE_TILT_RATE, NORM_G
REAL( KIND = EXTRA_LONG_REAL) :: LOCATIONS(1)
REAL( KIND = EXTRA_LONG_REAL), ALLOCATABLE :: GAUSS(:), SPEC_U_TOR(:), SPEC_U_POL(:), G(:,:)
REAL( KIND = EXTRA_LONG_REAL), ALLOCATABLE :: B_R_DOT(:,:), EVEC_TOR(:), EVEC_POL(:)
REAL( KIND = EXTRA_LONG_REAL), ALLOCATABLE :: GAUSS_N(:), SPEC_P(:)
REAL( KIND = 8) :: alpha, M1, M2, t0, tf, dt, ETA, RHS
INTEGER :: N_IT, HARMONIC_TYPE
REAL( KIND = 8) :: D_LOC, KE_TOR, KE_POL

REAL( KIND = EXTRA_LONG_REAL) :: OBS_COLAT, OBS_LONG, PHI_OBS
INTEGER :: RESTRICTION
INTEGER :: MINIM_FLAG, OUTPUT_RATE
CHARACTER(300) :: FILENAME
REAL( KIND = 8) :: SCALE_FACTOR, TARGET_RMS
REAL( KIND = 8) :: FACTOR

! GENERAL PARAMETERS
PRINT*, 'ENTER MAX DEGREE FOR U'
READ*, LMAX_U

PRINT*, 'ENTER MAX DEGREE FOR B'
READ*, LMAX_B_OBS

PRINT*, 'ENTER FILENAME FOR U'
READ*, FILENAME_U

PRINT*, 'ENTER FILENAME FOR B'
READ*, FILENAME_B

PRINT*, 'FLOW INITIAL CONDITION? (0: FROM FILE; 1: OPTIMAL FLOW)'
READ*, FLAG_U_INIT

! TIMESTEPPING PARAMETERS
PRINT*, 'MAGNETIC DIFFUSITVITY (in m^2/s)' 
READ*, ETA 

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

! OPTIMISATION STEP PARAMETERS
PRINT*, 'ENTER COLAT/LONG (IN DEGREES) FOR OBSERVATION POINT'
READ*, OBS_COLAT, OBS_LONG

PHI_OBS = OBS_LONG * Pi / 180.0_8

PRINT*, 'ENTER TARGET RMS VALUE FOR FLOW'
READ*, TARGET_RMS

PRINT*, 'ENTER SCALE FACTOR FOR ARROWS'
READ*, SCALE_FACTOR

PRINT*, 'RESTRICTION OF FLOW? 0: NONE, 1: POL ONLY, 2: TOR ONLY, 3: COLUMNAR'
READ*, RESTRICTION

PRINT*, 'ENTER MINIMIZATION FLAG. 0: DIPOLE TILT; 1: INTENSITY AT OBSERVATION POINT; 2: INTENSITY CHANGE AT INTENSITY MINIMUM (SAA)'
READ*, MINIM_FLAG

PRINT*, 'ENTER OUTPUT RATE (IN NUMBER OF TIMESTEPS). FIRST AND LAST TIMESTEP ARE ALWAYS GOING TO BE OUTPUT'
READ*, OUTPUT_RATE

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
  CLOSE(11)

IF( L .LT. LMAX_B_OBS ) THEN
  WRITE(6,'(A)') '********************'
  WRITE(6,'(A)')  'FIELD MODEL FILE CONTAINS TOO FEW COEFFICIENTS.'
  WRITE(6,'(A,I4,A,I4)') 'LMAX REQUIRED: ', LMAX_B_OBS,' LMAX IN FILE: ',L
  WRITE(6,'(A)')  'FIELD MODEL IS PADDED WITH ZEROS.'
  WRITE(6,'(A)') '********************'
ENDIF

ALLOCATE( SPEC_U_TOR(1: LMAX_U * (LMAX_U + 2) ), SPEC_U_POL(1: LMAX_U * (LMAX_U + 2) )   )

! Initialise flow field
! in case I want the initial condition be given in a file  
IF ( FLAG_U_INIT.EQ.0 ) THEN 

  ! Read in coefficients for flow.
  PRINT*, 'Read in coefficients for flow'

    OPEN(11, FILE = FILENAME_U, STATUS = 'OLD', FORM = 'FORMATTED', &
			      IOSTAT = IOS, ACTION = 'READ')
      IF( IOS .NE. 0) THEN
      PRINT*, 'ERROR IN OPENING FILE ', FILENAME_U
      STOP
      ENDIF

  NP = LMAX_U * (LMAX_U +2)
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
  SPEC_U_TOR(:) = SPEC_U_TOR(:) * CMB_RADIUS / 1000.0_LONG_REAL
  SPEC_U_POL(:) = SPEC_U_POL(:) * CMB_RADIUS / 1000.0_LONG_REAL
  
ENDIF 
  
! Transform ETA so that everything is in Km, yr units
ETA = ETA * (3600*24*365) / (1000.0_LONG_REAL)**2

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
DEALLOCATE( TEMP_THETA_TRANSFORM_DALF, TEMP_THETA_TRANSFORM_ALF )

DO I = 1, NUM_LEGENDRE_SV
  LEGENDRE_INV(1:NTHETA_GRID,I) = ALF(1:NTHETA_GRID,I) * LEGENDRE_GAUW(1:NTHETA_GRID) 
ENDDO

ALLOCATE( NONLIN_BR_CMB(NTHETA_GRID, 0:NPHI_GRID-1), SPEC(1:LMAX_SV * (LMAX_SV + 2)), SPEC_P(1:LMAX_SV * (LMAX_SV + 2)) )

ALLOCATE( B_R_DOT(NTHETA_GRID, 0:NPHI_GRID-1) )

ALLOCATE( B_R(1:NTHETA_GRID, 0:NPHI_GRID-1), GRAD_H_B_R(1:NTHETA_GRID,0:NPHI_GRID-1,2) )

ALLOCATE( U_H(1:NTHETA_GRID, 0:NPHI_GRID-1,2), DIV_H_U_H(1:NTHETA_GRID, 0:NPHI_GRID-1) )

ALLOCATE( SPEC_L(1: LMAX_SV))
ALLOCATE( FLOW_L_SPEC(1:LMAX_U) )

! Calculate Br on the CMB, on (theta,phi) grid
CALL EVALUATE_B_R_GRID( NTHETA_GRID, NPHI_GRID, GAUSS(:), HARMONICS, LMAX_B_OBS, LMAX_SV, B_R, GRAD_H_B_R, ALF, DALF, ONE_DIV_SIN_THETA )

  
SPEC_L(:) = 0.0_8


OPEN(12, FILE = 'MF_COEFFS.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')
OPEN(13, FILE = 'SV_COEFFS.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')
IF ( FLAG_U_INIT.EQ.1 ) THEN
  OPEN(14, FILE = 'U_TOR_COEFFS.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')
  OPEN(15, FILE = 'U_POL_COEFFS.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')
ENDIF
! SV coefficients do not include diffusion

! prepare format for gauss coefficients files
WRITE(STRNB, '(I5)') LMAX_B_OBS * (LMAX_B_OBS + 2)
FORM_B = '(' // TRIM(STRNB) // 'ES23.15)'

WRITE(STRNSV, '(I5)') LMAX_SV * (LMAX_SV + 2)
FORM_SV = '(' // TRIM(STRNSV) // 'ES23.15)'

IF ( FLAG_U_INIT.EQ.1 ) THEN
  WRITE(STRNU, '(I5)')  LMAX_U * (LMAX_U + 2)
  FORM_U = '(' // TRIM(STRNU) // 'ES23.15)'
ENDIF

PRINT*, 'Calculating integration coefficients'
! number of iterations
N_IT = 1
! Time-independent matrices
DO J=1, LMAX_B_OBS * (LMAX_B_OBS + 2)
    ! M1 = 1/dt
    M1 = 1/dt
    ! M2 = ETA*l*(l+1)/c^2
    M2 =  ETA * (HARMONICS(J)%L + 1) * HARMONICS(J)%L / CMB_RADIUS**(2)
ENDDO


ALLOCATE(  ALF_LOCATION(1,1: NUM_LEGENDRE_SV) )
ALLOCATE(  DALF_LOCATION(1,1: NUM_LEGENDRE_SV)  )
ALLOCATE( B_R_SURF(1:NTHETA_GRID, 0:NPHI_GRID-1), B_T_SURF(1:NTHETA_GRID, 0:NPHI_GRID-1), B_P_SURF(1:NTHETA_GRID, 0:NPHI_GRID-1), F_SURF(1:NTHETA_GRID, 0:NPHI_GRID-1) )
ALLOCATE( EVEC_TOR(1: LMAX_U * (LMAX_U + 2) ), EVEC_POL(1: LMAX_U * (LMAX_U + 2) ) )
ALLOCATE( G(1: LMAX_U * (LMAX_U + 2),1:2 ) )

! ***********************************
! 	big loop should start here
! ***********************************

PRINT*, ''
PRINT*, 'BEGIN OF TIME INTEGRATION'
PRINT*, ''

PRINT*,t0 + dt*(N_IT-1)
PRINT*,tf
DO

  IF ( t0 + dt*(N_IT-1) .GT. tf) EXIT
  
  PRINT*, '*********************************'
  PRINT*, 'TIME: ', t0 + dt*(N_IT-1)
  
  PRINT*, 'time integration'
  
  SPEC(:) = 0.0
  
  !! OPTIMISATION !!
  IF ( FLAG_U_INIT.EQ.1 .OR. N_IT.GT.1) THEN
    ! *********************************************************************** !
    ! THIS SHOULD ALL BE IN A SUBROUTINE, BUT IT DID NOT WORK WHEN I TRIED IT
    SPEC_U_TOR(:) = 0.0
    SPEC_U_POL(:) = 0.0
    G(:,:) = 0.0
    ! IF I WANT TO OPTIMIZE THE INTENSITY OF THE SAA I NEED TO OVERWRITE SOME THINGS 
    IF( MINIM_FLAG.EQ.2) THEN
    ! find intensity minimum location and define it as LOCATION
      CALL EVALUATE_B_SURF_GRID(NTHETA_GRID, NPHI_GRID, GAUSS, HARMONICS, LMAX_B_OBS, LMAX_SV, B_R_SURF, B_T_SURF, B_P_SURF, ALF, DALF, ONE_DIV_SIN_THETA )
      F_SURF = SQRT( B_R_SURF**2 + B_T_SURF**2 + B_P_SURF**2 )
      MIN_LOC = MINLOC(F_SURF)
      OBS_COLAT = 180.0_LONG_REAL - MIN_LOC(1) * 180.0_LONG_REAL / REAL(NTHETA_GRID, KIND = LONG_REAL) ! COLAT = 0 is the North pole
      OBS_LONG  = MIN_LOC(2) * 360.0_LONG_REAL / REAL(NPHI_GRID, KIND = LONG_REAL)
      PHI_OBS = OBS_LONG * Pi / 180.0_8
    ENDIF

    COSTHETA_LOC = COS(  OBS_COLAT * Pi / 180.0_8 )
    SINTHETA_LOC = SIN(  OBS_COLAT * Pi / 180.0_8 )
    LOCATIONS(1) = REAL( COSTHETA_LOC, KIND = EXTRA_LONG_REAL)
    
    ALLOCATE(  TEMP_THETA_TRANSFORM_DALF(1,1: NUM_LEGENDRE_SV) )    
    ALLOCATE(  TEMP_THETA_TRANSFORM_ALF(1,1: NUM_LEGENDRE_SV)  )
    CALL GET_LEGENDRE_FUNCTIONS( LOCATIONS, LMAX_SV, LMAX_SV, TEMP_THETA_TRANSFORM_ALF, TEMP_THETA_TRANSFORM_DALF)
    ALF_LOCATION = REAL(  TEMP_THETA_TRANSFORM_ALF, KIND = 8)
    DALF_LOCATION = REAL( TEMP_THETA_TRANSFORM_DALF, KIND = 8)
    DEALLOCATE( TEMP_THETA_TRANSFORM_DALF, TEMP_THETA_TRANSFORM_ALF )


    ! Loop over each component of u and transform
    DO HARMONIC_TYPE = 1, 2 ! Looks like it's toroidal first and poloidal afters
    DO HARMONIC = 1, LMAX_U * (LMAX_U + 2)

      ! Evaluate u_h and its horizontal divergence on (theta, phi) grid
      CALL EVALUATE_U_H_SINGLE_HARMONIC( NTHETA_GRID, NPHI_GRID, HARMONIC, HARMONIC_TYPE,  U_H, DIV_H_U_H, ALF, DALF, HARMONICS, LMAX_SV, ONE_DIV_SIN_THETA, COSTHETA_GRID  )

      ! Transform to spherical harmonic decomposition
      DO I_THETA = 1, NTHETA_GRID
      DO I_PHI = 0, NPHI_GRID-1
	B_R_DOT(I_THETA, I_PHI) = - U_H(I_THETA, I_PHI,1) * GRAD_H_B_R(I_THETA, I_PHI,1) - U_H(I_THETA, I_PHI,2) * GRAD_H_B_R(I_THETA, I_PHI,2) - B_R(I_THETA, I_PHI) * DIV_H_U_H(I_THETA, I_PHI)
      ENDDO
      ENDDO

      CALL REAL_2_SPEC( B_R_DOT, SPEC, LMAX_SV, LMAX_B_OBS + HARMONICS(HARMONIC)%L, HARMONICS, NPHI_GRID, NTHETA_GRID , LEGENDRE_INV )
      ! Extrapolate to observation point and evaluate the entry of matrix G.
    
	B_DOT_LOCATION(1:3) = 0.0_EXTRA_LONG_REAL
	B_LOCATION(1:3) = 0.0_EXTRA_LONG_REAL

      DO I = 1, LMAX_SV * (LMAX_SV + 2) ! 
	! evaluate B_r  and del_H^2 B_r at Earth surface
	IF( HARMONICS(I)%SINCOS .EQ. SINE_HARMONIC ) THEN
	  PHI_DEP = SIN(  HARMONICS(I)%M * PHI_OBS )
	  DERV_PHI_DEP = HARMONICS(I)%M * COS(  HARMONICS(I)%M * PHI_OBS )
	ELSE
	  PHI_DEP = COS(  HARMONICS(I)%M * PHI_OBS )
	  DERV_PHI_DEP = -HARMONICS(I)%M * SIN(  HARMONICS(I)%M * PHI_OBS )
	ENDIF

	INDEX_PLM = HARMONICS(I)%M * LMAX_SV +  (HARMONICS(I)%M * (3-HARMONICS(I)%M) ) /2 + 1 + (HARMONICS(I)%L - HARMONICS(I)%M)
	B_DOT_LOCATION(1) = B_DOT_LOCATION(1) + SPEC(I) * (CMB_RADIUS / ES_RADIUS )**(HARMONICS(I)%L+2) * ALF_LOCATION(1, INDEX_PLM) * PHI_DEP    !/ (HARMONICS(I)%L + 1). No extra factors as they cancel out.
	B_DOT_LOCATION(2) = B_DOT_LOCATION(2) - SPEC(I) * (CMB_RADIUS / ES_RADIUS )**(HARMONICS(I)%L+2) * DALF_LOCATION(1, INDEX_PLM) * PHI_DEP  / (HARMONICS(I)%L + 1)
	B_DOT_LOCATION(3) = B_DOT_LOCATION(3) - SPEC(I) * (CMB_RADIUS / ES_RADIUS )**(HARMONICS(I)%L+2) * ALF_LOCATION(1, INDEX_PLM) * DERV_PHI_DEP / SINTHETA_LOC / (HARMONICS(I)%L + 1)

	IF( HARMONICS(I)%L .LE. LMAX_B_OBS ) THEN !SV goes to higher L than B[obs], so exclude if index goes higher than required.
	  B_LOCATION(1) = B_LOCATION(1) + GAUSS(I) *  ALF_LOCATION(1, INDEX_PLM) * PHI_DEP * REAL(HARMONICS(I)%L + 1, KIND = LONG_REAL)
	  B_LOCATION(2) = B_LOCATION(2) - GAUSS(I) * DALF_LOCATION(1, INDEX_PLM) * PHI_DEP
	  B_LOCATION(3) = B_LOCATION(3) - GAUSS(I) * ALF_LOCATION(1, INDEX_PLM) * DERV_PHI_DEP / SINTHETA_LOC
	ENDIF
      ENDDO ! End evaluation of B_r
      ! ***********************************************
      ! ***********************************************
      ! Change the  definition of G
      ! ***********************************************
      ! ***********************************************

      ! if we want to optimize for the rate of change of declination:  
      D_LOC = atan2(B_LOCATION(3),-B_LOCATION(2))
      !    G( HARMONIC, HARMONIC_TYPE ) = (B_DOT_LOCATION(2) * B_LOCATION(3) - B_DOT_LOCATION(3) * B_LOCATION(2)) / ((B_LOCATION(2) * B_LOCATION(2)) + (B_LOCATION(3) * B_LOCATION(3)))

      F_LOC = SQRT( SUM(B_LOCATION**2) )
      NORM_G = (GAUSS(1)**2 + GAUSS(2)**2 + GAUSS(3)**2) * SQRT( GAUSS(2)**2 + GAUSS(3)**2 )

	  IF( MINIM_FLAG.EQ.0) THEN
      ! if we want to optimize for the rate of change of dipole tilt:      
	    G( HARMONIC, HARMONIC_TYPE ) = 0.5 * (CMB_RADIUS / ES_RADIUS)**3 * (SPEC(2) * GAUSS(1) * GAUSS(2) + SPEC(3) * GAUSS(1) * GAUSS(3) &
	    - SPEC(1) * (GAUSS(2)**2 +GAUSS(3)**2) )  / NORM_G
	  ELSE IF( (MINIM_FLAG.EQ.1) .OR. (MINIM_FLAG.EQ.2) ) THEN
      ! if we want to optimize for the rate of change of intensity:
	    G( HARMONIC, HARMONIC_TYPE ) = SUM( B_DOT_LOCATION(:) * B_LOCATION(:) ) / F_LOC
	  ELSE IF( (MINIM_FLAG.EQ.3) ) THEN
      ! if we want to optimize for the rate of change of declination at LOCATION:
	    G( HARMONIC, HARMONIC_TYPE ) = (B_DOT_LOCATION(2) * B_LOCATION(3) - B_DOT_LOCATION(3) * B_LOCATION(2)) / ((B_LOCATION(2) * B_LOCATION(2)) + (B_LOCATION(3) * B_LOCATION(3)))	
	  ELSE IF( (MINIM_FLAG.EQ.4) ) THEN
	  ! if we want to optimise the rate of change of g10
	    G( HARMONIC, HARMONIC_TYPE ) = 0.5 * (CMB_RADIUS / ES_RADIUS)**3 * SPEC(1) !* (3 / 4.0_LONG_REAL / Pi)
	  ELSE IF( (MINIM_FLAG.EQ.5) ) THEN
	  ! if we want to optimise the rate of change of inclination at LOCATION
	  !  G( HARMONIC, HARMONIC_TYPE ) = ( B_DOT_LOCATION(2) * B_LOCATION(2) * B_LOCATION(1) + B_DOT_LOCATION(3) * B_LOCATION(3) * B_LOCATION(1) &
	  !  - B_DOT_LOCATION(1) * (B_LOCATION(2)**2 + B_LOCATION(3)**2)  ) &
	  !  / ( F_LOC**2 * SQRT( B_LOCATION(2)**2 + B_LOCATION(3)**2 ) )
	  G( HARMONIC, HARMONIC_TYPE ) = ( B_DOT_LOCATION(2) * B_LOCATION(2) * B_LOCATION(1) + B_DOT_LOCATION(3) * B_LOCATION(3) * B_LOCATION(1) &
	    - B_DOT_LOCATION(1) * (B_LOCATION(2)**2 + B_LOCATION(3)**2)  ) &
	    / ( F_LOC**2 * SQRT( B_LOCATION(2)**2 + B_LOCATION(3)**2 ) )
	  ELSE IF( (MINIM_FLAG.EQ.6) ) THEN
	  ! OPTIMIZE G11
	  G( HARMONIC, HARMONIC_TYPE ) = 0.5 * (CMB_RADIUS / ES_RADIUS)**3 * SPEC(2) !* (3 / 4.0_LONG_REAL / Pi)
	  ENDIF
    ENDDO !end loop over harmonics
    ENDDO !and toroidal / poloidal types.
    
!    PRINT*,'MAXVAL(B_R_DOT)',MAXVAL(B_R_DOT)! same
!    PRINT*,'MAXVAL(SPEC)',MAXVAL(SPEC)
!    PRINT*,'G( 3, 1 )',G( 3, 1 )
!    PRINT*,'MAXVAL(U_H)',MAXVAL(U_H)! same
!    PRINT*,'B_DOT_LOCATION',B_DOT_LOCATION(1),B_DOT_LOCATION(2),B_DOT_LOCATION(3)
!    PRINT*,'B_LOCATION',B_LOCATION(1),B_LOCATION(2),B_LOCATION(3)! same
!    PRINT*,'NTHETA_GRID',NTHETA_GRID! same
!    PRINT*,'LEGENDRE_INV(1:3,1)',LEGENDRE_INV(1:3,1)
!    PRINT*,'SIZEOF(LEGENDRE_INV)',SIZEOF(LEGENDRE_INV)
!    PRINT*,'ALF(1:NTHETA_GRID,10)',ALF(1:NTHETA_GRID,10)
!    PRINT*,'SIZEOF(ALF)',SIZEOF(ALF)
!    PRINT*,'LEGENDRE_GAUW(1:NTHETA_GRID)',LEGENDRE_GAUW(1:NTHETA_GRID)
!    PRINT*,'SIZEOF(LEGENDRE_GAUW)',SIZEOF(LEGENDRE_GAUW)
!    PRINT*,'COSTHETA_GRID',COSTHETA_GRID ! same


    IF( RESTRICTION .eq. 1) THEN  !restrict to poloidal only
      G(:, TOROIDAL_FLOW ) = 0.0_LONG_REAL
    ELSE IF( RESTRICTION .eq. 2) THEN  !restrict to toroidal only
      G(:, POLOIDAL_FLOW ) = 0.0_LONG_REAL
    ENDIF

    IF( RESTRICTION .eq. 3) THEN  !restrict to COLUMNAR FLOW
      DO HARMONIC_TYPE = 1, 2
      DO HARMONIC = 1, LMAX_U * (LMAX_U + 2)
	IF ( HARMONIC_TYPE .eq. TOROIDAL_FLOW ) THEN ! coefficients are zero for l-m even
	  IF ( MOD(HARMONICS(HARMONIC)%L - HARMONICS(HARMONIC)%M,2) .eq. 0 ) G(HARMONIC, TOROIDAL_FLOW ) = 0.0_LONG_REAL
	ELSE IF ( HARMONIC_TYPE .eq. POLOIDAL_FLOW ) THEN ! zero for l-m odd
	  IF ( MOD(HARMONICS(HARMONIC)%L - HARMONICS(HARMONIC)%M,2) .eq. 1 ) G(HARMONIC, POLOIDAL_FLOW ) = 0.0_LONG_REAL
	ENDIF
      ENDDO
      ENDDO
    ENDIF


    KE_TOR = 0.0_LONG_REAL
    KE_POL = 0.0_LONG_REAL
    ! EVEC is the q in Phil's 2014 paper, and NORM_FACTOR is the diagonal matrix E element 
    DO I=1, LMAX_U * (LMAX_U + 2)
      NORM_FACTOR = HARMONICS(I)%L * (HARMONICS(I)%L + 1) * 4.0_LONG_REAL * Pi / REAL( 2 * HARMONICS(I)%L + 1, KIND = LONG_REAL )
      EVEC_TOR(I) = G(I, TOROIDAL_FLOW) / NORM_FACTOR
      EVEC_POL(I) = G(I, POLOIDAL_FLOW) / NORM_FACTOR
    ENDDO

    ! The flow (defined by the coefficients EVEC) is given in SI units of m/s.
    ! normalise: flow coeffs should describe a flow in km/yr
    EVEC_TOR(:) = EVEC_TOR(:) * 3600.0_LONG_REAL * 24.0_LONG_REAL * 365.0_LONG_REAL / 1000.0_LONG_REAL
    EVEC_POL(:) = EVEC_POL(:) * 3600.0_LONG_REAL * 24.0_LONG_REAL * 365.0_LONG_REAL / 1000.0_LONG_REAL

    DO I=1, LMAX_U * (LMAX_U + 2)
      NORM_FACTOR = HARMONICS(I)%L * (HARMONICS(I)%L + 1) * 4.0_LONG_REAL * Pi / REAL( 2 * HARMONICS(I)%L + 1, KIND = LONG_REAL )
      KE_TOR = KE_TOR + EVEC_TOR(I)**2 * NORM_FACTOR
      KE_POL = KE_POL + EVEC_POL(I)**2 * NORM_FACTOR
    ENDDO
    KE = KE_TOR + KE_POL

    ! KE is \int u^2 dS over the surface of a unit sphere. Convert to dimensional rms (in km/yr)
    TOTAL_RMS = SQRT( KE  / (4.0_LONG_REAL * Pi) )
    
    ! Scale solution by the target rms value in km/yr

    EVEC_TOR = EVEC_TOR * TARGET_RMS / TOTAL_RMS
    EVEC_POL = EVEC_POL * TARGET_RMS / TOTAL_RMS
    
    ! THE FOLLOWING IS UNNEDED, BUT I COPIED FROM THE OPTIMISATION FILE AND IT'S CONSISTENT WITH THE REST
    SPEC_U_TOR = EVEC_TOR  / (CMB_RADIUS/1000.0_LONG_REAL)
    SPEC_U_POL = EVEC_POL  / (CMB_RADIUS/1000.0_LONG_REAL)
    
    ! SKIPPING THE OUTPUT PART FROM THE OPTIMISATION SCRIPT
          
    ! *********************************************************************** !
    SPEC_U_TOR(:) = (-1)**REVERSE * SPEC_U_TOR(:) * CMB_RADIUS / 1000.0_LONG_REAL
    SPEC_U_POL(:) = (-1)**REVERSE * SPEC_U_POL(:) * CMB_RADIUS / 1000.0_LONG_REAL
    ! so spec is evec
  ENDIF ! end optimisation block

  !! CALCULATE INDUCTION EQUATION COEFFICIENTS !!
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
  NONLIN_BR_CMB(:,:) = NONLIN_BR_CMB(:,:) * 1000.0_LONG_REAL
  ! Calculate SH coefficients of the nonlinear term
  SPEC(:) = 0.0_8
  CALL REAL_2_SPEC( NONLIN_BR_CMB, SPEC, LMAX_SV, LMAX_SV, HARMONICS, NPHI_GRID, NTHETA_GRID , LEGENDRE_INV )
  ! Transform to Gauss-type coefficients:
  DO J=1, LMAX_SV * (LMAX_SV + 2)
    SPEC(J) = SPEC(J) * (CMB_RADIUS / ES_RADIUS )**(HARMONICS(J)%L+2) / REAL( HARMONICS(J)%L + 1, KIND = LONG_REAL)
  ENDDO

  !! OUTPUT CURRENT COEFFICIENTS !!
  IF ( (N_IT.EQ.1) .OR. (t0 + (N_IT-1)*dt.EQ.tf) .OR. (MOD(N_IT,OUTPUT_RATE).EQ.0) ) THEN
    PRINT*, 'Write Gauss coefficients'
    WRITE(12,'(F11.4)',advance="no")  t0 + (N_IT-1)*dt
    WRITE(12,FORM_B)  (GAUSS(I), I=1,LMAX_B_OBS * (LMAX_B_OBS + 2))

    PRINT*, 'Write SV coefficients'
    WRITE(13,'(F11.4)',advance="no")  t0 + (N_IT-1)*dt
    WRITE(13,FORM_SV)  (SPEC(I), I=1,LMAX_SV * (LMAX_SV + 2))

    IF ( FLAG_U_INIT.EQ.1 ) THEN  
      PRINT*, 'Write U coefficients'
      WRITE(14,'(F11.4)',advance="no")  t0 + (N_IT-1)*dt
      WRITE(14,FORM_U)  (SPEC_U_TOR(I)/ (CMB_RADIUS/1000.0_LONG_REAL), I=1, LMAX_U * (LMAX_U + 2))
      WRITE(15,'(F11.4)',advance="no")  t0 + (N_IT-1)*dt
      WRITE(15,FORM_U)  (SPEC_U_POL(I)/ (CMB_RADIUS/1000.0_LONG_REAL), I=1, LMAX_U * (LMAX_U + 2))
    ENDIF  
  ENDIF
  
  !! TIMESTEPPING !!
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


  PRINT*, 'update for the next timestep'
  SPEC_P = SPEC
  GAUSS = GAUSS_N
  N_IT = N_IT +1

ENDDO ! end time loop



! OUTPUTS and WRAP UP

CLOSE(12)
CLOSE(13)
CLOSE(14)
CLOSE(15)

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

       WRITE(6,'(A,x,x,F9.4)') 'RMS FLOW IS (KM/YR)', SQRT( SUM(FLOW_L_SPEC) / (4.0_8 * Pi) ) 
       WRITE(6,'(A,x,x,ES11.4)') 'RMS FLOW IS (M/S)', SQRT( SUM(FLOW_L_SPEC) / (4.0_8 * Pi) ) * 1000 / (3600 * 24 * 365)




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
     
END PROGRAM TIMESTEPPING_INDUCTION_OPT