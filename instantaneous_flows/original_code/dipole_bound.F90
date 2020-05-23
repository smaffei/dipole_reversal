PROGRAM DIPOLE_BOUND
USE SUBS
IMPLICIT NONE

! Computes flows that optimise the value of d g_1^0/dt
! Reads in a file of Gauss coefficients defining the magnetic field (units nT).

! Outputs a file of Gauss coefficients for the flow on the unit sphere, units yr^(-1).

TYPE (HARMONIC_STRUCTURE), ALLOCATABLE, DIMENSION(:) :: HARMONICS

REAL( KIND = 8) ::  THETA, PHI, PHI_OBS, VALUE(2)
INTEGER :: L, M, SINCOS, POLTOR, LMAX_B_OBS, LMAX_U, IOS, LMAX_FILE, NP, i, HARMONIC, LMAX_SV, NTHETA_GRID, NPHI_GRID, I_PHI, I_THETA, NUM_LEGENDRE_SV, INDEX_PLM, J, HARMONIC_TYPE, RESTRICTION, INDEX
CHARACTER(300) :: FILENAME, JUNK
CHARACTER(1) :: STRING_HARM_TYPE2


REAL( KIND = 8), ALLOCATABLE, DIMENSION(:) :: TOR_U, POL_U, POL_B_DOT, GAUSS, LEGENDRE_GAUW, COSTHETA_GRID
REAL( KIND = EXTRA_LONG_REAL), ALLOCATABLE, DIMENSION(:,:)  ::  TEMP_THETA_TRANSFORM_DALF, TEMP_THETA_TRANSFORM_ALF
REAL( KIND = 8), ALLOCATABLE :: GRAD_H_B_R(:,:,:), B_R(:,:), DIV_H_U_H(:,:), U_H(:,:,:),ALF_LOCATION(:,:), DALF_LOCATION(:,:), &
                 ALF(:,:), DALF(:,:), B_R_DOT(:,:), SPEC(:), LEGENDRE_INV(:,:), ONE_DIV_SIN_THETA(:), EVEC_TOR(:), EVEC_POL(:), FLOW_L_SPEC(:), G(:,:)
REAL( KIND = 8) :: COSTHETA_LOC, SINTHETA_LOC, PHI_DEP, DERV_PHI_DEP, B_DOT_LOCATION(3), F_LOC, MAX_DIPOLE_GROWTH_RATE, NORM_FACTOR, KE, TOTAL_RMS, SCALE_FACTOR, TARGET_RMS, KE_TOR, KE_POL


PRINT*, 'ENTER MAX DEGREE FOR U'
READ*, LMAX_U

PRINT*, 'ENTER MAX DEGREE FOR B_OBS'
READ*, LMAX_B_OBS

PRINT*, 'ENTER FILENAME FOR B_OBS: FORMAT SHOULD BE SINGLE LINE OF HEADER THEN GLM/HLM IN TABULAR FORMAT'
READ*, FILENAME

PRINT*, 'ENTER TARGET RMS VALUE FOR FLOW'
READ*, TARGET_RMS

PRINT*, 'RESTRICTION OF FLOW? 0: NONE, 1: POL ONLY, 2: TOR ONLY'
READ*, RESTRICTION

! Read in Gauss coefficients for field.

 OPEN(11, FILE = FILENAME, STATUS = 'OLD', FORM = 'FORMATTED', & 
                             IOSTAT = IOS, ACTION = 'READ')
    IF( IOS .NE. 0) THEN
    PRINT*, 'ERROR IN OPENING FILE ', FILENAME
    STOP
    ENDIF
   

     NP = LMAX_B_OBS*(LMAX_B_OBS +2)
     ALLOCATE( GAUSS(1:NP)  )
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
      WRITE(6,'(A,I3,A,I3)') 'LMAX REQUIRED: ', LMAX_B_OBS,' LMAX IN FILE: ',L
      WRITE(6,'(A)')  'FIELD MODEL IS PADDED WITH ZEROS.'
      WRITE(6,'(A)') '********************'
      ENDIF
  


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
      
! Get the Associated Legendre Functions (ALF) at the theta grid points, with their theta derivatives (DALF), of the ALF's up to degree LMAX_SV
   CALL GET_LEGENDRE_FUNCTIONS( REAL(COSTHETA_GRID, KIND = EXTRA_LONG_REAL), LMAX_SV, LMAX_SV, TEMP_THETA_TRANSFORM_ALF, TEMP_THETA_TRANSFORM_DALF) 
   ALF = REAL(  TEMP_THETA_TRANSFORM_ALF, KIND = 8)
   DALF = REAL( TEMP_THETA_TRANSFORM_DALF, KIND = 8)
   ONE_DIV_SIN_THETA(:) = 1.0_LONG_REAL / SQRT( 1.0_LONG_REAL - COSTHETA_GRID(:)**2 )

   DO I = 1, NUM_LEGENDRE_SV
   LEGENDRE_INV(1:NTHETA_GRID,I) = ALF(1:NTHETA_GRID,I) * LEGENDRE_GAUW(1:NTHETA_GRID) 
   ENDDO
   
    ALLOCATE( B_R_DOT(NTHETA_GRID, 0:NPHI_GRID-1), SPEC(1:LMAX_SV * (LMAX_SV + 2)) )
    ALLOCATE( B_R(1:NTHETA_GRID, 0:NPHI_GRID-1), GRAD_H_B_R(1:NTHETA_GRID,0:NPHI_GRID-1,2) )

    CALL EVALUATE_B_R_GRID( NTHETA_GRID, NPHI_GRID, GAUSS, HARMONICS, LMAX_B_OBS, LMAX_SV, B_R, GRAD_H_B_R, ALF, DALF, ONE_DIV_SIN_THETA )
    ALLOCATE( U_H(1:NTHETA_GRID, 0:NPHI_GRID-1,2), DIV_H_U_H(1:NTHETA_GRID, 0:NPHI_GRID-1) )

 
! Loop over each component of u and transform
     ALLOCATE( G(1: LMAX_U * (LMAX_U + 2),1:2 ) )
     DO HARMONIC_TYPE = 1, 2 
     DO HARMONIC = 1, LMAX_U * (LMAX_U + 2)

! Evaluate u_h and its horizontal divergence on (theta, phi) grid
   CALL EVALUATE_U_H_SINGLE_HARMONIC( NTHETA_GRID, NPHI_GRID, HARMONIC, HARMONIC_TYPE,  U_H, DIV_H_U_H, ALF, DALF, HARMONICS, LMAX_SV, ONE_DIV_SIN_THETA, COSTHETA_GRID  )    

! Transform to spherical harmonic decomposition
   DO I_THETA = 1, NTHETA_GRID
   DO I_PHI = 0, NPHI_GRID-1
   B_R_DOT(I_THETA, I_PHI) = - U_H(I_THETA, I_PHI,1) * GRAD_H_B_R(I_THETA, I_PHI,1) - U_H(I_THETA, I_PHI,2) * GRAD_H_B_R(I_THETA, I_PHI,2) - B_R(I_THETA, I_PHI) * DIV_H_U_H(I_THETA, I_PHI)
   ENDDO
   ENDDO

   !B_R_DOT = B_R
   CALL REAL_2_SPEC( B_R_DOT, SPEC, LMAX_SV, LMAX_B_OBS + HARMONICS(HARMONIC)%L, HARMONICS, NPHI_GRID, NTHETA_GRID , LEGENDRE_INV ) 

! The array SPEC contains the spherical harmonic transform of B_R. 
! B_r^dot  = \sum SPEC(i) Y_i = \sum (l+1) (a/r)^(l+2) g_l^m Y_l^m + h_l^m Y_l^m, so the coefficient g^dot[1^0] is then SPEC(1)/( 2 * (a/r)^3 )

    G( HARMONIC, HARMONIC_TYPE ) = 0.5_8 * SPEC(1) * (CMB_RADIUS / ES_RADIUS)**3
  !  PRINT*, G(HARMONIC, HARMONIC_TYPE)

     ENDDO !end loop over harmonics
     ENDDO !and toroidal / poloidal types.

     IF( RESTRICTION .eq. 1) THEN  !restrict to poloidal only
     G(:, TOROIDAL_FLOW ) = 0.0_LONG_REAL
     ELSE IF( RESTRICTION .eq. 2) THEN  !restrict to toroidal only
     G(:, POLOIDAL_FLOW ) = 0.0_LONG_REAL
     ENDIF

  
     ALLOCATE( EVEC_TOR(1: LMAX_U * (LMAX_U + 2) ), EVEC_POL(1: LMAX_U * (LMAX_U + 2) ) )
     KE_TOR = 0.0_LONG_REAL
     KE_POL = 0.0_LONG_REAL

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

! Units here:   flow is given in km/yr, B in nT, B_dot in nT/yr.
! Therefore, since we want B_dot to be computed using a flow given in m/yr (since radius r is measured in metres in the induction equation), multiply by 1000 to convert km -> m.

      MAX_DIPOLE_GROWTH_RATE = SUM(G(:, TOROIDAL_FLOW) * EVEC_TOR(:)) * 1000.0_LONG_REAL + SUM(G(:, POLOIDAL_FLOW) * EVEC_POL(:)) * 1000.0_LONG_REAL

      PRINT*, 'MAX RATE OF CHANGE OF G10 IS (nT/yr) ', MAX_DIPOLE_GROWTH_RATE
      WRITE(6,'(A,F8.2,A,F8.2,A)') 'FLOW IS ',KE_POL/KE * 100,' % POLOIDAL AND ',KE_TOR/KE * 100, '% TOROIDAL IN KE'
! write to disk, and convert gauss coeffs describing a flow in km/year to one describing a flow in yr^(-1) on the unit sphere.
       OPEN(22, FILE = 'OPTIMAL_FLOW.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')
       WRITE(22,'(A,ES10.2,A,I3,A,I3)') 'OPTIMAL FLOW WITH TARGET RMS ', TARGET_RMS, '     LMAX_U ', LMAX_U, ' LMAX_B', LMAX_B_OBS
       DO J=1, LMAX_U * (LMAX_U + 2)
       IF( HARMONICS(J)%M .NE. 0 .AND. HARMONICS(J)%SINCOS .EQ. COSINE_HARMONIC)  THEN
       WRITE(22,'(I3,x,I3,4(x,ES15.5,x))') HARMONICS(J)%L, HARMONICS(J)%M, EVEC_POL(J)  / (CMB_RADIUS/1000.0_LONG_REAL), EVEC_POL(J+1)  / (CMB_RADIUS/1000.0_LONG_REAL), &
                                           EVEC_TOR(J)  / (CMB_RADIUS/1000.0_LONG_REAL), EVEC_TOR(J+1)  / (CMB_RADIUS/1000.0_LONG_REAL)  
       ENDIF

       IF( HARMONICS(J)%M .EQ. 0 .AND. HARMONICS(J)%SINCOS .EQ. COSINE_HARMONIC)  THEN
       WRITE(22,'(I3,x,I3,4(x,ES15.5,x))') HARMONICS(J)%L, HARMONICS(J)%M, EVEC_POL(J)/ (CMB_RADIUS/1000.0_LONG_REAL), 0.0, EVEC_TOR(J)  / (CMB_RADIUS/1000.0_LONG_REAL), 0.0_LONG_REAL 
       ENDIF
       ENDDO

       CLOSE(22)



! Calculate spectrum of flow
       OPEN(23, FILE = 'FLOW_L_SPEC.DAT', STATUS = 'REPLACE', FORM = 'FORMATTED')
       ALLOCATE( FLOW_L_SPEC(1:LMAX_U) )
       FLOW_L_SPEC(:) = 0.0_LONG_REAL
       DO I=1, LMAX_U * (LMAX_U + 2)
       NORM_FACTOR = HARMONICS(I)%L * (HARMONICS(I)%L + 1) * 4.0_LONG_REAL * Pi / REAL( 2 * HARMONICS(I)%L + 1, KIND = LONG_REAL )
       FLOW_L_SPEC(HARMONICS(I)%L) = FLOW_L_SPEC(HARMONICS(I)%L) + ( EVEC_POL(I)**2 + EVEC_TOR(I)**2 ) * NORM_FACTOR
       ENDDO

       WRITE(23,*) 'KE SPECTRUM PER DEGREE L IN (KM/YR)^2: \int u^2 dOmega / (4 Pi)'
       DO I = 1, LMAX_U
       WRITE(23,*) I,  FLOW_L_SPEC(I)/(4.0_LONG_REAL * Pi)
       ENDDO
       CLOSE(23)
       WRITE(6,'(A,x,x,F8.3)') 'RMS FLOW IS (KM/YR)', SQRT( SUM(FLOW_L_SPEC) / (4.0_8 * Pi) ) 


! Write flow for Python
CALL WRITE_U_CENTRED( GMT_THETA, GMT_PHI, EVEC_TOR, EVEC_POL, HARMONICS, LMAX_SV)
CALL WRITE_U_RANDOM( EVEC_TOR, EVEC_POL, HARMONICS, LMAX_SV)

! Calculate the intensity change on the Earth's surface and output in a file that GMT can read
        
   CALL EVALUATE_U_H( NTHETA_GRID, NPHI_GRID, EVEC_TOR, EVEC_POL, U_H, DIV_H_U_H, ALF, DALF, HARMONICS, LMAX_U, LMAX_SV, ONE_DIV_SIN_THETA  )    
   CALL EVALUATE_B_R_GRID( NTHETA_GRID, NPHI_GRID, GAUSS, HARMONICS, LMAX_B_OBS, LMAX_SV, B_R, GRAD_H_B_R, ALF, DALF, ONE_DIV_SIN_THETA )

! Transform 
   DO I_THETA = 1, NTHETA_GRID
   DO I_PHI = 0, NPHI_GRID-1
   B_R_DOT(I_THETA, I_PHI) = - U_H(I_THETA, I_PHI,1) * GRAD_H_B_R(I_THETA, I_PHI,1) - U_H(I_THETA, I_PHI,2) * GRAD_H_B_R(I_THETA, I_PHI,2) - B_R(I_THETA, I_PHI) * DIV_H_U_H(I_THETA, I_PHI)
   ENDDO
   ENDDO
 
! B_R_DOT is produced from a flow in units of km/yr. Need to multiply by 1000 to get a magnetic field in units of nT/yr.
     B_R_DOT(:,:) = B_R_DOT(:,:) * 1000.0_LONG_REAL


   CALL REAL_2_SPEC( B_R_DOT, SPEC, LMAX_SV, LMAX_SV, HARMONICS, NPHI_GRID, NTHETA_GRID, LEGENDRE_INV ) 
  
! Write F and Fdot for output to GMT
! Convert SPEC coefficents to a Gauss-coefficient representation.
       DO J=1, LMAX_SV * (LMAX_SV + 2)
       SPEC(J) = SPEC(J) * (CMB_RADIUS / ES_RADIUS )**(HARMONICS(J)%L+2) / REAL( HARMONICS(J)%L + 1, KIND = LONG_REAL)
       ENDDO

   CALL WRITE_BR_BR_DOT (SPEC, GAUSS, LMAX_B_OBS, LMAX_SV, HARMONICS )


       CALL DFFTW_DESTROY_PLAN(PLAN_FFT_R2HC)
       DEALLOCATE(FFT_TRANSFORM_ARRAY)    
STOP
    


 
END PROGRAM DIPOLE_BOUND
