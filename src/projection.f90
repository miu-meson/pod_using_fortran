program PROJECT_POD_SIROVICH
!!   This routine takes a set of POD modes computed using extract_pod_sirovich
!! 	and performs a projection of primitive fields to obtain the corresponding
!! 	modal coefficients. The routines also computes different global quantities
!!		from the primitive fields and using low order reconstructions using 3, 6,
!!		9 and up to 12 POD modes is available.
!!
!!   Since the snapshots do not contain any information on the ghost cells,
!!   boundary conditions must be applied by hand. By default, conditions
!!   correspond to the RBNT case: iso-thermal top and bottom plates, adiabatic
!!   side-walls and no-slip boundary conditions everywhere.
!!
!!   If the code is compile with a _STRESSFREE directive, boundary conditions
!!   correspond to the RBST case, which is like the one above, but a stress-free
!!   boundary condition is imposed on the top boundary
!!
!!   Input files: The same input file used in `extract_pod_sirovich`
!!
!!		Output files:
!!		`proj-coefficients0-F.asc` 	The modal amplitudes
!!		`proj-comparison-F.asc` 		The global quantities from primitive fields
!!		`proj-recomparisonN-F.asc` 	The global quantities from reconstructions
!!
!!		For an example, see `make clean`  `make try`  `make proj`
  use :: DECLARATIONS
  use :: INPUT
  use :: OUTPUT

  implicit none

  integer :: I, J, R, S, M, COUNTER, MODE, NMODES, IOS, K
  real (kind=DP) :: N, A, B
  real (kind=DP), allocatable, dimension (:) :: X
  real (kind=DP), allocatable, dimension (:) :: Y
  real (kind=DP), allocatable, dimension (:) :: DX
  real (kind=DP), allocatable, dimension (:) :: DY
  real (kind=DP), allocatable, dimension (:, :) :: DV
  real (kind=DP), allocatable, dimension (:, :) :: T
  real (kind=DP), allocatable, dimension (:, :) :: U
  real (kind=DP), allocatable, dimension (:, :) :: V
  real (kind=DP), allocatable, dimension (:, :) :: YREF
  real (kind=DP), allocatable, dimension (:, :) :: PC0
  real (kind=DP), allocatable, dimension (:, :) :: PC1
  real (kind=DP), allocatable, dimension (:, :) :: PC2
  real (kind=DP), allocatable, dimension (:, :) :: PC3
  real (kind=DP), allocatable, dimension (:, :) :: EVECTORS_T
  real (kind=DP), allocatable, dimension (:, :) :: EVECTORS_U
  real (kind=DP), allocatable, dimension (:, :) :: EVECTORS_V
  real (kind=DP), allocatable, dimension (:, :) :: EVECTORS_Y
  real (kind=DP), allocatable, dimension (:, :) :: ENERGIES
  real (kind=DP), allocatable, dimension (:, :) :: RENERGIES
  real (kind=DP), allocatable, dimension (:, :) :: RENERGIES3
  real (kind=DP), allocatable, dimension (:, :) :: RENERGIES6
  real (kind=DP), allocatable, dimension (:, :) :: RENERGIES9
  real (kind=DP), allocatable, dimension (:, :) :: RENERGIES12
  real (kind=DP), allocatable, dimension (:, :, :) :: ORTHO
  character (len=1) :: FORMULATION

  integer :: ISETS, NSETS, ISNAPS, NSNAPS, NSKIP
  integer, allocatable, dimension (:) :: FIRST, LAST
  character (len=6) :: SWEEP, SUFFIX
  character (len=5) :: PREFIX
  character (len=300), allocatable, dimension (:) :: DIRLIST
  character (len=400) :: I_NAME, O_NAME

  real (kind=SP) :: N_FLOAT
  real (kind=SP), allocatable, dimension (:) :: X_FLOAT
  real (kind=SP), allocatable, dimension (:) :: Y_FLOAT
  real (kind=SP), allocatable, dimension (:, :) :: T_FLOAT
  real (kind=SP), allocatable, dimension (:, :) :: U_FLOAT
  real (kind=SP), allocatable, dimension (:, :) :: V_FLOAT
  real (kind=SP), allocatable, dimension (:, :) :: YREF_FLOAT


#define LTRIM(x) trim(adjustl(x))

! A. Read Inputs :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  read (5, '(A)') PREFIX
  read (5, '(A)') SUFFIX
  read (5, *) A, B
  read (5, *) NMODES, NMODES
  read (5, *) NSETS
  read (5, *) NSNAPS
  read (5, *) NSKIP
  allocate (DIRLIST(NSETS), FIRST(NSETS), LAST(NSETS))
  do ISETS = 1, NSETS
    read (5, *) DIRLIST(ISETS)
    FIRST(:) = 1
    LAST(:) = NSNAPS*NSKIP
    print *, ' POD projection of the complete series'
  end do


  print *, ' Read Existing POD Modes'
  M = 514
  N = M - 2

  allocate (X(M))
  allocate (Y(M))
  allocate (DX(M))
  allocate (DY(M))

  allocate (T(M,M))
  allocate (U(M,M))
  allocate (V(M,M))
  allocate (DV(M,M))
  allocate (YREF(M,M))

  allocate (EVECTORS_T(M*M,NMODES))
  allocate (EVECTORS_U(M*M,NMODES))
  allocate (EVECTORS_V(M*M,NMODES))
  allocate (EVECTORS_Y(M*M,NMODES))
  allocate (PC0(M*M,NMODES))
  allocate (PC1(M*M,NMODES))
  allocate (PC2(M*M,NMODES))
  allocate (PC3(M*M,NMODES))

  FORMULATION = 'F'
  allocate (ENERGIES(7,NSETS*NSNAPS))
  allocate (RENERGIES(7,NSETS*NSNAPS))
  allocate (RENERGIES3(7,NSETS*NSNAPS))
  allocate (RENERGIES6(7,NSETS*NSNAPS))
  allocate (RENERGIES9(7,NSETS*NSNAPS))
  allocate (RENERGIES12(7,NSETS*NSNAPS))

  open (unit=9, file=PREFIX//'-eigenvectors-'//FORMULATION//'-T-'//SUFFIX// &
  '.asc', action='read')
  do I = 1, M*M
    read (9, *) J, EVECTORS_T(I, 1:NMODES)
  end do
  close (unit=9)

  open (unit=9, file=PREFIX//'-eigenvectors-'//FORMULATION//'-U-'//SUFFIX// &
  '.asc', action='read')
  do I = 1, M*M
    read (9, *) J, EVECTORS_U(I, 1:NMODES)
  end do
  close (unit=9)

  open (unit=9, file=PREFIX//'-eigenvectors-'//FORMULATION//'-V-'//SUFFIX// &
  '.asc', action='read')
  do I = 1, M*M
    read (9, *) J, EVECTORS_V(I, 1:NMODES)
  end do
  close (unit=9)

  open (unit=9, file=PREFIX//'-eigenvectors-'//FORMULATION//'-Y-'//SUFFIX// &
  '.asc', action='read')
  do I = 1, M*M
    read (9, *) J, EVECTORS_Y(I, 1:NMODES)
  end do
  close (unit=9)


  print *, ' Define metrics...'
  DX(:) = 1.0_DP/N
  DY(:) = 1.0_DP/N
  do R = 1, M
    Y(R) = (R-1/2)*1.0_DP/N - 0.5_DP
    X(R) = (R-1/2)*1.0_DP/N - 0.5_DP
    do S = 1, M
      DV(R, S) = DX(R)*DY(S)
    end do
  end do

  allocate (ORTHO(NMODES,NMODES,3))
  print *, ' Verify Orthogonality...'
  ORTHO(:, :, :) = 0.
  do concurrent ( I = 1:MIN(12, NMODES), J = 1:MIN(12, NMODES) )
    do R = 1, M
      do S = 1, M
        ORTHO(I, J, 1) = ORTHO(I, J, 1) + DV(R, S)*&
            (A*EVECTORS_T((R-1)*M+S,I))*(A*EVECTORS_T((R-1)*M+S,J))
        ORTHO(I, J, 2) = ORTHO(I, J, 2) + DV(R, S)*&
            (B*EVECTORS_U((R-1)*M+S,I))*(B*EVECTORS_U((R-1)*M+S,J))
        ORTHO(I, J, 3) = ORTHO(I, J, 3) + DV(R, S)*&
            (B*EVECTORS_V((R-1)*M+S,I))*(B*EVECTORS_V((R-1)*M+S,J))
      end do
    end do
  end do

  do R = 1, 3
    print *, ' '
    print *, ' Variable No.', R
    do I = 1, NMODES
      write (*, '(12(f12.7,1X))') ORTHO(I, :, R)
    end do
  end do

  print *, ' '
  print *, ' Combined Orthogonality ... (2+3) '
  do I = 1, NMODES
    write (*, '(12(f12.7,1X))') ORTHO(I, :, 2) + ORTHO(I, :, 3)
  end do

  print *, ' '
  print *, ' Combined Orthogonality ... '
  do I = 1, NMODES
    write (*, '(12(f12.7,1X))') ORTHO(I, :, 1) + ORTHO(I, :, 2) + ORTHO(I, :, &
      3)
  end do


  ENERGIES(:, :) = 0.0_DP
  RENERGIES(:, :) = 0.0_DP
  RENERGIES3(:, :) = 0.0_DP
  RENERGIES6(:, :) = 0.0_DP
  RENERGIES9(:, :) = 0.0_DP
  RENERGIES12(:, :) = 0.0_DP

  COUNTER = 0
  do ISETS = 1, NSETS
    do ISNAPS = FIRST(ISETS), LAST(ISETS), NSKIP
      write (SWEEP, '(I6.6)') MIN(ISNAPS, 960)

      write (I_NAME, '(A160)') LTRIM(DIRLIST(ISETS)) // 'field-T-' // SWEEP // &
        '.bin'
      call INPUT_MATRIX(LTRIM(I_NAME), N_FLOAT, T_FLOAT, X_FLOAT, &
        Y_FLOAT, IOS, 0)

      write (I_NAME, '(A160)') LTRIM(DIRLIST(ISETS)) // 'field-U-' // SWEEP // &
        '.bin'
      call INPUT_MATRIX(LTRIM(I_NAME), N_FLOAT, U_FLOAT, X_FLOAT, &
        Y_FLOAT, IOS, 0)

      write (I_NAME, '(A160)') LTRIM(DIRLIST(ISETS)) // 'field-V-' // SWEEP // &
        '.bin'
      call INPUT_MATRIX(LTRIM(I_NAME), N_FLOAT, V_FLOAT, X_FLOAT, &
        Y_FLOAT, IOS, 0)


      call LORENZ_STATE(T_FLOAT, 1.0, M, YREF_FLOAT)
#ifndef _STRESSFREE
      call BOUNDARY_CONDITIONS_GHOST('N', M, T_FLOAT, U_FLOAT, V_FLOAT, &
        YREF_FLOAT)
#elif
      call BOUNDARY_CONDITIONS_GHOST('S', M, T_FLOAT, U_FLOAT, V_FLOAT, &
        YREF_FLOAT)
#endif

      T_FLOAT(:, :) = T_FLOAT(:, :) - 0.5
      YREF_FLOAT(:, :) = YREF_FLOAT(:, :) - 0.5

      COUNTER = COUNTER + 1
      do concurrent ( R = 1:M, S = 1:M )
        ENERGIES(1, COUNTER) = ENERGIES(1, COUNTER) + DV(R, S)*(T_FLOAT(R,S)**2.0_DP)
        ENERGIES(2, COUNTER) = ENERGIES(2, COUNTER) + DV(R, S)*(U_FLOAT(R,S)**2.0_DP)
        ENERGIES(3, COUNTER) = ENERGIES(3, COUNTER) + DV(R, S)*(V_FLOAT(R,S)**2.0_DP)
        ENERGIES(4, COUNTER) = ENERGIES(4, COUNTER) - DV(R, S)*Y(S)*T_FLOAT(R, S)
        ENERGIES(5, COUNTER) = ENERGIES(5, COUNTER) + DV(R, S)*T_FLOAT(R, S)
        ENERGIES(6, COUNTER) = ENERGIES(6, COUNTER) - DV(R, S)*YREF_FLOAT(R, S)*T_FLOAT(R, S)
        ENERGIES(7, COUNTER) = ENERGIES(7, COUNTER) + DV(R, S)*( &
          (A*T_FLOAT(R,S))**2.0_DP + &
          (B*U_FLOAT(R,S))**2.0_DP + &
          (B*V_FLOAT(R,S))**2.0_DP)
      end do

      do MODE = 1, MIN(12, NMODES)
        do concurrent ( R = 1:M, S = 1:M )
          PC1(COUNTER, MODE) = PC1(COUNTER, MODE) + &
            DV(R, S)*EVECTORS_T((R-1)*M+S, MODE)*T_FLOAT(R, S)

          PC2(COUNTER, MODE) = PC2(COUNTER, MODE) + &
            DV(R, S)*EVECTORS_U((R-1)*M+S, MODE)*U_FLOAT(R, S)

          PC3(COUNTER, MODE) = PC3(COUNTER, MODE) + &
            DV(R, S)*EVECTORS_V((R-1)*M+S, MODE)*V_FLOAT(R, S)

        end do
      end do

      PC0(COUNTER, :) = A*PC1(COUNTER, :) + B*PC2(COUNTER, :) + &
        B*PC3(COUNTER, :)

! **************************** 3 MODES ****************************************
      T(:, :) = 0.0_DP
      U(:, :) = 0.0_DP
      V(:, :) = 0.0_DP
      YREF(:, :) = 0.0_DP
      do MODE = 1, MIN(3, NMODES)
        do concurrent ( R = 1:M, S = 1:M )
          T(R, S) = T(R, S) + PC0(COUNTER, MODE)*EVECTORS_T((R-1)*M+S, MODE)
          U(R, S) = U(R, S) + PC0(COUNTER, MODE)*EVECTORS_U((R-1)*M+S, MODE)
          V(R, S) = V(R, S) + PC0(COUNTER, MODE)*EVECTORS_V((R-1)*M+S, MODE)
          YREF(R, S) = YREF(R, S) + PC0(COUNTER, MODE)*EVECTORS_Y((R-1)*M+S, &
              MODE)
        end do
      end do

      do concurrent ( R = 1:M, S = 1:M )
        RENERGIES3(1, COUNTER) = RENERGIES3(1, COUNTER) + DV(R, S)*(T(R,S)**2.0_DP)
        RENERGIES3(2, COUNTER) = RENERGIES3(2, COUNTER) + DV(R, S)*(U(R,S)**2.0_DP)
        RENERGIES3(3, COUNTER) = RENERGIES3(3, COUNTER) + DV(R, S)*(V(R,S)**2.0_DP)
        RENERGIES3(4, COUNTER) = RENERGIES3(4, COUNTER) - DV(R, S)*Y(S)*T(R, S)
        RENERGIES3(5, COUNTER) = RENERGIES3(5, COUNTER) + DV(R, S)*T(R, S)
        RENERGIES3(6, COUNTER) = RENERGIES3(6, COUNTER) - DV(R, S)*YREF(R, S)*T(R, S)
        RENERGIES3(7, COUNTER) = RENERGIES3(7, COUNTER) + DV(R, S)*(&
          (A*T(R,S))**2.0_DP+(B*U(R,S))**2.0_DP+(B*V(R,S))**2.0_DP)
      end do
! **************************** 6 MODES ****************************************
      T(:, :) = 0.0_DP
      U(:, :) = 0.0_DP
      V(:, :) = 0.0_DP
      YREF(:, :) = 0.0_DP
      do MODE = 1, MIN(6, NMODES)
        do concurrent ( R = 1:M, S = 1:M )
          T(R, S) = T(R, S) + PC0(COUNTER, MODE)*EVECTORS_T((R-1)*M+S, MODE)
          U(R, S) = U(R, S) + PC0(COUNTER, MODE)*EVECTORS_U((R-1)*M+S, MODE)
          V(R, S) = V(R, S) + PC0(COUNTER, MODE)*EVECTORS_V((R-1)*M+S, MODE)
          YREF(R, S) = YREF(R, S) + PC0(COUNTER, MODE)*EVECTORS_Y((R-1)*M+S, &
              MODE)
        end do
      end do

      do concurrent ( R = 1:M, S = 1:M )
        RENERGIES6(1, COUNTER) = RENERGIES6(1, COUNTER) + DV(R, S)*(T(R,S)**2.0_DP)
        RENERGIES6(2, COUNTER) = RENERGIES6(2, COUNTER) + DV(R, S)*(U(R,S)**2.0_DP)
        RENERGIES6(3, COUNTER) = RENERGIES6(3, COUNTER) + DV(R, S)*(V(R,S)**2.0_DP)
        RENERGIES6(4, COUNTER) = RENERGIES6(4, COUNTER) - DV(R, S)*Y(S)*T(R, S)
        RENERGIES6(5, COUNTER) = RENERGIES6(5, COUNTER) + DV(R, S)*T(R, S)
        RENERGIES6(6, COUNTER) = RENERGIES6(6, COUNTER) - DV(R, S)*YREF(R, S)*T(R, S)
        RENERGIES6(7, COUNTER) = RENERGIES6(7, COUNTER) + DV(R, S)*(&
          (A*T(R,S))**2.0_DP+(B*U(R,S))**2.0_DP+(B*V(R,S))**2.0_DP)
      end do
! **************************** 9 MODES ****************************************
      T(:, :) = 0.0_DP
      U(:, :) = 0.0_DP
      V(:, :) = 0.0_DP
      YREF(:, :) = 0.0_DP
      do MODE = 1, MIN(9, NMODES)
        do concurrent ( R = 1:M, S = 1:M )
          T(R, S) = T(R, S) + PC0(COUNTER, MODE)*EVECTORS_T((R-1)*M+S, MODE)
          U(R, S) = U(R, S) + PC0(COUNTER, MODE)*EVECTORS_U((R-1)*M+S, MODE)
          V(R, S) = V(R, S) + PC0(COUNTER, MODE)*EVECTORS_V((R-1)*M+S, MODE)
          YREF(R, S) = YREF(R, S) + PC0(COUNTER, MODE)*EVECTORS_Y((R-1)*M+S, &
              MODE)
        end do
      end do

      do concurrent ( R = 1:M, S = 1:M )
        RENERGIES9(1, COUNTER) = RENERGIES9(1, COUNTER) + DV(R, S)*(T(R,S)**2.0_DP)
        RENERGIES9(2, COUNTER) = RENERGIES9(2, COUNTER) + DV(R, S)*(U(R,S)**2.0_DP)
        RENERGIES9(3, COUNTER) = RENERGIES9(3, COUNTER) + DV(R, S)*(V(R,S)**2.0_DP)
        RENERGIES9(4, COUNTER) = RENERGIES9(4, COUNTER) - DV(R, S)*Y(S)*T(R, S)
        RENERGIES9(5, COUNTER) = RENERGIES9(5, COUNTER) + DV(R, S)*T(R, S)
        RENERGIES9(6, COUNTER) = RENERGIES9(6, COUNTER) - DV(R, S)*YREF(R, S)*T(R, S)
        RENERGIES9(7, COUNTER) = RENERGIES9(7, COUNTER) + DV(R, S)*(&
          (A*T(R,S))**2.0_DP+(B*U(R,S))**2.0_DP+(B*V(R,S))**2.0_DP)
      end do
! **************************** 12 MODES ****************************************
      T(:, :) = 0.0_DP
      U(:, :) = 0.0_DP
      V(:, :) = 0.0_DP
      YREF(:, :) = 0.0_DP
      do MODE = 1, MIN(12, NMODES)
        do concurrent ( R = 1:M, S = 1:M )
          T(R, S) = T(R, S) + PC0(COUNTER, MODE)*EVECTORS_T((R-1)*M+S, MODE)
          U(R, S) = U(R, S) + PC0(COUNTER, MODE)*EVECTORS_U((R-1)*M+S, MODE)
          V(R, S) = V(R, S) + PC0(COUNTER, MODE)*EVECTORS_V((R-1)*M+S, MODE)
          YREF(R, S) = YREF(R, S) + PC0(COUNTER, MODE)*EVECTORS_Y((R-1)*M+S, &
              MODE)
        end do
      end do

      do concurrent ( R = 1:M, S = 1:M )
        RENERGIES12(1, COUNTER) = RENERGIES12(1, COUNTER) + DV(R, S)*(T(R,S)**2.0_DP)
        RENERGIES12(2, COUNTER) = RENERGIES12(2, COUNTER) + DV(R, S)*(U(R,S)**2.0_DP)
        RENERGIES12(3, COUNTER) = RENERGIES12(3, COUNTER) + DV(R, S)*(V(R,S)**2.0_DP)
        RENERGIES12(4, COUNTER) = RENERGIES12(4, COUNTER) - DV(R, S)*Y(S)*T(R, S)
        RENERGIES12(5, COUNTER) = RENERGIES12(5, COUNTER) + DV(R, S)*T(R, S)
        RENERGIES12(6, COUNTER) = RENERGIES12(6, COUNTER) - DV(R, S)*YREF(R, S)*T(R, S)
        RENERGIES12(7, COUNTER) = RENERGIES12(7, COUNTER) + DV(R, S)*(&
          (A*T(R,S))**2.0_DP+(B*U(R,S))**2.0_DP+(B*V(R,S))**2.0_DP)
      end do
    end do
  end do

  open (unit=11, file='proj-comparison-'//FORMULATION//'.asc', action='write')
  open (unit=12+0, file='proj-coefficients0-'//FORMULATION//'.asc', action='write')
  open (unit=12+1, file='proj-coefficients1-'//FORMULATION//'.asc', action='write')
  open (unit=12+2, file='proj-coefficients2-'//FORMULATION//'.asc', action='write')
  open (unit=12+3, file='proj-coefficients3-'//FORMULATION//'.asc', action='write')
  open (unit=16+3, file='proj-recomparison3-'//FORMULATION//'.asc', action='write')
  open (unit=16+6, file='proj-recomparison6-'//FORMULATION//'.asc', action='write')
  open (unit=16+9, file='proj-recomparison9-'//FORMULATION//'.asc', action='write')
  open (unit=16+12, file='proj-recomparison12-'//FORMULATION//'.asc', action='write')

  do R = 1, NSETS*NSNAPS
    write (11, '(I6,8(E15.7))') R, ENERGIES(1:7, R)
    write (12+0, '(I6,1X,14(e15.7,1X))') R, PC0(R, 1:MIN(12,NMODES))
    write (12+1, '(I6,1X,14(e15.7,1X))') R, PC1(R, 1:MIN(12,NMODES))
    write (12+2, '(I6,1X,14(e15.7,1X))') R, PC2(R, 1:MIN(12,NMODES))
    write (12+3, '(I6,1X,14(e15.7,1X))') R, PC3(R, 1:MIN(12,NMODES))
    write (16+3, '(I6,8(E15.7))') R, RENERGIES3(1:7, R)
    write (16+6, '(I6,8(E15.7))') R, RENERGIES6(1:7, R)
    write (16+9, '(I6,8(E15.7))') R, RENERGIES9(1:7, R)
    write (16+12, '(I6,8(E15.7))') R, RENERGIES12(1:7, R)
  end do

  close (11)
  close (12)
  close (13)
  close (14)
  close (15)
  close (16+3)
  close (16+6)
  close (16+9)
  close (16+12)

end program
