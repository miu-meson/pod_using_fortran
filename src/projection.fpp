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
#:for FIELD in ['X', 'Y', 'DX', 'DY']
  real (kind=DP), allocatable, dimension (:) :: ${FIELD}$
#:endfor
#:for FIELD in ['DV', 'T', 'U', 'V', 'YREF']
  real (kind=DP), allocatable, dimension (:, :) :: ${FIELD}$
#:endfor
#:for I in range(4)
  real (kind=DP), allocatable, dimension (:, :) :: PC${I}$
#:endfor
#:for FIELD in ['T', 'U', 'V', 'Y']
  real (kind=DP), allocatable, dimension (:, :) :: EVECTORS_${FIELD}$
#:endfor
  real (kind=DP), allocatable, dimension (:, :) :: ENERGIES
  real (kind=DP), allocatable, dimension (:, :) :: RENERGIES
#:for I in range(3,15,3)
  real (kind=DP), allocatable, dimension (:, :) :: RENERGIES${I}$
#:endfor
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
#:for FIELD in ['T', 'U', 'V', 'YREF']
  real (kind=SP), allocatable, dimension (:, :) :: ${FIELD}$_FLOAT
#:endfor


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

#:for FIELD in ['X', 'Y', 'DX', 'DY']
  allocate (${FIELD}$(M))
#:endfor

#:for FIELD in ['T', 'U', 'V', 'DV', 'YREF']
  allocate (${FIELD}$(M,M))
#:endfor

#:for FIELD in ['T', 'U', 'V', 'Y']
  allocate (EVECTORS_${FIELD}$(M*M,NMODES))
#:endfor
#:for I in range(4)
  allocate (PC${I}$(M*M,NMODES))
#:endfor

  FORMULATION = 'F'
  allocate (ENERGIES(7,NSETS*NSNAPS))
  allocate (RENERGIES(7,NSETS*NSNAPS))
#:for I in range(3,15,3)
  allocate (RENERGIES${I}$(7,NSETS*NSNAPS))
#:endfor

#:for FIELD in ['T', 'U', 'V', 'Y']
  open (unit=9, file=PREFIX//'-eigenvectors-'//FORMULATION//'-${FIELD}$-'//SUFFIX// &
  '.asc', action='read')
  do I = 1, M*M
    read (9, *) J, EVECTORS_${FIELD}$(I, 1:NMODES)
  end do
  close (unit=9)

#:endfor

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
#:for I in range(3,15,3)
  RENERGIES${I}$(:, :) = 0.0_DP
#:endfor

  COUNTER = 0
  do ISETS = 1, NSETS
    do ISNAPS = FIRST(ISETS), LAST(ISETS), NSKIP
      write (SWEEP, '(I6.6)') MIN(ISNAPS, 960)

#:for FIELD in ['T', 'U', 'V']
      write (I_NAME, '(A160)') LTRIM(DIRLIST(ISETS)) // 'field-${FIELD}$-' // SWEEP // &
        '.bin'
      call INPUT_MATRIX(LTRIM(I_NAME), N_FLOAT, ${FIELD}$_FLOAT, X_FLOAT, &
        Y_FLOAT, IOS, 0)

#:endfor

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
#:for I, FIELD in enumerate(['T', 'U', 'V'])
        ENERGIES(${I+1}$, COUNTER) = ENERGIES(${I+1}$, COUNTER) + DV(R, S)*(${FIELD}$_FLOAT(R,S)**2.0_DP)
#:endfor
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
#:for I, FIELD in enumerate(['T', 'U', 'V'])
          PC${I+1}$(COUNTER, MODE) = PC${I+1}$(COUNTER, MODE) + &
            DV(R, S)*EVECTORS_${FIELD}$((R-1)*M+S, MODE)*${FIELD}$_FLOAT(R, S)

#:endfor
        end do
      end do

      PC0(COUNTER, :) = A*PC1(COUNTER, :) + B*PC2(COUNTER, :) + &
        B*PC3(COUNTER, :)

#:for J in range(3,15,3)
! **************************** ${J}$ MODES ****************************************
#:for I, FIELD in enumerate(['T', 'U', 'V', 'YREF'])
      ${FIELD}$(:, :) = 0.0_DP
#:endfor
      do MODE = 1, MIN(${J}$, NMODES)
        do concurrent ( R = 1:M, S = 1:M )
#:for I, FIELD in enumerate(['T', 'U', 'V'])
          ${FIELD}$(R, S) = ${FIELD}$(R, S) + PC0(COUNTER, MODE)*EVECTORS_${FIELD}$((R-1)*M+S, MODE)
#:endfor
          YREF(R, S) = YREF(R, S) + PC0(COUNTER, MODE)*EVECTORS_Y((R-1)*M+S, &
              MODE)
        end do
      end do

      do concurrent ( R = 1:M, S = 1:M )
#:for I, FIELD in enumerate(['T', 'U', 'V'])
        RENERGIES${J}$(${I+1}$, COUNTER) = RENERGIES${J}$(${I+1}$, COUNTER) + DV(R, S)*(${FIELD}$(R,S)**2.0_DP)
#:endfor
        RENERGIES${J}$(4, COUNTER) = RENERGIES${J}$(4, COUNTER) - DV(R, S)*Y(S)*T(R, S)
        RENERGIES${J}$(5, COUNTER) = RENERGIES${J}$(5, COUNTER) + DV(R, S)*T(R, S)
        RENERGIES${J}$(6, COUNTER) = RENERGIES${J}$(6, COUNTER) - DV(R, S)*YREF(R, S)*T(R, S)
        RENERGIES${J}$(7, COUNTER) = RENERGIES${J}$(7, COUNTER) + DV(R, S)*(&
          (A*T(R,S))**2.0_DP+(B*U(R,S))**2.0_DP+(B*V(R,S))**2.0_DP)
      end do
#:endfor
    end do
  end do

  open (unit=11, file='proj-comparison-'//FORMULATION//'.asc', action='write')
#:for J in range(4)
  open (unit=12+${J}$, file='proj-coefficients${J}$-'//FORMULATION//'.asc', action='write')
#:endfor
#:for J in range(3,15,3)
  open (unit=16+${J}$, file='proj-recomparison${J}$-'//FORMULATION//'.asc', action='write')
#:endfor

  do R = 1, NSETS*NSNAPS
    write (11, '(I6,8(E15.7))') R, ENERGIES(1:7, R)
#:for J in range(4)
    write (12+${J}$, '(I6,1X,14(e15.7,1X))') R, PC${J}$(R, 1:MIN(12,NMODES))
#:endfor
#:for J in range(3,15,3)
    write (16+${J}$, '(I6,8(E15.7))') R, RENERGIES${J}$(1:7, R)
#:endfor
  end do

#:for J in range(11,16)
  close (${J}$)
#:endfor
#:for J in range(3,15,3)
  close (16+${J}$)
#:endfor

end program
