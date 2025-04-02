program EXTRACT_POD_SIROVICH
!!   This routine takes a sequence of snapshots from 2-D DNS in Basilisk
!!   and extract the POD modes using Sirovich's method and a joint
!!   velocity-temperature formulation with scale factors a and b to be defined
!!   in the input file.
!!   If a=1, b=0 we obtain the temperature formulation
!!   If a=0, b=1 we obtain the velocity
!!   If a=b=1 we obtain the formulation used in the article
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
!!   Input files: The input file contais the following
!!   A prefix (6 characters)
!!   A suffix (6 characters)
!!   Scale parameters a and b
!!   The number of POD modes to compute and to write
!!   Number of total sets (our data may be contained inside different folders)
!!   Number of snapshots per set
!!   Skip every n snapshots
!!   and the Path to snapshots (one per set)
!!   See the example input.dat file inside the tst folder
!!   Also see "make try"
!!
!!   Output files:
!!   prefix-coefficient-F-suffix.asc     - Modal amplitudes
!!   prefix-eigenvalue-F-suffix.asc      - Eigenvalues
!!   prefix-eigenvector-F-suffix.asc     - Spatial eigenfunctions
!!   prefix-eigenvectors-F-T-suffix.asc  - Educed Modes (for projection)
!!   prefix-T-000001-suffix.bin          - Educed Modes (for visualisation)

  use :: DECLARATIONS
  use :: INPUT
  use :: OUTPUT
  use :: PIERCE

  implicit none

  integer :: I, J, R, S, M, COUNTER, MSNAPS, MODE, ICOVCOR
  integer :: ISETS, ISNAPS, IOS, IMINNXNT, IMAXNXNT, JASCENDING, JDESCENDING
  integer :: NSETS, NSNAPS, NSKIP
  integer, allocatable, dimension (:) :: FIRST, LAST
  character (len=300), allocatable, dimension (:) :: DIRLIST
  character (len=6) :: SWEEP, SUFFIX
  character (len=5) :: PREFIX
  character (len=400) :: I_NAME, O_NAME

  logical :: DOSWITCHED
  integer :: ORDEROFS, NMODES, NMODES_SAVE

  real (kind=SP) :: N

#:for I in ['X', 'Y', 'DX', 'DY']
  real (kind=SP), allocatable, dimension (:) :: ${I}$
#:endfor

#:for I in ['T', 'U', 'V', 'YREF', 'PSI', 'DV', 'ENERGIES']
  real (kind=SP), allocatable, dimension (:, :) :: ${I}$
#:endfor

  real (kind=DP) :: SCALE_FACTOR, A, B
  real (kind=DP), dimension (4) :: MEAN_ENERGIES

#:for I in ['SPACKED', 'EVALS', 'VOLUME', 'HEIGHT']
  real (kind=DP), allocatable, dimension (:) :: ${I}$
#:endfor

#:for I in ['TEMPERATURE', 'XVELOCITY', 'YVELOCITY', 'REFERENCE_HEIGHT', 'ORTHO', 'EVECS', 'PRINPC']
  real (kind=DP), allocatable, dimension (:, :) :: ${I}$
#:endfor

#:for I in ['TENTPC', 'EVECTORS']
#:for J in ['T', 'U', 'V', 'YREF']
  real (kind=DP), allocatable, dimension (:, :) :: ${I}$_${J}$
#:endfor
#:endfor

  real (kind=DP), allocatable, dimension (:, :, :) :: DATASET


! Preamble. Define Macros
#define LTRIM(x) trim(adjustl(x))

! Default values
  A = 1.0_DP
  B = 1.0_DP

! A. Read Inputs :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  read (5, '(A)') PREFIX
  read (5, '(A)') SUFFIX
  read (5, *) A, B
  read (5, *) NMODES, NMODES_SAVE
  read (5, *) NSETS
  read (5, *) NSNAPS
  read (5, *) NSKIP
  allocate (DIRLIST(NSETS), FIRST(NSETS), LAST(NSETS))


  do ISETS = 1, NSETS
    read (5, *) DIRLIST(ISETS)
#ifdef _REVERSAL
    read (5, *) FIRST(ISETS)
    read (5, *) LAST(ISETS)
    print *, ' POD decomposition of reversals only '
#else
    FIRST(:) = 1
    LAST(:) = NSNAPS*NSKIP
    print *, ' POD decomposition of the complete series'
#endif
  end do

! A.1. Obtain time-averaged fields ---------------------------------------------
  do ISETS = 1, NSETS
    print ('(1X,A,3X,I8,3X,A)'), ' Set No.: ', ISETS, DIRLIST(ISETS)
    write (I_NAME, '(A160)') LTRIM(DIRLIST(ISETS)) // 'field-T-000001.bin'
    call INPUT_MATRIX(I_NAME, N, T, X, Y, IOS, 0)
    M = INT(N) + 2
  end do

  allocate (DX(M), DY(M), DV(M,M))
  DX(:) = 1.0_DP/N
  DY(:) = 1.0_DP/N
  do R = 1, M
    Y(R) = (R-1/2)*1.d0/N - 0.5
    do S = 1, M
      DV(R, S) = DX(R)*DY(S)
    end do
  end do

  print *, ' POD decomposition using complete fields '
  print *, ' End of part A.1 '


! A.2. Load Snapshots --------------------------------------------------------
#ifdef _REVERSAL
  COUNTER = 0
  do concurrent ( ISETS = 1:NSETS, ISNAPS = FIRST(ISETS):LAST(ISETS):NSKIP)
    COUNTER = COUNTER + 1
  end do
  print *, ' POD decomposition of reversals only '
  print *, ' Total number of Snapshots ', COUNTER
  do ISETS = 1, NSETS
    print *, ' Reversal No. ', ISETS, ' between ', FIRST(ISETS), ' - ', &
      LAST(ISETS)
  end do
#else
  COUNTER = NSETS*NSNAPS
  print *, ' POD decomposition complete fields '
  print *, ' Total number of Snapshots ', COUNTER
#endif
  MSNAPS = COUNTER

  IMINNXNT = MIN(M*M, MSNAPS)
  IMAXNXNT = MAX(M*M, MSNAPS)


  allocate (ENERGIES(4,MSNAPS))
  allocate (VOLUME(M*M))
  allocate (HEIGHT(M*M))
#:for FIELD in ['TEMPERATURE', 'XVELOCITY', 'YVELOCITY', 'REFERENCE_HEIGHT']
  allocate (${FIELD}$(M*M,MSNAPS))
#:endfor

  COUNTER = 0
#:for FIELD in ['TEMPERATURE', 'XVELOCITY', 'YVELOCITY', 'ENERGIES']
  ${FIELD}$(:, :) = 0.0_DP
#:endfor
  MEAN_ENERGIES(:) = 0.0_DP

  open (unit=10, file='time-energies-'//SUFFIX//'.asc', action='write')

  do ISETS = 1, NSETS
    do ISNAPS = FIRST(ISETS), LAST(ISETS), NSKIP
      write (SWEEP, '(I6.6)') MIN(ISNAPS, 960)

#:for FIELD in ['T', 'U', 'V']
      write (I_NAME, '(A160)') LTRIM(DIRLIST(ISETS)) // 'field-${FIELD}$-' // SWEEP // '.bin'
      call INPUT_MATRIX(I_NAME, N, ${FIELD}$, X, Y, IOS, 0)

#:endfor

      call LORENZ_STATE(T, 1.0, M, YREF)
#ifndef _STRESSFREE
      call BOUNDARY_CONDITIONS_GHOST('N', M, T, U, V, YREF)
#elif
      call BOUNDARY_CONDITIONS_GHOST('S', M, T, U, V, YREF)
#endif

      T(:, :) = T(:, :) - 0.5
      YREF(:, :) = YREF(:, :) - 0.5

      COUNTER = COUNTER + 1
      do concurrent ( R = 1:M, S = 1:M )
        HEIGHT((R-1)*M+S) = Y(S)
        VOLUME((R-1)*M+S) = DV(R, S)

#:for SHORT, LONG in zip(['T', 'U', 'V', 'YREF'], ['TEMPERATURE', 'XVELOCITY', 'YVELOCITY', 'REFERENCE_HEIGHT'])
        ${LONG}$((R-1)*M+S, COUNTER) = ${SHORT}$(R, S)
#:endfor

#:for I, FIELD in enumerate(['T', 'U', 'V'])
        ENERGIES(${I+1}$, COUNTER) = ENERGIES(${I+1}$, COUNTER) + DV(R, S)*(${FIELD}$(R,S)**2.0_DP)
#:endfor
        ENERGIES(4, COUNTER) = ENERGIES(4, COUNTER) - DV(R, S)*(Y(S)*T(R,S))
      end do
      write (10, *) COUNTER, ENERGIES(1:4, COUNTER)
    end do
  end do
  close (unit=10)

  print *, ' Finished Loading Snapshots '
  print *, ' End of part A.2 '

! A.3. Evaluate Total energies for posterior verification

#:for I in range(4)
  MEAN_ENERGIES(${I+1}$) = SUM(ENERGIES(${I+1}$, :))/REAL(MSNAPS)
#:endfor

  print *, ' Total energies in Sampled Snapshots'
  print *, ' Temperature : ', MEAN_ENERGIES(1)
  print *, ' U-Velocity  : ', MEAN_ENERGIES(2)
  print *, ' V-Velocity  : ', MEAN_ENERGIES(3)
  print *, ' Kinetic     : ', SUM(MEAN_ENERGIES(2:3))
  print *, ' Total (u+T) : ', SUM(MEAN_ENERGIES(1:3))
  print *, ' Potential   : ', MEAN_ENERGIES(4)
  print *, ' Total (K+T) : ', SUM(MEAN_ENERGIES(2:4))
  print *, ' End of part A.3 '
  print *, ' End of part A '


! B. Assembly Matrix K :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!	------------------------------------------------------------------------------
!	Figure out whether we should do the 'regular' EOF or the one with switched X
!	and T axes.  The correlation matrix for the regular method is size (nx,nx) and
!	for the  switched method it is size (nt,nt); choose based on which of these is
!	smaller. See D. Pierce
!	------------------------------------------------------------------------------

  print *, ' B. Assembly Matrix K '
  print *, ' Joint velocity temperature formulation '
  print *, ' Psi(x,t)= (a*theta(x,t), b*u(x,t), b*v(x,t)) '
  print *, ' with a=', A, ' and b=', B

  print *, ' Figuring switched or not'
  DOSWITCHED = (M*M>MSNAPS)
  if (DOSWITCHED) then
    ORDEROFS = MSNAPS
    print *, ' deof: Working in switched mode'
  else
    ORDEROFS = M*M
    print *, ' deof: Working in unswitched mode'
  end if
  if (ORDEROFS>IMINNXNT) then
    write (0, *) ' Error!  EOF routine must be supplied '
    write (0, *) ' with enough workspace; passed parameter'
    write (0, *) ' iminnxnt must be at least ', ORDEROFS
    write (0, *) ' Passed value was iminnxnt=', IMINNXNT
    stop 'eof'
  end if
  if (NMODES>ORDEROFS) then
    write (0, *) ' Error! EOF routine called requesting more'
    write (0, *) ' modes than exist!  Request=', NMODES, ' exist=', ORDEROFS
    stop 'eof'
  end if

!	------------------------------------------------------------------------------
! Form the covariance or correlation matrix, put it  into 's'.  Note that 's' is
! always symmetric -- the correlation between X and Y is the same as  between Y
! and X -- so use packed storage for 's'. The packed storage scheme we use is
! the same as LAPACK uses so we can pass 's' directly to the solver routine: the
! matrix is lower triangular, and s(i,j) = spacked(i+(j-1)*(2*n-j)/2).
! See D. Pierce
!	------------------------------------------------------------------------------
  ICOVCOR = 0
  print *, ' Assembly Covariance Matrix: Joint Formulation'
  allocate (SPACKED(IMINNXNT*(IMINNXNT+1)/2), DATASET(M*M,MSNAPS,3))
  DATASET(:, :, 1) = A*TEMPERATURE(:, :)
  DATASET(:, :, 2) = B*XVELOCITY(:, :)
  DATASET(:, :, 3) = B*YVELOCITY(:, :)
  call DEOFCOVCOR(DATASET, ICOVCOR, 0, SPACKED, DOSWITCHED, VOLUME)
  deallocate (DATASET)
  print *, ' End of part B '

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! C. Solve Eigenvalue Problem ::::::::::::::::::::::::::::::::::::::::::::::::

  print *, ' C. Solve Eigenvalue Problem: Joint Formulation'
  allocate (EVALS(IMINNXNT), EVECS(IMINNXNT,NMODES))
  EVALS(:) = 0.0_DP
  EVECS(:, :) = 0.0_DP

  call SOLVE_EIGENVALUES(SPACKED, EVALS, EVECS, ORDEROFS, NMODES, IMINNXNT)
  deallocate (SPACKED)

  print *, ' Sum of ', NMODES, ' eigenvalues: ', SUM(EVALS(1:NMODES)), &
    SUM(EVALS(1:NMODES))/REAL(MSNAPS), SUM(MEAN_ENERGIES(1:3))
  print *, ' Sum of ', NMODES_SAVE, 'eigenvalues: ', &
    SUM(EVALS(NMODES_SAVE:NMODES)), SUM(EVALS(NMODES_SAVE:NMODES))/ &
    REAL(MSNAPS)

  open (unit=9, file=PREFIX//'-eigenvector-F-'//SUFFIX//'.asc', &
    action='write')
  open (unit=10, file=PREFIX//'-eigenvalue-F-'//SUFFIX//'.asc', &
    action='write')
  do J = 1, IMINNXNT
    write (9, *) J, EVECS(J, 1:NMODES)
  end do
  write (10, '(1(e14.7))') EVALS
  close (unit=9)
  close (unit=10)
  print *, ' End of part C '

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! D. Obtain Orthogonal Base ::::::::::::::::::::::::::::::::::::::::::::::::::
!	------------------------------------------------------------------------------
! Compute the tentative principal components; they are 'tentative' because they
! might be the PCs  of the switched data.  Put them into 'tentpcs', which is of
! size (nt,nmodes) if we are doing regular EOFs and (nx,nmodes) if we are doing
! switched EOFs.  These PCs come out in order corresponding to the order of
! 'evecs', which is in ASCENDING order.
!	------------------------------------------------------------------------------

  print *, ' D. Obtain educed modes '
  allocate (DATASET(M*M,MSNAPS,1))
#:for I, FIELD in enumerate(['T', 'U', 'V', 'YREF'])
  allocate (TENTPC_${FIELD}$(IMAXNXNT, NMODES))
#:endfor

#:for SHORT, LONG in zip(['T', 'U', 'V', 'YREF'], ['TEMPERATURE', 'XVELOCITY', 'YVELOCITY', 'REFERENCE_HEIGHT'])
  DATASET(:, :, 1) = ${LONG}$(:, :)
  call DEOFPCS(DATASET, EVECS, DOSWITCHED, TENTPC_${SHORT}$)

#:endfor
  deallocate (DATASET)
  print *, ' ... '
  print *, ' ... '
  print *, ' ... '
  print *, ' Completed educed modes '
  print *, ' End of part D '

!	----------------------------------------------------------------------------
!	Now we have all the pieces to assemble our final  result.  How we actually
!	assemble them depends on whether we are doing switched or unswitched EOFs
!	(except for the eigenVALUES, which are the same either way).
!	----------------------------------------------------------------------------
  allocate (PRINPC(MSNAPS,NMODES))
#:for I, FIELD in enumerate(['T', 'U', 'V', 'YREF'])
  allocate (EVECTORS_${FIELD}$(M*M,NMODES))
#:endfor

  if (DOSWITCHED) then
! ------------------------------------------------------------------------------
! In this case we must switch the principal components and the eigenvectors,
! applying the proper normalization.  First get the unswitched eigenvectors, which
! are the switched (tentative) principal components divided by the square root of
! the appropriate eigenvalue.  Recall that the LAPACK values are in ASCENDING
! order while we want them in DESCENDING order; do the switch in this loop.
! ------------------------------------------------------------------------------
    do JASCENDING = 1, NMODES
      JDESCENDING = NMODES - JASCENDING + 1
      SCALE_FACTOR = 1.d0/SQRT(EVALS(JASCENDING))
      do I = 1, M*M
#:for I, FIELD in enumerate(['T', 'U', 'V', 'YREF'])
        EVECTORS_${FIELD}$(I, JDESCENDING) = TENTPC_${FIELD}$(I, JASCENDING)* &
          SCALE_FACTOR

#:endfor
      end do
    end do
!	------------------------------------------------------------------------------
! Next get unswitched principal components, which are the switched eigenvectors
! multiplied by the appropriate eigenvalues.
!	------------------------------------------------------------------------------
    do JASCENDING = 1, NMODES
      JDESCENDING = NMODES - JASCENDING + 1
      SCALE_FACTOR = SQRT(EVALS(JASCENDING))
      do I = 1, MSNAPS
        PRINPC(I, JDESCENDING) = EVECS(I, JASCENDING)*SCALE_FACTOR
      end do
    end do
  end if

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! E. Normalize :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!	----------------------------------------------------------------------------
! Scale the eigenvectors to have a magnitude of 1; scale the corresponding
! principal components to reproduce the original data.
!	----------------------------------------------------------------------------

  print *, ' E. Normalization of the Eigenvectors and Eigenvalues '
  do MODE = 1, NMODES_SAVE
! --------------------------------------------------------------------------
! Get the normalization factor for each mode
! --------------------------------------------------------------------------
    SCALE_FACTOR = 0
    do I = 1, M*M
      SCALE_FACTOR = SCALE_FACTOR + VOLUME(I)*(&
        (A*EVECTORS_T(I, MODE))**2.0_DP + &
        (B*EVECTORS_U(I, MODE))**2.0_DP + &
        (B*EVECTORS_V(I, MODE))**2.0_DP)
    end do
    SCALE_FACTOR = SQRT(SCALE_FACTOR)
!	--------------------------------------------------------------------------
!	Normalize the eigenvectors
! --------------------------------------------------------------------------
    do I = 1, M*M
#:for I, FIELD in enumerate(['T', 'U', 'V', 'YREF'])
      EVECTORS_${FIELD}$(I, MODE) = EVECTORS_${FIELD}$(I, MODE)/SCALE_FACTOR
#:endfor
    end do
!	--------------------------------------------------------------------------
!	Normalize the principal components
! --------------------------------------------------------------------------
    do I = 1, MSNAPS
      PRINPC(I, MODE) = PRINPC(I, MODE)*SCALE_FACTOR
    end do
  end do


! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! F. Verify Orthogonality ::::::::::::::::::::::::::::::::::::::::::::::::::::
  allocate (ORTHO(NMODES_SAVE,NMODES_SAVE))
  ORTHO(:, :) = 0.
  do I = 1, NMODES_SAVE
    do J = 1, NMODES_SAVE
      do R = 1, M*M

        ORTHO(I, J) = ORTHO(I, J) + VOLUME(I)*(&
          (A*EVECTORS_T(R, I))*(A*EVECTORS_T(R,J)) + &
          (B*EVECTORS_U(R, I))*(B*EVECTORS_U(R,J)) + &
          (B*EVECTORS_V(R, I))*(B*EVECTORS_V(R,J)))

      end do
    end do
  end do
  print *, ' ... '
  print *, ' ... '
  print *, ' ... '
  print *, ' End of part E '


  print *, ' H. Verify Orthogonality'
  do I = 1, NMODES_SAVE
    write (*, '(12(e15.7,1X))') ORTHO(I, :)
  end do
  print *, ' End of part H '


! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! G. Wrap-up :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  print *, ' G. Wrap-up '

  print *, ' Save spatial eigenfunctions in ASCII format (double precision) '
#:for I, FIELD in enumerate(['T', 'U', 'V', 'YREF'])
  open (unit=9, file=PREFIX//'-eigenvectors-F-${FIELD}$-'//SUFFIX//'.asc', &
  action='write')
  do J = 1, M*M
    write (9, '(I6,14(e15.7))') J, EVECTORS_${FIELD}$(J, 1:MIN(NMODES,NMODES_SAVE))
  end do
  close (unit=9)

#:endfor

  print *, ' Save modal coefficients in ASCII format '
  open (unit=9, file=PREFIX//'-coefficient-F-'//SUFFIX//'.asc', &
    action='write')
  do J = 1, MSNAPS
    write (9, '(I6,14(e15.7))') J, PRINPC(J, 1:MIN(NMODES,NMODES_SAVE))
  end do
  close (unit=9)


  print *, ' Save spatial eigenfunctions in gnuplot compatible &
    &format (single precision) '
  allocate (PSI(M,M))
  do MODE = 1, MIN(NMODES, NMODES_SAVE)
    do concurrent ( R = 1:M, S = 1:M )
#:for I, FIELD in enumerate(['T', 'U', 'V'])
      ${FIELD}$(R, S) = EVECTORS_${FIELD}$((R-1)*M+S, MODE)
#:endfor
    end do

! **********************************************************************
! **********************************************************************
! Create pgm files from input fields for quick verification
    if (MODE==1) then
      call STREAMFUNCTION(U, V, M, PSI)
#:for I, FIELD in enumerate(['T', 'U', 'V', 'PSI'])
      call WRITEPPM2MATRIX(${FIELD}$, 'field-${FIELD}$')
#:endfor
    end if

#:for I, FIELD in enumerate(['T', 'U', 'V'])
    write (O_NAME, '(a,i6.6,a)') PREFIX // '-${FIELD}$-', MODE, &
      '-' // SUFFIX // '.bin'
    call OUTPUT_MATRIX(O_NAME, M, ${FIELD}$, X, Y)
#:endfor
  end do
  print *, ' End of part G '
  print *, ' End of Program '

end program
