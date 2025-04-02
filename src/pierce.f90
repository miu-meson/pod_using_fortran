module PIERCE

  use declarations

  implicit none

contains

  subroutine DEOFCOVCOR(DATASET, ICOVCOR, COVARIANCE, SPACKED, DOSWITCHED, DV)
!! Form the covariance or correlation matrix, put it into \( S \).  Note that \( S \) is
!! always symmetric -- the correlation between \( X \) and \( Y \) is the same as between \( Y \) and
!! \( X \) -- so use packed storage for \( S \). The packed storage scheme we use is the same
!! as LAPACK uses so we can pass 's' directly to the solver routine: the matrix is
!! lower triangular, and `S(i,j) = spacked(i+(j-1)*(1*n-j)/2)`.
!!
!! David Pierce <br>
!! Scripps Institution of Oceanography <br>
!! Climate Research Division <br>
!! dpierce@ucsd.edu <br>
!! Jan 29, 1996 <br>
!!
    real (kind=DP), intent(in) :: DATASET(:, :, :)
!!  The basic dataset array.  THIS MUST BE ANOMALIES.
    real (kind=DP), intent(in) :: DV(:)
!!  Control volume.
    integer, intent(in) :: ICOVCOR
!!  If .eq. covariance, then calculate the covariance array;
!!  otherwise, calculate the correlation array. [INTEGER]
    integer, intent(in) :: COVARIANCE
!!  Values to indicate each of these options. [INTEGER]
    logical, intent(in) :: DOSWITCHED
!!  if .TRUE., then calculate the 'switched' array (which is of size (nt,nt));
!!  if .FALSE., then calculate the normal array of size (nx,nx). [LOGICAL]
    real (kind=DP), intent(out) :: SPACKED(:)
!!  the covariance or correlation array.  This is in packed form corresponding
!!  to LAPACK's lower triangular form.

    integer :: I, J, K, R, NX, NT, NF
    real (kind=DP) :: FACTOR(3), SCALE

    NX = ubound( DATASET, 1 )
    NT = ubound( DATASET, 2 )
    NF = ubound( DATASET, 3 )

    if (NX<=1) then
      write (0, *) 'deofcovcor.F: error: nx too small!! nx=', NX
      stop 0
    end if

    if (ICOVCOR==COVARIANCE) then
      print *, ' Using covariance matrix'
    else
      print *, ' Using correlation matrix'
    end if

    if (DOSWITCHED) then
      do J = 1, NT
        do I = J, NT
          FACTOR(:) = 0.0
          do concurrent (K = 1:NX, R = 1:NF)

            FACTOR(1) = FACTOR(1) + DATASET(K, I, R)*DATASET(K, J, R)*DV(K)
            FACTOR(2) = FACTOR(2) + DATASET(K, I, R)*DATASET(K, I, R)*DV(K)
            FACTOR(3) = FACTOR(3) + DATASET(K, J, R)*DATASET(K, J, R)*DV(K)

          end do

          if (ICOVCOR==COVARIANCE) then
            SCALE = 1.0
          else
            SCALE = 1.0/(SQRT(FACTOR(2))*SQRT(FACTOR(3)))
          end if

          SPACKED(I+(J-1)*(2*NT-J)/2) = FACTOR(1)*SCALE
        end do
      end do
    end if

  end subroutine

  subroutine SOLVE_EIGENVALUES(SPACKED, EVALS, EVECS, ORDEROFS, NMODES, &
    IMINNXNT)
!! Now call the LAPACK solver to get the eigenvalues and eigenvectors.  The
!! eigenvalues express the amount of variance explained by the various modes, so
!! choose to return those 'nmodes' modes which explain the most variance.
!!
!! Remember that the calculated eigenvectors may not be the ones we really want
!! if we are doing switched X and T axes.  However the eigenvalues are the same
!! either way. Note that the LAPACK routine returns
!! the eigenvalues (and corresponding eigenvectors) in ASCENDING order, but this
!! routine returns them in DESCENDING order;  this will be switched in the final
!! assembly phase, below.

    integer, intent(in) :: ORDEROFS, NMODES, IMINNXNT
    real (kind=DP), allocatable, dimension (:), intent(in) :: SPACKED
!! the covariance or correlation array.  This is in packed form corresponding
!! to LAPACK's lower triangular form.
    real (kind=DP), allocatable, dimension (:), intent(inout) :: EVALS
!! the calculated eigenvalues
    real (kind=DP), allocatable, dimension (:, :), intent(inout) :: EVECS
!! the calculated eigenvectors

    integer :: I, J, N, M, IL, IU, LDZ, INFO
    integer, allocatable, dimension (:) :: ILAWORK, IFAIL
    character (len=1) :: JOBZ, RANGE, UPLO
    real (kind=DP) :: VL, VU, ABSTOL
    real (kind=DP), allocatable, dimension (:) :: RLAWORK

    allocate (RLAWORK(8*IMINNXNT), ILAWORK(5*IMINNXNT), IFAIL(IMINNXNT))

    ! Options for Lapack routine
    
    JOBZ = 'V'          ! Both eigenvalues and eigenvectors
    RANGE = 'I'         ! Specify range of eigenvalues to get.
    UPLO = 'L'          ! 'spacked' has lower triangular part of S
    N = ORDEROFS
    IL = N - NMODES + 1 ! Smallest eigenvalue to get
    IU = N              ! Largest eigenvalue to get
    LDZ = N
    ABSTOL = 0.0        ! See LAPACK documentation

    EVALS(:) = 0.0_DP
    EVECS(:, :) = 0.0_DP

    print *, 'about to call dspevx, spacked=', SPACKED(1:3)
    call DSPEVX(JOBZ, RANGE, UPLO, N, SPACKED, VL, VU, IL, IU, ABSTOL, M, &
      EVALS, EVECS, LDZ, RLAWORK, ILAWORK, IFAIL, INFO)

    if (INFO/=0) then
      if (INFO<0) then
        write (0, *) 'LAPACK error: argument ', -INFO, ' had illegal value'
        stop 'eof pierce.f90 line 136'
      else
        write (0, *) 'LAPACK error: ', INFO, &
          'eigenvectors failed to converge!'
        write (0, *) 'Consult the LAPACK docs!'
        stop 'eof pierce.f90 line 141'
      end if
    end if

    do I = 1, NMODES
      if (EVALS(I)<=0.0) then
        write (0, *) 'Error! LAPACK routine returned'
        write (0, *) 'eigenvalue <= 0!! ', I, EVALS(I)
        do J = 1, NMODES
          print *, J, EVALS(J)
        end do
        write (0, *) 'Note: This often means you are asking'
        write (0, *) 'for more modes than the dataset supports.'
        write (0, *) 'Try reducing the number of requested'
        write (0, *) 'modes.'
        stop 'eof line 156'
      end if
    end do

  end subroutine

  subroutine DEOFPCS(DATASET, EVECS, DOSWITCHED, TENTPCS)
!! Compute the tentative principal components; they are 'tentative' because they
!! might be the PCs of the switched data.  Put them into 'tentpcs', which is of
!! size (nt,nmodes) if we are doing regular EOFs and (nx,nmodes) if we are doing
!! switched EOFs.  These PCs come out in order corresponding to the order of
!! 'evecs', which is in ASCENDING order.
!!
!! David W. Pierce <br>
!! Scripps Institution of Oceanography <br>
!! Climate Research Division <br>
!! dpierce@ucsd.edu <br>
!! Jan 29, 1996 <br>
!!
    real (kind=DP), intent(in) :: DATASET(:, :, :)
!!   data(nx,nt): The input data.  THESE MUST BE ANOMALIES.
!!   nmodes: # of modes to calculate.
!!   iminnxnt: min(nx,nt)
    real (kind=DP), intent(in) :: EVECS(:, :)
!!   evecs(iminnxnt,nmodes): the eigenvectors (which might be switched).
    logical, intent(in) :: DOSWITCHED
!!   doswitched: if .TRUE., then we are doing switched (space,time) calculation;
!!     otherwise, regular (time,space) calculation.
    real (kind=DP), intent(out) :: TENTPCS(:, :)
!!   tentpcs(imaxnxnt,nmodes): the tentative	(possibly switched) principal components.
!! ------------------------------------------------------------------------------

    integer :: I, J, K, R, NX, NT, NF, IMINNXNT, IMAXNXNT, NMODES
    real (kind=DP) :: VAR

    NX = ubound( DATASET, 1 )
    NT = ubound( DATASET, 2 )
    NF = ubound( DATASET, 3 )
    IMINNXNT = ubound( EVECS, 1 )
    IMAXNXNT = ubound( TENTPCS, 1 )
    NMODES = ubound( TENTPCS, 2 )

    if (DOSWITCHED) then
      do concurrent ( I = 1:NX, J = 1:NMODES)
        VAR = 0.0_DP
        do concurrent ( K = 1:NT, R = 1:NF )
          VAR = VAR + DATASET(I, K, R)*EVECS(K, J)
        end do
        TENTPCS(I, J) = VAR
      end do
    end if

  end subroutine

end module
