module INPUT

  use declarations

  implicit none!

contains

  !> Reads a gnuplot-compatible binary matrix file
  subroutine INPUT_MATRIX(NAME, N, F, X, Y, IOS, CCW)

    character (len=300), intent (in) :: NAME
    real (kind=SP), intent (out), allocatable, dimension (:) :: X, Y
    real (kind=SP), intent (out), allocatable, dimension (:, :) :: F
    real (kind=SP), intent (out) :: N
    integer, optional :: CCW

    integer, parameter :: NUNIT = 2
    integer :: IOS
    integer :: I, J, M

    print *, 'name=', TRIM(ADJUSTL(NAME))
    open (unit=NUNIT, file=TRIM(ADJUSTL(NAME)), iostat=IOS, status='old', &
      access='stream')
    if (IOS/=0) print *, 'Error in input_matrix()'

    read (NUNIT, pos=1) N

    M = INT(N) + 2

    if (.not. ALLOCATED(X)) allocate (X(M))
    if (.not. ALLOCATED(Y)) allocate (Y(M))
    if (.not. ALLOCATED(F)) allocate (F(M,M))

    F(1, :) = 0.0
    F(M, :) = 0.0
    F(:, 1) = 0.0
    F(:, M) = 0.0

    Y(1) = -0.5
    do J = 2, M - 1
      read (NUNIT) Y(J)
    end do
    Y(M) = 0.5

    if ( present( CCW )) then
      if (CCW==1) then
        X(1) = -0.5
        do I = 2, M - 1
          read (NUNIT) X(I)
          do J = 2, M - 1
            read (NUNIT) F(I, J)
          end do
        end do
        X(M) = 0.5
      else
        X(M) = 0.5
        do I = M - 1, 2, -1
          read (NUNIT) X(I)
          do J = 2, M - 1
            read (NUNIT) F(I, J)
          end do
        end do
        X(1) = -0.5
      end if
    else
      X(1) = -0.5
      do I = 2, M - 1
        read (NUNIT) X(I)
        do J = 2, M - 1
          read (NUNIT) F(I, J)
        end do
      end do
      X(M) = 0.5
    end if

    print *, 'n=', INT(N)
    close (NUNIT)

  end subroutine

  !> Updates the boundary conditions for the RBNT configuration
  subroutine BOUNDARY_CONDITIONS(M, T, U, V, YREF)

    integer, intent (in) :: M
    real (kind=SP), intent (inout), allocatable, dimension (:, :), &
      optional :: T, U, V, YREF

!1. Temperature BC
    if (PRESENT(T)) then
      T(:, 1) = 1.0 ! Bottom
      T(:, M) = 0.0 ! Top
      T(1, :) = T(2, :) ! Left
      T(M, :) = T(M-1, :) ! Right
    end if
!2. X-Velocity BC: no slip
    if (PRESENT(U)) then
      U(:, 1) = 0.0 ! Bottom
      U(:, M) = 0.0 ! Top
      U(1, :) = 0.0 ! Left
      U(M, :) = 0.0 ! Right
    end if
!3. Y-Velocity BC:  no slip
    if (PRESENT(V)) then
      V(:, 1) = 0.0 ! Bottom
      V(:, M) = 0.0 ! Top
      V(1, :) = 0.0 ! Left
      V(M, :) = 0.0 ! Right
    end if
!4. Reference State BC
    if (PRESENT(YREF)) then
      YREF(:, 1) = 1.0 ! Bottom
      YREF(:, M) = 0.0 ! Top
      YREF(1, :) = YREF(2, :) ! Left
      YREF(M, :) = YREF(M-1, :) ! Right
    end if
  end subroutine

  subroutine BOUNDARY_CONDITIONS_GHOST(BC, M, T, U, V, YREF)

    character, intent (in) :: BC
    integer, intent (in) :: M
    real (kind=SP) :: S
    real (kind=SP), intent (inout), dimension (:, :), &
      optional :: T, U, V, YREF

    S = -1.
    if (BC=='S') then
      S = 1.
    end if

!1. Temperature BC
    if (PRESENT(T)) then
      T(:, 1) = 2.0 - T(:, 2) ! Bottom
      T(:, M) = 0.0 - T(:, M-1) ! Top
      T(1, :) = T(2, :) ! Left
      T(M, :) = T(M-1, :) ! Right
    end if
!2. X-Velocity BC: no slip
    if (PRESENT(U)) then
      U(:, 1) = -U(:, 2) ! Bottom
      U(:, M) = S*U(:, M-1) ! Top
      U(1, :) = -U(2, :) ! Left
      U(M, :) = -U(M-1, :) ! Right
    end if
!3. Y-Velocity BC:  no slip
    if (PRESENT(V)) then
      V(:, 1) = -V(:, 2) ! Bottom
      V(:, M) = -V(:, M-1) ! Top
      V(1, :) = -V(2, :) ! Left
      V(M, :) = -V(M-1, :) ! Right
    end if
!4. Reference State BC
    if (PRESENT(YREF)) then
      YREF(:, 1) = 2.0 - YREF(:, 2) ! Bottom
      YREF(:, M) = 0.0 - YREF(:, M-1) ! Top
      YREF(1, :) = YREF(2, :) ! Left
      YREF(M, :) = YREF(M-1, :) ! Right
    end if
  end subroutine

  !> Updates the boundary conditions for the RBST configuration
  subroutine BOUNDARY_CONDITIONS_RBST(M, T, U, V, YREF)
    integer, intent (in) :: M
    real (kind=SP), intent (inout), allocatable, dimension (:, :), &
      optional :: T, U, V, YREF

!1. Temperature BC
    if (PRESENT(T)) then
      T(:, 1) = T(:, 2) ! Bottom
      T(:, M) = T(:, M-1) ! Top
      T(1, :) = T(2, :) ! Left
      T(M, :) = T(M-1, :) ! Right
    end if

!2. X-Velocity BC:
    if (PRESENT(U)) then
      U(:, 1) = 0.0 ! Bottom
      U(:, M) = U(:, M-1) ! Top (free-slip)
      U(1, :) = 0.0 ! Left
      U(M, :) = 0.0 ! Right
    end if
!3. Y-Velocity BC:  no slip
    if (PRESENT(V)) then
      V(:, 1) = 0.0 ! Bottom
      V(:, M) = 0.0 ! Top
      V(1, :) = 0.0 ! Left
      V(M, :) = 0.0 ! Right
    end if
!4. Reference State BC
    if (PRESENT(YREF)) then
      YREF(:, 1) = 1.0 ! Bottom
      YREF(:, M) = 0.0 ! Top
      YREF(1, :) = YREF(2, :) ! Left
      YREF(M, :) = YREF(M-1, :) ! Right
    end if
  end subroutine

  subroutine BOUNDARY_CONDITIONS_INTERPOL(T, U, V, M)
  !> Updates the boundary conditions for the no-slip condition
    integer, intent (in) :: M
    real (kind=SP), intent (inout), allocatable, dimension (:, :) :: T, U, V

!1. Temperature BC
    T(:, 1) = 1.0 ! Bottom
    T(:, M) = 0.0 ! Top
    T(1, :) = T(2, :) ! Left
    T(M, :) = T(M-1, :) ! Right

!2. X-Velocity BC: no slip
    U(:, 1) = 0.0 ! Bottom
    U(:, M) = 0.0 ! Top
    U(1, :) = 0.0 ! Left
    U(M, :) = 0.0 ! Right

!3. Y-Velocity BC:  no slip
    V(:, 1) = 0.0 ! Bottom
    V(:, M) = 0.0 ! Top
    V(1, :) = 0.0 ! Left
    V(M, :) = 0.0 ! Right

  end subroutine

end module
