module OUTPUT

  use DECLARATIONS

contains

  subroutine WRITEPPM2MATRIX(M, TEXT)

    character (len=*) :: TEXT
    integer :: COLS, ROWS, I, J
    real (kind=SP) :: MAXVALUE, MINVALUE
    real (kind=SP), allocatable :: M(:, :)

! Open File
    open (unit=100, file=TRIM(TEXT)//'.pgm', status='unknown')

! Write Header and ppm file type
    write (100, '( A )') 'P2'
    write (100, '( A )') '# PPM Type 2 File (generated with fortran)'

! Write Image Size
    COLS = SIZE(M, 2)
    ROWS = SIZE(M, 1)
    write (100, '( I5,1X,I5)') COLS, ROWS

! Write Maximum Value
    MAXVALUE = MAX(MAXVAL(M), 0.5)
    MINVALUE = MIN(MINVAL(M), -0.5)
    write (100, '(I5)') 1024

! Write Image
    do I = 1, ROWS
      do J = 1, COLS
        write (100, '( I5,1X )', advance='no') NINT(1024.*(M(I, &
          J)-MINVALUE)/(MAXVALUE-MINVALUE))
      end do
      write (100, *) ! Endline
    end do
    close (100)

  end subroutine

  pure subroutine LINT(XA, YA, N, X, Y)
! Given  the  arrays xa(1:n) and ya(1:n) of length n,  which  tabulate  a
! function and  given  a value  of x, this routine  returns a linear
! interpolated  value y
    integer, intent (in) :: N
    real (kind=SP), intent (in) :: XA(N), YA(N), X
    real (kind=SP), intent (out) :: Y

    integer :: K, KHI, KLO
    real (kind=SP) :: H, A, B

    KLO = 1
    KHI = N
    do while (KHI-KLO>1)
      K = (KHI+KLO)/2
      if (XA(K)>X) then
        KHI = K
      else
        KLO = K
      end if
    end do
    H = XA(KHI) - XA(KLO)
!if (h.eq.0.) pause
    A = (XA(KHI)-X)/H
    B = (X-XA(KLO))/H
    Y = A*YA(KLO) + B*YA(KHI)

  end subroutine

  subroutine LORENZ_STATE(T, THOT, M, YREF)

    integer, intent (in) :: M
    real (kind=SP), intent (in) :: THOT
    real (kind=SP), intent (in), allocatable, dimension (:, :) :: T

    integer :: I, J, BIN
    integer, parameter :: NBIN = 512
    real (kind=SP), parameter :: TCOLD = 0.
    real (kind=SP), dimension (NBIN) :: PDF, CDF, THETA
    real (kind=SP), intent (out), allocatable, dimension (:, :) :: YREF

    THETA(1) = TCOLD
    do I = 2, NBIN - 1
      THETA(I) = (FLOAT(I)-0.5)*(THOT-TCOLD)/FLOAT(NBIN-2) + TCOLD
    end do
    THETA(NBIN) = THOT

    PDF(:) = 0.
    do I = 2, M - 1
      do J = 2, M - 1
        BIN = NINT((T(I,J)-TCOLD)/(THOT-TCOLD)*(NBIN-1)) + 1
        BIN = MAX(BIN, 1)
        BIN = MIN(BIN, NBIN-1)
        PDF(BIN) = PDF(BIN) + 1.
      end do
    end do

    CDF(:) = 0.
    do I = 2, NBIN - 1
      CDF(I) = CDF(I-1) + PDF(I)/(M-2)/(M-2)
    end do
    CDF(NBIN) = 1.

    allocate (YREF(M,M))
    YREF(:, :) = 0.0_DP
    do I = 2, M - 1
      do J = 2, M - 1
        call LINT(THETA, CDF, NBIN, T(I,J), YREF(I,J))
      end do
    end do

  end subroutine

  subroutine STREAMFUNCTION(U, V, M, PSI)

    integer, intent (in) :: M
    real (kind=SP), intent (in), allocatable, dimension (:, :) :: U, V
    real (kind=SP), intent (out), allocatable, dimension (:, :) :: PSI

    integer :: I, J
    real (kind=SP) :: DX, DY

    allocate (PSI(M,M))
    PSI(:, :) = 0.
    DX = 0.5/FLOAT(M-2)
    DY = 1./FLOAT(M-2)
    do J = 2, M - 1
      if ((J==2) .or. (J==M+1)) DY = 0.5/(M-2)
      do I = 2, M - 1
        if ((I==2) .or. (I==M+1)) DX = 0.5/(M-2)
        PSI(I+1, J) = PSI(I, J) - 0.5*(V(I,J)+V(I+1,J))/DX
      end do
      PSI(2, J+1) = PSI(2, J) + 0.5*(U(2,J+1)+U(2,J))/DY
    end do

  end subroutine

  subroutine GRADIENT(F, X, M, GX)
! This routine provides a second order finite differences approximation of the
! gradient in the x-direction. Using CDS with an irregular grid, the first
! derivative can be approximated as :
!
!     - h² f(x-g) + (h² - g²) f(x) + g²  f(x+h)
! f'(x) = ----------------------------------------
!         (g²h + gh²)
!
! while using FDS, it can be approximated as :
!
!      (h²-(g+h)²) f(x) + (g+h)² f(x+h) - h² f(x+h+g)
! f'(x) = -----------------------------------------------
!        h(g+h)² - h²(g+h)

    integer, intent (in) :: M
    real (kind=SP), intent (in), dimension (M) :: X
    real (kind=SP), intent (in), dimension (M) :: F
    real (kind=DP), dimension (M) :: GX

    integer :: I
    real (kind=SP) :: G, H
    real (kind=DP) :: D
    real (kind=DP), dimension (3) :: C


    H = X(2) - X(1)
    G = X(3) - X(2)
    D = H*((G+H)**2.d0) - (H**2.d0)*(G+H)
    C(1) = (H**2.d0) - ((G+H)**2.d0)
    C(2) = ((G+H)**2.d0)
    C(3) = -(H**2.d0)
    GX(1) = SUM(F(1:3)*C(1:3))/D


    H = X(M) - X(M-1)
    G = X(M-1) - X(M-2)
    D = H*((G+H)**2.d0) - (H**2.d0)*(G+H)
    C(1) = (H**2.d0) - ((G+H)**2.d0)
    C(2) = ((G+H)**2.d0)
    C(3) = -(H**2.d0)
    GX(M) = SUM(-F(M:M-2:-1)*C(1:3))/D

    do I = 2, M - 1
      H = X(I+1) - X(I)
      G = X(I) - X(I-1)
      D = (G**2.d0)*H + G*(H**2.d0)
      C(1) = -H**2.d0
      C(2) = H**2.d0 - G**2.d0
      C(3) = G**2.d0
      GX(I) = SUM(F(I-1:I+1)*C(1:3))/D
    end do

  end subroutine

  subroutine INTERPOL(FILENAME, PCT, U, V, T, X, Y, CCW)

    use :: INPUT, only: INPUT_MATRIX, BOUNDARY_CONDITIONS_INTERPOL

    integer :: M, IOS
    real (kind=SP) :: N, PCT
    real (kind=SP), allocatable, dimension (:) :: X, Y
    real (kind=SP), allocatable, dimension (:, :) :: U1, U2, V1, V2, T1, T2
    real (kind=SP), allocatable, dimension (:, :) :: U, V, T
    character (len=300), dimension (6) :: FILENAME
    integer :: CCW

    call INPUT_MATRIX(FILENAME(1), N, T1, X, Y, IOS, CCW)
    call INPUT_MATRIX(FILENAME(2), N, U1, X, Y, IOS, CCW)
    call INPUT_MATRIX(FILENAME(3), N, V1, X, Y, IOS, CCW)
    call INPUT_MATRIX(FILENAME(4), N, T2, X, Y, IOS, CCW)
    call INPUT_MATRIX(FILENAME(5), N, U2, X, Y, IOS, CCW)
    call INPUT_MATRIX(FILENAME(6), N, V2, X, Y, IOS, CCW)

    M = INT(N) + 2

    if (.not. ALLOCATED(T)) allocate (T(M,M))
    T(:, :) = (1.-PCT)*T1(:, :) + PCT*T2(:, :)

    if (.not. ALLOCATED(U)) allocate (U(M,M))
    U(:, :) = (1.-PCT)*U1(:, :) + PCT*U2(:, :)

    if (.not. ALLOCATED(V)) allocate (V(M,M))
    V(:, :) = (1.-PCT)*V1(:, :) + PCT*V2(:, :)


    call BOUNDARY_CONDITIONS_INTERPOL(T, U, V, M)

  end subroutine

  subroutine OUTPUT_MATRIX(NAME, M, F, X, Y)

    character (len=300), intent (in) :: NAME
    real (kind=SP), intent (in), allocatable, dimension (:) :: X, Y
    real (kind=SP), intent (in), allocatable, dimension (:, :) :: F
    integer, intent (in) :: M

    integer, parameter :: NUNIT = 80

    integer :: I, J, RECLEN
    real (kind=SP) :: N

    N = FLOAT(M) - 2.

    inquire (iolength=RECLEN) N
    open (unit=NUNIT, file=TRIM(NAME), form='UNFORMATTED', access='STREAM')
    write (unit=NUNIT) N
    do J = 2, M - 1
      write (unit=NUNIT) Y(J)
    end do
    do I = 2, M - 1
      write (unit=NUNIT) X(I)
      do J = 2, M - 1
        write (unit=NUNIT) F(I, J)
      end do
    end do
    close (unit=NUNIT)

  end subroutine


  subroutine INPUT_MATRIX_VTK(NAME, NSNAP)

    character (len=300), intent (in) :: NAME
    integer, intent(in) :: NSNAP

    integer, parameter :: NUNIT = 2
    integer :: IOS

    character (len=300), dimension (3) :: FILENAMEOUT
    character (len=50) :: HEADER
    integer :: NPOINTS, N, I, J
    real (kind=SP), allocatable, dimension (:) :: X, Y
    real (kind=SP), allocatable, dimension (:, :) :: T, YREF, U, V

    print *, 'name=', TRIM(NAME)
    open (unit=NUNIT, file=TRIM(NAME), iostat=IOS, status='old')
    if (IOS/=0) print *, 'Error in input_matrix()'

    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a11,i3,i3,i3)') HEADER, N
    print *, N
    read (NUNIT, '(a7,i6)') HEADER, NPOINTS
    print *, NPOINTS

    allocate (X(N), Y(N), T(N,N), YREF(N,N), U(N,N), V(N,N))

    do I = 1, N
      do J = 1, N
        read (NUNIT, *) X(I), Y(J)
      end do
    end do

    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    do I = 1, N
      do J = 1, N
        read (NUNIT, *) T(I, J)
      end do
    end do

    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    do I = 1, N
      do J = 1, N
        read (NUNIT, *) YREF(I, J)
      end do
    end do

    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    do I = 1, N
      do J = 1, N
        read (NUNIT, *) U(I, J)
      end do
    end do

    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    read (NUNIT, '(a)') HEADER
    print *, HEADER
    do I = 1, N
      do J = 1, N
        read (NUNIT, *) V(I, J)
      end do
    end do


    write (FILENAMEOUT(1), '(a,i6.6,a)') 'field-T-', NSNAP, '.bin'
    write (FILENAMEOUT(2), '(a,i6.6,a)') 'field-U-', NSNAP, '.bin'
    write (FILENAMEOUT(3), '(a,i6.6,a)') 'field-V-', NSNAP, '.bin'

    call OUTPUT_MATRIX(FILENAMEOUT(1), N, T, X, Y)
    call OUTPUT_MATRIX(FILENAMEOUT(2), N, U, X, Y)
    call OUTPUT_MATRIX(FILENAMEOUT(3), N, V, X, Y)

    close (NUNIT)

  end subroutine

end module
