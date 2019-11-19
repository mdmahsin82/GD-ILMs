MODULE RDistributionss
implicit none
! Inspired from: http://www.johndcook.com/julia_rng.html
! Original author in julia : John D Cook
! coded : Sukhbinder in fortran
! Date : 28th Feb 2012
!
!		rand_uniform(a,b)   						function
!		rand_normal(mean,stdev)   					function
!		rand_exponential(mean)   					function
!		rand_gamma(shape, SCALE)   					function
!		rand_chi_square(dof)   						function
!		rand_inverse_gamma(shape, SCALE)   			function
!		rand_weibull(shape, SCALE)   				function
!		rand_cauchy(median, SCALE)   				function
!		rand_student_t(dof)   						function
!		rand_laplace(mean, SCALE)   				function
!		rand_log_normal(mu, sigma)   				function
!		rand_beta(a, b)   							function
!		m7_cumnor ( arg, presult, ccum ) 			subroutine
!		init_random_seed()    						subroutine
!		dnorm(z,mu,varin)   						function
!		density_gamma(xx,a,b)   					function
!		density_exp(xx,a)  							function
!		r8_gamma_log ( x )   						function
!		r8_gamma_pdf ( alph , bet , rval )   		function
!		half_normal_pdf ( x, a, b, pdf ) 			subroutine
!		beta_pdf ( x, a, b, pdf )  					subroutine
!		binomial_coef ( n, k, cnk )  				subroutine
!		binomial_pdf ( x, a, b, pdf )  				subroutine
!		gamma_log ( x )   							function
!		beta ( a, b )   							function
!
! Non uniform random Number Generators in Fortran
!
      DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462
      CONTAINS



      FUNCTION rand_uniform(a,b) RESULT(c)
       DOUBLE PRECISION :: a,b,c,temp
       CALL RANDOM_NUMBER(temp)
       c= a+temp*(b-a)
      END FUNCTION

!
! Random Sample from normal (Gaussian) distribution
!
      FUNCTION rand_normal(mean,stdev) RESULT(c)
       DOUBLE PRECISION :: mean,stdev,c,temp(2),r,theta
      IF(stdev <= 0.0d0) THEN

        WRITE(*,*) "Standard Deviation must be +ve"
      ELSE
        CALL RANDOM_NUMBER(temp)
        r=(-2.0d0*log(temp(1)))**0.5
        theta = 2.0d0*PI*temp(2)
        c= mean+stdev*r*sin(theta)
      END IF
      END FUNCTION


      
SUBROUTINE random_mvnorm(n, h, d, f, first, x, ier)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! N.B. An extra argument, ier, has been added to Dagpunar's routine

!     SUBROUTINE GENERATES AN N VARIATE RANDOM NORMAL
!     VECTOR USING A CHOLESKY DECOMPOSITION.

! ARGUMENTS:
!        N = NUMBER OF VARIATES IN VECTOR
!           (INPUT,INTEGER >= 1)
!     H(J) = J'TH ELEMENT OF VECTOR OF MEANS
!           (INPUT,REAL)
!     X(J) = J'TH ELEMENT OF DELIVERED VECTOR
!           (OUTPUT,REAL)
!
!    D(J*(J-1)/2+I) = (I,J)'TH ELEMENT OF VARIANCE MATRIX (J> = I)
!            (INPUT,REAL)
!    F((J-1)*(2*N-J)/2+I) = (I,J)'TH ELEMENT OF LOWER TRIANGULAR
!           DECOMPOSITION OF VARIANCE MATRIX (J <= I)
!            (OUTPUT,REAL)

!    FIRST = .TRUE. IF THIS IS THE FIRST CALL OF THE ROUTINE
!    OR IF THE DISTRIBUTION HAS CHANGED SINCE THE LAST CALL OF THE ROUTINE.
!    OTHERWISE SET TO .FALSE.
!            (INPUT,LOGICAL)

!    ier = 1 if the input covariance matrix is not +ve definite
!        = 0 otherwise

INTEGER, INTENT(IN)   :: n
double precision, INTENT(IN)      :: h(:), d(:)   ! d(n*(n+1)/2)
double precision, INTENT(IN OUT)     :: f(:)         ! f(n*(n+1)/2)
double precision, INTENT(OUT)     :: x(:)
LOGICAL, INTENT(IN)   :: first
INTEGER, INTENT(OUT)  :: ier

!     Local variables
INTEGER       :: j, i, m
double precision          :: y, v
INTEGER, SAVE :: n2

IF (n < 1) THEN
  WRITE(*, *) 'SIZE OF VECTOR IS NON POSITIVE'
  STOP
END IF

ier = 0
IF (first) THEN                        ! Initialization, if necessary
  n2 = 2*n
  IF (d(1) < 0.0d0) THEN
    ier = 1
    RETURN
  END IF

  f(1) = SQRT(d(1))
  y = 1.0d0/f(1)
  DO j = 2,n
    f(j) = d(1+j*(j-1)/2) * y
  END DO

  DO i = 2,n
    v = d(i*(i-1)/2+i)
    DO m = 1,i-1
      v = v - f((m-1)*(n2-m)/2+i)**2
    END DO

    IF (v < 0.0d0) THEN
      ier = 1
      RETURN
    END IF

    v = SQRT(v)
    y = 1.0d0/v
    f((i-1)*(n2-i)/2+i) = v
    DO j = i+1,n
      v = d(j*(j-1)/2+i)
      DO m = 1,i-1
        v = v - f((m-1)*(n2-m)/2+i)*f((m-1)*(n2-m)/2 + j)
      END DO ! m = 1,i-1
      f((i-1)*(n2-i)/2 + j) = v*y
    END DO ! j = i+1,n
  END DO ! i = 2,n
END IF

x(1:n) = h(1:n)
DO j = 1,n
  y = rand_normal(0.0d0,1.0d0)
!  y = random_normal()
  DO i = j,n
    x(i) = x(i) + f((j-1)*(n2-j)/2 + i) * y
  END DO ! i = j,n
END DO ! j = 1,n

RETURN
END SUBROUTINE random_mvnorm

!
! Random smaple from an exponential distribution
!
      FUNCTION rand_exponential(mean) RESULT(c)
      DOUBLE PRECISION :: mean,c,temp
      IF (mean <= 0.0d0) THEN

        WRITE(*,*) "mean must be positive"
      ELSE
       CALL RANDOM_NUMBER(temp)
       c=-(1/mean)*log(temp)
      END IF
      END FUNCTION

!
! Return a random sample from a gamma distribution
!
      RECURSIVE FUNCTION rand_gamma(shape, SCALE) RESULT(ans)
      DOUBLE PRECISION SHAPE,scale,u,w,d,c,x,xsq,g,ans,v
      IF (shape <= 0.0d0) THEN

        WRITE(*,*) "Shape PARAMETER must be positive"
      END IF
      IF (scale <= 0.0d0) THEN

        WRITE(*,*) "Scale PARAMETER must be positive"
      END IF
!
! ## Implementation based on "A Simple Method for Generating Gamma Variables"
! ## by George Marsaglia and Wai Wan Tsang.
! ## ACM Transactions on Mathematical Software

! ## Vol 26, No 3, September 2000, pages 363-372.
!
      IF (shape >= 1.0d0) THEN
        d = SHAPE - (1.0d0/3.0d0)
        c = 1.0d0/((9.0d0 * d)**0.5)
        DO while (.true.)
            x = rand_normal(0.0d0, 1.0d0)
            v = 1.0 + c*x
            DO while (v <= 0.0d0)
                x = rand_normal(0.0d0, 1.0d0)
                v = 1.0d0 + c*x
            END DO

            v = v*v*v
            CALL RANDOM_NUMBER(u)
            xsq = x*x
            IF ((u < 1.0d0 -.0331d0*xsq*xsq) .OR.  &
              (log(u) < 0.5d0*xsq + d*(1.0d0 - v + log(v))) )then
                ans=scale*d*v
                RETURN
            END IF

        END DO
      ELSE
        g = rand_gamma(shape+1.0d0, 1.0d0)
        CALL RANDOM_NUMBER(w)
        ans=scale*g*(w**(1.0d0/shape))
        RETURN
      END IF

      END FUNCTION
!
! ## return a random sample from a chi square distribution
! ## with the specified degrees of freedom
!
      FUNCTION rand_chi_square(dof) RESULT(ans)
      DOUBLE PRECISION ans,dof
         ans=rand_gamma(0.5d0, 2.0d0*dof)
      END FUNCTION

!
! ## return a random sample from an inverse gamma random variable
!
      FUNCTION rand_inverse_gamma(shape, SCALE) RESULT(ans)
      DOUBLE PRECISION SHAPE,scale,ans

! ## If X is gamma(shape, scale) then
! ## 1/Y is inverse gamma(shape, 1/scale)
      ans= 1.0d0 / rand_gamma(shape, 1.0d0 / SCALE)
      END FUNCTION
!
!## return a sample from a Weibull distribution
!

      FUNCTION rand_weibull(shape, SCALE) RESULT(ans)
      DOUBLE PRECISION SHAPE,scale,temp,ans
      IF (shape <= 0.0d0) THEN

        WRITE(*,*) "Shape PARAMETER must be positive"
      END IF
      IF (scale <= 0.0d0) THEN

        WRITE(*,*) "Scale PARAMETER must be positive"
      END IF
      CALL RANDOM_NUMBER(temp)
      ans= SCALE * (-log(temp))**(1.0 / SHAPE)
      END FUNCTION

!
!## return a random sample from a Cauchy distribution
!
      FUNCTION rand_cauchy(median, SCALE) RESULT(ans)
      DOUBLE PRECISION ans,median,scale,p

      IF (scale <= 0.0d0) THEN

        WRITE(*,*) "Scale PARAMETER must be positive"
      END IF
      CALL RANDOM_NUMBER(p)
      ans = median + SCALE*tan(PI*(p - 0.5))
      END FUNCTION

!
!## return a random sample from a Student t distribution
!
      FUNCTION rand_student_t(dof) RESULT(ans)
      DOUBLE PRECISION ans,dof,y1,y2
      IF (dof <= 0.d0) THEN

        WRITE(*,*) "Degrees of freedom must be positive"
      END IF
!
! ## See Seminumerical Algorithms by Knuth
      y1 = rand_normal(0.0d0, 1.0d0)
      y2 = rand_chi_square(dof)
      ans= y1 / (y2 / DOf)**0.50d0
!

      END FUNCTION

!
!## return a random sample from a Laplace distribution
!## The Laplace distribution is also known as the double exponential distribution.
!
      FUNCTION rand_laplace(mean, SCALE)  RESULT(ans)
      DOUBLE PRECISION ans,mean,scale,u
      IF (scale <= 0.0d0) THEN

        WRITE(*,*) "Scale PARAMETER must be positive"
      END IF
      CALL RANDOM_NUMBER(u)
      IF (u < 0.5d0) THEN

        ans = mean + SCALE*log(2.0*u)
      ELSE
        ans = mean - SCALE*log(2*(1-u))
      END IF

      END FUNCTION

!
! ## return a random sample from a log-normal distribution
!
      FUNCTION rand_log_normal(mu, sigma) RESULT(ans)
      DOUBLE PRECISION ans,mu,sigma
        ans= EXP(rand_normal(mu, sigma))
      END FUNCTION

!
! ## return a random sample from a beta distribution
!
      FUNCTION rand_beta(a, b) RESULT(ans)
      DOUBLE PRECISION :: a,b,ans,u,v
      IF ((a <= 0.0d0) .OR. (b <= 0.0d0)) THEN

        WRITE(*,*) "Beta PARAMETERs must be positive"
      END IF

! ## There are more efficient methods for generating beta samples.
! ## However such methods are a little more efficient and much more complicated.
! ## For an explanation of why the following method works, see
! ## http://www.johndcook.com/distribution_chart.html#gamma_beta

       u = rand_gamma(a, 1.0d0)
       v = rand_gamma(b, 1.0d0)
       ans = u / (u + v)
      END FUNCTION


! Normal CDF 
SUBROUTINE m7_cumnor ( arg, presult, ccum )
  !
  !*******************************************************************************
  !
  !! CUMNOR computes the cumulative normal distribution.
  !
  !
  !     the integral from -infinity to x of
  !          (1/sqrt(2*pi)) exp(-u*u/2) du
  !
  !  Author:
  !  -------
  !  Original source:
  !
  !    W. J. Cody    Mathematics and Computer Science Division
  !                  Argonne National Laboratory
  !                  Argonne, IL 60439
  !
  !    DCDFLIB is attributed to Barry Brown, James Lovato, and Kathy Russell
  !            bwb@odin.mda.uth.tmc.edu.
  !
  !    Adopted to ECHAM/M7:
  !
  !    Philip Stier  (MPI-MET)                    2001
  !    Luis Kornblueh (MPI-MET)                   2006
  !       (bugfixes for portability)
  !
  !  Reference:
  !  ----------
  !
  !    W D Cody, 
  !    "ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special 
  !    Function Routines and Test Drivers"
  !    ACM Transactions on Mathematical Software,
  !    Volume 19, 1993, pages 22-32.
  !
  !  Parameters:
  !
  !     ARG --> Upper limit of integration.
  !                                        X is double precision
  !
  !     presult <-- Cumulative normal distribution.
  !                                        presult is double precision
  !
  !     CCUM <-- Complement of Cumulative normal distribution.
  !                                        CCUM is double precision
  !
  !
  ! Original Comments:
  !
  !
  ! This function evaluates the normal distribution function:
  !
  !                              / x
  !                     1       |       -t*t/2
  !          P(x) = ----------- |      e       dt
  !                 sqrt(2 pi)  |
  !                             /-oo
  !
  !   The main computation evaluates near-minimax approximations
  !   derived from those in "Rational Chebyshev approximations for
  !   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
  !   This transportable program uses rational functions that
  !   theoretically approximate the normal distribution function to
  !   at least 18 significant decimal digits.  The accuracy achieved
  !   depends on the arithmetic system, the compiler, the intrinsic
  !   functions, and proper selection of the machine-dependent
  !   constants.
  !
  !  Explanation of machine-dependent constants.
  !
  !   MIN   = smallest machine representable number.
  !
  !   EPS   = argument below which anorm(x) may be represented by
  !           0.5  and above which  x*x  will not underflow.
  !           A conservative value is the largest machine number X
  !           such that   1.0 + X = 1.0   to machine precision.
  !
  !  Error returns
  !
  !  The program returns  ANORM = 0     for  ARG .LE. XLOW.
  !
  !  Author: 
  !
  !    W. J. Cody
  !    Mathematics and Computer Science Division
  !    Argonne National Laboratory
  !    Argonne, IL 60439
  !
  !  Latest modification: March 15, 1992
  !
  USE mo_kind, ONLY: dp
  !
  IMPLICIT NONE
  !
  real(dp), PARAMETER, DIMENSION ( 5 ) :: a = (/ &
       2.2352520354606839287e00_dp, &
       1.6102823106855587881e02_dp, &
       1.0676894854603709582e03_dp, &
       1.8154981253343561249e04_dp, &
       6.5682337918207449113e-2_dp /)
  real(dp) :: arg
  real(dp), PARAMETER, DIMENSION ( 4 ) :: b = (/ &
       4.7202581904688241870e01_dp, &
       9.7609855173777669322e02_dp, &
       1.0260932208618978205e04_dp, &
       4.5507789335026729956e04_dp /)
  real(dp), PARAMETER, DIMENSION ( 9 ) :: c = (/ &
       3.9894151208813466764e-1_dp, &
       8.8831497943883759412e00_dp, &
       9.3506656132177855979e01_dp, &
       5.9727027639480026226e02_dp, &
       2.4945375852903726711e03_dp, &
       6.8481904505362823326e03_dp, &
       1.1602651437647350124e04_dp, &
       9.8427148383839780218e03_dp, &
       1.0765576773720192317e-8_dp /)
  real(dp) :: ccum
  real(dp), PARAMETER, DIMENSION ( 8 ) :: d = (/ &
       2.2266688044328115691e01_dp, &
       2.3538790178262499861e02_dp, &
       1.5193775994075548050e03_dp, &
       6.4855582982667607550e03_dp, &
       1.8615571640885098091e04_dp, &
       3.4900952721145977266e04_dp, &
       3.8912003286093271411e04_dp, &
       1.9685429676859990727e04_dp /)
  real(dp) :: del
!@@@ real(dp) :: dpmpar
  real(dp) :: eps
  INTEGER :: i
  real(dp) :: zmin
  real(dp), PARAMETER, DIMENSION ( 6 ) :: p = (/ &
       2.1589853405795699e-1_dp, &
       1.274011611602473639e-1_dp, &
       2.2235277870649807e-2_dp, &
       1.421619193227893466e-3_dp, &
       2.9112874951168792e-5_dp, &
       2.307344176494017303e-2_dp /)
  real(dp), PARAMETER, DIMENSION ( 5 ) :: q = (/ &
       1.28426009614491121e00_dp, &
       4.68238212480865118e-1_dp, &
       6.59881378689285515e-2_dp, &
       3.78239633202758244e-3_dp, &
       7.29751555083966205e-5_dp /)
  real(dp) :: presult
  real(dp), PARAMETER :: root32 = 5.656854248_dp
  real(dp), PARAMETER :: sixten = 16.0_dp
  real(dp) :: temp
  real(dp), PARAMETER :: sqrpi = 3.9894228040143267794e-1_dp
  real(dp), PARAMETER :: thrsh = 0.66291_dp
  real(dp) :: x
  real(dp) :: xden
  real(dp) :: xnum
  real(dp) :: y
  real(dp) :: xsq
  !
  !  Machine dependent constants
  !
  eps = EPSILON ( 1.0_dp ) * 0.5_dp
  !
  !@@@ Simplified calculation of the smallest machine representable number
  !    (Higher accuracy than needed!)
  !
  !@@@ min = dpmpar(2)

  zmin = EPSILON ( 1.0_dp )

  x = arg
  y = ABS ( x )

  IF ( y <= thrsh ) THEN
     !
     !  Evaluate  anorm  for  |X| <= 0.66291
     !
     IF ( y > eps ) THEN
        xsq = x * x
     ELSE
        xsq = 0.0_dp
     END IF

     xnum = a(5) * xsq
     xden = xsq
     DO i = 1, 3
        xnum = ( xnum + a(i) ) * xsq
        xden = ( xden + b(i) ) * xsq
     END DO
     presult = x * ( xnum + a(4) ) / ( xden + b(4) )
     temp = presult
     presult = 0.5_dp + temp
     ccum = 0.5_dp - temp
     !
     !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
     !
  ELSE IF ( y <= root32 ) THEN

     xnum = c(9) * y
     xden = y
!CDIR UNROLL=7
     DO i = 1, 7
        xnum = ( xnum + c(i) ) * y
        xden = ( xden + d(i) ) * y
     END DO
     presult = ( xnum + c(8) ) / ( xden + d(8) )
     xsq = AINT ( y * sixten ) / sixten
     del = ( y - xsq ) * ( y + xsq )
     presult = EXP(-xsq*xsq*0.5_dp) * EXP(-del*0.5_dp) * presult
     ccum = 1.0_dp - presult

     IF ( x > 0.0_dp ) THEN
        temp = presult
        presult = ccum
        ccum = temp
     END IF
     !
     !  Evaluate  anorm  for |X| > sqrt(32).
     !
  ELSE

     presult = 0.0_dp
     xsq = 1.0_dp / ( x * x )
     xnum = p(6) * xsq
     xden = xsq
     DO i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
     END DO

     presult = xsq * ( xnum + p(5) ) / ( xden + q(5) )
     presult = ( sqrpi - presult ) / y
     xsq = AINT ( x * sixten ) / sixten
     del = ( x - xsq ) * ( x + xsq )
     presult = EXP ( - xsq * xsq * 0.5_dp ) * EXP ( - del * 0.5_dp ) * presult
     ccum = 1.0_dp - presult  

     IF ( x > 0.0_dp ) THEN
        temp = presult
        presult = ccum
        ccum = temp
     END IF

  END IF

  IF ( presult < zmin ) THEN
     presult = 0.0_dp
  END IF

  IF ( ccum < zmin ) THEN
     ccum = 0.0_dp
  END IF

END SUBROUTINE m7_cumnor

subroutine init_random_seed()
            use iso_fortran_env, only: int64
!            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid, getpid
            integer(int64) :: t

          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if


            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for double precision work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg

end subroutine init_random_seed

double precision function dnorm(z,mu,varin) result(normpdf)
implicit none
double precision,intent(in):: z,mu,varin
double precision,parameter:: logroot2pi = 0.9189385
double precision::arg,var
var=varin
var = var**2
arg = -0.5*((z-mu)**2)/var
arg = -logroot2pi - 0.5*dlog(var) + arg
normpdf = dexp(arg)
end function dnorm

! Gamma density 
function density_gamma(xx,a,bb) result(dn)
implicit none
double precision:: xx,a,bb,dn
    dn = (xx**(a-1)) * dexp(- (xx*bb))!((a-1)* dlog(xx)) - (b*xx)
end function

! Exponential density 
function density_exp(xx,aa) result(ddn)
implicit none
double precision:: xx,aa,ddn
ddn = (aa) * dexp(-(xx*aa))
end function

! Beta distribution
subroutine beta_pdf(x, a, b, pdf)
implicit none
double precision, intent(in) :: x , a, b
double precision, intent(out) :: pdf
double precision :: y

if ( x <= 0.0d0 .and. x >= 1.0d0 ) then
    pdf = 0.0d00
  else
  y = x
  pdf = (y**(a - 1))* ((1 - y)**(b - 1))
end if

end subroutine beta_pdf

! Inverse-gamma density

subroutine inverse_gamma_pdf(x, a, b, pdf)
implicit none
double precision, intent(in) :: x, a, b
double precision, intent(out) :: pdf
double precision :: y

if ( x <= 0.0d0 ) then
    pdf = 0.0d00
  else
  y = x
  pdf = (y**(- a - 1))* dexp(-b/y)
end if

end subroutine inverse_gamma_pdf





subroutine half_normal_pdf ( x, b, pdf )
!*****************************************************************************80
!
!! HALF_NORMAL_PDF evaluates the Half Normal PDF.
!
!  Discussion:
!
!    PDF(A,B;X) =
!      sqrt ( 2 / PI ) * ( 1 / B ) * exp ( - 0.5D+00 * ( ( X - A ) / B )^2 )
!
!    for A <= X
!
!    The Half Normal PDF is a special case of both the Chi PDF and the
!    Folded Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, double precision ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, double precision ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none
  double precision,intent(in)::b,x
  double precision,intent(out)::pdf
  double precision::y

  double precision, parameter :: pi = 3.141592653589793D+00
  if ( x <= 0.0d0 ) then
    pdf = 0.0d00
  else
    y = x
        pdf = dsqrt(2.0d0/(b*b*pi))* dexp(- (y*y)/(2*b*b) )
  end if

end subroutine half_normal_pdf

subroutine normal_pdf ( x, a, b, pdf )

  implicit none

  double precision:: a,b,pdf,x,y
  double precision, parameter :: pi = 3.141592653589793D+00

  y = ( x - a ) / b

  pdf = exp ( - 0.5D+00 * y * y )  / ( b * sqrt ( 2.0D+00 * pi ) )

end subroutine normal_pdf

function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end function

subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x^(n-1) + a(n) * x^(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X^N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mag
  character plus_minus
  integer ( kind = 4 ) r8poly_degree
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( a(n) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n) )

  if ( 2 <= n ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n
  else if ( n == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n - 1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end subroutine r8poly_print

function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end FUNCTION

function r8poly_value_horner ( m, c, x ) result(value)

!*****************************************************************************80
!
!! R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
!
!  Discussion:
!
!    The polynomial 
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the value X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the degree.
!
!    Input, real ( kind = 8 ) C(0:M), the polynomial coefficients.  
!    C(I) is the coefficient of X^I.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE_HORNER, the polynomial value.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(0:m)
  integer ( kind = 4 ) i
  !real ( kind = 8 ):: r8poly_value
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = c(m)
  do i = m - 1, 0, -1
    value = value * x + c(i)
  end do

!  r8poly_value_horner = value

!  return
end function r8poly_value_horner

subroutine normal_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
  real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
  real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
  real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
  real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
  real ( kind = 8 ), parameter :: b1 = 3.8052D-08
  real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
  real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
  real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
  real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
  real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
  real ( kind = 8 ),intent(out):: cdf
  real ( kind = 8 ) q
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28D+00 ) then

    y = 0.5D+00 * x * x

    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7D+00 ) then

    y = 0.5D+00 * x * x

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0D+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

!  return
end subroutine normal_01_cdf


subroutine normal_01_cdf_inv ( p, x )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10^16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2015
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an
!    "infinite" value will be returned.
!
!    Output, real ( kind = 8 ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = 8 ), parameter :: const1 = 0.180625D+00
  real ( kind = 8 ), parameter :: const2 = 1.6D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
!  real ( kind = 8 ) r8poly_value_horner
  real ( kind = 8 ), parameter :: split1 = 0.425D+00
  real ( kind = 8 ), parameter :: split2 = 5.0D+00
  real ( kind = 8 ),intent(out):: x

  if ( p <= 0.0D+00 ) then
    x = - huge ( x )
    return
  end if

  if ( 1.0D+00 <= p ) then
    x = huge ( x )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    x = q * r8poly_value_horner ( 7, a, r ) /&
        &  r8poly_value_horner ( 7, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then

      x = huge ( x )

    else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value_horner ( 7, c, r ) /&
          &  r8poly_value_horner ( 7, d, r )

      else

        r = r - split2
        x = r8poly_value_horner ( 7, e, r )/ &
          & r8poly_value_horner ( 7, f, r )

      end if

    end if

    if ( q < 0.0D+00 ) then
      x = -x
    end if

  end if

!  return
end subroutine normal_01_cdf_inv

subroutine truncated_normal_ab_sample ( mu, sigma, a, b, seed, x )

!*****************************************************************************80
!
!! TRUNCATED_NORMAL_AB_SAMPLE samples the truncated Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) MU, SIGMA, the mean and standard deviation of the
!    parent Normal distribution.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper truncation limits.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha_cdf
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) beta_cdf
  real ( kind = 8 ) mu
!  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) sigma
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u
  real ( kind = 8 ),intent(out):: x
  real ( kind = 8 ) xi
  real ( kind = 8 ) xi_cdf

  alpha = ( a - mu ) / sigma
  beta = ( b - mu ) / sigma

  call normal_01_cdf ( alpha, alpha_cdf )
  call normal_01_cdf ( beta, beta_cdf )

  !u = r8_uniform_01 ( seed )
  call random_number(u)
  xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf )
  call normal_01_cdf_inv ( xi_cdf, xi )

  x = mu + sigma * xi

!  return
end subroutine truncated_normal_ab_sample




END MODULE RDistributionss
