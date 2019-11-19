MODULE RDistributions
!implicit none
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
            integer :: i, n, un, istat, dt(8), pid
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

END MODULE RDistributions