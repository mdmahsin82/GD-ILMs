!******************************************************************
! Random-walk MCMC iteration for region-restricted GD-ILMs disease under SIR compartmental framework
! Author: Md Mahsin (md.mahsin@ucalgary.ca)
! Date:  April 20, 2017
! Modified: May 28, 2018
! Mcmc updates of parameters
! run :: ifort ncarmcmc.f90 ditribution.f90 knd.f90 random.f90 -O2 
! True value : alpha = 0.30, alpha1 = 0.40, beta = 4.0, lambda = 0.8, sigma = 0.6, infectious period = 3
! ******************************************************************

program carepidemic
use RDistributionss
use mo_kind
!use random
USE ieee_arithmetic

implicit none
! n define the number of observations and w is 
! the number of iterations of MCMC
! total local geography areas in Calgary metropolitan is 16 

integer, parameter :: a = 16, n = 1570, w = 310000, mm = 136
integer :: T, i, j, k, m, l, it
integer :: acc1, acc2, acc3, acc4, acc10
double precision, dimension(n) :: x, y, pop
double precision, dimension(w) :: alpha, beta, lambda, sigma, tau, alpha1
integer, dimension(a) :: tallyphi
double precision, dimension(w, a) :: phit, psit
double precision, dimension(a) :: phi=0.10d0, phisv, prophi, psi
double precision, dimension(a) :: phinew 
double precision, dimension(a, 1) :: tphinew 
double precision, dimension(a, a):: Wstar
double precision, dimension(a) :: Wstarval
double precision, dimension(a, a) :: idty
integer, dimension(a, 8) :: neighbor
integer, dimension(n) :: rownumber,  area, inftime
integer :: abs1, abs2, abs3, abs4, abs10
double precision :: start, finish
double precision, parameter :: sigalpha0 = 10.0d0, sigbeta0 = 10.0d0 ! alpha
double precision, parameter :: sigalpha1 = 10.0d0, sigbeta1 = 10.0d0 !alpha1
double precision, parameter :: sigalpha2 = 10.0d0, sigbeta2 = 10.0d0 ! beta
double precision, parameter ::  a1 = 0.50d0,  b1 = 0.50d0, a2 = 0.50d0, b2 = 0.50d0, a3 = 3.0d0, b3 = 0.50d0
double precision ::  cdf1, ccdf1, num1, den1, pr1, accp1, u1, prioralpha1, prioralpha2
double precision :: cdf2, ccdf2, num2, den2, pr2, accp2, u2, priorbeta1, priorbeta2
double precision :: y33, y22, num3, den3, pr3, pr33, accp3, u3, plambda1, plambda2, f11
double precision :: y44, num4, y11, den4, pr4, accp4, u4, psigma1, psigma2
double precision :: num5, den5, pr5, accp5, u5, phiprior2, phiprior1
double precision :: u10, den10, num10, pr10, accp10, prioralpha21, prioralpha20
double precision :: u ,  v, phis, Inf, phimean, phivar, avgphi
double precision, dimension (a, a) :: propqlambda, slambda
double precision, dimension (a) :: proppre 

double precision, dimension (a, a) :: qlambda, plambda
double precision, dimension (a) :: lpre, prolamb, proslamb
double precision :: detp, ppre, detlambda, propdetlambda, propsigma2, propsigma1, priordenom 

double precision :: galpha = 0.5d0, gbeta = 0.5d0, lalpha = 4.0d0, lbeta = 2.0d0
double precision :: taurate, taushape, tauscale, phisd, clambda
double precision, dimension(mm) :: blocksig
double precision, dimension(mm) :: aower
double precision :: rlikn, rlikd

IF (ieee_support_inf(Inf)) THEN
    Inf = ieee_value(Inf,  ieee_positive_inf)
  END IF
! ************************************
! Fixed effects parameters
! Save the output file inside the Loop
! Fixed effects parameters
open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter1.txt")
!open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter2.txt")
!open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter3.txt")
!open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter4.txt")
!open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter5.txt")
!open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter6.txt")
!open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter7.txt")
!open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter8.txt")
!open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter9.txt")
!open(unit=10, file = "/home/md.mahsin/nsimulation/Highcov1/parameter10.txt")

! Random effects parameter

open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi1.txt")
!open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi2.txt")
!open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi3.txt")
!open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi4.txt")
!open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi5.txt")
!open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi6.txt")
!open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi7.txt")
!open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi8.txt")
!open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi9.txt")
!open(unit=100, file="/home/md.mahsin/nsimulation/Highcov1/phi10.txt")

! ************************************
! input file of neighboring counties
! ************************************
open(50, file = "/home/md.mahsin/nsimulation/Highcov1/neib.txt")
!************************
do i = 1, a
read(50, *) neighbor(i, 1:8)
end do
close(50)
! **********************************
! input file precision matrix for LCAR

! input file precision matrix for LCAR
open(80, file = "/home/md.mahsin/nsimulation/Highcov1/Q_rho.txt" )
!**************************
do i = 1, a
read(80, *) Wstar(i, 1:a)
end do 
close(80)


! input file for identity matrix

open(8000, file = "/home/md.mahsin/nsimulation/Highcov1/identity.txt" )
!****************************
do i = 1, a
read(8000, *) idty(i, 1:a)
end do
close(8000)

! **********************************
! Simulated epidemics data 
! **********************************
open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim1.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim2.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim3.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim4.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim5.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim6.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim7.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim8.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim9.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim10.txt" )
!*************************************
do i = 1, n
read(500, *) rownumber(i), x(i), y(i), area(i), inftime(i), pop(i)
end do
close(500)
! *************************************
alpha(1) = 0.50d0
alpha1(1) = 0.50d0
beta(1) = 5.50d0
lambda(1) = 0.50d0
tau(1) = 2.5d0
phit(1, 1:a) = phi 
psit(1, 1:a) = 1.0d0 
acc1 = 0
acc2 = 0
acc3 = 0
acc4 = 0
acc10 = 0
! ********************************

call init_random_seed()
call cpu_time (start)

! *******************************
! start MCMC 

! *****************************

do i = 1, w

! updating alpha parameter
! proposal density
abs1  = 0
do while (abs1 .eq. 0)
y11 = rand_normal(alpha(i), 0.05d0)
	if (y11 .gt. 0.0d0) then
	abs1 = 1
	else
	abs1 = 0
	endif
end do

! prior distribution; half normal prior
call half_normal_pdf(y11, sigalpha0, prioralpha2)
call half_normal_pdf(alpha(i), sigalpha0, prioralpha1)


! posterior distribution for alpha
num1 = likelihood(y11, alpha1(i), beta(i), phit(i,:), 100) + log(prioralpha2)
den1 = likelihood(alpha(i), alpha1(i), beta(i), phit(i, :), 100) + log(prioralpha1)
pr1 = exp(num1 - den1)

accp1 = min(1.0d0, pr1)


call random_number(u1)

if ( u1 .le. accp1)then
alpha(i+1) = y11
acc1 = acc1 + 1
else
alpha(i+1) = alpha(i)
acc1 = acc1
end if





! updating alpha1 parameter
! proposal density
abs10  = 0
do while (abs10 .eq. 0)
f11 = rand_normal(alpha1(i), 0.05d0)
!y11 = alpha(i) + rand_uniform(a1, b1)
	if (f11 .gt. 0.0d0) then
	abs10 = 1
	else
	abs10 = 0
	endif
end do

! prior distribution; half normal prior
call half_normal_pdf(f11, sigalpha1, prioralpha21)
call half_normal_pdf(alpha1(i), sigalpha1, prioralpha20)

! posterior distribution for alpha
num10 = likelihood(alpha(i + 1), f11, beta(i), phit(i,:), 100) + log(prioralpha21)
den10 = likelihood(alpha(i + 1), alpha1(i), beta(i), phit(i, :), 100) + log(prioralpha20)

pr10 = exp(num10 - den10)

accp10 = min(1.0d0, pr10)


call random_number(u10)

if ( u10 .le. accp10)then
alpha1(i+1) = f11
acc10 = acc10 + 1
else
alpha1(i+1) = alpha1(i)
acc10 = acc10
end if






! updating beta (delta in the article) parameter

! proposal density
abs2  = 0
do while (abs2 .eq. 0)
       y22 = rand_normal(beta(i), 0.10d0)
	if (y22 .gt. 0.0d0) then
	abs2 = 1
	else
	abs2 = 0
	endif
end do

! prior distribution; half normal

call half_normal_pdf (y22, sigalpha2, priorbeta2)
call half_normal_pdf (beta(i), sigalpha2, priorbeta1)


! posterior distribution for beta 

num2 = likelihood(alpha(i + 1), alpha1(i + 1), y22, phit(i, :), 100) + log(priorbeta2)
den2 = likelihood(alpha(i + 1), alpha1(i + 1), beta(i), phit(i, :), 100) + log(priorbeta1)

pr2 = exp(num2 - den2)

accp2 = min(1.0d0, pr2)


call random_number(u2)

if ( u2 .le. accp2)then
beta(i+1) = y22
acc2 = acc2 + 1
else
beta(i+1) = beta(i)
acc2 = acc2
end if





! updating the lambda parameter 
! proposal density

abs3  = 0
do while (abs3 .eq. 0)
call truncated_normal_ab_sample(lambda(i), 0.10d0, 0.0d0, 1.0d0, 123,  y33)
	if (y33 .gt. 0.0d0 .and. y33 .lt. 1.0d0 ) then
	abs3 = 1
	else
	abs3 = 0
	endif
end do
! multivariate LCAR
qlambda = (lambda(i)*Wstar + (1 - lambda(i))*idty)*tau(i)

do l = 1, a
lpre(l) = sum(qlambda(l, :)* phit(i, :))
end do


propqlambda = (y33*Wstar + (1 - y33)*idty)*tau(i)
do l = 1, a
proppre(l) = sum(propqlambda(l, :)*phit(i,:))
end do
!!

! posterior distribution for lambda

num3 = 0.5*(log(det(propqlambda, a))) - 0.5*dot_product(phit(i, :), proppre ) + (lalpha - 1)*log(y33) + (lbeta - 1)*log(1 - y33)

den3 = 0.5*(log(det(qlambda, a)))  -  0.5*dot_product(phit(i, :), lpre) + (lalpha - 1)*log(lambda(i)) + (lbeta - 1)*log(lambda(i))


pr3 = exp(num3 - den3)

accp3 = min(1.0d0, pr3)

call random_number(u3)

if (u3 .le. accp3)then
lambda(i+1) = y33
acc3 = acc3 + 1
else
lambda(i+1) = lambda(i)
acc3 = acc3
end if





! updating sigma parameter (Gibbs sampler)


slambda = lambda(i + 1)*Wstar + (1 - lambda(i + 1))*idty
do l = 1, a
proslamb(l) = sum(slambda(l, :)*phit(i, :))
end do
taurate = 0.5*dot_product(phit(i, :), proslamb) + gbeta
taushape = dble(a)/2 + galpha
tauscale = 1/taurate
tau(i+1) = rand_gamma(taushape, tauscale)
acc4 = acc4 + 1
sigma(i + 1) = 1/sqrt(tau(i + 1))






! updating the spatial random effects (phi) parameters


do j = 1, a
 m = 0; u = 0.0d0
do k = 1, 8
 if(neighbor (j, k) .gt. 0) then
u = u + phi(neighbor(j,k))
m = m + 1
end if
  end do
  priordenom = 1 - lambda(i + 1) + lambda(i + 1)*dble(m)
  phimean = lambda(i + 1)*u/priordenom
  phivar = 1.0d0/(tau(i + 1)*priordenom)

! Proposal distribution
phis = rand_normal(phi(j), 0.10d0)
phisv = phi
phisv(j) = phis


! prior distribution (LCAR)
num5 = likelihood(alpha(i + 1), alpha1(i + 1), beta(i + 1), phisv, 100)  - (0.5/phivar)*((phis - phimean)**2) 
den5 = likelihood(alpha(i + 1), alpha1(i + 1), beta(i + 1), phi, 100)  - (0.5/phivar)*((phi(j) - phimean)**2) 

pr5 = exp(num5 - den5)
accp5 = min(1.0d0, pr5)

call random_number(u5)

if (u5 .le.  accp5) then
phi(j) = phis

else
phi(j) = phi(j)
end if
if(mod(i, 10) == 0) write(100,*) phit(i,j) ! saving every 10th iteration
end do


phit(i + 1, 1:a) = phi(1:a) 

print 75, i +1, alpha(i+ 1), acc1, alpha1(i + 1), acc10, beta(i+1), acc2, lambda(i+1), acc3, sigma(i+1), acc4, likelihood(alpha(i + 1), alpha1(i + 1), beta(i + 1), phit(i + 1, :), 100) 
75 format(i6, 1x, f14.8, 1x, i6, 1x, f14.8, 1x, i6, f14.8, 1x, i6, 1x, f14.8, 1x, i6, 1x, f14.8, 1x, i6, 1x,  f24.8)

! write the output inside the loop

! saving every 10th iteration of the fixed effects parameters
if (mod(i, 10) == 0) write(10, 200) alpha(i + 1), alpha1(i + 1), beta(i + 1), lambda(i + 1), sigma(i + 1), likelihood(alpha(i + 1), alpha1(i + 1), beta(i + 1), phit(i + 1, :), 100)
200 format(f14.8, 1x, f14.8, 1x, f14.8, 1x, f14.8, 1x, f14.8, 1x, f24.8)

end do  ! w loop

! *****************************
! END: main MCMC
! *****************************

! output MCMV iterates

! **************************************
call cpu_time (finish)

close(unit = 100)

close(unit = 10)
!close(unit = 300)
!END program carepidemic   ! End of the program
! *****************************************************
! Compute the likelihood function 

! *******************************

contains

function likelihood(alpha, alpha1, beta, phi, ind) result(likk)

implicit none

integer, parameter :: a = 16, n = 1570
integer :: t, i, j, k, m ,  infperiod, ltime, htime
integer, intent(in) :: ind
double precision,  dimension(a) :: phi
double precision,  dimension(n) :: x, y, pop
integer, dimension(n) :: area, inftime, rownumber
double precision, intent(in) ::  alpha, beta, alpha1
integer,  allocatable, dimension(:) :: sareas, sareai
integer, dimension(a) :: lga, lga1
integer, dimension(a, 8) :: neighbor
double precision, dimension(n , 6) :: xx
double precision, allocatable, dimension (:, :) :: newdatas, newdatai
double precision :: likk
double precision, dimension(a) :: lik
double precision ::  dist1, dist2,  epsilon, sprob, inprob, likli
epsilon = 0.0001d0  

! **********************************
open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim1.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim2.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim3.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim4.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim5.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim6.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim7.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim8.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim9.txt" )
!open(500, file = "/home/md.mahsin/nsimulation/Highcov1/nsim10.txt" )
!*********************************
do i = 1, n
read(500, *) rownumber(i), x(i), y(i), area(i), inftime(i), pop(i)
xx(i, 1) = dble(rownumber(i))
xx(i, 2) = x(i)
xx(i, 3) = y(i)
xx(i, 4) = dble(area(i))
xx(i, 5) = dble(inftime(i))
xx(i, 6) = pop(i)
end do
close(500)
! *************************************
! input file of neighboring counties
! ************************************
open(50, file = "/home/md.mahsin/nsimulation/Highcov1/neib.txt")
!***************************
do i = 1, a
read(50, *) neighbor(i, 1:8)
end do
close(50)
!************************
! likelihood computation for each area

DO k = 1, a  !! area loop
	lga(k) = size(pack(xx(:, 1), xx(:, 4) .eq. k))
	!lga1(k) = size(pack(xx(:, 1), xx(:, 4) .eq. k))
	lga1(k) = size(pack(xx(:, 1), xx(:, 4) .eq. k .or. xx(:, 4) .eq. neighbor(k, 1) &
      .or. xx(:, 4)  .eq. neighbor(k, 2) .or. xx(:, 4) .eq. neighbor(k, 3) &
      .or. xx(:, 4) .eq. neighbor(k, 4) .or. xx(:, 4) .eq. neighbor(k, 5) & 
      .or. xx(:, 4) .eq. neighbor(k, 6) .or. xx(:, 4) .eq.  neighbor(k, 7) &
      .or. xx(:, 4) .eq. neighbor(k, 8)) )

	allocate(sareas(lga(k)))
	allocate(sareai(lga1(k)))
	allocate(newdatas(lga(k), 6))
	allocate(newdatai(lga1(k), 6))
	sareas = pack(xx(:, 1), xx(:, 4) .eq. k)
	sareai = pack(xx(:, 1), xx(:, 4) .eq. k .or. xx(:, 4) .eq. neighbor(k, 1) & 
	.or. xx(:, 4) .eq. neighbor(k, 2) .or. xx(:, 4) .eq. neighbor(k, 3) &
	.or. xx(:, 4) .eq. neighbor(k, 4) .or.  xx(:, 4) .eq. neighbor(k, 5) &
	.or. xx(:, 4) .eq. neighbor(k, 6) .or. xx(:, 4) .eq. neighbor(k, 7) &
	.or. xx(:, 4) .eq. neighbor(k, 8))
	newdatas = xx(sareas, 1:6)
	newdatai = xx(sareai, 1:6)
	ltime = minval(newdatas(: , 5))
        htime = maxval(newdatas(: , 5)) + 3  - 1 
	
        likli = 0.0d0
	DO t = ltime, htime
	!! susceptible part likelihood

	sprob = 0.0d0
	DO i = 1, lga(k)
	IF (newdatas(i , 5) == 0 .or. newdatas(i , 5) .gt. (t + 1) ) THEN 
	dist1 = 0.0d0
		DO j = 1, lga1(k)
		IF (newdatai(j,  5) .ne. 0) THEN
		IF (newdatai(j , 5) .le. t  .and. (newdatai(j , 5) + 3) .gt. t) THEN
dist1 = dist1 + ((SQRT((newdatas(i, 2) - newdatai(j, 2))**2 + (newdatas(i, 3) - newdatai(j, 3))**2))/500)**(-beta)
ENDIF
ENDIF
ENDDO
sprob = sprob - (exp(alpha + alpha1*newdatas(i, 6) + phi(k)) *dist1 + epsilon)
ENDIF
ENDDO

!! Infectious part likelihood

inprob = 0.0d0

DO i = 1, lga(k)
IF (newdatas(i , 5) ==  (t + 1) ) THEN 
dist2 = 0.0d0
DO j = 1, lga1(k)
IF (newdatai(j, 5) .ne. 0) THEN
IF (newdatai(j , 5) .le. t .and. (newdatai(j, 5) + 3) .gt. t) THEN 
dist2 = dist2 + ((SQRT((newdatas(i, 2) - newdatai(j, 2))**2 + (newdatas(i, 3) - newdatai(j, 3))**2))/500)**(-beta)
END IF
ENDIF
END DO
inprob = inprob + log(1 - exp(-(exp(alpha + alpha1*newdatas(i, 6) + phi(k)) *dist2 + epsilon)))
END IF
END DO ! end of infectious part


likli = likli + sprob + inprob
END DO ! end of likelihood for time loop
lik(k) = likli
deallocate(sareas)
deallocate(sareai)
deallocate(newdatas)
deallocate(newdatai)
END DO ! end area loop


if (ind == 100) then
likk = sum(lik)
else
DO m = 1, a
if (ind == m) then
likk = lik(m)
endif
END DO
end if


END function likelihood


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Determinant of a matrix (Numerical Recipes)     	  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function det(tab,n)

integer:: n
double precision,intent(in)::tab(n,n)
double precision:: d ,det
double precision::tabbis(n,n)
integer::indx(n)
integer::j

tabbis=tab
call LUdecomp(tabbis,n,indx,d)
do j=1,n
d = d*tabbis(j,j)
end do
det = d

end function det

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LU decomposition of a matrix	(Numerical Recipes)	  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LUdecomp(tab, n, indx,d)

integer,parameter:: nmax=100
double precision,parameter:: tiny=1.E-20
integer, intent(in)::n
double precision, intent(inout)::tab(n, n)
double precision, intent(out):: d
integer, intent(out)::indx(n)
double precision vv(n), sum, dum,aamax
integer ::imax,i,j,ii,k

d=1.
do i=1, n
	aamax=0.
	do j=1, n
		if(abs(tab(i, j))>aamax) aamax=abs(tab(i, j))
	end do
	if(aamax==0.)then
	 stop !pause !'ERRONEOUS RESULT DETECTED'
	else
	vv(i)=1./aamax
	end if
end do

do j=1, n
	do i=1, j-1
		sum=tab(i, j)
		do k=1, i-1
			sum=sum-tab(i, k)*tab(k, j)
		end do
		tab(i, j)=sum
	end do
	aamax=0.
	do i=j, n
		sum=tab(i, j)
		do k=1, j-1
			sum=sum-tab(i, k)*tab(k, j)
		end do
		tab(i, j)=sum
		dum=vv(i)*abs(sum)
		if(dum>=aamax) then
			imax=i
			aamax=dum
		end if
	end do
	if(j.NE.imax) then
		do k=1, n
			dum=tab(imax, k)
			tab(imax, k)=tab(j, k)
			tab(j, k)=dum
		end do
		d=-d
		vv(imax)=vv(j)
	end if
	indx(j)=imax
	if(tab(j, j)==0) tab(j, j)=tiny
	if(j.NE.n) then
		dum=1./tab(j, j)
		do i=j+1, n
			tab(i, j)=tab(i, j)*dum
		end do
	end if
end do
return

end subroutine LUdecomp



END program carepidemic   ! End of the program

! ************************************************











































