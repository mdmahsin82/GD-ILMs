! gfortran nmahsin-cont-SEIR.f90
program datsimSEIR
USE ieee_arithmetic
implicit none

integer, parameter:: n = 1570, m =16, nsuspar = 2, observednum = 1, narea= 16
integer, dimension(n):: arealevel
double precision, dimension(n):: x, y, effect
double precision, dimension(m):: phi
integer, dimension(m, m):: nei
double precision, dimension(n,n):: dd
double precision:: spark, deltain1, deltain2, deltanr1, deltanr2, tmax, kernelpar
double precision, dimension(nsuspar):: suspar
double precision, dimension(observednum, 9):: observedepi
double precision, dimension(n,9):: epidat
integer::i, ninfected, temp
double precision::start,finish,Inf,likk
double precision, dimension(n) :: covarea

IF (ieee_support_inf(Inf)) THEN
  Inf = ieee_value(Inf,  ieee_positive_inf)
END IF

open(232,  file = "/Users/mahsin/Documents/OneDriveM/Research/SEIR/usim28-mahsin-SEIR.txt" )

call calgaryned(x, y, arealevel, covarea)
call nphic(phi)
call neighb(nei)

tmax = 1000.0d0
spark = 0.0d0
deltain1 = 100.0d0
deltain2 = 100.0d0
deltanr1 = 12.0d0
deltanr2 = 36.0d0

suspar = (/0.3d0, 0.3d0/)
kernelpar = 2.0d0

observedepi(1, :) =(/ (0.0d0, i= 1, 9) /)
!observedepi(1,:) = (/ 189.0d0 , 3.0d0, 3.0d0, 0.0d0 /)
temp = 0!523782

call cpu_time(start)


call datasimulationsinr(n, observednum, observedepi, tmax, temp, suspar, &
& kernelpar, spark, deltain1, deltain2, deltanr1, deltanr2, &
& x, y, phi, covarea, arealevel, nei, epidat)

call cpu_time(finish)
print '("Time = ",f14.8," seconds.")',finish-start

ninfected = count(epidat(:,2) .ne. inf)
print*, ninfected

do i = 1,ninfected
!print*, i, epidat(i,:)
end do

do i = 1,n
write(232,*) epidat(i,:)
end do


contains


subroutine calgaryned(x,y,arealevel,effect)
implicit none
integer,parameter:: n=1570
integer:: m
integer, intent(out), dimension(n):: arealevel
integer, dimension(n):: inft,id
double precision,intent(out), dimension(n):: x, y, effect
open(555, file = "/Users/mahsin/Documents/OneDriveM/Research/SEIR/fclaimdata3.txt")
do m = 1 , n
read(555,*) id(m), x(m), y(m), arealevel(m), inft(m), effect(m)
end do
close(555)
end subroutine calgaryned

subroutine nphic(phi)
implicit none
integer,parameter:: n=16
integer:: m
double precision,intent(out), dimension(n):: phi
open(5556, file = "/Users/mahsin/Documents/OneDriveM/Research/SEIR/nphic.txt")
do m = 1 , n
read(5556,*) phi(m)
end do
close(5556)
end subroutine nphic

subroutine neighb(nei)
implicit none
integer,parameter:: n=16
integer:: m
integer,intent(out), dimension(n,n):: nei
open(5556, file = "/Users/mahsin/Documents/OneDriveM/Research/SEIR/nei.txt")
do m = 1 , n
read(5556,*) nei(m,:)
end do
close(5556)
end subroutine neighb



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 			EPIDEMIC SIMULATION subroutine		 	 	 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine datasimulationsinr(n, observednum, observedepi, tmax, temp, suspar, &
& kernelpar, spark, deltain1, deltain2, deltanr1, deltanr2, &
& x, y, phi, covarea, arealevel, neighb, epidat)

USE ieee_arithmetic

implicit none

integer, intent(in) :: n, observednum, temp   			! integers
double precision, intent(in), dimension(16) :: phi
double precision, intent(in), dimension(n) :: covarea
integer, intent(in), dimension(n) :: arealevel      						! susceptibility covariates
integer, intent(in), dimension(16, 16) :: neighb  							! transmissibility covariates
double precision, intent(in), dimension(n) :: x, y                    		! network & distance matrices
double precision, intent(in), dimension(2) :: suspar        			! susceptibility parameters
double precision, intent(in) :: spark, tmax                     			! spark& notification effect& max infec. time
double precision, intent(in) :: deltain1, deltain2, deltanr1, deltanr2  	! Parameters of the infectious period distribution
double precision, intent(in) :: kernelpar               					! parameter of the kernel function
double precision, intent(in), dimension(observednum, 9) :: observedepi 		! observed epidemic to start
double precision, dimension(observednum, 9) :: observedepi1 				! observed epidemic to start
double precision, dimension(n, 9), intent(out) :: epidat               		! OUTPUT
integer :: nnn1                                                           	! # of infected by the end of epidemic
integer, dimension(n, 2) :: xx                               				! Auxiliary variable
double precision, dimension(1, 2) :: ts                               		! OUTPUT from the rate subroutine
double precision :: t0                                              		! current infection time during the simulation
double precision :: Inf, u                                             		! defining Infinity
integer :: ctr, i, j, sdg, mg
integer, allocatable, dimension(:) :: mmg

if (temp .ne. 0) then
	call initrandomseedsinr(temp)
else
	call initrandomseed2()
end if

IF (ieee_support_inf(Inf)) THEN
	Inf = ieee_value(Inf,  ieee_positive_inf)
END IF


! defining auxiliary variable (xx) for the status of each individual:
! 0 : susceptible
! 1 : exposed
! 2 : infectious
! 3 : removed

	xx          = 0
	xx(:, 1)     = (/(j, j=1, n)/)
	epidat      = 0.0d0
	observedepi1 = observedepi
! initial observed epidemic:
	if (observednum .eq. 1) then
		if (observedepi1(1,1) .eq. 0) then
			call random_number(u)
			observedepi1(1,1) = int(u*n) + 1
		end if
	end if

	do j = 1,  observednum
		if (observedepi1(j, 2) .eq. 0.0d0) then
			epidat(j, 1)  = observedepi1(j, 1)
			epidat(j, 6)  = observedepi1(j, 6)
			epidat(j, 5)  = randgamma22(deltain1, 1.0d0/deltain2)
			epidat(j, 4)  = epidat(j, 5) + epidat(j, 6)
			epidat(j, 3)  = randgamma22(deltanr1, 1.0d0/deltanr2)
			epidat(j, 2)  = epidat(j, 4) + epidat(j, 3)
			epidat(1, 7)  = dble(arealevel(int(epidat(j, 1))))
			epidat(1, 8)  = covarea(int(epidat(j, 1)))
			epidat(1, 9)  = dble(phi(arealevel(int(epidat(j, 1)))))
		else
			epidat(j, :)  = observedepi1(j, :)
		end if
	end do

! the current infection time to start from:
	t0 = epidat(observednum, 6)
	xx(int(epidat(1:observednum, 1)), 2) = 1

	mg  = size(pack(int(epidat(1:(observednum-1), 1)), epidat(1:(observednum-1), 4) .lt. epidat(observednum, 6) ))
	if (mg .gt. 0) then
		allocate(mmg(mg))
		mmg = pack(int(epidat(1:(observednum-1), 1)), epidat(1:(observednum-1), 4) .lt. &
					& epidat(observednum, 6) )
		xx(mmg, 2) = 2
		deallocate(mmg)
	end if

	mg  = size(pack(int(epidat(1:(observednum-1), 1)), epidat(1:(observednum-1), 2) .lt. &
			& epidat(observednum, 6) ))
	if (mg .gt. 0) then
		allocate(mmg(mg))
		mmg = pack(int(epidat(1:(observednum-1), 1)), epidat(1:(observednum-1), 2) .lt. &
					 & epidat(observednum, 6) )
		xx(mmg, 2) = 3
		deallocate(mmg)
	end if

! starting simulating the epidemic
	ctr = observednum
	do while( (ctr .le. n) )
		ctr = ctr + 1
! to provide the next infected individual with minmum waiting time to infection:
		call rateSINR(n, suspar, kernelpar, spark, xx, x, y, phi, covarea, arealevel, neighb, ts)

		! to judge stopping the epidemic or keep generating:
		if ( (ts(1, 2) .ne. Inf) .and. (ts(1, 1) .ne. 0.0d0) ) then
			ts = ts
		else
			where(xx(:, 2) .eq. 1) xx(:, 2) = 2
			exit
		end if

!making sure there is still infectious individuals that can transmit the disease

		sdg = 0
		do i = 1,  (ctr-1)
			if ( (epidat(i, 2) .gt. (ts(1, 2)+t0)) .and. &
          & neighb(arealevel(int(epidat(i, 1))),arealevel(int(ts(1,1)))).eq. 1 ) then
				sdg = sdg +1
			else
				sdg = sdg
			end if
		end do

! assigning infection time,  incubation period,  notification time,  delay period and
! removal time for the newly infected:

		if (sdg .eq. 0 ) then
			where(xx(:, 2) .eq. 1) xx(:, 2) = 2
			exit
		else
			epidat(ctr, 6) = ts(1, 2) + t0
			epidat(ctr, 5) = randgamma22(deltain1, 1.0d0/deltain2)
			epidat(ctr, 4) = epidat(ctr, 5) + epidat(ctr, 6)
			epidat(ctr, 3) = randgamma22(deltanr1, 1.0d0/deltanr2)
			epidat(ctr, 2) = epidat(ctr, 3) + epidat(ctr, 4)
			epidat(ctr, 1) = ts(1, 1)
			t0 = epidat(ctr, 6)
			xx(int(epidat(ctr, 1)), 2) = 1
			epidat(ctr, 7)  = dble(arealevel(int(epidat(ctr, 1))))
			epidat(ctr, 8)  = covarea(int(epidat(ctr, 1)))
			epidat(ctr, 9)  = dble(phi(arealevel(int(epidat(ctr, 1)))))
		end if

		if ( (epidat(ctr, 4) .gt. tmax) ) then
			epidat(ctr, 2) = 0.0d0
			epidat(ctr, 3) = 0.0d0
			epidat(ctr, 4) = 0.0d0
			epidat(ctr, 5) = 0.0d0
			epidat(ctr, 6) = 0.0d0
			epidat(ctr, 7) = 0.0d0
			epidat(ctr, 8) = 0.0d0
			epidat(ctr, 9) = 0.0d0
			exit
		end if

! update the auxiliary variable of the status of individuals:
		mg  = size(pack(int(epidat(1:ctr-1, 1)), epidat(1:ctr-1, 4) .lt. epidat(ctr, 6) ))
		if (mg .gt. 0) then
			allocate(mmg(mg))
			mmg = pack(int(epidat(1:ctr-1, 1)), epidat(1:ctr-1, 4) .lt. epidat(ctr, 6) )
			xx(mmg, 2) = 2
			deallocate(mmg)
		end if

		mg  = size(pack(int(epidat(1:ctr-1, 1)), epidat(1:ctr-1, 2) .lt. epidat(ctr, 6) ))
		if (mg .gt. 0) then
			allocate(mmg(mg))
			mmg = pack(int(epidat(1:ctr-1, 1)), epidat(1:ctr-1, 2) .lt. epidat(ctr, 6) )
			xx(mmg, 2) = 3
			deallocate(mmg)
		end if

	end do

! assigning infinity values for those uninfected by the end of the epidemic

	nnn1 = count(epidat(:, 2)  .ne. 0.0d0)
	do i = (nnn1+1),  n
		do j = 1, n
			if (all(int(epidat(1:(i-1), 1)) .ne. j)) then
				epidat(i, 1) = dble(j)
			end if
		end do
		epidat(i, 2) = Inf
		epidat(i, 3) = 0.0d0
		epidat(i, 4) = Inf
		epidat(i, 5) = 0.0d0
		epidat(i, 6) = Inf
        epidat(i, 7)  = dble(arealevel(int(epidat(i, 1))))
        epidat(i, 8)  = covarea(int(epidat(i, 1)))
        epidat(i, 9)  = dble(phi(arealevel(int(epidat(i, 1)))))
	end do



end subroutine datasimulationsinr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 			INFECTIVITY rateSINR subroutine 		 	 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rateSINR(n, suspar, kernelpar, spark, xx, x, y, phi, covarea, arealevel, neighb, mms)

USE ieee_arithmetic

implicit none

integer ::i, mg, mg1, mg2, j, m

integer, intent(in) :: n                   						   !integers
double precision, intent(in), dimension(16) :: phi
double precision, intent(in), dimension(n) :: covarea
integer, intent(in), dimension(n) :: arealevel                     ! susceptibility covariates
integer, intent(in), dimension(16, 16) :: neighb                   ! transmissibility covariates
integer, intent(in), dimension(n, 2) :: xx
double precision, intent(in), dimension(n) :: x,y                  ! network and distance matrices
double precision, intent(in), dimension(nsuspar) :: suspar         ! susceptibility parameters
double precision, intent(in) :: spark		                       ! spark & notification effec parameters
double precision, intent(in) :: kernelpar            			   ! parameters of the kernel function
double precision :: Inf, disdis, disdis1
double precision, dimension(1, 2) :: mms

integer, allocatable, dimension(:) :: mmg, mmg1, mmg2
double precision, allocatable, dimension(:, :) :: rr

IF (ieee_support_inf(Inf)) THEN
Inf = ieee_value(Inf,  ieee_positive_inf)
END IF


! Calculating the infectivity rateSINR of distance-based ILM with spark term with
! "powerlaw" kernel.

! defining infectious individuals:
        mg  = size( pack(xx(:, 1), xx(:, 2).eq. 1.0d0 ) )
        allocate(mmg(mg))
        mmg = pack(xx(:, 1), xx(:, 2).eq. 1.0d0 )

! defining susceptible individuals:
        mg1 = size( pack(xx(:, 1), xx(:, 2).eq. 0.0d0 ) )
        allocate(mmg1(mg1))
        mmg1 = pack(xx(:, 1), xx(:, 2).eq. 0.0d0 )

! defining notified individuals:
        mg2 = size( pack(xx(:, 1), xx(:, 2).eq. 2.0d0 ) )
        allocate(mmg2(mg2))
        mmg2 = pack(xx(:, 1), xx(:, 2).eq. 2.0d0 )

! declaring a variable with size equal to the number of susceptible individuals
! the first column is the id number of the susceptible individuals
        allocate(rr(mg1, 2))
        rr(:, 1) = mmg1

! start calculating the infectivity rate for each susceptible individuals:
	do i = 1, mg1

		disdis = 0.0d0
		do j = 1,  mg
		  if(neighb(arealevel(mmg1(i)),arealevel(mmg(j))).eq. 1)then
			disdis = disdis + ((SQRT((((x(mmg1(i)) - x(mmg(j)))**2) + &
			& ((y(mmg1(i)) - y(mmg(j)))**2)))/dble(500))**(-kernelpar))
		  else
			disdis = disdis
		  end if
		end do

		disdis1 = 0.0d0
		do j = 1,  mg2
		  if(neighb(arealevel(mmg1(i)),arealevel(mmg2(j))).eq. 1)then
			disdis1 = disdis1 + ((SQRT((((x(mmg1(i)) - x(mmg2(j)))**2) + &
			& ((y(mmg1(i)) - y(mmg2(j)))**2)))/dble(500))**(-kernelpar))
		  else
			disdis1 = disdis1
		  end if
		end do

        rr(i, 2) = (exp(suspar(1)+(suspar(2)*covarea(mmg1(i)))+phi(arealevel(mmg1(i))))*(disdis + disdis1)) + spark

	end do

! assigning waiting time to infection for each susceptible individual:
	do i = 1, mg1
		if (rr(i, 2) .eq. 0.0d0) then
			rr(i, 2) = Inf
		else
			rr(i, 2) = randgamma22(1.0d0, 1.0d0/rr(i, 2))
		end if
	end do

! choose the one with minmum waiting time as the newly infected individual:
	if (all(rr(:, 2) .eq. Inf).eqv. .true.) then
		mms(1, 1)  = 0.0d0
		mms(1, 2)  = Inf
	else
		mms(1, 1:2)  = rr(int(minloc(rr(:, 2), 1)), :)
	end if

	deallocate(rr)
	deallocate(mmg2)
	deallocate(mmg1)
	deallocate(mmg)

end subroutine rateSINR

    subroutine initrandomseedsinr(temp)
    implicit none
    integer :: n
    integer, intent(in):: temp
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))
    seed = temp
    call random_seed(PUT = seed)
    deallocate(seed)

    end subroutine initrandomseedsinr

    subroutine initrandomseed2()
!The seed for the random number generation method random_number() has been reset

    implicit none

    integer :: i
    integer :: n
    integer :: clock
    integer,  dimension(:),  allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(COUNT=clock)

        seed = clock + 37 * (/ (i - 1,  i = 1,  n) /)
        call random_seed(PUT = seed)

        deallocate(seed)
    end subroutine initrandomseed2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Generating random variables for diffierent distributions !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!####################  NORMAL distribution ######################

    FUNCTION randnormal22(mean, stdev) RESULT(c)

    implicit none

    double precision :: mean, stdev, c, temp(2), r, theta
    double precision,  PARAMETER :: PI=3.141592653589793238462d0

        CALL RANDOM_NUMBER(temp)
        r=(-2.0d0*log(temp(1)))**0.5d0
        theta = 2.0d0*PI*temp(2)
        c= mean+stdev*r*sin(theta)

    END FUNCTION randnormal22

!#################### GAMMA distribution ######################

    RECURSIVE FUNCTION randgamma22(shape,  SCALE) RESULT(ans)
    double precision :: SHAPE, scale, u, w, d, c, x, xsq, g, ans, v

! DESCRIPTION: Implementation based on "A Simple Method for Generating Gamma Variables"
! by George Marsaglia and Wai Wan Tsang.
! ACM Transactions on Mathematical Software and released in public domain.
! ## Vol 26,  No 3,  September 2000,  pages 363-372.

        IF (shape >= 1.0d0) THEN
            d = SHAPE - (1.0d0/3.0d0)
            c = 1.0d0/((9.0d0 * d)**0.5)
            DO while (.true.)
                x = randnormal22(0.0d0,  1.0d0)
                v = 1.0 + c*x
                DO while (v <= 0.0d0)
                    x = randnormal22(0.0d0,  1.0d0)
                    v = 1.0d0 + c*x
                END DO
                v = v*v*v
                CALL RANDOM_NUMBER(u)
                xsq = x*x
                IF ((u < 1.0d0 -.0331d0*xsq*xsq) .OR.  &
                (log(u) < 0.5d0*xsq + d*(1.0d0 - v + log(v))) ) then
                    ans=scale*d*v
                    RETURN
                END IF

            END DO
        ELSE
            g = randgamma22(shape+1.0d0,  1.0d0)
            CALL RANDOM_NUMBER(w)
            ans=scale*g*(w**(1.0d0/shape))
            RETURN
        END IF

    END FUNCTION randgamma22

end program datsimSEIR
