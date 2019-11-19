! Simulation of epidemics based on discrete time region-retricted and global GD-ILMs under SIR compartmental framework 
! Author: Md Mahsin (md.mahsin@ucalgary.ca)
! Last update: May 18, 2018 
! How to run it: gfortran ndata_sim.f90 distributions.f95 kind.f95 random.f90
! True parameter: alpha = 0.30, alpha1 = 0.40, beta = 4.0, lambda = 0.80, sigma = 0.60


program dsimulation
use RDistributions
use mo_kind
use random
USE ieee_arithmetic


implicit none

! total area is 16 for Calgary Metropolitan 

integer, parameter :: a = 16, n = 1570
integer :: T, i, j, k, m , tmax, infperiod
double precision, dimension(a) :: phi
double precision, dimension(n) :: x, y, pop
integer, dimension(n) :: area, inftime
double precision :: alpha, alpha1, beta, epsilon, dist, prob, u
integer, allocatable, dimension(:) :: sareas, sareai
integer, dimension(a) :: lga, lga1
integer, dimension(a, 8) :: neighbor
double precision, dimension(n , 6) :: xx
double precision, allocatable, dimension (:, :) :: newdatas, newdatai

!! Set parameter values
tmax = 20
alpha = 0.30d0
beta = 4.0d0
alpha1 = 0.40d0
epsilon = 0.0d0
infperiod = 3

!call car_data(x, y, area, inftime)
!call car(phi)

call init_random_seed()
!*****************************
! Starting with 9 infected individuals 
open(10, file = "/Users/mahsin/Documents/OneDriveM/Research/nsimulation/calgarynda.txt")
! Only one initial infected
!open(10, file = "/Users/mahsin/Documents/OneDriveM/Research/nsimulation/calgarynnda.txt")
do i = 1, n
read(10, *) x(i), y(i), area(i), inftime(i), pop(i)

xx(i, 1) = dble(i)
xx(i, 2) = x(i)
xx(i, 3) = y(i)
xx(i, 4) = dble(area(i))
xx(i, 5) = dble(inftime(i))
xx(i, 6) = dble(pop(i))
end do
close(10)

!**********************
! Spatial random effects for strong correlation
!open(100, file = "/Users/mahsin/Documents/OneDriveM/Research/nsimulation/nphic.txt")

! Spatial random effects for medium correlation
!open(100, file = "/Users/mahsin/Documents/OneDriveM/Research/nsimulation/mphic.txt")

! Spatial random effects for weak correlation
open(100, file = "/Users/mahsin/Documents/OneDriveM/Research/nsimulation/lphic.txt")
do k = 1, a
read(100, *) phi(k)
end do
close(100)
!******************************
! ************************************
! First-order neighborhood structure
open(50, file = "/Users/mahsin/Documents/OneDriveM/Research/nsimulation/neib.txt")
do i = 1, a
read(50, *) neighbor(i, 1:8)
end do
close(50)
!*************************

!! Generate the simulated epidemics from region-restricted and global models

DO T = 1, tmax
	DO k = 1, a  !! area loop
	lga(k) = size(pack(xx(:, 1), xx(:, 4) .eq. k))
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
		DO i = 1, lga(k)  !! loop for suscetibles
		IF (newdatas(i , 5) == 0 ) THEN
		dist = 0.0
			DO j = 1, lga1(k) !! loop for infectious individuals for region-restricted model
			!DO j = 1, n !! loop for infectious individuals for global model
			IF (newdatai(j, 5) .NE. 0 ) THEN
			!IF (xx(j, 5) .NE. 0 ) THEN !! 
			IF (newdatai(j, 5) .LE. T .AND. (newdatai(j, 5) + infperiod) .GT. T ) THEN
			!IF (xx(j, 5) .LE. T .AND. (xx(j, 5) + infperiod) .GT. T ) THEN
			dist = dist + ((SQRT((newdatas(i, 2) - newdatai(j, 2))**2 + (newdatas(i, 3) - newdatai(j, 3))**2))/500)**(-beta)
			!dist = dist + ((SQRT((newdatas(i, 2) - xx(j, 2))**2 + (newdatas(i, 3) - xx(j, 3))**2))/500)**(-beta)
			ENDIF
			ENDIF
			ENDDO
		prob = 1 - exp(-(exp(alpha + alpha1*newdatas(i, 6) + phi(k)) *dist + epsilon))
		
		call random_number(u)
		IF(prob .GE. u)THEN
		newdatas(i, 5)= T + 1
		ENDIF
		ENDIF
		ENDDO
	
	
	xx(sareas, :) = newdatas(:, 1:6)
	deallocate(sareas)
	deallocate(sareai)
	deallocate(newdatas)
	deallocate(newdatai)

ENDDO
ENDDO


!!! Writing and saving the data

open(200, file = "/Users/mahsin/Documents/OneDriveM/Research/nsimulation/usim1.txt" )
DO m = 1, n
write(200, *) xx(m, 1), xx(m, 2), xx(m, 3), xx(m, 4), xx(m, 5), xx(m, 6)
ENDDO
close(200)

end program dsimulation







 





