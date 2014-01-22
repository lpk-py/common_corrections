PROGRAM TEST

!=========================================
! Program name: TEST
! Goal: run and test BLNK ellipticity
!       correction code
!
! Dependencies:
! elcordir.tbl  ellip.f  libsun.f  libtau.f
! ttlim.inc 
!=========================================

IMPLICIT NONE

character(8) :: phase="Pdiff"
real :: degrad, edist, edepth, elat, ecolat 
real :: eazim_deg, eazim, azim, tcor, zs
real, dimension(2) :: usrc
logical :: abrt

!degree ---> radian
degrad = 45.0/atan(1.0)

! open ellipticity correction file
open(21,file='elcordir.tbl',access='direct',form='formatted',recl=80)

elat = 0.
!ecolat: co-latitude of the souece (in radians) ---> (90.0-latitude)*pi/180.
ecolat = (90.0 - elat)/degrad
call ellref(ecolat)

!depth of event
zs = 10.0
call depset(zs,usrc)
edepth = zs

!epicentral distance to station (in degrees)
edist = 120.0 
!azimuth from source to station (in radians)
eazim_deg = 90.0
azim = eazim_deg/degrad

WRITE(*,*) "phase: ", phase
WRITE(*,*) "edist: ", edist
WRITE(*,*) "depth: ", edepth
WRITE(*,*) "ecolat: ", ecolat
WRITE(*,*) "azim: ", azim

call ellcor(phase, edist, edepth, ecolat, azim, tcor, abrt)

WRITE(*,*) "corrected time: ", tcor

END PROGRAM TEST
