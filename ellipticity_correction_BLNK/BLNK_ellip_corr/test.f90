PROGRAM TEST

character(8) :: phase="Pdiff"
real :: edist, edepth, elat, ecolat, eazim_deg, eazim, azim, tcor
logical :: abrt
dimension usrc(2)

!degree ---> radian
degrad = 45.0/atan(1.0)

! open ellipticity correction file
open(21,file='elcordir.tbl',access='direct',form='formatted',recl=80)

!ecolat: co-latitude of the souece (in radians) ---> (90.0-latitude)*pi/180.
elat = 0.
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
