SUBROUTINE ellip_blnk(phase_rd, del_rd, sdep_rd, scolat_rd, az_rd, tcor)

IMPLICIT NONE

character(8), INTENT(IN) :: phase_rd
real, INTENT(IN) :: del_rd, sdep_rd, scolat_rd, az_rd
real, INTENT(OUT) :: tcor
real :: degrad, edist, eazim_deg, azim, edepth, elat, ecolat, zs
real, dimension(2) :: usrc
logical :: abrt

!degree ---> radian
degrad = 45.0/atan(1.0)

! open ellipticity correction file
open(21,file='elcordir.tbl',access='direct',form='formatted',recl=80)

!epicentral distance to station (in degrees)
edist = del_rd

!azimuth from source to station (in radians)
eazim_deg = az_rd
azim = eazim_deg/degrad

!depth of event
zs = sdep_rd
call depset(zs,usrc)
edepth = zs

!ecolat: co-latitude of the souece (in radians) ---> (90.0-latitude)*pi/180.
ecolat = scolat_rd/degrad
call ellref(ecolat)

call ellcor(phase_rd, edist, edepth, ecolat, azim, tcor, abrt)
!WRITE(*,*) "corrected time: ", tcor

END SUBROUTINE ellip_blnk
