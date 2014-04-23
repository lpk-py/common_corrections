PROGRAM elcorr_file
! ==============================================
! GOAL: Apply ellipticity correction on a file
!
! INPUT:
! datar: file which will be read with the following format:
! dataw: file which will be written similar to datar
! with two new colomns at the end (tt-tcor, tcor)
! ==============================================

IMPLICIT NONE

character(8) :: phase
character(30) :: datar, dataw
integer :: ios
real :: del, sdep, scolat, slat, az, tt
real :: tcor
logical :: rythr

print *,'Give file to be read:'
read(5,fmt='(a)') datar
inquire(file=datar,exist=rythr)

if(.not. rythr) then
    print *,'ERROR: cannot find: ',datar
    stop
else
    print *,'Reading...',datar
endif  

print *,'Give file to be written:'
read(5,fmt='(a)') dataw

open(3, file=datar)
open(4, file=dataw)


DO

read(3,*,iostat=ios) phase,del,sdep,slat,az,tt


IF (ios<0) THEN
EXIT
END IF

scolat = 90. - slat
call ellip_blnk(phase, del, sdep, scolat, az, tcor)

WRITE(4,fmt=44) &
    phase,del,sdep,slat,az,tt,&
    & tt-tcor,tcor
44 format(a,',',f10.6,',',f10.5,',',f10.6,',',f10.6,',',f10.6,',',f10.6,',',f10.6)

END DO

END PROGRAM elcorr_file
