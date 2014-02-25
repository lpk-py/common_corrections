      ! compile: g7 raytable raytracesubs azdel ellip

      ! WARNING - this version may not be up to date with the latest
      ! versions of the subroutines. Check the call statements.

      program raytable

      ! traces rays - was used to test raydata subroutines but may be
      ! used to compute travel time tables etc or GMT plot files.

      ! output:
      ! raytable.xxxx

      ! the program is most efficient if the data are
      ! sorted in segments of lines with equal source depths

      dimension rseg(20),ktseg(20),kdwn(20)                     ! ray segments
      dimension ystart(4),y(4,9000),ytbl(4,4000),jsg(9000)  ! ray vector

      logical abrt
      character*72 fname,tablef,chray,dataf,datacomment
      character*8 stationcode,phase

      ! background model
      common /mod/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/

      jdebug=1

      ! set some constants
      h=15.     ! integration step size (must be equal that in mkrtable)
      itmax=20  ! max nr of iterations for ray "shooting"

      print *,'Give model file:'
      read(5,fmt='(a)') fname
      call model(fname)

      write(6,*) 'Radius of Earth model:',r(n)
      write(6,*) 'Radius of CMB:',r(noc)
      write(6,*) 'Radius of ICB:',r(nic)
      write(6,*) 'Radius of Moho:',r(nmoh)

      ! open the raw data file for input
      print *,'Give ray file with extra header or data file name:'
      read(5,fmt='(a)') dataf
      open(3,file=dataf)

      read(3,*)         ! skip first line (allows input files from raydata.f)

      ! read ray info from the raw data file
      call rdray(3,0,chray,phase,rseg,ktseg,kdwn,nseg)
      open(10,file='tdelta.xy')
      open(11,file='pdelta.xy')
      open(12,file='rdelta.xy')
      open(14,file='adelta.xy')

      ! compute table of i vs phi (angle vs. distance)
      ang1=90.0+(kdwn(1)-0.5)*0.2
      ang2=179.99
      if(kdwn(1).eq.0) then
        ang2=ang1
        ang1=0.
      endif  
      rseg(1)=r(n)-sdep
      rsrc=rseg(1)
      call velo1(rsrc,iqsrc,1-kdwn(1),ktseg(1),csrc,dc,d2c)
      call mkrtable(rseg,ktseg,kdwn,nseg,ang1,ang2,ytbl,ntbl)
      do i=1,ntbl
        delta=ytbl(2,i)
        ang=ytbl(1,i)/r2d
        p=rsrc*sin(ang)/csrc
        t=ytbl(3,i)
        rturn=ytbl(4,i)
        write(10,99) delta,t
        write(11,99) delta,p
        write(12,99) delta,rturn
        write(14,99) delta,ytbl(1,i)
99      format(2f10.2)        
      enddo

      end
