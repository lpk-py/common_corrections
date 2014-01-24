      ! compile: gfortran -fdefault-real-8 -o raydyntrace raydyntrace.f 
      !          raytracesubs.f azdel.f refl.f90

      program raydyntrace

      ! This is version 2.1, 17 feb 2013

      ! Testing and documentation:


      ! Y. Tian, S.-H. Hung, G. Nolet, R. Montelli, F.A. Dahlen, Dynamic ray 
      ! tracing and traveltime corrections for global seismic tomography, 
      ! J. Comp. Phys., 226:672-687, 2007.

      ! Y. Tian, R. Montelli, G. Nolet, F.A. Dahlen, Computing traveltime and 
      ! amplitude sensitivity kernels in finite-frequency tomography, 
      ! J. Comp. Phys., 226:2271-2288, 2007.

      ! See file Raydyntrace.manual in directory Documentation for users instructions

      ! Background reading:
      ! Nolet, G., A Breviary of Seismic Tomography, Cambridge Univ. Press, 2008.
      ! ISBN 978-0-521-88244-6

      !
      ! August 2006: removed relev from elevation corrections when kd<5
      !              added stop for sources above local surface
      !              corrected ptlat/lon for surface source in PPP etc.
      !              changed azimuth in givloc call for deltas > 180
      ! Sept 2006: Depth of source shallower than CRUST2.0 surface is replaced 
      !            with CRUST2.0 elev; a list of shallow sources is in src.shallow
      !            Output warning message of imprecise hessian against analytical expression
      ! Dec 2007: Replaced Kennetts ellipticity routine with getelcor. 
      ! March 2008: added table option and fixed bug for ellipticity corrections of 
      !           ghost phases like pP etc (in raytracesubs.f)
      ! March 2009: added option to compute tables with p,T,Delta and Rrs only.
      ! Nov 2009: added travel time to ray nodes in all output (i.e. also if option no ,
      !           table), changed name of common /mod/ to /modl/, closed CRUST2.0 files
      !           in rcrust2
      ! Jan 2013: fixed a number of issues with crustal corrections in oceanic
      !           areas, including island stations not recognized in CRUST2.0
      !           Removed list of shallow sources is in src.shallow - sources
      !           in the air are now recogizable by a telev>0 in first ray
      !           segment.

      

      ! input if tableflag >0 or <0:
      ! [1] model file name
      ! [2] maxarr, tdiffmax  (nr of arrivals, largest differential time)
      ! [3] output file ident
      ! [4] tableflag
      ! [5] source radius (6371 - source depth)
      ! [6] ray definition file
      !     ... if tablflag <0 continue with more ray def files, if any
      ! [7] stop
       
      ! input if tableflag=0:
      ! [1] model file name
      ! [2] maxarr, tdiffmax  (nr of arrivals, largest differential time)
      ! [3] data file ident
      ! [4] tableflag (=0)
      ! make sure the CRUST2.0 files are available in the directory
      ! the data file has lines with:
      ! idate,slat,slon,sdep,rlat,rlon,relev,stationcode,kpole
      ! where idate is the event date (eg 081223), slat,slon,sdep are
      ! source latitude, longitude (in deg) and depth (in km) and
      ! rlat,rlon,relev the receiver lon, lat and elevation (km). kpole
      ! should be 1 if the distance exceeds 180 deg, otherwise 0.


      ! output:
      ! raydyntrace.xxxx
      ! if model input is MKS (m/s etc), length unit is converted to km.

      ! if tableflag>0, the output file raydyntrace.xxxx has a number of rays, with at each
      ! ray node (nomenclature follows Nolets book):
      ! r (km), i (rad), Delta (rad), T (s), c (km/s), q0=1/Qs, h11, h22 (s/km^2), Rxs c(km), s (km).
      ! and a file out.raydyntrace.xxxx with one entry for each ray.

      parameter (ityp=360,nla=90,nlo=180)                        ! for rcrust2 arrays
      parameter (NFR=10)                                        ! nr of freq bands

      dimension rseg(20),ktseg(20),kdwn(20)                     ! ray segments
      parameter (LRAY=20000)
      dimension ystart(4),y(4,LRAY),ytbl(4,4000),jsg(LRAY)  ! ray vector
      dimension a1(10),a2(10),d1(10),d2(10),trtime(10),indx(10)   ! ray times/dist
      dimension hmf(2,LRAY),hmb(2,LRAY),rayvel(LRAY),rayq(LRAY) ! Hessian
      dimension ypq(10,LRAY)                                    ! idem
      dimension rcr2(16),vcr2(16)                               ! CRUST2.0
      dimension rref(50),vref(50)                               ! Reference crust
      complex zz,rps(3,3),tps(3,3)              ! refl/transm coeff

      integer tableflag

      logical rythr
      character*72 fname,tablef,chray
      character*30 dataf
      character*16 stationcode
      character*8 phase
      character*2 ctype(ityp),types,atype(nlo)
      character*1 yn,star
      character*3 comp

      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      ! the following are filled by subroutine rcrust2:
      common/crust2/amapvp(8,nlo,nla),amaprho(8,nlo,nla),
     &              amapvs(8,nlo,nla),amapthi(7,nlo,nla),
     &              amapele(nlo,nla)                    ! model CRUST2.0
      common/glbgrd/types(nlo,nla)                      ! model CRUST2.0

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/

      jdebug=0

      ! intialize CMB and ICB to 0
      nic=0
      noc=0

      ! set some constants
      h=10.     ! integration step size (must be equal that in mkrtable)
      itmax=20  ! max nr of iterations for ray shooting

      print *,'Give model file:'
      read(5,fmt='(a)') fname
      inquire(file=fname,exist=rythr)
      if(.not. rythr) then
        print *,'ERROR: cannot find: ',fname
        stop
      endif  
      call model(fname)
!     open(15,file='spline.'//fname)
!     write(15,*) '   i         r        Vp     coefficients'
!     do i=1,n
!       write(15,fmt='(i5,2f10.3,3e14.5)') i,r(i),vp(i),qvp(1,i),
!    &        qvp(2,i),qvp(3,i)
!     enddo
!     write(15,*) '   i         r        Vs     coefficients'
!     do i=1,n
!       write(15,fmt='(i5,2f10.3,3e14.5)') i,r(i),vs(i),qvs(1,i),
!    &        qvs(2,i),qvs(3,i)
!     enddo
!     write(15,*) '   i         r         Qs    coefficients'
!     do i=1,n
!       write(15,fmt='(i5,2f10.3,3e14.5)') i,r(i),Qs(i),qqs(1,i),
!    &        qqs(2,i),qqs(3,i)
!     enddo
!     close(15)

      write(6,*) 'Radius of Earth model:',r(n)
      if(abs(r(n)-6371.).gt.0.1) write(6,*) 'WARNING! Model radius',
     & ' deviates from 6371 - crustal corrections are not OK'
      write(6,*) 'Radius of CMB:',r(noc)
      write(6,*) 'Radius of ICB:',r(nic)
      write(6,*) 'Radius of Moho:',r(nmoh)

      print *,'Give max number of arrivals to include in computations'
      print *,'and their max difference (in sec) w.r.t. first arrival:'
      read(5,*) maxarr,tdifmax

      ! open the input file
      print *,'Give input file name IN THIS DIRECTORY (eg Pshallow),'
      print *,'or give output file ident if tableflag != 0:'
      read(5,fmt='(a)') dataf

      ! open the big output file
      tablef='raydyntrace.'//dataf
      print *,'Output is to file ',tablef
      open(2,file=tablef)

      ! open the table output file
      tablef='out.raydyntrace.'//dataf
      print *,'Table output is to file ',tablef
      open(8,file=tablef)

      ! get tableflag. If tableflag=0, the program expects a file with all
      ! source-station pairs and will converge to the exact ray for each pair
      ! and output the dynamic raytracing variables needed to apply eq 108
      ! in Dahlen et al. GJI 141:157-174, 2000 or any of the equations for
      ! various kernels in Nolets Breviary, e.g. 7.28 for traveltime, 8.7
      ! for amplitude focusing, 8.14 for Q.
      ! If tableflag>0, the program spits out a table and leaves it to the
      ! user to interpolate this table to find the necessary ingredients for
      ! kernel computation for various earthquake depth/distances. 
      ! If tableflag<0, the program just computes travel times tables for
      ! multiple ray definition files.

      print *,'Give tableflag (<0 for T-X files only,',
     &      ' >0 to create traveltime table, else 0):'
      read(5,*,iostat=ios) tableflag
      if(ios.ne.0) then
        print *,'Warning: tableflag missing in inputfile. 0 assumed'
        tableflag=0
      endif  

      if(tableflag.ne.0) then
        print *,'Give source radius in km (overrides ray definition):'
        read(5,*) srcrad
        rsrc=srcrad
      endif  

      print *,'Give file with ray definition:'
      read(5,'(a)') fname
      open(4,file=fname)
      iout=2
      if(tableflag.lt.0) iout=0
      if(jdebug>0) write(13,*) 'calling rdray'
      call rdray(4,iout,chray,phase,rseg,ktseg,kdwn,nseg)
      if(tableflag.ne.0) rseg(1)=srcrad
      iqsrc=n           ! initial guess for source layer (at surface)

      if(tableflag.gt.0) then   ! just make table, then quit

        rewind 8
        write(8,'(a,6(/,a))') 'Significance of terms in the table:',
     &    'Rez = real product of refl/transmission coefficients',
     &    'Imz = imaginary product of refl/transmission coefficients',
     &    'Ux = horizontal surface amplification',
     &    'Uz = vertical surface amplification',
     &    'phase = phase for Uz. If /= 0, phase for Ux is (phase-pi/2)',
     &    'Mas = Maslov index.'
        write(8,fmt='(3a)')'     #     N         p      rmin         i', 
     &   '     Delta         T        t*         Rxs',
     &   '     Rez     Imz      Ux      Uz   phase  Mas'
        write(8,fmt='(17x,2a)')'s/rad        km       deg       deg',
     &    '       sec       sec          km'   

        call mkrtable2(rseg,ktseg,kdwn,nseg,tableflag,
     &  hmf,hmb,ypq,rayvel,rayq)

        stop 'End of raydyntrace - it did the tableflag>0 option'

      else if(tableflag.lt.0) then

        do while (fname(1:4).ne.'stop')
          ang1=90.1
          ang2=179.99
          if(kdwn(1).eq.0) then   ! ray starts upward
            ang1=0.1
            ang2=89.99
          endif  
          write(8,fmt='(a)') fname
          open(12,file='XT.'//phase )
          call velo1(rsrc,iqsrc,1-kdwn(1),ktseg(1),csrc,dc,d2c)
          if(jdebug>0) write(13,*) 'rscr, csrc=',rsrc,csrc
          call mkrtable2(rseg,ktseg,kdwn,nseg,tableflag,
     &         hmf,hmb,ypq,rayvel,rayq)
          close(12)

          ! read new ray info 
          print *,'Give file with ray definition or type <stop>:'
          read(5,'(a)',iostat=ios) fname
          if(ios.eq.0.and.fname(1:4).ne.'stop') then
            close(4)
            open(4,file=fname)
            iout=2
            if(tableflag.lt.0) iout=0
            if(jdebug>0) write(13,*) 'calling rdray'
            call rdray(4,iout,chray,phase,rseg,ktseg,kdwn,nseg)
            rseg(1)=srcrad
            iqsrc=n           ! initial guess for source layer (at surface)
          else
            fname='stop'
          endif  
        enddo

        write(8,fmt='(2a)')'     #     N         p      rmin         i', 
     &        '     Delta         T         Rrs'
        write(8,fmt='(17x,2a)')'s/rad        km       deg       deg',
     &    '       sec          km'   
        stop 'End of raydyntrace - it did the tableflag<0 option'

      endif  

      write(8,5)
5     format('    Date N Receiver    Rlat    Rlon    Slat    Slon',
     &    '  Sdep     Time  Elcor     t*       p')

      ! read CRUST2.0 data
      call rcrust2

      inquire(file=dataf,exist=rythr)
      if(.not. rythr) then
        print *,'ERROR: cannot find: ',dataf
        stop
      endif  
      open(3,file=dataf)

      sdepold=-9999.
      ndata=0
      itm=time()        ! g77 - may need change on other compilers
      write(6,*) 'Now processing data:'
      inline=1

10    read(3,*,end=900,iostat=ios) idate,slat,slon,sdep,rlat,rlon,relev,
     &     stationcode,kpole
      if(ios.ne.0) then
        write(6,*) 'Error reading input file line',inline
        stop 'input file error in header line'
      endif  
      inline=inline+1

      if(jdebug.gt.0) write(13,*) 'datum: ',idate,' ',stationcode,' ',
     &slat,slon,sdep,rlat,rlon,relev

      ! recompute table of i vs phi (angle vs. distance) if new depth
      if(abs(sdep-sdepold).gt.0.01) then
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
        if(jdebug.gt.0) write(13,*) 'Ray table',ang1,ang2,ntbl,sdep,
     &          sdepold
        sdepold=sdep
      endif

      ! compute epicentral distance del. Note stadis needs geographical
      ! coordinates
      scolat=90.-slat
      rcolat=90.-rlat
      call stadis(scolat,slon,rcolat,rlon,del,az,t0,p0,t2,p1)
      del=360*((kpole+1)/2)+2*(0.5-mod(kpole,2))*del
      if(jdebug.gt.0) write(13,*) 'delta, src azimuth:',del,az
      rdel=del/r2d
      raz=az/r2d

      scolatr=geocen((90.-slat)/r2d)
      rcolatr=geocen((90.-rlat)/r2d)    ! needed for getelcor
      slonr=slon/r2d
      rlonr=rlon/r2d

      ! find all ray arrivals at distance del within tdifmax of first
      call findarr(del,tdifmax,ytbl,ntbl,a1,a2,d1,d2,trtime,narr)
      if(jdebug.eq.1) then
        write(13,*) 'findarr found the following brackets, del=',del
        do i=1,narr
          write(13,fmt='(i2,4f8.2)') i,a1(i),a2(i),d1(i),d2(i)
        enddo
      endif  

      narr=min(maxarr,narr)     ! limit number of arrivals
      ! add del,az and narr to the data line
      ndata=ndata+1
      itm=time()                ! works with g77 but is not standard

      ystart(1)=rseg(1)         ! radius (km)
      ystart(3)=0.              ! epicentral distance (radians)
      ystart(4)=0.              ! time (sec)

      ! converge to each arrival in the list
      iconv=0                   ! counts actual arrivals
      do iar=1,narr
        iter=0
        nlegs=0                 ! some dummy values if no convergence
        trtime(iar)=0.
        tstar=0
        qray=0
        p=0
        if(jdebug.eq.1) then
          write(13,*) 'Iterating for iar=',iar
          write(13,*) 'iter        a1        a2        d1        d2'
        endif  
        if(abs(d2(iar)-d1(iar)).lt.0.01) then
          anew=0.5*(a1(iar)+a2(iar))
          goto 80
        endif  
40      anew=(del-d1(iar))*(a2(iar)-a1(iar))/(d2(iar)-d1(iar))+a1(iar)
        if(anew.ge.a2(iar).or.anew.le.a1(iar).or.mod(iter,2).eq.0)
     &       anew=0.5*(a1(iar)+a2(iar))         ! bisect
        ystart(2)=anew/r2d      ! angle i with vertical (radians,0=down)
        iter=iter+1
        call tracer(ystart,rseg,ktseg,kdwn,nseg,y,jsg,rmin,tstar,h,nray)
        if(jdebug.gt.0.and.nray.le.0) 
     &    write(13,*) 'WARNING: tracer returned nray=',nray
        if(nray.le.0) goto 100
        qray=1.0e20
        if(tstar.gt.0.) qray=y(4,nray)/tstar
        aold=anew
        delnew=y(3,nray)*r2d
        if(jdebug.gt.0)
     &    write(13,43) iter,a1(iar),a2(iar),d1(iar),d2(iar),anew,delnew
43        format(i5,6f10.3)     
        if(abs(del-delnew).lt.0.01) goto 80
        if(iter.gt.itmax) then
          nray=0
          goto 100      ! lack of convergence
        endif  

        if((delnew-del)*(d1(iar)-del).gt.0.) then
          d1(iar)=delnew
          a1(iar)=anew
        else
          d2(iar)=delnew
          a2(iar)=anew
        endif  
        if((d1(iar)-del)*(d2(iar)-del).gt.0.) then
          nray=0
          goto 100   ! no convergence
        endif
        if(abs(d2(iar)-d1(iar)).lt.0.01) goto 80
        goto 40

        ! converged. Now do dynamic ray tracing
        ! anew is re-computed (but not used...)
80      if(abs(d1(iar)-d2(iar)).gt.0.)
     &   anew=(del-d1(iar))*(a2(iar)-a1(iar))/(d2(iar)-d1(iar))+a1(iar)

        
        ! compute ellipticity correction 
        ! (ADD ecorr to spherical value to get elliptical Earth time)
        !if(phase/="Pdiff" .and. del.le.95.0) then
        call getelcor(scolatr,slonr,rcolatr,rlonr,y,nray,jsg,rseg,ktseg,
     &      kdwn,p,rayvel,rayq,ecorr)
        WRITE(*,*) "Ellipticity Correction [getelcor]: ", ecorr
        
        !else
        ! compute ellipticity correction using BLN Kennet method
        call ellip_blnk(phase, del, sdep, scolat, az, tcor)
        WRITE(*,*) "Ellipticity Correction [BLNK-met]: ", tcor
        !ecorr = tcor
        !end if
        WRITE(*,*) "Difference:", abs(tcor-ecorr)/tcor*100.0, "%"
        WRITE(*,*) '==================='

        call hessian(y,ypq,nray,jsg,rseg,ktseg,kdwn,hmf,hmb,p,rayvel,
     &      rayq,Rxs,Mxs)
        ! correct travel time for dDelta
        trtime(iar)=y(4,nray)+p*(rdel-y(3,nray)) 
        if(jdebug.gt.0) write(13,*) 'Converged a,t=',anew,trtime(iar)

        iconv=iconv+1

        if(iconv.gt.1) goto 100

        ! compute the crustal corrections for each segment

        ! Changes made in Jan 2013 w.r.t. original code:
        ! terrain elevation correction at surface reflection is now
        !   incorporated into the CRUST2.0 crustal correction. Only
        !   at begin or endpoint may telev be not zero, and reflects the
        !   difference between CRUST2.0 topography and source or station 
        !   elevation (source elevation handles case that source is 
        !   *above* CRUST2.0 topo!)
        ! The delay for each segment is normally printed with the segment
        ! end point, except for a ray segment downward from the source.
        ! If the source is above the topography f CRUST2.0, the ray
        ! is computed from CRUST2.0 topo but telev>0.

        icor=0
        iray=1
        ptlat=slat
        ptlon=slon          ! geographical location of the point
        rtarget=rsrc
        kseg=0
        do while (kseg < nseg)

          ptlat0=ptlat          ! coordinates of start point
          ptlon0=ptlon
          tau=0.
          telev=0.


          kseg=kseg+1
          kd=kdwn(kseg)
          kdn=kd
          if(kseg<nseg) kdn=kdwn(kseg+1)
          ktype=ktseg(kseg)     ! is segment P or S?
          rstart=rseg(kseg)
          rtarget=rseg(kseg+1)  ! r at end of segment

          if(jdebug>0) write(13,*) 'segment ',kseg,kd,rstart,rtarget

          ! find ray index at end of segment
          if(kdn.eq.2) then     ! if heading towards turning point
            do while(y(1,iray)>y(1,iray+1)-0.01.and.iray.lt.nray-1)
              iray=iray+1
            enddo  
            rturn=y(1,iray)
          else  
            do while(abs(y(1,iray)-rtarget).gt.0.01.and.iray.lt.nray)
              iray=iray+1
            enddo  
          endif

          delx=y(3,iray)*r2d    ! distance to epicentre for segment end

          ! adjust azimuth for end point if needed
          if((mod(del,360.)-180.)*(mod(delx,360.)-180.).ge.0.) then
            azx=az
          else
            azx=mod(az+180.,360.)
          endif

          ! calculate coordinates of endpoint:
          if(iray.eq.nray) then
            ptlat=rlat
            ptlon=rlon
          else
            call givloc(scolat,slon,delx,azx,ptlat,ptlon)
          endif  
          if(jdebug.eq.1) write(13,*) 'iray,delx,ptlat,ptlon=',
     &       iray,delx,ptlat,ptlon

          ! ignore cases without crustal corrections; 

          if(kd.eq.0) then
            ! get CRUST2.0 at source
            call getcr2(ptlat0,ptlon0,ktype,vcr2,rcr2,ncr2,ecr2,rmoho,
     &          vsed)
            if(rstart.gt.rcr2(1)) then          ! if source is in the air/water
              sini=p*vsed/rstart
              cosi=sqrt(max(0.,1.0-sini*sini))
              telev=(rcr2(1)-rstart)*cosi/vsed
              rstart=rcr2(1)
            endif
            goto 90           ! for p or s wait until reflecting
          endif
          if(kd.eq.2) rstart=rturn      ! show turning r in output
          if(kd.eq.2.or.kd.eq.4) goto 90

          ! assume reflections below 6360 to be
          ! below the crust and not in need of corrections 
          if(kd.eq.3.and.rstart.lt.6360.0) goto 90

          if(kd.eq.1) then         ! case downgoing ray from source            
            ! get CRUST2.0 at source
            call getcr2(ptlat0,ptlon0,ktype,vcr2,rcr2,ncr2,ecr2,rmoho,
     &          vsed)
            if(rstart.gt.rcr2(1)) then          ! if source is in the air/water
              sini=p*vsed/rstart
              cosi=sqrt(max(0.,1.0-sini*sini))
              telev=(rstart-rcr2(1))*cosi/vsed
              rstart=rcr2(1)
            endif
            call inttau(p,rcr2,rstart,vcr2,ncr2,taucr2,tausrc2)
            tau=taucr2-tausrc2
            ! subtract reference
            call getref(ktype,vref,rref,nref,rmoho)
            call inttau(p,rref,rsrc,vref,nref,tauref,tausrcr)
            tau=tau-(tauref-tausrcr)
            if(jdebug>0) write(13,*) 'case kd=1 ',taucr2,tauref,tausrc2,
     &        ptlat0,ptlon0,ktype
            goto 90

          else if(kd.eq.5) then         ! case end of ray  

            ktype=ktseg(kseg-1)
            ! get CRUST2.0 at end of segment 
            call getcr2(ptlat,ptlon,ktype,vcr2,rcr2,ncr2,ecr2,rmoho,
     &          vsed)
            ! for station elevation correction assume basalt velocity
            ! in oceanic regions (assuming volcanic island)
            if(rcr2(1).le.6370.) vsed=9.3-2.8*min(2,ktype)
            telev=cosi*(relev-rcr2(1)+6371.0)/vsed
            call inttau(p,rcr2,rsrc,vcr2,ncr2,tau,tausrc2)
            ! subtract reference
            call getref(ktype,vref,rref,nref,rmoho)
            call inttau(p,rref,rsrc,vref,nref,tauref,tausrcr)
            tau=tau-tauref
            if(kseg.eq.2) tau=tausrc2-tausrcr
            if(jdebug>0) write(13,*) 'case kd=5 ',tau,tauref,tausrc2,
     &        ptlat,ptlon,ktype
            goto 90

          else if(kd.eq.3) then

            ! get CRUST2.0 at reflection point and do incoming segment
            ktype=ktseg(kseg-1)         ! incoming type P or S
            call getcr2(ptlat0,ptlon0,ktype,vcr2,rcr2,ncr2,ecr2,rmoho,
     &          vsed)
            call inttau(p,rcr2,rsrc,vcr2,ncr2,tau,tausrc2)
            ! subtract reference
            call getref(ktype,vref,rref,nref,rmoho)
            call inttau(p,rref,rsrc,vref,nref,tauref,tausrcr)
            tau=tau-tauref
            if(kseg.eq.2) tau=tausrc2-tausrcr
            if(jdebug>0) write(13,*) 'case kd=3a',tau,tauref,tausrc2,
     &        ptlat0,ptlon0,ktype

            ! now do downgoing leg
            ktype1=ktype
            ktype=ktseg(kseg)            ! is segment P or S?
            if(ktype.ne.ktype1) then
              call getcr2(ptlat0,ptlon0,ktype,vcr2,rcr2,ncr2,ecr2,rmoho,
     &          vsed)
              call getref(ktype,vref,rref,nref,rmoho)
              call inttau(p,rref,rsrc,vref,nref,tauref,tausrcr)
            endif
            call inttau(p,rcr2,rsrc,vcr2,ncr2,tau1,tausrc2)
            ! subtract reference
            tau=tau+tau1-tauref         ! add to upcoming leg tau
            if(jdebug>0) write(13,*) 'case kd=3b',tau1,tauref,tausrc2,
     &        ptlat0,ptlon0,ktype

            goto 90

          else
            print *,'ERROR in crustal correction'
            stop 'fatal bug'
          endif  
            

90        write(2,95) kseg,kd,ptlat0,ptlon0,rstart,tau,telev
95        format(2i3,2f8.2,f8.1,2f8.3)
          if(jdebug.eq.1)
     &      write(13,*) 'Final crustal correction:',tau,' Elev:',telev

        enddo         ! end of loop over segments for first arriving ray

100     write(2,102) idate,slat,slon,sdep,rlat,rlon,relev,stationcode
102     format(i8,2f9.3,f7.1,2f9.3,f7.3,1x,a)
        write(2,105) iar,nray,nlegs,trtime(iar),ecorr,tstar,qray,p
105     format(i2,2i5,3f10.2,2e12.3)
        if(jdebug>0) write(2,106) 'r','i','phi','t','c','1/Qs','h11',
     &      'h22','Rxs','Rxr','L','ReZ','ImZ','uzph'
106     format(6a10,4a14,a10,2a9,a12)
        write(8,107) idate,iar,stationcode,rlat,rlon,slat,slon,sdep,
     &      trtime(iar),ecorr,tstar,p
107     format(i8,i2,1x,a8,4f8.2,f6.1,f10.2,f9.2,f9.3,f8.1)
        iq=n
        y(1,nray+1)=0.          ! avoids error in if statement at the end
        y(3,nray+1)=0.          ! avoids error in if statement at the end
        jsg(nray+1)=jsg(nray)
        uzph=0.
        do i=1,nray
          zz=(1.,0.)            ! refl/transm coefficient
          pkm=p/y(1,i)          ! slowness in s/km
          k=jsg(i)              ! k<0 if ray travels upward
          k1=abs(jsg(i+1))
          ktrans=1
          if(sign(1,jsg(i)).ne.sign(1,jsg(i+1))) ktrans=0  ! reflection
          kk=abs(k)
          ktype=ktseg(kk)
          ktype1=ktseg(k1)
          if(y(1,i).eq.y(1,i+1).and.y(3,i).eq.y(3,i+1)) then
            if(ktrans.eq.0) then   ! reflection
              rr=y(1,i)  
              if(k.gt.0) then              ! down->up
                call velo2(rr,iq,1,vp1,vs1,rho1)
                call velo2(rr,iq,0,vp2,vs2,rho2)
                call refl(rps,vp1,vp2,vs1,vs2,rho1,rho2,pkm)
                zz=rps(ktype,ktype1)
              else                         ! up > down  
                call velo2(rr,iq,0,vp1,vs1,rho1)
                call velo2(rr,iq,1,vp2,vs2,rho2)
                call refl(rps,vp1,vp2,vs1,vs2,rho1,rho2,pkm)
                zz=rps(ktype,ktype1)
              endif
            endif
            if(y(1,i).eq.y(1,i+1).and.ktrans.eq.1) then   ! transmission
              rr=y(1,i)  
              if(k.gt.0) then              ! down->down
                call velo2(rr,iq,1,vp1,vs1,rho1)
                call velo2(rr,iq,0,vp2,vs2,rho2)
                call trans(tps,vp1,vp2,vs1,vs2,rho1,rho2,pkm)
                zz=tps(ktype,ktype1)
              else                         ! up -> up
                call velo2(rr,iq,0,vp1,vs1,rho1)
                call velo2(rr,iq,1,vp2,vs2,rho2)
                call trans(tps,vp1,vp2,vs1,vs2,rho1,rho2,pkm)
                zz=tps(ktype,ktype1)
              endif
            endif  
          endif  
          if(i.eq.nray) then       ! at surface zz is surface Ux,Uz factor
            call surface(pkm,ktype,ux,uz,uzph)
            zz=cmplx(ux,uz)        ! only a trick to get it printed
          endif

          c=abs(rayvel(i))
          q0=rayq(i)
          h11=hmf(1,i)+hmb(1,i)
          h22=hmf(2,i)+hmb(2,i)
          if(i.gt.1) then
            dr=y(1,i)-y(1,i-1)
            rddel=y(1,i)*(y(3,i)-y(3,i-1))
            ds=sqrt(dr*dr+rddel*rddel)    ! ray step length
            raylen=raylen+ds
          else
            raylen=0.
          endif  
          Rxs=sqrt(abs(ypq(3,i)*ypq(4,i)))/abs(rayvel(1))      ! (A56)
          Rxr=sqrt(abs(ypq(7,i)*ypq(8,i)))/abs(rayvel(nray))      ! (A57)
          ! write 
          ! y(1-4,i), i=1,...,nray: r,angle,delta,time of ray
          if(tableflag.ge.0)
     &      write(2,120) (y(j,i),j=1,4),c,q0,h11,h22,Rxs,Rxr,raylen,
     &      zz,uzph
120       format(f10.3,2f10.5,f10.3,f10.5,f10.6,4e14.6,f10.1,2f9.3,
     &            g12.3)
        enddo

      enddo     ! end of do loop over iar

      goto 10

900   print *, 'End of data file reached'
      print *, 'Total number of raypaths:',ndata
      print *
      print *,'***************************************************'
      print *,'*Warning: output for option tableflag=0 now also  *'
      print *,'*includes partial travel times; format has changed*'
      print *,'***************************************************'
      print *
      close(7)

      end
      
      subroutine findarr(del,tdifmax,ytbl,ntbl,a1,a2,d1,d2,trtime,narr)

      ! find all arrivals at distance del and sort them, first
      ! arrival first

      ! input:
      ! del=epicentral distance in degrees
      ! tdifmax=max time differential for later arrivals
      ! ytbl(1,i)=ray angles, (2,i)=distances (deg), (3,i)=times (s)
      ! ntbl=number of rays in table ytbl
      ! ray angles can be ascending or descending but must be monotonic

      ! output
      ! narr=number of arrivals, then for each one of these in arrays:
      ! a1,a2=bracket for ray angle (0<a1<a2<90 or 90<a1<a2<180) 
      ! d1,d2=deltas belong to a1,a2 (degr)
      
      dimension a1(10),a2(10),d1(10),d2(10),trtime(10)
      dimension ytbl(4,4000),indx(10)

      jdebug=0

      if(jdebug.gt.0) write(13,*) 'subroutine findarr,del=',del
      if(jdebug.gt.0) write(13,*) '#      i        Delta(deg)   T(s)'

      narr=0
      if(jdebug.gt.0)
     &  write(13,*) 1,ytbl(1,1),ytbl(2,1),ytbl(3,1)
      do i=2,ntbl
        if(jdebug.gt.0)
     &    write(13,*) i,ytbl(1,i),ytbl(2,i),ytbl(3,i)
        if((del-ytbl(2,i-1))*(del-ytbl(2,i)).lt.0.0
     &        .or.del.eq.ytbl(2,i-1)) then
          narr=min(10,narr+1)
          a1(narr)=ytbl(1,i-1)
          a2(narr)=ytbl(1,i)
          d1(narr)=ytbl(2,i-1)
          d2(narr)=ytbl(2,i)
          dd=d2(narr)-d1(narr)
          trtime(narr)=ytbl(3,i-1)
          if(dd.ne.0.) trtime(narr)=ytbl(3,i-1)+(del-d1(narr))*
     &          (ytbl(3,i)-ytbl(3,i-1))/dd
          if(jdebug.gt.0)
     &      write(13,*) 'bracket:',narr,dd,trtime(narr)
          if(narr.ge.10) goto 20
        endif
      enddo

      ! swap if a1>a2
      do i=1,narr
        if(a1(i).gt.a2(i)) then
          aa=a2(i)
          a2(i)=a1(i)
          a1(i)=aa
          aa=d2(i)
          d2(i)=d1(i)
          d1(i)=aa
        endif
      enddo  

20    if(narr.le.1) return

c     write(13,*) 'sort:'
c     do i=1,narr
c       write(13,9) i,trtime(i),a1(i),a2(i),d1(i),d2(i)
c     enddo
      call indexx(narr,trtime,indx)
      if(jdebug.eq.1) write(13,*) 'indx=',(indx(i),i=1,narr)
      call sort(narr,a1,indx)
      call sort(narr,a2,indx)
      call sort(narr,d1,indx)
      call sort(narr,d2,indx)
      call sort(narr,trtime,indx)

c     write(13,*) 'gives:'
c     do i=1,narr
c       write(13,9) i,trtime(i),a1(i),a2(i),d1(i),d2(i)
c     enddo

      nn=narr
      do i=1,nn
        if(trtime(i)-trtime(1).le.tdifmax) narr=i
      enddo  

c     write(13,*) 'tdifmax',tdifmax,' gives narr=',narr
9     format(i3,5f10.3)

      return
      end

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END

      SUBROUTINE SORT(N,RA,INDX)
      DIMENSION RA(10),WKSP(10),INDX(10)
      DO 11 J=1,N
        WKSP(J)=RA(J)
11    CONTINUE
      DO 12 J=1,N
        RA(J)=WKSP(INDX(J))
12    CONTINUE
      RETURN
      END

      subroutine surface(p,ktype,ux,uz,uzph)

      ! computes surface displacement factors for incoming P
      ! (ktype=1), SV (=2) or SH (=3) waves
      ! See Aki&Richards Problem 5.6
      ! Temporary fix: if oceanic, seabottom is assumed to be at 
      ! surface (NEEDS TO BE IMPROVED)

      ! input: slowness p in s/km, ktype
      ! output amplitudes ux,uz, uzph is phase of uz in case of evanescent
      !   P reflection from SV
      ! Note uzph is the phase for Uz. The phase for Ux is (uzph-pi/2), if
      ! P is evanescent (i.e. if uzph /= 0), else 0.
      ! sign convention for phase is as in Aki and Richards


      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      jdebug=0

      if(ktype.eq.3) then               ! SH
        ux=2.0
        uz=0.
        return
      endif  

      a=vp(nsl)
      b=vs(nsl)
      a2=a*a
      b2=b*b
      p2=p*p
      sini=p*a
      cosi=sqrt(abs(1.0-sini*sini))
      uzph=0.
      sinj=p*b
      cosj=sqrt(1.0-sinj*sinj)
      b2p=(1./b2-2*p2)
      if(sini.le.1.0) then
        Dre=b2p*b2p+4*p2*cosi*cosj/(a*b)
        Dim=0.
      else
        Dre=b2p*b2p
        Dim=4*p2*cosi*cosj/(a*b)
        uzph=atan2(Dre,Dim)            ! phase for Uz
      endif  

      if(jdebug.gt.0) write(13,*) 'surf:',p,a,b,sini,cosi,Dre,Dim,uzph

      if(ktype.eq.1) then               ! P
        ux=4*p*cosi*cosj/(Dre*b**3)
        uz=-2*cosi*b2p/(Dre*b2)
        return
      else                              ! SV
        if(sini.le.1.0) then
          ux=2*cosj*b2p/(Dre*b2)
          uz=4*p*cosi*cosj/(a*b2*Dre)
        else
          D=sqrt(Dre*Dre+Dim*Dim)
          ux=2*cosj*b2p/(D*b2)
          uz=4*p*cosi*cosj/(a*b2*D)
        endif
        return
      endif

      end
