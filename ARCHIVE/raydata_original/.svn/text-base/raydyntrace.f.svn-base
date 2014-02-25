      ! compile: g77 -o raydyntrace raydyntrace.f raytracesubs.f azdel.f 

      program raydyntrace

      ! Testing and documentation:
      ! Tian, Y., et al., Dynamic ray tracing and travel time corrections
      ! for global seismic tomography, subm. to XX, 2007.

      ! Background reading:
      ! Nolet, G., A Breviary for Seismic Tomography, Cambridge Univ. Press, accepted
      ! for publication, 2007.

      !
      ! August 2006: removed relev from elevation corrections when kd<5
      !              added stop for sources above local surface
      !              corrected ptlat/lon for surface source in PPP etc.
      !              changed azimuth in givloc call for delta's > 180
      ! Sept 2006: Depth of source shallower than CRUST2.0 surface is replaced 
      !            with CRUST2.0 elev; a list of shallow sources is in src.shallow
      !            Output warning message of imprecise hessian against analytical expression
      ! Dec 2007: Replaced Kennett's ellipticity routine with getelcor. 
      

      ! input:
      ! [1] a ray definition file (see rdray)
      ! [2] a file with source and receiver coordinates, each line has:
      ! idate,slat,slon,sdepth,rlat,rlon,relev,stationcode,kpole
      ! (lat, lon in degrees, depth and elevation in km)
      ! [3] a velocity model for the Earth (see subroutine model)
      ! [4] model CRUST2.0 (see subroutine rcrust)

      ! output:
      ! raydyntrace.xxxx

      ! the program is most efficient if the data are
      ! sorted in segments of lines with equal source depths

      parameter (ityp=360,nla=90,nlo=180)                        ! for rcrust2 arrays
      parameter (NFR=10)                                        ! nr of freq bands

      dimension rseg(20),ktseg(20),kdwn(20)                     ! ray segments
      dimension ystart(4),y(4,9000),ytbl(4,4000),jsg(9000)  ! ray vector
      dimension a1(10),a2(10),d1(10),d2(10),trtime(10),indx(10)   ! ray times/dist
      dimension hmf(2,9000),hmb(2,9000),rayvel(9000),rayq(9000) ! Hessian
      dimension ypq(10,9000)                                    ! idem
      dimension rcr2(16),vcr2(16)                               ! CRUST2.0
      dimension rref(50),vref(50)                               ! Reference crust

      logical rythr
      character*72 fname,tablef,chray
      character*30 dataf
      character*16 stationcode
      character*8 phase
      character*2 ctype(ityp),types,atype(nlo)
      character*1 yn,star
      character*3 comp

      ! background model
      common /mod/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &     Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      ! the following are filled by subroutine rcrust2:
      common/crust2/amapvp(8,nlo,nla),amaprho(8,nlo,nla),
     &              amapvs(8,nlo,nla),amapthi(7,nlo,nla),
     &              amapele(nlo,nla)                    ! model CRUST2.0
      common/glbgrd/types(nlo,nla)                      ! model CRUST2.0

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/

      jdebug=0

      ! set some constants
      h=15.     ! integration step size (must be equal that in mkrtable)
      itmax=20  ! max nr of iterations for ray "shooting"

      ! read CRUST2.0 data
      call rcrust2

      print *,'Give model file:'
      read(5,fmt='(a)') fname
      inquire(file=fname,exist=rythr)
      if(.not. rythr) then
        print *,'ERROR: cannot find: ',fname
        stop
      endif  
      call model(fname)

      write(6,*) 'Radius of Earth model:',r(n)
      if(abs(r(n)-6371.).gt.0.1) write(6,*) 'WARNING! Model radius',
     & ' deviates from 6371 - crustal corrections are not OK'
      write(6,*) 'Radius of CMB:',r(noc)
      write(6,*) 'Radius of ICB:',r(nic)
      write(6,*) 'Radius of Moho:',r(nmoh)

      print *,'Give max number of arrivals to include in computations'
      print *,'and their max difference (in sec) w.r.t. first arrival:'
      read *,maxarr,tdifmax

      ! open the input file
      print *,'Give input file name IN THIS DIRECTORY (eg Pshallow):'
      read(5,fmt='(a)') dataf
      inquire(file=dataf,exist=rythr)
      if(.not. rythr) then
        print *,'ERROR: cannot find: ',dataf
        stop
      endif  
      open(3,file=dataf)

      ! open the big output file
      tablef='raydyntrace.'//dataf
      print *,'Output is to file ',tablef
      open(2,file=tablef)

      ! open the diagnostic output file
      tablef='out.raydyntrace.'//dataf
      print *,'Diagnostic output is to file ',tablef
      open(8,file=tablef)
      write(8,5)
5     format('    Date N Receiver    Rlat    Rlon    Slat    Slon',
     &    '  Sdep    Time Elcor    t*       p')


      ! read ray info from the raw data file, write it to the processed file
      print *,'Give file with ray definition:'
      read(5,'(a)') fname
      open(4,file=fname)
      call rdray(4,2,chray,phase,rseg,ktseg,kdwn,nseg)

      sdepold=-9999.
      ndata=0
      itm=time()        ! g77 - may need change on other compilers
      write(6,*) 'Now processing data:'
      open(7,file='src.shallow.'//dataf) ! sources shallower than CRUST2.0 surface
      write(7,*) 'Sources shallower than CRUST2.0 surface:'
      write(7,*) '  idate   srclat   srclon    srcR CRUST2.0R    dt0',
     & ' kd'
      jshallow=0 ! jshallow counts sources above CRUST2.0 surface
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
        call getelcor(scolatr,slonr,rcolatr,rlonr,y,nray,jsg,rseg,ktseg,
     &      kdwn,p,rayvel,rayq,ecorr)
        call hessian(y,ypq,nray,jsg,rseg,ktseg,kdwn,hmf,hmb,p,rayvel,
     &      rayq)
        ! correct travel time for dDelta
        trtime(iar)=y(4,nray)+p*(rdel-y(3,nray)) 
        if(jdebug.gt.0) write(13,*) 'Converged a,t=',anew,trtime(iar)

        iconv=iconv+1

        if(iconv.gt.1) goto 100

        ! compute the crustal corrections at each reflection point
        icor=0
        iray=1
        ptlat=slat
        ptlon=slon          ! geographical location of the point
        rtarget=rsrc
        jup=0               ! sense switch for surface phases like pP
        if(kdwn(1).eq.0) jup=1
        do kseg=1,nseg
          telev=0.               ! station elevation correction
          kd=kdwn(kseg)
          if(kseg.gt.2) jup=0
          rtarget=rseg(kseg)
          if(kd.eq.5) iray=nray
          if(jdebug.eq.1) write(13,*) 'kseg,kd,jup,rtarget:',kseg,
     &               kd,jup,rtarget
          tau=0.
          if(kd.eq.2.or.kd.eq.4) goto 90         
          do while(abs(y(1,iray)-rtarget).gt.0.01.and.iray.lt.nray)
            iray=iray+1
          enddo  
          delx=y(3,iray)*r2d  ! distance to epicentre for segment end
          ! calculate azimuth of reflection or end point
          if((mod(del,360.)-180.)*(mod(delx,360.)-180.).ge.0.) then
            azx=az
          else
            azx=mod(az+180.,360.)
          endif
          ! calculate coordinates of reflecion or endpoint:
          if(iray.eq.nray) then
            ptlat=rlat
            ptlon=rlon
          else
            call givloc(scolat,slon,delx,azx,ptlat,ptlon)
          endif  
          if(jdebug.eq.1) write(13,*) 'iray,delx,ptlat,ptlon=',
     &       iray,delx,ptlat,ptlon
          iray=min(nray,iray+2)  ! avoid problem in next section for PPP
          ! get CRUST2.0 at start of segment kseg
          ktype=ktseg(kseg)            ! is segment P or S?
          ktype1=ktype
          if(kd.eq.5) ktype=ktseg(kseg-1)
          if(kseg.gt.1) ktype1=ktseg(kseg-1)  ! and incoming?
          ! do departing ray at start of segment
          call getcr2(ptlat,ptlon,ktype,vcr2,rcr2,ncr2,ecr2,rmoho,
     &          vsed)
          vsurf=vsed
          ! for islands take Vp=6.5, Vs=3.7 in elevation correction:
          if(rcr2(1).le.6370.) vsurf=9.3-2.8*ktype  ! basalt velocity
          ! if source is above CRUST2.0 surface
          ! replace source radius with CRUST2.0 surface radius
          ! this assumes r(n)=6371, otherwise the crustal correction
          ! makes no sense.
          if(kseg.eq.1.and.rsrc.gt.rcr2(1)) then
            vsurfref=vp(n)
            if(ktype.eq.2) vsurfref=vs(n)
            cosi=sqrt(1.0-(vsurfref*p/6371.0)**2)
            dt0=(rsrc-rcr2(1))*cosi/vsurfref ! origin time correction
            if(kd.eq.0) dt0=-dt0        ! ray is up from source
            write(7,50) idate,slat,slon,rsrc,rcr2(1),dt0,kd
            rsrc=rcr2(1) 
            jshallow=jshallow+1
50          format(i8,2f9.2,f8.1,f10.1,f7.3,i3)
          endif
          call inttau(p,rcr2,rsrc,vcr2,ncr2,taucr2,tausrc2)
          tau=0                       ! intercept time
          ! Note: station elevation relev is wrt sealevel (r=6371), 
          ! not surface radius rcr2(1) of CRUST2.0

          ! find CR2 tau at segment start for kd=0,1 or 3 (i.e. if
          ! it starts from source or reflects at surface)
          if(kd.eq.0) then            ! if ray is up from source
            tau=tausrc2
            cosi=sqrt(1.0-(vsurf*p/6371.0)**2)
            telev=cosi*(-rcr2(1)+6371.0)/vsurf  
            if(jdebug.gt.0) write(13,*) 'elev,rcr2,telev=',
     &         -rcr2(1)+6371.,rcr2(1),telev
          else if(kd.eq.1) then       ! down from source  
            tau=taucr2-tausrc2
          ! the following assumes surface reflection if r>6360
          else if(kd.eq.3.and.rtarget.gt.6360.) then 
            tau=taucr2
            cosi=sqrt(1.0-(vsurf*p/6371.0)**2)
            telev=cosi*(-rcr2(1)+6371.0)/vsurf
            if(jdebug.gt.0) write(13,*) 'elev,rcr2,telev=',
     &         -rcr2(1)+6371.,rcr2(1),telev
          endif  
          if(jdebug.eq.1.and.kd.lt.5) 
     &          write(13,*) 'Outgoing tau,tel=',tau,telev

          ! at this point tau is the tau travelled in CR2 from start
          ! of the segment
          ! now add contribution from ray *to* start of segment in
          ! case kd=3 or 5 (i.e. surface reflection or end point)
          if(ktype.ne.ktype1) then      ! if wave type differs
            call getcr2(ptlat,ptlon,ktype1,vcr2,rcr2,ncr2,ecr2,rmoho,
     &          vsed)
            vsurf=vsed
            if(rcr2(1).le.6370.) vsurf=9.3-2.8*ktype1  ! basalt velocity 
            call inttau(p,rcr2,rsrc,vcr2,ncr2,taucr2,tausrc2)
          endif
          ! avoid double counting of pP if jup=1
          if(kd.eq.3.and.rtarget.gt.6360.and.jup.eq.0) then
            tau=tau+taucr2
            cosi=sqrt(1.0-(vsurf*p/6371.0)**2)
            telev=telev+cosi*(-rcr2(1)+6371.0)/vsurf
            if(jdebug.gt.0) write(13,*) 'elev,rcr2,telev=',
     &         -rcr2(1)+6371.,rcr2(1),telev
          else if(kd.eq.5.and.jup.eq.0) then
            tau=taucr2
            cosi=sqrt(1.0-(vsurf*p/6371.0)**2)
            telev=cosi*(relev-rcr2(1)+6371.0)/vsurf
            if(jdebug.gt.0) write(13,*) 'elev,rcr2,telev=',
     &         relev-rcr2(1)+6371.,rcr2(1),telev
          endif  
          if(jdebug.eq.1.and.kseg.gt.1) 
     &          write(13,*) 'Incoming tau,telev=',taucr2,telev

          ! subtract the reference model tau
          ktype=ktseg(kseg)          ! is departing ray P or S?
          if(kd.eq.5) ktype=ktseg(kseg-1)
          if(kseg.gt.1) ktype1=ktseg(kseg-1)  ! and incoming?
          ! do departing first
          call getref(ktype,vref,rref,nref,rmoho)
          call inttau(p,rref,rsrc,vref,nref,tauref,tausrcr)
          if(kd.eq.0) then            ! if ray is up from source
            tau=tau-tausrcr
          else if(kd.eq.1) then       ! down from source  
            tau=tau-(tauref-tausrcr)
          ! the following assumes surface reflection if r>6360
          else if(kd.eq.3.and.rtarget.gt.6360.) then 
            tau=tau-tauref
          endif  
          if(jdebug.eq.1.and.kd.lt.5) 
     &            write(13,*) 'Outgoing ref tau=',tauref,tausrcr
          ! now add incoming ray
          if(ktype.ne.ktype1) then    ! if wavetype differs
            call getref(ktype1,vref,rref,nref,rmoho)
            call inttau(p,rref,rsrc,vref,nref,tauref,tausrcr)
          endif  
          if(kd.eq.3.and.rtarget.gt.6360.and.jup.eq.0) then
            tau=tau-tauref
          else if(kd.eq.5.and.jup.eq.0) then
            tau=tau-tauref
          endif  
          if(jdebug.eq.1.and.kseg.gt.1) 
     &          write(13,*) 'Incoming ref tau=',tauref

90        if(kd.ne.2.and.kd.ne.4) then
            write(2,95) kseg,kd,ptlat,ptlon,rtarget,tau,telev
95          format(2i3,2f8.2,f8.1,2f8.3)
          else
            write(2,96) kseg,kd,0,0,rtarget,0,0
96          format(2i3,2i8,f8.1,2i8)
          endif  
          if(jdebug.eq.1)
     &      write(13,*) 'Final crustal correction:',tau,' Elev:',telev

        enddo         ! end of loop over segments for first arriving ray

100     write(2,102) idate,slat,slon,sdep,rlat,rlon,relev,stationcode
102     format(i8,2f9.3,f7.1,2f9.3,f7.3,1x,a)
        write(2,105) iar,nray,nlegs,trtime(iar),ecorr,tstar,qray,p
105     format(i2,2i5,3f10.2,2e12.3)
        write(8,107) idate,iar,stationcode,rlat,rlon,slat,slon,sdep,
     &      trtime(iar),ecorr,tstar,p
107     format(i8,i2,1x,a8,4f8.2,f6.1,f8.1,2f6.1,f8.1)
        iq=n
        do i=1,nray
          c=rayvel(i)
          q0=rayq(i)
          h11=hmf(1,i)+hmb(1,i)
          h22=hmf(2,i)+hmb(2,i)
          write(2,120) (y(j,i),j=1,3),c,q0,h11,h22
120       format(f10.3,3f10.5,f10.6,4e14.6)
        enddo

      enddo     ! end of do loop over iar

      goto 10

900   print *, 'End of data file reached'
      print *, 'Total number of raypaths:',ndata
      close(7)
      if(jshallow.gt.0) then
        print *, 'WARNING: The source depths in ',jshallow,' paths with'
        print *, 'sources shallower than CRUST2.0 surface were replaced'
        print *, 'with source depth at CRUST2.0 surface.'
        print *, 'See file src.shallow.'//dataf
      endif

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
      ! narr=number of arrivals, then for each one:
      ! a1,a2=bracket for ray angle (a1<a2) in degrees
      ! d1,d2=delta's belong to a1,a2 (degr)
      
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
