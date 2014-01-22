      ! Version 1.2.1 July 2010

      ! testing done on a limited set of phases

      ! Known issues:
      ! Bug: crustal correction not correct for pS and similar phases if
      ! p reflects at an oceanic surface.
      ! Bug: crustal corrections not correct for OBS on ocean bottom
      ! Todo: crustal correction not correct for stations situated on
      ! an ice layer
      ! Crustal corrections may become degenerate if the ray does not
      ! turn below the deepest Moho in Crust2.0. This is not a bug but
      ! reflects the limitations of linearized crustal corrections.

      ! List of changes:
      ! June 06 - warnings for ray length changed to 8998
      ! July 06 - changed mkraytable to avoid rays trapped in LVZ
      ! Aug 06 - fixed bug in trapezoidal rule in inttau for sources
      !          not on a model interface
      !          fixed nleg bug in getlegs for reflected phases
      !          fixed y() index bug in getlegs for reflected phases
      !          adjusted criterion to recognize reflections in getlegs
      !          fixed bug in continuity of reflected phases at CMB
      ! Sep 06 - in hessian:
      !          Adjusted criterion to avoid taking a turning point as reflection point
      !          Treat the region near source/receiver as a homogeneous
      !                sphere to avoid oscillation of hessian elements.
      ! Nov 06 - Output warning message of imprecise hessian against analytical expression
      !          Add precision test of A53 on backward H11
      ! April 08 - corrected getelcor for ghost phases like pP. Version 1.0 would
      !          still go wrong for layers thinner than the change in r for an
      !          elliptical Earth. Improved spacing of rays in the table if tableflag=1
      !          and fixed problem with rays near 180 degrees.
      ! Dec 08 - tables can now be generated for *any* source depth. Output
      !          for tableflag>0 option expanded to include Rxr as well as
      !          Rxs, and useful information repeated in header
      ! March 09 - negative tableflag option added for computation of
      !          travel time curves only with limited output
      ! Nov 09   changed unit numbers for reading crust2.0 files
      ! June 11 - changed ray length to LRAY and model size to 999
      ! Jan 12 - added computation of reflection/transmission coefficients
      ! Feb 13 - corrected maslov Mxs for forward ray


      subroutine hessian(y,ypq,nray,jsg,rseg,ktseg,kdwn,hmf,hmb,p,
     &      rayvel,rayq,Rxs,Mxs)

      ! calculate forward and backward hessian matrix (hmf,hmb)
      ! equation numbers refer to Dahlen et al, GJI 141:157-174,2000

      ! input:
      ! y(i,j) - ray geometry from subroutine tracer, j=1,nray
      !        y1=r, y2=angle i, y3=distance delta, y4=time
      ! nray - number of ray points
      ! jsg(i) - segment number belonging to ray node i
      ! rseg,ktseg,kdwn - ray description (see subroutine rdray)

      ! output
      ! hmf - M11 and M22 for forward Hessian (A47)
      ! hmb - M11 and M22 for backward Hessian (A48)
      !  (e.g. hmf(2,j) is M22 forward in node j).
      ! slowness p (s/rad)
      ! rayvel - array with model velocity at ray node, <0 if Vs
      ! rayq(i) - 1/Qs at node i
      ! Rxs - geometrical spreading for forward ray
      ! Mxs - maslov index for forward ray

c     double precision ypq,ynew
      parameter (LRAY=20000)
      dimension y(4,LRAY),jsg(LRAY),rayvel(LRAY),dv(LRAY),rayq(LRAY)
      dimension hmf(2,LRAY),hmb(2,LRAY),ypq(10,LRAY)
      dimension rseg(20),ktseg(20),kdwn(20),ynew(10)
      integer nray

      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/,eps/0.0/
      ! TODO: safer to set eps=1e-4 etc, but gives problem when the turning point
      ! is very close to the 1-D model node
      data r0/20./  ! radius of homogeneous sphere, the same r0 as in raymatrix
      data err0/1.0/  ! error tolerance of hessian against analytical expression

      jdebug=0     ! 0 no debug; 1 check accuracy against analytical expression
                   ! 2 for more debug output
      mray=nray

      ! set initial conditions

      ! forward elements of the partition matrices P and Q
      ypq(1,1)=1.0   ! P1 (A35)
      ypq(2,1)=1.0   ! P2 (A35)
      ypq(3,1)=0.0   ! Q1 (A35)
      ypq(4,1)=0.0   ! Q2 (A35)

      !backward
      ypq(5,1)=0.0   ! P~1 (A37)
      ypq(6,1)=0.0   ! P~2 (A37)
      ypq(7,1)=1.0   ! Q~1 (A37)
      ypq(8,1)=1.0   ! Q~2 (A37)

      ! ray radius and angle with vertical
      ypq(9,1)=y(1,1)  ! r
      ypq(10,1)=y(2,1) ! i (0 for upward)

      iq=n
      iseg=1
      ktype=ktseg(1)

      ! find source velocity (assume source in the entering layer
      ! if at a discontinuity)
      call velo1(rseg(1),iq,1-kdwn(1),ktseg(1),csrc,dc,d2c)
      p=rseg(1)*sin(y(2,1))/csrc     ! slowness

      ! temp debug
      ! if(p>476.) jdebug=2

      p2=p*p
      rayvel(1)=csrc
      if(ktseg(1).gt.1) rayvel(1)=-rayvel(1)
      rayq(1)=99999.
      if(Qs(iq).ne.0.)rayq(1)=1.0/Qs(iq)
      dv(1)=dc
      if(jdebug.gt.1) write(13,*) 'hessian p,csrc,i0=',p,csrc,y(2,1)

      if(jdebug.gt.1) then
      write(13,98)
98    format(5x,'i',3x,'jsg',4x,'kd',5x,'phi',9x,'r',6x,'c',3x,'ypq1-8')
      write(13,99) 1,jsg(1),kdwn(iseg),y(3,1)*r2d,y(1,1),rayvel(1),
     &       (ypq(j,1),j=1,8)
      endif

      ! integrate with phi as independent variable, unless a
      ! discontinuity is crossed (when we apply boundary conditions)
      Rxs2old=1.0                       ! Spreading for forward ray
      Mxs=0                             ! Maslov for forward ray
      do i=2,mray                       ! step from i-1 to i
        jdown=1
        if(jsg(i).lt.0) jdown=0
        dphi=y(3,i)-y(3,i-1)
        rr=y(1,i-1)
        dr=y(1,i)-rr
        call findiq(iq,jdown,rr)
        drnode=abs(rr-r(iq))    ! distance from node
        
        ! check if node i-1 and i are at same radius of a discontinuity
        if(abs(dr).le.eps.and.drnode.le.eps) then    
          rr2=rr*rr
          csi=cos(y(2,i-1))
          if(jdebug>1) write(13,*)'Discont i,jsgs=',i,jsg(i-1),jsg(i)
          ! c is velocity of incoming ray, c1 of outgoing ray at
          ! this discontinuity
          kupr=jdown    ! case of transmitted ray
          if(jsg(i)*jsg(i-1).lt.0) kupr=1-jdown       ! reflection
          call velo1(rr,iq,kupr,ktype,c,dc,d2c)
          if(jdebug.gt.1)
     &    write(13,*) i-1,' Incoming rr,c=',rr,c,dc,' kupr=',kupr

          ktype1=ktseg(abs(jsg(i)))     ! type may change at point i
          if(jsg(i)*jsg(i-1).gt.0) kupr=1-kupr        ! transmission
          call velo1(rr,iq,kupr,ktype1,c1,dc1,d2c)

          rayvel(i)=c1
          if(ktype1.gt.1) rayvel(i)=-rayvel(i)
          rayq(i)=99999.
          if(Qs(iq).ne.0.) rayq(i)=1.0/Qs(iq)
          dv(i)=dc1
          sini=p*c1/rr
          if(jdebug.gt.1)
     &    write(13,*) i,' Outgoing rr,c1=',rr,c1,dc1,' kupr=',kupr

          ! The following uses (A32)-(A34). To re-derive subsitute
          ! (A34) first
          if(jdebug.gt.1)
     &      write(13,*) '(refl/transmit, kdwn,sini=',kdwn(abs(jsg(i))),
     &          sini,')'
          csi1=cos(y(2,i))
          ratcs=csi1/csi
          ypq(1,i)=ypq(1,i-1)/ratcs
     &         +(p2*dc/rr2-1./(rr*c))*ypq(3,i-1)/csi1
     &         -(p2*dc1/rr2-1./(rr*c1))*ypq(3,i-1)/csi
          ypq(2,i)=ypq(2,i-1)+ypq(4,i-1)*(csi1/c1-csi/c)/rr   ! (A33)
          ypq(3,i)=ypq(3,i-1)*ratcs                           ! (A32)
          ypq(4,i)=ypq(4,i-1)                                 ! (A33)
          ypq(5,i)=ypq(5,i-1)/ratcs
     &         +(p2*dc/rr2-1./(rr*c))*ypq(7,i-1)/csi1
     &         -(p2*dc1/rr2-1./(rr*c1))*ypq(7,i-1)/csi
          ypq(6,i)=ypq(6,i-1)+ypq(8,i-1)*(csi1/c1-csi/c)/rr   ! (A33)
          ypq(7,i)=ypq(7,i-1)*ratcs                           ! (A32)
          ypq(8,i)=ypq(8,i-1)                                 ! (A33)
          ypq(9,i)=y(1,i)
          ypq(10,i)=y(2,i)
          Rxs2=ypq(3,i)*ypq(4,i)                ! (A56)
          if(jdebug.gt.1) then
            write(13,*) 'terms:',ypq(1,i-1)/ratcs,
     &           (p2*dc/rr2-1./(rr*c))*ypq(3,i-1)/csi1,
     &           -(p2*dc1/rr2-1./(rr*c1))*ypq(3,i-1)/csi
            a32=csi*ypq(1,i-1)+(p2*dc/rr2-1./(rr*c))*ypq(3,i-1)
            a33=ypq(2,i-1)-csi*ypq(4,i-1)/(rr*c)
            a34=ypq(3,i-1)/csi
            write(13,'(i5,a,3e15.6)') i-1,' Incoming a32-4:',a32,a33,a34
            a32=csi1*ypq(1,i)+(p2*dc1/rr2-1./(rr*c1))*ypq(3,i)
            a33=ypq(2,i)-csi1*ypq(4,i)/(rr*c1)
            a34=ypq(3,i)/csi1
            write(13,'(i5,a,3e15.6)') i,' Outgoing a32-4:',a32,a33,a34
            write(13,*) 'terms:',ypq(5,i-1)/ratcs,
     &           (p2*dc/rr2-1./(rr*c))*ypq(7,i-1)/csi1,
     &           -(p2*dc1/rr2-1./(rr*c1))*ypq(7,i-1)/csi
            a32=csi*ypq(5,i-1)+(p2*dc/rr2-1./(rr*c))*ypq(7,i-1)
            a33=ypq(6,i-1)-csi*ypq(8,i-1)/(rr*c)
            a34=ypq(7,i-1)/csi
            write(13,'(i5,a,3e15.6)') i-1,' Incoming a32bw:',a32,a33,a34
            a32=csi1*ypq(5,i)+(p2*dc1/rr2-1./(rr*c1))*ypq(7,i)
            a33=ypq(6,i)-csi1*ypq(8,i)/(rr*c1)
            a34=ypq(7,i)/csi1
            write(13,'(i5,a,3e15.6)') i,' Outgoing a32bw:',a32,a33,a34
            write(13,99)i,jsg(i),kdwn(iseg),y(3,i)*r2d,y(1,i),rayvel(i),
     &         (ypq(j,i),j=1,8)
          endif
        else            ! normal propagation (no discontinuity)
          rr=y(1,i)
          ktype1=ktseg(abs(jsg(i)))     ! type may change at point i
          call velo1(rr,iq,jdown,ktype1,c1,dc1,d2c)
          rayvel(i)=c1
          if(ktype1.gt.1) rayvel(i)=-rayvel(i)
          rayq(i)=99999.
          if(Qs(iq).ne.0.) rayq(i)=1.0/Qs(iq)
          dv(i)=dc1                     ! store first derivative of c at node i
          call rk4pq(ypq(1,i-1),iq,dphi,ktype,jdown,ynew,p)
          do j=1,10
            ypq(j,i)=ynew(j)
          enddo  
          q2th=rr*rseg(1)*sin(y(3,i))/p         ! (A49) safeguard precision loss
          Rxs2=ypq(3,i)*q2th                    ! (A56) argument with sign
          if(Rxs2*Rxs2old<0.) Mxs=Mxs+1         ! ignore sign change at discon

          ! debug
          if(jdebug.gt.1)
     &    write(13,99) i,jsg(i),kdwn(iseg),y(3,i)*r2d,y(1,i),rayvel(i),
     &       (ypq(j,i),j=1,8),Rxs2,Mxs
99        format(3i6,f8.2,f10.2,f7.3,9e12.3,i3)

          ypq(9,i)=y(1,i)
          ypq(10,i)=y(2,i)

        endif

        Rxs2old=Rxs2
        iseg=abs(jsg(i))           ! segment number for next step
        ktype=ktseg(iseg)          ! wave type for next step
        rstop=rseg(iseg)           ! r at end of segment

      end do

      q1tr=ypq(7,mray)
      q2tr=ypq(8,mray)
      q1fr=ypq(3,mray)
      q2fr=ypq(4,mray)
      if(jdebug.gt.1) write(13,*) 'qifr etc.:',q1fr,q2fr,q1tr,q2tr,mray

      ! constants used in the following do loop of raynodes
      y3d=y(3,mray)
      sid=sin(y3d)
      rcosis=y(1,1)*cos(y(2,1))
      rcosisb=-y(1,mray)*cos(y(2,mray))
      xs=y(1,1)*cos(y(3,1))    ! source coordinates in equatorial plane
      ys=y(1,1)*sin(y(3,1))    ! x=r*cos phi,y=r*sin phi
      xr=y(1,mray)*cos(y(3,mray))   ! receiver coordinates in equatorial plane
      yr=y(1,mray)*sin(y(3,mray))    ! x=r*cos phi,y=r*sin phi

      if(jdebug>1) then
        write(13,'(12x,a7,4a7,3(11x,a4))') 'vel','A49P','A49Q','A50',
     &  'A51','A53f','A53b','A54'
        write(13,'(a19,4a7,3(11x,a4))') 'should be:','%','%','%',
     &  '%','cont','cont','c+/-'
      endif

      Rxs=sqrt(abs(Rxs2))/csrc                          ! (A56)

      dr=1.0
      do k=1,mray

c replace ypq(5-8,k) as backward P1,P2,Q1,Q2 see eqs (A43) to (A46)
        ypq(5,k)=ypq(1,k)*q1tr-ypq(5,k)*q1fr         ! P3
        ypq(6,k)=ypq(2,k)*q2tr-ypq(6,k)*q2fr         ! P4
        ypq(7,k)=ypq(3,k)*q1tr-ypq(7,k)*q1fr         ! -Q3
        ypq(8,k)=ypq(4,k)*q2tr-ypq(8,k)*q2fr         ! -Q4

c forward hessian M11=P1/Q1 and M22=P2/Q2  
        if (k.eq.1.or.ypq(3,k).eq.0.0) then          
          hmf(1,k)=sign(999.,ypq(1,k))
        else
          hmf(1,k)=ypq(1,k)/ypq(3,k)        ! M'1=P1/Q1
        endif

        if (k.eq.1.or.ypq(4,k).eq.0.0) then
          hmf(2,k)=sign(999.,ypq(2,k))
        else
          hmf(2,k)=ypq(2,k)/ypq(4,k)       ! M'2=P2/Q2
        endif

        ! get distance between ray node and source
        x0=y(1,k)*cos(y(3,k))    ! ray node coordinates
        y0=y(1,k)*sin(y(3,k))
        diss=sqrt((x0-xs)**2+(y0-ys)**2)
        
        ! discontinuity?
        if(k>1) dr=abs(y(1,k)-y(1,k-1))

        ! near source. 
        if(k.ne.1.and.diss.lt.r0) then  
          ss=0.  ! ray length from source
          vela=0.  ! average velocity
          do kk=2,k
            ds=abs((y(1,kk-1)-y(1,kk))/cos(y(2,kk-1)))
            ss=ss+ds
            vela=vela+abs((rayvel(kk)+rayvel(kk-1))/2.*ds)
          enddo
          vela=abs(vela/ss)
          hmflim=1./ss/vela
          if(abs(hmf(1,k)).gt.hmflim) hmf(1,k)=hmflim  ! Q1~0 gives trouble
          if(abs(hmf(2,k)).gt.hmflim) hmf(2,k)=hmflim  ! Q2~0 gives trouble
        endif

c backward hessian M11=P3/Q3 and M22=P4/Q4 
        if (k.eq.mray.or.ypq(7,k).eq.0.0) then
          hmb(1,k)=sign(999.,ypq(5,k))
        else
          hmb(1,k)=-ypq(5,k)/ypq(7,k)       ! M"1=P3/Q3
        endif

        if (k.eq.mray.or.ypq(8,k).eq.0.0) then
          hmb(2,k)=sign(999.0,ypq(6,k))
        else
          hmb(2,k)=-ypq(6,k)/ypq(8,k)       ! M"2=P4/Q4
        endif                  

!       if(jdebug>1) then
!         write(13,*) k,hmf(1,k),hmf(2,k),hmb(1,k),hmb(2,k)
!       endif  

        ! get distance between ray node and receiver
        disr=sqrt((x0-xr)**2+(y0-yr)**2)

        ! near receiver
        if(k.ne.mray.and.disr.lt.r0) then  
          sr=0.  ! ray length from receiver
          vela=0.
          do kk=k,mray-1
            ds=abs((y(1,kk)-y(1,kk+1))/cos(y(2,kk+1)))
            sr=sr+ds
            vela=vela+abs((rayvel(kk)+rayvel(kk+1))/2.*ds)
          enddo
          vela=abs(vela/sr)
          hmblim=1./sr/vela
          if(abs(hmb(1,k)).gt.hmblim) hmb(1,k)=hmblim  ! Q3~0 gives trouble
          if(abs(hmb(2,k)).gt.hmblim) hmb(2,k)=hmblim  ! Q4~0 gives trouble
        endif

        ! check for accuracy against analytical expressions
        ! for nodes where Q1,Q2,Q1~,Q2~ are not zero
        if(jdebug.gt.0 .and. k.gt.1.and.k.lt.mray .and.
     &    ypq(3,k).ne.0.0.and.ypq(4,k).ne.0.0 .and.
     &    ypq(7,k).ne.0.0.and.ypq(8,k).ne.0.0) then  
          sip=sin(y(3,k))               ! sin(phi)
          sipi=sin(y(3,k)+y(2,k))       ! sin(phi+i)
          vk=abs(rayvel(k))
          rr=y(1,k)
          sii=sin(y(2,k))
          csi=cos(y(2,k))
          tai=sii/csi
          cti=csi/sii
          a53=rr*csi*ypq(1,k)+(p*p*dv(k)/rr-1.0/vk)*ypq(3,k) ! (A53)
          a53b=-rr*csi*ypq(5,k)-(p*p*dv(k)/rr-1.0/vk)*ypq(7,k) ! (A53) backward
          ! analytical expressions for M2,P2,Q2:
          am2=sipi/(vk*rr*sip)                               ! (A50) 
          ap2=rseg(1)*sipi/(p*vk)                            ! (A49) 
          aq2=rr*rseg(1)*sip/p                               ! (A49) 
          err=(hmf(2,k)-am2)/am2*100.
          sidp=sin(y3d-y(3,k))
          summ=p*sid/(rr*rr*sip*sidp)                        ! (A51) 
          am2b=summ-am2
          h22k=hmf(2,k)+hmb(2,k)
          h22k1=hmf(2,k-1)+hmb(2,k-1)
          csi1=cos(y(2,k-1))
          h11k=csi*csi*(hmf(1,k)+hmb(1,k))
          h11k1=csi1*csi1*(hmf(1,k-1)+hmb(1,k-1))
          er49P=min(99.,max(-99.,100.*(ap2-ypq(2,k))/ap2))
          er49Q=min(99.,max(-99.,100.*(aq2-ypq(4,k))/aq2))
          er50=min(99.,max(-99.,100.*(am2-hmf(2,k))/am2))
          er51=min(99.,max(-99.,100.*(h22k-summ)/summ))
          write(13,'(i5,f7.1,5f7.2,3e15.6)') k,y(1,k),vk,er49P,
     &          er49Q,er50,er51,a53,a53b,h11k
        endif

      end do

      return
      end 

      subroutine tracer(ystart,rseg,ktseg,kdwn,nseg,y,jsg,rmin,
     &      tstar,h,nray)

      ! traces ray specified by starting conditions in ystart and
      ! reflection/transmission/turning behaviour in rseg etc.

      ! input
      ! ystart(1)=r, ystart(2)=angle i, ystart(3)=delta, ystart(4)=time
      !        of the first ray point (usually 3 and 4 are zero)
      !        convention: angle i=0 for vertical up, 180 for down.
      ! rseg,ktseg,kdwn, for i=1,...,nseg describe nature of ray (see
      !        main program)
      ! h=integration step in km for Runge Kutta algorithm

      ! output
      ! y(1-4,i), i=1,...,nray: r,angle,delta,time of ray
      ! jsg(i), i=1,...,nray: segment number of node i, with - sign if
      !          the ray travels in the upward direction
      ! nray=number of ray nodes in y (<= 0 if error, eg when
      !      ray depth exceeds turning point depth in rseg)
      ! rmin=minimum radius reached by the ray
      ! tstar=traveltime/qray

      ! "error" messages:
      ! nray= 0: ray turns before reaching deeper segment
      !      -1: idem, but at discontinuity (not segment boundary)
      !      -2: reaches deepest segment level before turning
      !      -3: downward transmitted ray becomes evanescent at disc
      !      -4: upward transmitted ray becomes evanescent at discont
      !      -5: upward transmitted ray reflects (not segment boundary)
      !      -6: ray get closer than 10 km to Earth's center

      ! debug output if jdebug=1

      dimension rseg(20),ktseg(20),kdwn(20)
      parameter (LRAY=20000)
      dimension ystart(4),y(4,LRAY),dydx(4),jsg(LRAY)
      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/
      data epsdr/0.1/           ! discontinuity tolerance

      jdebug=0

      if(jdebug.eq.1) then
        write(13,*) 'call to tracer with ystart=',ystart
        write(13,*) 'r,k,kdwn=',rseg(1),ktseg(1),kdwn(1),' h=',h
      endif  
      do i=1,4
        y(i,1)=ystart(i)
      enddo
      iq=n
      call velo1(rseg(1),iq,1-kdwn(1),ktseg(1),csrc,dc,d2c)
      p=rseg(1)*sin(ystart(2))/csrc     ! slowness
      if(jdebug.gt.0) write(13,*) 'p,csrc=',p,csrc
c     modified statement to get reduced debug output for one ray:
      ! if(p>477.) jdebug=1      
      iseg=1
      nray=1
      iqstp=n
      jdown=kdwn(1)           ! 1 if ray heads downwards, 0 up
      rmin=y(1,1)
      jsg(1)=1
      if(jdown.eq.0) jsg(1)=-1
      epstp=0.                ! makes sure we pass a discontinuity

      if(jdebug.gt.0)
     &write(13,98)
98    format(2x,'nray',9x,'r',9x,'i',7x,'del',9x,'t jd',8x,'dr',5x,'iq')
      if(jdebug.gt.0)
     &write(13,99) nray,y(1,nray),y(2,nray)*r2d,y(3,nray)*r2d,
     &        y(4,nray),jdown,0.,iq,jsg(nray)

      ! the following do loop fills y for node nray+1:
      do while(iseg.lt.nseg)
        ! don't trespass rstop as determined by rstp:
        rstart=y(1,nray)-jdown*epstp
        call rstp(rstart,iqstp,jdown,rseg,iseg,rstop,jsegsw)   
        ! step h in the ray direction
        if(jdebug.gt.0)
     &  write(13,*) 'call to rk4l,iq,rstop,h=',iq,rstop,h
        call rk4l(y(1,nray),iq,h,ktseg(iseg),jdown,y(1,nray+1)) 
        dr=y(1,nray+1)-rstop    ! -distance to rstop
        if(jdown.eq.1) dr=-dr
        jsg(nray+1)=iseg
        if(jdown.eq.0) jsg(nray+1)=-iseg
        if(jdebug.gt.0)
     &  write(13,99) nray+1,y(1,nray+1),y(2,nray+1)*r2d,y(3,nray+1)*r2d,
     &        y(4,nray+1),jdown,dr,iq,jsg(nray)
99      format(i6,4f10.3,i3,f10.3,2i7)

c----------------------------------------------------------------
c
c  IF THE RAY HITS A DISCONTINUITY
c
c----------------------------------------------------------------

        if(dr.gt.-epsdr) then       ! near or past rstop, redo to rstop exactly
          epstp=epsdr               ! makes sure it doesn't get stuck here    
          hr=rstop-y(1,nray)        ! adjusted step size    
          if(jdebug.gt.0)
     &    write(13,*) 'call to rk4r,iqstp,rstop,hr=',iqstp,rstop,hr
          ! step hr in r direction
          call rk4r(y(1,nray),iq,hr,ktseg(iseg),jdown,y(1,nray+1)) 
          if(jdebug.gt.0)
     &    write(13,99) nray+1,y(1,nray+1),y(2,nray+1)*r2d,
     &          y(3,nray+1)*r2d,y(4,nray+1),jdown,hr,iq,jsg(nray)

          if(y(1,nray).lt.1.0) then     ! avoid division by r=0
            nray=-6
            if(jdebug.gt.0) write(13,*) 'Ray hits Earth center'
            return
          endif  
          ! apply continuity conditions

c----------------------------------------------------------------
c         AND THIS IS NOT THE END OF A SEGMENT
c----------------------------------------------------------------
          if(jsegsw.eq.0.) then         ! no segment change

            ! we filled y already to nray+1, at one side of the discnty
            nray=nray+2
            if(nray.gt.LRAY-2) then
              print *,'Ray too long for slowness ',p
              print *,'ystart(2)=',ystart(2),' rad'
              print *,'Last distance ',y(3,nray-1)*r2d
              stop 'A nray>LRAY-2'
            endif
            y(1,nray)=y(1,nray-1)       ! r continuous
            y(3,nray)=y(3,nray-1)       ! delta continuous
            y(4,nray)=y(4,nray-1)       ! time continuous
            jsg(nray)=iseg
            if(jdown.eq.0) jsg(nray)=-iseg
            call findiq(iq,1-jdown,rstop)  ! find lower node if downwards
            if(jdown.eq.1) then
              sini=p*vps(iq,ktseg(iseg))/rstop
              if(jdebug.gt.0) write(13,*) 'at nray,iq=',nray,iq,
     &              ' sini=',sini,rstop,vps(iq,ktseg(iseg)),p
              if(sini.ge.1) then        ! turning point
                iseg=iseg+1
                if(kdwn(iseg).ne.2) then        ! turning unexpected
                  nray=-1
                  return
                else
                  y(2,nray)=pi-y(2,nray-1)
                  jdown=0
                  jsg(nray)=-abs(jsg(nray))
                endif  
              else
                y(2,nray)=pi-asin(sini)   ! i=0 for vertical ray UP
              endif  
            else        ! case of upward ray
              sini=p*vps(iq,ktseg(iseg))/rstop
              if(abs(sini).ge.1.0) then   ! reflects back down!
                nray=-5
                if(jdebug.gt.0) write(13,*) 'p,v,rstop=',p,vps(iq,
     &                ktseg(iseg)),rstop,' nray=-5'
                return
              else  
                y(2,nray)=asin(sini)
              endif  
            endif
            if(jdebug.gt.0)
     &      write(13,99) nray,y(1,nray),y(2,nray)*r2d,
     &          y(3,nray)*r2d,y(4,nray),jdown,dr,iq,jsg(nray)

c-----------------------------------------------------------
c           BUT IF THIS DISCONTUITY IS ALSO END OF SEGMENT
c-----------------------------------------------------------

          else                   ! careful, next segment reached

            iseg=iseg+1
            if(jdebug.gt.0)
     &      write(13,*) 'new segment',iseg,rseg(iseg),kdwn(iseg),
     &                  ktseg(iseg)
            if(kdwn(iseg).eq.5) then            ! end of ray
              nray=nray+1
              jsg(nray)=iseg
              if(jdown.eq.0) jsg(nray)=-iseg
              if(nray.gt.LRAY-2) then
                print *,'ystart(2)=',ystart(2),' rad'
                print *,'Ray too long for slowness ',p
                print *,'Last distance ',y(3,nray-1)*r2d
                stop 'B nray>LRAY-2'
              endif
              if(jdebug.gt.0)
     &        write(13,*) 'end of ray nrry=',nray
              call gettstar(y,nray,jsg,ktseg,tstar,qray)
              return
            else if(kdwn(iseg).eq.2) then       ! turning point
              nray=-2
              if(jdebug.gt.0)
     &        write(13,*) 'turning point at disc, nray=-2'
              return
            else if(kdwn(iseg).eq.3) then       ! reflection
              nray=nray+2
              if(nray.gt.LRAY-2) then
                print *,'ystart(2)=',ystart(2),' rad'
                print *,'Ray too long for slowness ',p
                print *,'Last distance ',y(3,nray-1)*r2d
                stop 'C nray>LRAY-2'
              endif
              y(1,nray)=y(1,nray-1)     ! r continuous
              y(3,nray)=y(3,nray-1)     ! delta continuous
              y(4,nray)=y(4,nray-1)     ! time continuous
              jdown=1-jdown             ! reverse direction
              jsg(nray)=iseg
              if(jdown.eq.0) jsg(nray)=-iseg
              if(ktseg(iseg-1).eq.ktseg(iseg)) then   ! if no P/S conversion
                y(2,nray)=pi-y(2,nray-1)
              else
                sinin=sin(y(2,nray-1))
                call findiq(iq,1-jdown,rstop)
                vratio=vps(iq,ktseg(iseg-1))/vps(iq,ktseg(iseg))
                sinout=sinin/vratio
                if(abs(sinout).ge.1.) then        ! evanescent wave
                  nray=-3
                  if(jdebug.gt.0)
     &            write(13,*) 'evanescent wave'
                  return
                endif  
                y(2,nray)=asin(sinout)
                if(jdown.eq.1) y(2,nray)=pi-y(2,nray)
              endif
              if(jdebug.gt.0)
     &        write(13,99) nray,y(1,nray),y(2,nray)*r2d,
     &          y(3,nray)*r2d,y(4,nray),jdown,dr,iq,jsg(nray)
            else                                ! transmission
              nray=nray+2
              if(nray.gt.LRAY-2) then
                print *,'ystart(2)=',ystart(2),' rad'
                print *,'Ray too long for slowness ',p
                print *,'Last distance ',y(3,nray-1)*r2d
                stop 'D nray>LRAY-2'
              endif
              y(1,nray)=y(1,nray-1)     ! r continuous
              y(3,nray)=y(3,nray-1)     ! delta continuous
              y(4,nray)=y(4,nray-1)     ! time continuous
              jsg(nray)=iseg
              if(jdown.eq.0) jsg(nray)=-iseg
              call findiq(iq,1-jdown,rstop)
              if(jdown.eq.1) then
                sini=p*vps(iq,ktseg(iseg))/rstop
                if(abs(sini).ge.1.) then        ! evanescent wave
                  nray=-3
                  if(jdebug.gt.0)
     &            write(13,*) 'evanescent wave'
                  return
                endif  
                y(2,nray)=pi-asin(sini) ! i=0 for vertical ray UP
              else      ! case upward traveling ray
                sini=p*vps(iq,ktseg(iseg))/rstop
                if(abs(sini).ge.1.) then        ! evanescent wave
                  nray=-4
                  if(jdebug.gt.0)
     &            write(13,*) 'evanescent wave'
                  return
                endif  
                y(2,nray)=asin(sini)
                if(jdebug.gt.0)
     &          write(13,99), nray,y(1,nray),y(2,nray)*r2d,
     &            y(3,nray)*r2d,y(4,nray),jdown,dr,iq,jsg(nray)
              endif             ! for (jdown.eq.1) fork
            endif               ! for (kdwn.eq.5) fork
          endif                 ! for (jsegsw.eq.0) fork


c-----------------------------------------------------------
c
c       NORMAL INTEGRATION STEP (NO DISCONTINUITY)
c
c-----------------------------------------------------------

        ! we did not cross rstop (see fork on dr.ge.0)
        else 
          nray=nray+1
          epstp=0.
          if(nray.gt.LRAY-2) then
            print *,'ystart(2)=',ystart(2),' rad'
            print *,'Ray too long for slowness ',p
            print *,'Last distance ',y(3,nray-1)*r2d
            stop 'E nray>LRAY-2'
          endif

          ! finally, check if we passed the turning point
          if(jdown.eq.1.and.y(2,nray).lt.halfpi) then
            jdown=0
            iseg=iseg+1
            if(jdebug.gt.0)
     &      write(13,*) 'TP: new segment',iseg,rseg(iseg),kdwn(iseg),
     &                  ktseg(iseg)
            if(kdwn(iseg).ne.2) then  ! ray was not supposed to turn!
              nray=0
              if(jdebug.gt.0)
     &        write(13,*) 'Ray was not supposed to turn, nray=0'
              return
            endif  
            if(kdwn(iseg).eq.2.and.rseg(iseg+1).lt.y(1,nray)) then
              nray=0
              if(jdebug.gt.0)
     &        write(13,*) 'Ray turns before reaching reflection, nray=0'         
              return
            endif
          endif                 ! for turning point fork
          jsg(nray)=iseg
          if(jdown.eq.0) jsg(nray)=-iseg
        endif                   ! for (dr.ge.0) fork
        rmin=min(rmin,y(1,nray))
      end do                    ! for (iseg.lt.nseg) loop

      print *, 'Hmmm - you should never have reached this line ....'
      stop 'software error in tracer'

      return
      end

      subroutine rdray(kin,kout,chray,phase,rseg,ktseg,kdwn,nseg)

      ! read ray definition

      ! input:
      ! kin - I/O unit number for input (5 for screen) [assumed open]
      ! kout - I/O unit number for output (0 for no outp) [assumed open]
      ! output:
      ! chray - ray identification (e.g. PKiKP) from first line of kin
      ! phase - (unused) ray ident for ellipticity corrections (max 8 char)
      ! rseg(i) - start radius of segment i
      ! ktseg(i) - wave type of segment i (1=P,2=SV,3=SH,last point of ray 
      !            may be 0, but will be set to 1 or 2 on output)
      ! kdwn(i) - ray direction (1=down,0=up,2=turning point,3=refl,
      !           4=transmission, 5=last point of ray) from node i
      ! nseg - number of segments

      ! the ray is defined by segments and their starting radii
      ! the first (source) one specifies: radius, up/down, P/S
      ! for each subsequent segment: radius,refl/trans/turn,P/S
      ! and for the last point of the ray: (radius,5,0))
      ! note that a turning point separates two different segments
      ! (the radius is then interpreted as the MINIMUM radius for 
      ! the turning point).

      ! If tableflag=0 in main program, the source depth is adjusted
      ! to the actual source depth in the data file.

      ! Input examples  from unit kin (without two header lines, segments only):
      ! P wave, source at 10 km depth:
      ! 6361 1 1
      ! 3480 2 1
      ! 6371 5 0
      ! PcS wave, source 400 km depth:
      ! 5971 1 1
      ! 3480 3 2
      ! 6371 5 0
      ! pPKP wave, 100 km depth:
      ! 6271 0 1
      ! 6371 3 1
      ! 3480 4 1
      ! 1221.5 2 1
      ! 3480 4 1
      ! 6371 5 0


      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh
      dimension rseg(20),ktseg(20),kdwn(20)
      character*72 chray
      character*2 wtyp(0:3)
      character*4 ktyp(0:5)
      character*8 phase

      data wtyp/' ','P','SV','SH'/ 
      data ktyp/'up','down','turn','refl','tran','end'/

      jdebug=0

      rcmb=-10.
      ricb=-10.
      if(noc.gt.0) rcmb=r(noc)
      if(nic.gt.0) rinc=r(nic)

      if(kin.eq.5) then
        print *,'Give ray description, then ellipticity ident, then'
        print *,'Input of ray segments (one per line):'
        print *,'Give radius (km), kdwn, wavetype (P=1,S=2),'
        print *,'where kdwn=1 for downward start, 0 for upward start,'
        print *,'2 for max turning depth, 3 for reflection, 4 for'
        print *,'transmission, 5 for end of ray'
      endif

      ifSH=0            ! (currently ignored, model is either SH or SV)
      ierror=0
      i=1
      read(kin,fmt='(a)') chray
      print *,'Ray ident: ',chray
      if(kout.ne.0) write(kout,fmt='(a)') chray
      read(kin,fmt='(a)') phase
      if(kout.ne.0) write(kout,fmt='(a)') phase
      read(kin,*) rseg(i),kdwn(i),ktseg(i)

      ! adjust to CMB or ICB if close (avoids having to write new defs)
      if(abs(rseg(i)-rcmb).lt.10.) rseg(i)=rcmb
      if(abs(rseg(i)-rinc).lt.10.) rseg(i)=rinc

      if(ktseg(i).eq.3) ifSH=1
      if(kout.ne.0) 
     &      write(kout,40) rseg(i),kdwn(i),ktseg(i),ktyp(kdwn(i)),
     &      wtyp(ktseg(i))
40    format(f8.1,2i5,2x,a4,2x,a2)
      if(ktseg(i).lt.1.or.ktseg(i).gt.3) stop 'wave type must be 1,2,3'
      if(kdwn(1).gt.1) then
        print *,'First ray point should have kdwn 1 or 0'
        if(kout.ne.0) write(kout,*) 
     &        'First ray point should have kdwn 1 or 0'
        ierror=ierror+1
      endif  
      rmin=rseg(1)
      rmax=rmin
      do while(kdwn(i).ne.5)
        i=i+1
        if(i.gt.20) stop 'Max segments is 20'
        read(kin,*) rseg(i),kdwn(i),ktseg(i)

        if(abs(rseg(i)-rcmb).lt.10.) rseg(i)=rcmb
        if(abs(rseg(i)-rinc).lt.10.) rseg(i)=rinc

        if(ktseg(i).eq.0) ktseg(i)=ktseg(i-1)
        if(ktseg(i).lt.1.or.ktseg(i).gt.3) stop 'type must be 1,2,3'
        if(ktseg(i).eq.3) ifSH=1
        if(kout.ne.0) 
     &        write(kout,40) rseg(i),kdwn(i),ktseg(i),ktyp(kdwn(i)),
     &        wtyp(ktseg(i))
        if(kdwn(i).lt.2) then
          print *,'kdwn < 2 only allowed for first segment'
          if(kout.ne.0) write(kout,*) 
     &          'kdwn < 2 only allowed for first segment'
          ierror=ierror+1
        endif  
        if(kdwn(i).eq.3.or.kdwn(i).eq.4) then
          call findiq(iq,1,rseg(i))     ! choose upper layer at discnty
          if(iq.lt.2.or.iq.gt.n) then
            print *,'Radius outside of model, segment',i
            if(kout.ne.0) write(kout,*) 'Radius outside of model'
            ierror=ierror+1
          else if(r(iq).lt.6370.0.and.abs(r(iq)-r(iq-1)).gt.0.01.or.
     &      abs(rseg(i)-r(iq)).gt.0.01) then
            print *,'Refl/transmission not at discontinuity, segment',i
            print *,'rseg=',rseg(i),'iq, r(iq)=',iq,r(iq)
            if(kout.ne.0) write(kout,*) 
     &            'Refl/transmission not at discontinuity'
            ierror=ierror+1
          endif
        endif  
            
        rmin=min(rmin,rseg(i))
        rmax=max(rmax,rseg(i))
        if(jdebug.eq.1) write(13,*) rseg(i),kdwn(i),ktseg(i),
     &      ktyp(kdwn(i)), wtyp(ktseg(i))
      end do
      if(jdebug.eq.1) write(13,*) 'nseg=',i
      nseg=i
      print *,'Rays between r=',rmin,rmax
      
      if(ifSH.eq.1) then       ! check if all are SH if one segment is
        do i=1,nseg
          if(ktseg(i).ne.3) then
            print *,'Layer ',i,' type=',ktseg(i),' but elsewhere 3 (SH)'
            ierror=ierror+1
          endif
        enddo
      endif  

      if(ierror.gt.0) then
        print *,'Stop because of',ierror,' errors in ray definition'
        stop
      endif  

      return
      end

      subroutine rdband(kin,kout,bndomega,dotm,nfreq,nband)

      ! reads spectral band info
      ! input: I/O unit numbers for input and output kin,kout
      !        (assumed open, no output if kout=0).
      ! output: bndomega=frequencies (rad/s), dotm=mdot(omega), nfreq=# of frequencs

      parameter (NFR=16)
      dimension bndomega(500,NFR),dotm(500,NFR),nfreq(NFR)

      read(kin,*) nband
      if(nband.gt.NFR) stop 'nband > NFR, increase dimensions'
      if(kout.ne.0) write(kout,fmt='(i5)') nband
      do ib=1,nband
        read(kin,*) nfreq(ib)
        if(nfreq(ib).gt.500) stop 'nfreq>500, increase dimensions'
        if(kout.ne.0) write(kout,fmt='(i5)') nfreq(ib)
        do j=1,nfreq(ib)
          read(kin,*) bndomega(j,ib),dotm(j,ib)
          if(kout.ne.0) 
     &          write(kout,fmt='(2f8.4)') bndomega(j,ib),dotm(j,ib)
        end do
      end do

      return
      end


      subroutine rstp(rstart,iq,jdown,rseg,iseg,rstop,jsegsw)   
      dimension rseg(20),ktseg(20),kdwn(20)
      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data eps/0.01/

      jdebug=0
      ! if(rstart<3480.) jdebug=1

      ! determines if discontinuity or segment will be crossed
      ! or if end of segment depth is reached

      ! input: rstart=current radius, 
      !        jdown=direction (1 down, else 0)
      !        iq=guess for layer number belonging to rstart
      !        rseg=segment start radii, iseg=current segment #
      ! output: rstop=r limit
      !        jsegsw=0 if simply a discontinuity (ie transmission)
      !        jsegsw=1 if the r limit is also a segment boundary

      rr=rstart
      call findiq(iq,1-jdown,rr)  ! r(iq) is directly below (or =) rr
      rs=rseg(iseg+1)           ! start r of next segment
      if(jdebug.gt.0) write(13,*) 'rstp called:',rstart,iq,jdown,rs
      i=iq
      if(jdebug.gt.0) write(13,*) 'rstp at i,r=',i,r(i-1),r(i),r(i+1)
      ! find first discontinuity before r(i) reaches rs:
      ! (step down or up depending on jdown)
      if(jdown.eq.1) then       ! downward ray
        do while((i.gt.1).and.r(i)-r(i-1).gt.eps.and.r(i)-rs.gt.0.)
          i=i-1
          if(jdebug.gt.0) write(13,*) i,r(i),' stepping down'
        enddo  
      else  
        i=min(n,iq+1)
        do while((i.lt.n).and.r(i+1)-r(i).gt.eps.and.rs-r(i).gt.eps)
          i=i+1
          if(jdebug.gt.0) write(13,*) i,r(i),' stepping up'
        enddo  
      endif  
      jsegsw=0
      rstop=r(i)
      if(jdown.eq.1.and.rstop.lt.rs) rstop=rs
      if(jdown.eq.0.and.rstop.gt.rs) rstop=rs
      if(abs(rstop-rs).lt.eps) jsegsw=1
      if(jdebug.gt.0) write(13,*) 'rstop =',rstop

      return
      end

      subroutine velo1(rr,iq,kupr,ktype,c,dc,d2c)

      ! input: 
      ! rr=radius 
      ! iq is guess for model node directly beneath or at rr
      !    (warning: iq is also output and may be changed, so do 
      !    never use a numerical value like <1> or a common block 
      !    variable like <n> or <noc> etc!)
      ! kupr=1 if the upper layer of a discontinuity must be selected,
      !      0 if the lower one. kupr is ignored if r is not on a
      !      discontinuity.
      ! ktype=1 for P, 2 for S wave
      ! common block /modl/ must be filled (see subroutine model.f)

      ! output: 
      ! iq= node # directly below or at rr 
      !     (is modified if wrong on input!)
      ! c=velocity at radius rr
      ! dc= dc/dr
      ! d2c = d2c/dr2 (second derivative)

      ! method: spline interpolation

      !        The ambiguity if rr is at a discontinuity is decided 
      !        by kupr:
      !        for example, if kdown=1 (downgoing ray) and the
      !        point of interest is at the end of the ray segment,
      !        you shall wish to choose the top of the layer (i.e.
      !        kupr=1=kdown). If it is at the start of the segment
      !        choose kupr=0=1-kdown. Similarly, if the ray goes
      !        up (kdown=0) then kupr=0=kdown at the end, and
      !        kupr=1=1-kdown at the start. Note the expressions
      !        in terms of kdown are the same for up and downgoing rays
      !        i.e.:
      !        at start call with          1-kdown
      !        and end call with           kdown


      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data eps/0.01/            ! precision 10 m

      jdebug=0

      if(jdebug.gt.0) write(13,*) 'velo1 called:',rr,iq,kupr,ktype

      ! find node immediately below r
      call findiq(iq,kupr,rr)

      t=rr-r(iq)                ! distance from node
      if(ktype.eq.1) then
        c=vp(iq)+t*(qvp(1,iq)+t*(qvp(2,iq)+t*qvp(3,iq)))

        ! debug: try linear interpolation
        dr=r(iq+1)-r(iq)
        if(dr>0.) c=vp(iq)+t*(vp(iq+1)-vp(iq))/dr

        dc=qvp(1,iq)+2.*qvp(2,iq)*t+3*qvp(3,iq)*t*t
        d2c=2*qvp(2,iq)+6*t*qvp(3,iq)
      else  
        c=vs(iq)+t*(qvs(1,iq)+t*(qvs(2,iq)+t*qvs(3,iq)))
        if(c.le.0.) then        ! avoid problem in tracer
          c=0.01
          dc=0
          d2c=0
          return
        endif
        dc=qvs(1,iq)+2.*qvs(2,iq)*t+3*qvs(3,iq)*t*t
        d2c=2*qvs(2,iq)+6*t*qvs(3,iq)
      endif  

      return
      end

      subroutine velo2(rr,iq,kupr,a,b,d)

      ! input: 
      ! rr=radius - must be lower node of discontinuity
      ! iq is guess for model node directly beneath or at rr
      !    (warning: iq is also output and may be changed, so do 
      !    never use a numerical value like <1> or a common block 
      !    variable like <n> or <noc> etc!)
      ! kupr=1 if the upper layer of a discontinuity must be selected,
      !      0 if the lower one. kupr is ignored if r is not on a
      !      discontinuity.
      ! common block /modl/ must be filled (see subroutine model.f)

      ! output: 
      ! iq= node # directly below or at rr 
      !     (is modified if wrong on input!)
      ! a,b,d=Vp,Vs,density


      !        The ambiguity if rr is at a discontinuity is decided 
      !        by kupr:
      !        for example, if kdown=1 (downgoing ray) and the
      !        point of interest is at the end of the ray segment,
      !        you shall wish to choose the top of the layer (i.e.
      !        kupr=1=kdown). If it is at the start of the segment
      !        choose kupr=0=1-kdown. Similarly, if the ray goes
      !        up (kdown=0) then kupr=0=kdown at the end, and
      !        kupr=1=1-kdown at the start. Note the expressions
      !        in terms of kdown are the same for up and downgoing rays
      !        i.e.:
      !        at start call with          1-kdown
      !        and end call with           kdown


      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      jdebug=0

      if(jdebug.gt.0) write(13,*) 'velo2 called:',rr,iq,kupr

      ! find node at r
      iqin=iq                   ! for debug only
      call findiq(iq,kupr,rr)

      t=rr-r(iq)                ! distance from node
      if(abs(t)>0.1) then          ! for debug only
        write(13,*) 'WARNING: velo2 called at rr=',rr,', t=',t
        write(13,*) 'iqin=',iqin,', out=',iq,' kupr=',kupr
        print *,'BUG in velo2, check fort.13'
        stop
      endif  

      a=vp(iq)
      b=vs(iq)
      d=rho(iq)

      if(jdebug.gt.0) write(13,*) 'velo2 output:',iq,a,b,d

      return
      end

      subroutine findiq(iq,kupr,rr)

      ! input: iq=guess for layer number
      !        kupr= handles discontinuity (see below)
      !        rr=radius
      ! output: iq= layer number at or directly below rr. 

      !        The ambiguity if rr is at a discontinuity is decided 
      !        by kupr:
      !        if rr is at a discontinuity and:
      !           kupr=1, iq is layer number of the upper layer
      !           kupr=0, iq is layer number of the lower layer
      !        for example, if kdown=1 (downgoing ray) and the
      !        point of interest is at the end of the ray segment,
      !        you shall wish to choose the top of the layer (i.e.
      !        kupr=1=kdown). If it is at the start of the segment
      !        choose kupr=0=1-kdown. Similarly, if the ray goes
      !        up (kdown=0) then kupr=0=kdown at the end, and
      !        kupr=1=1-kdown at the start. Note the expressions
      !        in terms of kdown are the same for up and downgoing rays
      !        i.e.:
      !        at start call with          1-kdown
      !        and end call with           kdown


      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data eps/0.01/

      jdebug=0

      if(jdebug.gt.0) write(13,*) 'findiq called: ',iq,kupr,rr

      ! find node immediately below rr
      iq=min(n,max(iq,1))
      if(jdebug.gt.0) write(13,*) 'corrected iq=',iq
      do while(iq.gt.1.and.r(iq).gt.rr+eps)
        if(jdebug>1) write(13,*) 'iq-1 for',iq,r(iq),'>',rr+eps
        iq=iq-1
      enddo
      do while(iq.lt.n.and.r(iq+1).lt.rr+eps)
        if(jdebug>1) write(13,*) 'iq+1 for',iq,r(iq+1),'<',rr+eps
        iq=iq+1
      enddo  

      if(jdebug.gt.0) write(13,*) 'final iq=',iq

      if(jdebug>0 .and. rr-r(iq).ge.eps) write(13,*) 'no disc iq=',iq
      if(rr-r(iq).ge.eps) return

      ! resolve discontinuity ambiguitiy
      if(jdebug>0) write(13,*) 'hit disc:',r(iq+1),r(iq)
      if(iq.lt.n.and.abs(r(iq+1)-r(iq)).lt.eps.and.kupr.eq.1) 
     &          iq=iq+1
      if(iq.gt.1.and.abs(r(iq)-r(iq-1)).lt.eps.and.kupr.eq.0)
     &          iq=iq-1 


      if(jdebug.gt.0) write(13,*) 'adjusted for discontinuity iq=',iq

      return
      end

      function vps(i,ktype)
      ! returns vp at node i if ktype=1, else vs
      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      if(ktype.eq.1) then
        vps=vp(i)
      else
        vps=vs(i)
      endif
      return
      end

      subroutine model(fname)

      ! Reads a model from file <fname>
      ! Model units are assumed to be MKS, unless the first density
      ! read is less than 100, in which case km and km/s units are
      ! assumed; final model is in km and km/s.
      ! This is the same file format as used for love.f and rayleigh.f

      ! Format of the model file [unit 1, opened/closed by the routine]:
      ! Line 1: ignored
      ! Line 2: ignored
      ! Line 3: n (# of nodes), nic (top node of inner core), noc (top
      !         node of outer core), nmoh (bottom[!] node of the crust)
      ! Lines 4ff: radius, density, Vp, Vs, Qs

      ! Anisotropy: the model file should either have SV or (more likely)
      ! SH. You'd have to run the program twice on different model files
      ! to get SV and SH travel times. Azimuthal anisotropy cannot be
      ! handled.

      ! The nodes are ordered with INCREASING radius
      ! Discontinuities are indicated by two subsequent nodes
      ! with the same radius. Nodes closer than 0.01 km are collapsed
      ! into a discontinuity.

      ! input: model file name fname (character*72)
      ! output: common block /modl/ is filled with model information
      !         including spline coefficients.

      character*72 fname
      dimension wrk(3,999)
      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      jdebug=0

      nmoh=0            ! safeguard against old model files
      open(1,file=fname)
      read(1,*) 
      read(1,*) 
      read(1,*) n,nic,noc,nmoh
      if(n.gt.999) stop 'Increase model dimension'
      if(nmoh.le.noc.or.nmoh.gt.n) stop 'model: Error in Moho index'
      mks=1
      read(1,*) r(1),rho(1),vp(1),vs(1)
      if(rho(1).lt.100.) mks=0
      backspace(1)
      do i=1,n
        read(1,*) r(i),rho(i),vp(i),vs(i),Qs(i)
        if(rho(i).le.0.) stop 'subroutine model: rho=0 not allowed'
        if(mks.eq.1) then
          r(i)=0.001*r(i)
          vp(i)=0.001*vp(i)
          vs(i)=0.001*vs(i)
          rho(i)=0.001*rho(i)
        endif
        if(jdebug.gt.0) write(13,*) i,r(i),vp(i),vs(i),rho(i),Qs(i)
      end do  
      close(1)

   60 nsl=n   
      if(vs(nsl).gt.0.) go to 70
   65 nsl=nsl-1
      if(vs(nsl).le.0.) go to 65
   70 nicp1=nic+1
      nocp1=noc+1
      nslp1=nsl+1
      if(r(nic).ne.r(nicp1)) stop 'subroutine model: error in nic'
      if(r(noc).ne.r(nocp1)) stop 'subroutine model: error in noc'
      if(r(nmoh).ne.r(nmoh-1)) stop 'subroutine model: error in nmoh'

!     if(nsl.ne.n) stop 'Oceanic reference model not allowed'

c*** spline ***
      rn=r(n)
!!    do 45 i=1,n
!!    if(i.gt.1.and.abs(r(i)-r(i-1)).lt.0.01) r(i)=r(i-1)
! Intel gives the error message: Subscript #1 of the array R has value 0 which is less than the lower bound of 1
! g77 doesn't have this problem
      do 45 i=2,n
      if(abs(r(i)-r(i-1)).lt.0.01) r(i)=r(i-1)
   45 continue
      call drspln(1,n,r,vp,qvp,wrk)
      call drspln(1,n,r,vs,qvs,wrk)
      call drspln(1,n,r,Qs,qqs,wrk)
      call drspln(1,n,r,rho,qrh,wrk)

      if(jdebug.eq.0) return

      write(13,*) ' i    r    qvp'
      do i=1,n
        write(13,*) i,r(i),qvp(1,i),qvp(2,i),qvp(3,i)
      enddo

! added 19/11/06 ----------------------
      if(jdebug.eq.0) return

      ! write GMT files for debugging
      open(41,file='vp.xy')
      open(42,file='dvp.xy')
      open(43,file='ddvp.xy')
      open(44,file='vs.xy')
      open(45,file='dvs.xy')
      open(46,file='ddvs.xy')

      iq=1
      do i=2,n
        step=0.1000001*(r(i)-r(i-1))
        step=max(1.0,step)              ! avoid step=0
        rr=r(i-1)
        do while(rr.lt.r(i))
          call velo1(rr,iq,0,1,c,dc,d2c)        ! P velocity
          write(41,50) rr,c
          write(42,51) rr,dc
          write(43,51) rr,d2c
          if(i.gt.noc) then
            call velo1(rr,iq,0,2,c,dc,d2c)      ! S velocity
            write(44,50) rr,c
            write(45,51) rr,dc
            write(46,51) rr,d2c
          endif
          rr=rr+step
        enddo
      enddo
      close(41)
      close(42)
      close(43)
      close(44)
      close(45)
      close(46)
50    format(f8.1,f8.3)
51    format(f8.1,e12.3)
!-----------------------------------------

      return
      end        

      subroutine drspln(i1,i2,x,y,q,f)

c   rspln computes cubic spline interpolation coefficients
c   for y(x) between grid points i1 and i2 saving them in q.  the
c   interpolation is continuous with continuous first and second
c   derivitives.  it agrees exactly with y at grid points and with the
c   three point first derivitives at both end points (i1 and i2).
c   x must be monotonic but if two successive values of x are equal
c   a discontinuity is assumed and seperate interpolation is done on
c   each strictly monotonic segment.  the arrays must be dimensioned at
c   least - x(i2), y(i2), q(3,i2), and f(3,i2).  f is working storage
c   for rspln.
c                                                     -rpb
      dimension x(i2),y(i2),q(3,i2),f(3,i2),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0./
      j1=i1+1
      y0=0.
c   bail out if there are less than two points total.
      if(i2-i1)13,17,8
 8    a0=x(j1-1)
c   search for discontinuities.
      do 3 i=j1,i2
      b0=a0
      a0=x(i)
      if(a0-b0)3,4,3
 3    continue
 17   j1=j1-1
      j2=i2-2
      go to 5
 4    j1=j1-1
      j2=i-3
c   see if there are enough points to interpolate (at least three).
 5    if(j2+1-j1)9,10,11
c   only two points.  use linear interpolation.
 10   j2=j2+2
      y0=(y(j2)-y(j1))/(x(j2)-x(j1))
      do 15 j=1,3
      q(j,j1)=yy(j)
 15   q(j,j2)=yy(j)
      go to 12
c   more than two points.  do spline interpolation.
 11   a0=0.
      h=x(j1+1)-x(j1)
      h2=x(j1+2)-x(j1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
c   calculate derivitive at near end.
      b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
      b1=b0
c   explicitly reduce banded matrix to an upper banded matrix.
      do 1 i=j1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.*a0
      h3a=2.*h-3.*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.*y0*ha)/(h*h3a)
      a0=q(3,i)
 1    b0=f(3,i)
c   take care of last two rows.
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
c   calculate derivitive at far end.
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
c   solve upper banded matrix by reverse iteration.
      do 2 j=j1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
 2    i=k
      q(1,i)=b1
c   fill in the last point with a linear extrapolation.
 9    j2=j2+2
      do 14 j=1,3
 14   q(j,j2)=yy(j)
c   see if this discontinuity is the last.
 12   if(j2-i2)6,13,13
c   no.  go back for more.
 6    j1=j2+2
      if(j1-i2)8,8,7
c   there is only one point left after the latest discontinuity.
 7    do 16 j=1,3
 16   q(j,i2)=yy(j)
c   fini.
 13   return
      end

c Eq numbers from Dahlen et al (GJI,2000,appendix)

      subroutine derivs(y,dydx,ir,kdis,kt,c,dc,d2c)

      ! computes derivatives of r,i,phi and t w.r.t. l

      dimension dydx(4),y(4)

      ! input:
      ! y(1)=r, y(2)=i
      ! ir: layer index
      ! kdis: 1 if ray approaches node ir from above, else 0.
      ! kt: wave type (1=P, 2=S)

      ! output
      ! dydx: derivatives of y w.r.t. ray distance dl
      ! c: velocity
      ! dc: dc/dr
      ! dc2: d^2c/dr^2

c calculate velocity and dv/dr at a given r
      call velo1(y(1),ir,kdis,kt,c,dc,d2c)

      sii=sin(y(2))
      coi=cos(y(2))
      ri=1.0/max(0.01,y(1))

      dydx(1)=coi                       ! dr/dl (A5)
      dydx(2)=sii*(dc/c-ri)             ! di/dl
      dydx(3)=sii*ri                    ! dphi/dl
      dydx(4)=1.0/c                     ! dT/dl

      return
      end

      subroutine derivr(y,dydx,ir,kdis,kt,c,dc,d2c)

      ! derivatives of r,i,phi and t w.r.t. r (not l!)

      dimension dydx(4),y(4)

      ! input:
      ! y(1)=r, y(2)=i, y(3)=phi y(4)=t
      ! ir: layer index
      ! kdis: 1 if ray approaches node ir from above, else 0.
      ! kt: wave type (1=P, 2=S)

      ! output
      ! dydx: derivatives of y w.r.t. ray distance dl
      ! c: velocity
      ! dc: dc/dr
      ! dc2: d^2c/dr^2

      ! debug - remove
c     if(y(1).lt.3484.0) then
c       print *, 'deriv is in OC!, y(1)=',y(1)
c     endif  

c calculate velocity and dv/dr at a given r
      rr=y(1)
      call velo1(rr,ir,kdis,kt,c,dc,d2c)

      tgi=tan(y(2))
      coi=cos(y(2))
      ri=1.0/max(0.01,y(1))

      dydx(1)=1.0                       ! dr/dr 
      dydx(2)=tgi*(dc/c-ri)             ! di/dr
      dydx(3)=tgi*ri                    ! dphi/dr
      dydx(4)=1.0/(c*coi)               ! dT/dr

      return
      end

      subroutine rk4l(y,iq,h,ktype,jdown,yout)

      ! integrates r,phi,i,t as function of l

      ! 4th order Runge Kutta solution to y'=f(x,y):
      ! y(n+1)=y(n)+(k1+2k2+2k3+k4)/6

      ! input:
      ! y(1)=r, y(2)=i, y(3)=phi, y(4)=t
      ! iq=guess for node # directly below y(1)
      ! h: step size (km)
      ! ktype=1 for P, 2 for S
      ! jdown=1 if step is decreasing r, 0 if upward ray

      ! output
      ! yout(1)-yout(4): next ray point
      ! iq is adjusted if needed

      dimension dydx(4),y(4),yout(4)
      dimension dym(4),dyt(4),yt(4)
      data halfpi/1.570796326795/

      jdebug = 0
      hh=h*0.5
      h6=h/6.0
      n=4

      kupr=jdown      ! choose proper layer if y(1) at discontinuity
      call derivs(y,dydx,iq,1-kupr,ktype,c,dc,d2c)
      if(jdebug.gt.0) then
        write(13,*) 'Enter rk4l, y,h=',y,h
        write(13,97)
        write(13,96) 'k1',y(1),y(2),dydx(1),dydx(2),dydx(3),dydx(4),c,dc
      endif
      do  i=1,n
        yt(i)=y(i)+hh*dydx(i)  
      end do

      call derivs(yt,dyt,iq,kupr,ktype,c,dc,d2c)
      if(jdebug.gt.0)
     &  write(13,96) 'k2',yt(1),yt(2),dyt(1),dyt(2),dyt(3),dyt(4),c,dc
      do i=1,n
        yt(i)=y(i)+hh*dyt(i)    
      enddo
      call derivs(yt,dym,iq,kupr,ktype,c,dc,d2c)
      if(jdebug.gt.0)
     &  write(13,96) 'k3',yt(1),yt(2),dym(1),dym(2),dym(3),dym(4),c,dc

      do i=1,n
        yt(i)=y(i)+h*dym(i)   
        dym(i)=dyt(i)+dym(i)
      end do
      call derivs(yt,dyt,iq,kupr,ktype,c,dc,d2c)
      if(jdebug.gt.0)
     &  write(13,96) 'k4',yt(1),yt(2),dyt(1),dyt(2),dyt(3),dyt(4),c,dc
96    format(a2,8f10.4)
97    format(11x,'r',9x,'i',5x,'dr/dl',5x,'di/dl',5x,'dD/dl',
     &  5x,'dt/dl',9x,'c',5x,'dc/dr')

      do  i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
      enddo

      if(jdebug>0) write(13,*) 'yout=',yout

      return
      end

      subroutine rk4r(y,iq,h,ktype,jdown,yout)

      ! integrates l,phi,i,t as function of r

      ! 4th order Runge Kutta solution to y'=f(x,y):
      ! y(n+1)=y(n)+(k1+2k2+2k3+k4)/6

      ! input (note y(1) is not r but l!):
      ! y(1)=l, y(2)=i, y(3)=phi, y(4)=t
      ! iq=guess for node # directly below y(1)
      ! h: step size (km)
      ! ktype=1 for P, 2 for S
      ! jdown=1 if step is decreasing r, 0 if upward ray

      ! output
      ! yout(1)-yout(4): next ray point
      ! iq is adjusted if needed

      dimension dydx(4),y(4),yout(4)
      dimension dym(4),dyt(4),yt(4)

      hh=h*0.5
      h6=h/6.0
      n=4

      kupr=jdown      ! choose proper layer if discontinuity
      call derivr(y,dydx,iq,1-kupr,ktype,c,dc,d2c)
      do  i=1,n
        yt(i)=y(i)+hh*dydx(i)  
      end do

      call derivr(yt,dyt,iq,kupr,ktype,c,dc,d2c)

      do i=1,n
        yt(i)=y(i)+hh*dyt(i)    
      enddo

      call derivr(yt,dym,iq,kupr,ktype,c,dc,d2c)

      do i=1,n
        yt(i)=y(i)+h*dym(i)   
        dym(i)=dyt(i)+dym(i)
      end do

      call derivr(yt,dyt,iq,kupr,ktype,c,dc,d2c)

      do  i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
      enddo

      return
      end

      subroutine mkrtable(rseg,ktseg,kdwn,nseg,a1,a2,ytbl,ntbl)

      ! traces rays and produces ray table ytbl

      ! input:
      ! ray geometry in rseg,ktseg,kdwn,nseg (see subroutine rdray)
      ! a1,a2 start/final ray path angle at source (in degrees,
      !    0<= a1 < a2 < 90 if upward, else 90 < a1 < a2 < 180)
      ! a1 and a2 override ray segment angles. Set (0,90) or (90,180)
      !    to let tracer decide what is the full sweep of angles
      ! model in common /modl/ (see subroutine model)

      ! output:
      ! ytbl(1,j)=i, ytbl(2,j)=Delta, ytbl(3,j)=T,  ytbl(4,j)=rturn
      ! for j=1,ntbl, with i and Delta in degrees, T in seconds
      ! a1,a2 - adjusted if out of bounds dictated by the ray

      ! shortcoming:
      ! p and P have to be treated with different ray definition
      ! files

      dimension rseg(20),ktseg(20),kdwn(20),ytbl(4,4000)
      parameter (LRAY=20000)
      dimension ystart(4),y(4,LRAY),jsg(LRAY)
      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/, eps/1.0/

      jdebug=0

      ! careful... step size h must be same in mkrtable2 and main program
      h=10.    ! integration step size (keep equal to h in calling program)
      if(a1.gt.a2) stop 'mkrtable: a1>a2'

      iqsrc=n
      ktsrc=ktseg(1)
      rsrc=rseg(1)
      call velo1(rseg(1),iqsrc,1-kdwn(1),ktsrc,csrc,dc,d2c) ! find csrc

      ! check max slowness to avoid that upward ray turns before surface
      ! TODO: this does not yet work if a mode conversion occurs
      ! before the surface and the converted wave has the larger
      ! angle with the vertical. 
c     vlim=vps(iqsrc,ktsrc)
      plim=0.999*rsrc/csrc
      do i=iqsrc+1,n
c       vlim=max(vlim,vps(i,ktsrc))
        vr=vps(i,ktsrc)
        plim=min(plim,r(i)/vr)
      enddo
c     plim=rsrc*0.99/vlim
      vsurf=vps(n,ktsrc)
c     alim=asin(vsurf*plim/r(n))
      alim=asin(csrc*plim/rsrc)
c     if(jdebug.eq.1) write(13,*) 'vlim=',vlim,'  alim=',alim*r2d
      if(jdebug.eq.1) write(13,*) 'plim=',plim,'  alim=',alim*r2d

      if(jdebug.gt.0) write(13,*) 'rsrc, csrc=',rsrc,csrc
      ra1=a1/r2d
      ra2=a2/r2d
      if(kdwn(1).eq.0) then
        amin=max(0.00001,ra1)   ! upward: 0<angle<halfpi
        amax=min(alim,ra2)
      else
        amin=max(pi-alim,ra1)    ! down: halfpi<angle<pi
        amax=min(pi-0.00001,ra2)
      endif  
      if(jdebug.eq.1) write(13,*) 'amin,amax=',amin*r2d,amax*r2d

      ! find slowness range for defined turning points.
      pturn=0.
      ptrns=(rsrc-eps)/csrc
      p=ptrns   ! (defines p for first debug output only)
      do i=1,nseg
        if(kdwn(i).eq.2) then   ! if this is a turning ray
          rturn=rseg(i)+eps
          kturn=ktseg(i)
          call findiq(iqmin,1,rturn)    ! choose upper layer at discnty
          cturn=vps(iqmin,kturn)
          p=rturn/cturn
          pturn=max(p,pturn)
          if(jdebug.eq.1) write(13,*) 'cturn,p,pturn=',cturn,p,pturn
          sintrn=pturn*csrc/rsrc
          if(abs(sintrn).gt.1.0) stop 'Unexpected isin>1'
          aturn=asin(sintrn)  ! angle for turning wave
          if(kdwn(1).eq.0) then
            amin=max(amin+0.00001,aturn)        ! up: raise minimum angle
          else
            amax=min(amax-0.00001,pi-aturn)     ! down: lower max angle
          endif  
          if(jdebug.eq.1) write(13,*)'aturn,amin,amax=',aturn,amin,amax
        else if(kdwn(i).eq.3.or.kdwn(i).eq.4) then    ! reflection/transmission
          rtrns=rseg(i)-eps     ! ray needs to reach at least here
          ktrns=ktseg(i-1)      ! to reflect at or transmit this layer 
          call findiq(iqmin,1,rseg(i))  ! choose upper layer at discnty
          ctrns=vps(iqmin,ktrns)
          p=rtrns/ctrns
          ptrns=min(p,ptrns)
          ktrns=ktseg(i)
          ctrns=vps(iqmin,ktrns)
          p=rtrns/ctrns
          ptrns=min(p,ptrns)
          if(jdebug.eq.1) write(13,*) 'crefl/tr,p,ptrns=',ctrns,p,ptrns
          sintrs=ptrns*csrc/rsrc
          atrns=asin(sintrs)  ! angle for first refl/transmitted wave
          if(kdwn(1).eq.0) then
            amax=min(amax,atrns+0.00001)        ! up: lower max angle
          else
            amin=max(amin,pi-atrns-0.00001)     ! down : raise min angle
          endif  
          if(jdebug.eq.1) write(13,*)'atrns,amin,amax=',atrns,amin,amax
        endif
        if(jdebug.eq.1) write(13,*) i,kdwn(i),p,amin,amax
      enddo  

      ra1=amin
      ra2=amax
      iqmax=n
      if(kdwn(1).eq.1) iqmax=iqsrc

      ddtarget=0.2                      ! target delta spacing (deg)
      drtarget=ddtarget/r2d

      ystart(1)=rseg(1)         ! source radius
      ystart(3)=0.              ! source delta (phi in Dahlen 2000)
      ystart(4)=0.              ! source travel time
      ntbl=0


      ! find dpdd at start of table
      dra=0.0025*min(1.0,ra2-ra1)
10    ystart(2)=ra1
      if(jdebug.gt.0) write(13,*) 'start mkrtable a=',ra1*r2d,' to ',
     & ra2*r2d
      if(jdebug.gt.0) write(13,*) 'ntbl angle       delta    time dra'
      call tracer(ystart,rseg,ktseg,kdwn,nseg,y,jsg,rmin,tstar,h,nray)
      if(nray.le.0) then
        ra1=ra1+dra
        if(jdebug.gt.0) write(13,*) 'tracer returns nray=',nray,' i=',
     &     ystart(2)*r2d
        if(ra1.gt.ra2) then
          print *,'last raytable try for i=',ystart(2)*r2d
          print *,'which exceeds ra2=',ra2*r2d,' degrees'
          print *,'tracer returns nray=',nray
          print *,'Check ray definition or try running with jdebug =1'
          stop 'failure to make ray table'
        endif  
        goto 10
      endif  
      ntbl=ntbl+1
      ytbl(1,ntbl)=ystart(2)*r2d        ! store i
      ytbl(2,ntbl)=y(3,nray)*r2d        ! delta
      ytbl(3,ntbl)=y(4,nray)            ! time
      ytbl(4,ntbl)=rmin                 ! turning point radius
      if(jdebug.gt.0)
     &  write(13,*) ntbl,ystart(2)*r2d,y(3,nray)*r2d,y(4,nray),dra*r2d
      alast=ra1
      afirst=ra1
      dlast=y(3,nray)
      dadd=dra/drtarget         ! default in case first ddelta is zero

      do while(ystart(2).lt.ra2-0.0001)
        ystart(2)=min(ra2,ystart(2)+dra)
        if(ystart(2)-alast.lt.0.0001) goto 20   ! don't get stuck
        call tracer(ystart,rseg,ktseg,kdwn,nseg,y,jsg,rmin,tstar,h,nray)
        if(nray.le.0) then      ! still prblematic at low angles, restart
          ra1=ystart(2)+dra
          ntbl=0
          if(jdebug.gt.0) write(13,*) 'unexpected nray=',nray,' i=',
     &     ystart(2)*r2d
          if(ra1.lt.ra2) goto 10
          if(jdebug.gt.0) write(13,*) 'Unexpected stop of mkrtable'
          return
        endif
        ntbl=ntbl+1
        if(ntbl.gt.4000) stop 'ntbl>4000'
        ytbl(1,ntbl)=ystart(2)*r2d        ! store i
        ytbl(2,ntbl)=y(3,nray)*r2d        ! delta
        ytbl(3,ntbl)=y(4,nray)            ! time
        ytbl(4,ntbl)=rmin                 ! turning point radius
        if(jdebug.gt.0)
     &    write(13,*) ntbl,ystart(2)*r2d,y(3,nray)*r2d,y(4,nray),dra*r2d
        ! step between 0.057 and 0.57 degree (0.001-0.01 rad)
        dra=min(0.01,max(0.001,abs(drtarget*dadd)))
        alast=ystart(2)
        dlast=y(3,nray)
      enddo
20    a1=min(afirst,alast)*r2d
      a2=max(afirst,alast)*r2d
      if(jdebug.gt.0) write(13,*) 'Final table between i=',ra1*r2d,
     &   ' and ',alast*r2d

      return

      end



      subroutine derivpq(kupr,ktype,iq,y,dydx,p)

      ! computes derivatives (Dahlen et al., GJI 141:157-174, 2000)

      ! input:
      ! kupr - 1 to select c at upper layer of discontinuity, else 0
      ! ktype - 1 for P, 2 for S
      ! iq - guess for node number directly below rr
      ! y - 10 dimensional vector of P,Q,r and i
      ! p - slowness (s/rad)

      ! vector y is (P1,P2,Q1,Q2,~P1,~P2,~Q1,~Q2,r,i) in
      ! Dahlen's notation

      ! output:
      ! dydx - derivatives of y w.r.t. phi (Dahlen eqs A30, A31)

c     double precision y,dydx
      dimension y(10),dydx(10)

      jdebug=0

      rr=y(9)
      call velo1(rr,iq,kupr,ktype,c,dc,d2c)
      ri=1.0/rr
      r2=rr*rr
      tai=tan(y(10))
      cti=1.0/tai
      cti2=cti*cti
      if(jdebug.gt.0) write(13,*) 'rr,i,d,dc,d2c=',y(9),y(10),c,dc,d2c
c     if(jdebug.gt.0) write(13,*) 'tai,cti,ct2=',tai,cti,cti2
c     if(jdebug.gt.0) write(13,90) '    y=',(y(i),i=1,8)
      c3=c*c*c
      dydx(1)=-y(3)*p*(d2c+dc*cti2*ri)/c        ! (A30)
      dydx(2)=-y(4)*rr*dc/(c3*p)                 ! (A31)
      dydx(3)=y(1)*r2/p                         ! (A30)
      dydx(4)=y(2)*r2/p                         ! (A31)
      dydx(5)=-y(7)*p*(d2c+dc*cti2*ri)/c        ! (A30)
      dydx(6)=-y(8)*rr*dc/(c3*p)                 ! (A31)
      dydx(7)=y(5)*r2/p                         ! (A30)
      dydx(8)=y(6)*r2/p                         ! (A31)
      dydx(9)=rr*cti                             ! (A6)
      dydx(10)=rr*dc/c-1.0                       ! (A6)

      if(jdebug.gt.0) write(13,90) 'deriv:',(dydx(i),i=1,8)
90    format(a6,8e12.3)

      return
      end

      subroutine rk4pq(y,iq,h,ktype,jdown,yout,p)

      ! integrates y w.r.t. phi

      ! 4th order Runge Kutta solution to y'=f(x,y):
      ! y(n+1)=y(n)+(k1+2k2+2k3+k4)/6

      ! input:
      ! y(1-10) (see subroutine hessian ypq)
      ! iq=guess for node # directly below y(1)
      ! h: step size (rad)
      ! ktype=1 for P, 2 for S
      ! jdown=1 if step is decreasing r, 0 if upward ray

      ! output
      ! yout(1)-yout(10): next ray point
      ! iq is adjusted if needed

c     double precision dydx,y,yout,dym,dyt,yt
      dimension dydx(10),y(10),yout(10)
      dimension dym(10),dyt(10),yt(10)
      data halfpi/1.570796326795/

      jdebug=0

      hh=h*0.5
      h6=h/6.0
      n=10

      kupr=jdown      ! choose proper layer if y(1) at discontinuity
      call derivpq(1-kupr,ktype,iq,y,dydx,p)
      do  i=1,n
        yt(i)=y(i)+hh*dydx(i)  
      end do

      call derivpq(kupr,ktype,iq,yt,dyt,p)
      do i=1,n
        yt(i)=y(i)+hh*dyt(i)    
      enddo
      call derivpq(kupr,ktype,iq,yt,dym,p)

      do i=1,n
        yt(i)=y(i)+h*dym(i)   
        dym(i)=dyt(i)+dym(i)
      end do
      call derivpq(kupr,ktype,iq,yt,dyt,p)

      do  i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
      enddo

      if(jdebug.eq.1) write(13,90) yout
90    format('yout=',10e12.3)

      return
      end

      subroutine getcr2(flat,flon,ktype,vcr2,rcr2,ncr2,ecr2,rmoho,vsed)

      ! returns CRUST2.0 model in a format usable for subroutine inttau.

      ! But getcr2 ignores ice and water layer: it assumes you are 
      ! recordeing on an ocean bottom seismometer if oceanic.

      ! input: latitude and longitude flat, flon in degrees, ktype (1=P,2=S)
      ! output: arrays vcr2 (velocity) and rcr2 (radius) with ncr2 elements,
      !         and elevation ecr2 in km, rmoho the radius of the 
      !         model bottom, only in case this exceeds the crustal
      !         thickness of the reference model is this the CR2 moho...
      !         vsed is the velocity in the hard sedimentary layer.
      !         rcrust is the radius of upper crust surface
      ! ASSUMES r(n)=6371

      parameter(ityp=360,nla=90,nlo=180)                        ! for rcrust2 arrays

      character*2 types
      dimension vcr2(16),rcr2(16)

      common/crust2/amapvp(8,nlo,nla),amaprho(8,nlo,nla),
     &              amapvs(8,nlo,nla),amapthi(7,nlo,nla),
     &              amapele(nlo,nla)                    ! model CRUST2.0
      common/glbgrd/types(nlo,nla)                      ! model CRUST2.0
      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      jdebug=0
      if(jdebug>0) write(13,*) 'debug getcr2'

      cola=90.-flat
      cola=min(cola,179.9)              ! fixes S Pole bug
      if(flon.gt.180.)flon=flon-360.
      dx=360.0/nlo
      ilat=min(90,int(cola/dx)+1)
      ilon=min(180,int((flon+180.)/dx)+1)
      ecr2=amapele(ilon,ilat)/1000.             ! converts m to km
      vsed=amapvp(4,ilon,ilat)          ! 4=hard sediment layer in CRUST2.0
      if(ktype.gt.1) vsed=amapvs(4,ilon,ilat)
      a=6371.0+ecr2                             ! r of last solid layer
      ncr2=0
      rr=a
      if(jdebug.eq.1) write(13,*) 'getcr2',flat,flon,ilat,ilon,ktype
      do i=3,7                  ! ignore ice and ocean layer 1,2
        if(amapthi(i,ilon,ilat).gt.0.) then
          ncr2=ncr2+1
          vcr2(ncr2)=amapvp(i,ilon,ilat)
          if(ktype.gt.1) vcr2(ncr2)=amapvs(i,ilon,ilat)
          rcr2(ncr2)=rr
          if(jdebug.eq.1) write(13,*)i,ncr2,vcr2(ncr2),rcr2(ncr2)
          rr=rr-amapthi(i,ilon,ilat)
          ncr2=ncr2+1
          vcr2(ncr2)=amapvp(i,ilon,ilat)
          if(ktype.gt.1) vcr2(ncr2)=amapvs(i,ilon,ilat)
          rcr2(ncr2)=rr
          if(jdebug.eq.1) write(13,*) i,ncr2,vcr2(ncr2),rcr2(ncr2)
        endif
      enddo
      ! add mantle layer if reference Moho is deeper
      if(rr.gt.r(nmoh)) then
        ncr2=ncr2+1
        rcr2(ncr2)=rr
        vcr2(ncr2)=vp(nmoh-1)           ! Reference model Vp
        if(ktype.gt.1) vcr2(ncr2)=vs(nmoh-1)
        ncr2=ncr2+1
        vcr2(ncr2)=vcr2(ncr2-1)
        rcr2(ncr2)=r(nmoh)
        if(jdebug.eq.1) then
          write(13,*) nmoh,ncr2-1,vcr2(ncr2),rr
          write(13,*) nmoh-1,ncr2,vcr2(ncr2),rcr2(ncr2)
          write(13,*) 'CR2 Moho is at',rr,' versus reference',r(nmoh)
          write(13,*) 'Added:',rr-r(nmoh),' km with v=',vcr2(ncr2)
        endif
        rr=r(nmoh)
      endif  
      rmoho=rr
      if(jdebug.eq.1) write(13,*) 'Moho at r=',rmoho,' z=',6371-rmoho
      return
      end

      subroutine getref(ktype,vref,rref,nref,rmoho)

      ! returns reference model in a format usable for subroutine inttau.
      ! input: rmoho (from getcr2), the moho radius in km, ktype (1=P,2=S).
      ! output: arrays vref (velocity) and rref (radius) with nref elements.

      dimension vref(50),rref(50)

      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      jdebug=0
      if(jdebug>0) write(13,*) 'debug subroutine getref'

      nref=1
      vref(1)=vp(n)
      if(ktype.gt.1) vref(1)=vs(n)
      rref(1)=r(n)
      i=n
      if(jdebug>0) write(13,*) 'nref    i    rref    vref ktype=',ktype
      if(jdebug>0) write(13,"(2i5,f8.1,f8.3)") nref,i,rref(1),vref(1)

      ! map reference model in reversed order until Moho is reached
      do while(r(i).gt.rmoho)
        nref=nref+1
        i=i-1
        if(nref.gt.50) stop 'getref: nref>50, increase dimensions'
        if(i.lt.1) stop 'getref: i<1; crustal thickness error?'
        vref(nref)=vp(i)
        if(ktype.gt.1) vref(nref)=vs(i)
        rref(nref)=r(i)
        if(jdebug>0)write(13,"(2i5,f8.1,f8.3)") nref,i,
     &        rref(nref),vref(nref)
      enddo

      if(abs(rref(nref)-rmoho).lt.0.01) return

      ! interpolate last model node at rmoho
      dr=r(i+1)-r(i)
      if(dr.le.0.) stop 'getref: dr<0'
      vref(nref)=vp(i)+(rmoho-r(i))*(vp(i+1)-vp(i))/dr
      if(ktype.gt.1) vref(nref)=vs(i)+(rmoho-r(i))*(vs(i+1)-vs(i))/dr
      rref(nref)=rmoho
      if(jdebug>0) write(13,"(i5,5x,f8.1,f8.3)") nref,rmoho,vref(nref)

      return
      end

      subroutine inttau(p,r,rsrc,v,n,tau,tausrc)

      ! integrates for the intercept time tau throughout all layers
      ! of the model given in input model
      ! input: slowness p (s/rad), radii r (km), velocities (km/s) 
      ! in n nodes, starting at surface
      ! output: intercept time tau (s) over the whole model and tausrc 
      ! from rsrc to the surface
      ! tausrc is only correct if r is in descending order!

      dimension r(50),v(50),f(50)

      jdebug=0
      if(jdebug.eq.1) write(13,98) p
98    format('inttau slowness ',f8.1,/,'  i       r       v',
     & '    sini      dr    dtau     tau')
99    format(i3,f8.1,5f8.3)
97    format(a3,f8.1,5f8.3)

      ! integrate tau = integral cos(i)/v dr using trapezium rule
      tau=0
      tausrc=0                 ! default if no source
      dtau=0

      ! compute the integrand
      do i=1,n
        sini=p*v(i)/r(i)
        if(sini.ge.1.0) then
          f(i)=0.
        else  
          f(i)=sqrt(1.0-sini*sini)/v(i)
        endif  
      enddo  

      ! trapezoidal quadrature
      if(jdebug.eq.1) write(13,99) 1,r(1),v(1),p*v(1)/r(1),0.,0.,0.
      do i=2,n
        dr=r(i)-r(i-1)
        dr=abs(dr)                              ! r is descending
        dtau=0.5*dr*(f(i)+f(i-1))
        tau=tau+dtau
        if(jdebug.eq.1) write(13,99) i,r(i),v(i),p*v(i)/r(i),dr,dtau,tau
        ! the following assumes descending r - are we past the source?
        if(tausrc.eq.0.0.and.r(i).le.rsrc) then
          tau0=tau-dtau         ! tau at r(i-1)
          dr=r(i)-r(i-1)
          if(dr.ne.0.) tausrc=tau0+dtau*(rsrc-r(i-1))/dr
          if(jdebug.eq.1) write(13,97) 'src',rsrc,v(i),p*v(i)/r(i),
     &          rsrc-r(i),0.,tausrc
        endif  
      enddo
      if(rsrc.ge.r(1)) tausrc=0.   ! correct if source above bathymetrie or topo
      if(rsrc.le.r(n)) tausrc=tau
      if(jdebug.eq.1) write(13,*) 'Exit inttau tau,tausrc=',tau,tausrc

      return
      end

      subroutine rcrust2

c read infos to use to compute scripps crustal corrections
c required by vpcrust2
c 360 key profile: ctype(ityp),fvel(ityp,8),fdep(ityp,8)
c types(ilat,ilon)

c Reads the CRUST2.0 model files CNtype2_key.txt, 'CNtype2.txt and CNelevatio2.txt.
c It first looks for these files in $ directory specified by environment
c variable CRUST20, if it cannot find them it looks in the current directory

c fdep contains thickness of the crustal layers, fdep(ityp,8)
c is the total thickness of the crustal model along that profile

c CRUST2.0 has the following layers:
c 1 ice
c 2 water
c 3 soft sediments
c 4 hard sediments
c 5 upper crust
c 6 middle crust
c 7 lower crust
c 8 mantle

c For example, ocean (Vp,Vs,rho,thickness):
c A9      normal oceanic 9 km seds.                                       
c 3.81    1.5     2.3     4.5     5.0     6.6     7.1     8.15
c 1.94    0       1.2     2.5     2.5     3.65    3.9     4.65
c 0.92    1.02    2.2     2.6     2.6     2.9     3.05    3.35
c 0       5       1.5     7.5     1.7     2.3     2.5     inf.    20.5

c Continent:
c D1      Platform 11.5 km seds.
c 3.81    1.5     2.5     4.9     6.2     6.6     7.3     8.2
c 1.94    0       1.2     2.8     3.6     3.7     4       4.7
c 0.92    1.02    2.1     2.6     2.8     2.9     3.1     3.4
c 0       0       1.5     10.     10.5    11      8       inf.    41

c But on output ocean/water layers are switched - this is a peculiarity
c of Gabi's routine. These layers are ignored in getcr2, so it makes
c no difference to us.

      parameter(ityp=360,nla=90,nlo=180)

      dimension fvel(ityp,8),fvels(ityp,8),frho(ityp,8),fthi(ityp,8)
      character dum*1,line*506
      character dum0*5,fname*72
      character*2 ctype(ityp),types,atype(nlo)
      logical ruthere

      common/crust2/amapvp(8,nlo,nla),amaprho(8,nlo,nla),
     &              amapvs(8,nlo,nla),amapthi(7,nlo,nla),
     &              amapele(nlo,nla)
      common/glbgrd/types(nlo,nla)

c open files with Crust2.0 model

      fname='/Users/guust/data/models/3d/Crust2.0/CNtype2_key.txt'
      write(6,'(2a)') 'Looking for ',fname
      inquire(file=fname,exist=ruthere)
      if(ruthere) then
        print *,'Opening CRUST2. files in directory '
        print *,'/Users/guust/data/models/3d/Crust2.0/'
        open(72,file=fname,status='old')
        fname='/Users/guust/data/models/3d/Crust2.0/CNtype2.txt'
        open(73,file=fname,status='old')
        fname='/Users/guust/data/models/3d/Crust2.0/CNelevatio2.txt'
        open(74,file=fname,status='old')
      else  
        inquire(file='CNtype2_key.txt',exist=ruthere)
        if(.not.ruthere) then
          print *,'Cannot find Crust2.0 files. Either move them to '
          print *,'the current directory or change fname=...'
          print *,'to the directory where these files are.'
          stop 'Error locating Crust2.0 files such as CNtype2_key.txt'
        endif  
        print *,'Opening CRUST2. files in current directory'
        open(72,file='CNtype2_key.txt',status='old')
        open(73,file='CNtype2.txt',status='old')
        open(74,file='CNelevatio2.txt',status='old')
      endif  

      write(*,*)'reading crust2.0'

c ---- read key profiles --- file CNtype2_key.txt

      read(72,890)dum
      do 101 i=1,ityp
         read(72,899)ctype(i)
         read(72,899)line
         read(line,*)(fvel(i,l),l=1,8)
         read(72,899)line
         read(line,*)(fvels(i,l),l=1,8)
         read(72,899)line
         read(line,*)(frho(i,l),l=1,8)
         read(72,899)line
         read(line,*)(fthi(i,l),l=1,7)
c flip layers (NEEDS MODIFICATION IF STATION ON ICE LAYER)
         aux=fvel(i,1)
         fvel(i,1)=fvel(i,2)
         fvel(i,2)=aux
         aux=fvels(i,1)
         fvels(i,1)=fvels(i,2)
         fvels(i,2)=aux
         aux=frho(i,1)
         frho(i,1)=frho(i,2)
         frho(i,2)=aux
         aux=fthi(i,1)
         fthi(i,1)=fthi(i,2)
         fthi(i,2)=aux
 101  continue

c ---- read in key per cell of global grid (2x2 degree)
c ---- file CNtype2.txt

      read(73,899)line
      read(line,*)flons
      read(74,899)line
      do j=1,nla
        read(73,901)ilat,atype
        read(74,*) ilat,(amapele(i,j),i=1,nlo)
c       write(99,*)ilat,(amapele(i,j),i=1,nlo)
        do i=1,nlo
          do l=1,ityp
            if(atype(i).eq.ctype(l))then
              types(i,j)=ctype(l)
              do k=1,8
                amapvp(k,i,j)=fvel(l,k)
                amapvs(k,i,j)=fvels(l,k)
                amaprho(k,i,j)=frho(l,k)
              enddo
              do k=1,7
                amapthi(k,i,j)=fthi(l,k)
              enddo
              goto 10
            endif
         enddo
  10   enddo

      enddo

 890  format(////a)
 899  format(a)
 901  format(i4,1x,180(2x,a2,1x))

      write(*,*)'done read crust2.0'

      ! close files...
      close(72)
      close(73)
      close(74)

      return
      end

      subroutine getlegs(y,rayvel,nray,legend,nlegs)


      ! extracts 'legs' of the ray: the segments in between
      ! the end- or reflection points.
      ! transmissions should in principle also be treated
      ! as 'legs', but the error is hopefully small. The
      ! commented lines are ok if you wish to include them.

      ! input:
      ! y(4,i), i=1,...,nray: ray vector
      !         (y1=radius, y2=delta)
      ! rayvel(i), seismic velocity at node i, <0 if Vs

      ! output:
      ! legend(i), i=1,nlegs: end node of each leg
      ! nlegs: number of legs

      parameter (LRAY=20000)
      dimension y(4,LRAY),legend(190),rayvel(LRAY)

      data eps/0.01/    ! must be compatible with eps in rstp
      data halfpi/1.570796327/

      do i=1,190
        legend(i)=0
      enddo  
      nlegs=0
      do i=2,nray-1

c commented lines in case transmissions need to be handled explicitly
c       if(abs(y(1,i)-y(1,i-1)).lt.eps.and.     ! refl/transm
c    &     abs(y(2,i)-y(2,i-1)).lt.0.001.and.
c    &     abs(rayvel(i)-rayvel(i-1)).gt.0.001) then

        ! recognize reflection point
        if(abs(y(1,i)-y(1,i-1)).lt.eps.and.            ! r remains same
     &     (y(2,i)-halfpi)*(y(2,i-1)-halfpi).lt.-eps) then  ! angle flip
          nlegs=nlegs+1
          legend(nlegs)=i-1
          if(nlegs.ge.190) stop 'getlegs: nlegs>190'
        endif

      enddo

      nlegs=nlegs+1
      legend(nlegs)=nray

      return
      end

      subroutine gettstar(y,nray,jsg,ktseg,tstar,qray)

      parameter (LRAY=20000)
      dimension y(4,LRAY),jsg(LRAY),ktseg(20)
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      ! computes tstar = \int dt/Q = T/Qray from y,jsg and ktseg 
      ! where y is computed by subroutine tracer

      jdebug=0

      iq=n
      tstar=0.
      if(jdebug.gt.0) write(13,*) 'i    rr       qq     dt      tstar'
      kk=1
      do i=2,nray
        dt=y(4,i)-y(4,i-1)       ! travel time in step
        rr=0.5*(y(1,i)+y(1,i-1)) ! average r of step
        call findiq(iq,kk,rr)     ! find nearest node below rr
        t=rr-r(iq)               ! distance to nearest node
        ktype=ktseg(abs(jsg(i-1)))     ! P or S
        qq=Qs(iq)+t*(qqs(1,iq)+t*(qqs(2,iq)+t*qqs(3,iq)))
        if(ktype.eq.1) qq=0.75*(vp(iq)/max(0.01,vs(iq)))**2*qq
        if(vs(iq)*qq.ne.0.) tstar=tstar+dt/qq
        if(jdebug.gt.0) write(13,*) i,rr,qq,dt,tstar,t,iq
      enddo

      qray=1.0e20
      if(tstar.gt.0.) qray=y(4,nray)/tstar

      return
      end

      subroutine getelcor(scolat,slon,rcolat,rlon,y,nray,jsg,rseg,ktseg,
     &      kdwn,p,rayvel,rayq,ecorr)

      ! calculate ellipticity correction

      ! input:
      ! scolat,slon - source colatitude, longitude (rad, geocentric)
      ! rcolat,rlon - receiver coordinates (rad)
      ! y(k,i) - ray geometry from subroutine tracer: k=r,rayangle,phi,t
      !          where r in km, rayangle and phi in rad (t not used).
      ! nray - number of ray nodes i=1,nray
      ! jsg(i) - segment belonging to ray node i, negative if ray upward
      ! rseg(i) - start radius of segment i (km)
      ! ktseg(i) - wave type of segment i (1=P,2 or 3=S)
      ! kdwn(i) - ray direction (1=down,0=up,2=turning point,3=refl,
      !           4=transmission, 5=last point of ray)

      ! output:
      ! ecorr - time to be added to spherical earth travel time to
      !         obtain time in elliptical earth
      ! p - ray parameter in rad/s
      ! rayvel - ray velocity at node (>0 if Vp, negative if Vs)
      ! rayq - array with 1/Qs at ray node

      parameter (LRAY=20000)
      dimension y(4,LRAY),jsg(LRAY),rayvel(LRAY),rayq(LRAY)
      dimension rseg(20),ktseg(20),kdwn(20)
      dimension e2g(3,3)

      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/, third/0.333333333333333/

      ! for extensive debug set jdebug=2, just tracking=1, none=0
      jdebug=0

      ! compute transformation matrix from equatorial->geocentric
      call euler1(scolat,slon,rcolat,rlon,e2g)
      
      ! anchor ellipticity shift eldr to 0 at elliptical surface
      epsa=geteps(r(n))*r(n)
      
      ! get c and ellipticity at source radius
      jdown=kdwn(1)
      jghost=0
      if(kdwn(1).eq.0) jghost=1                 ! identifies pP,pS,sP or sS leg
      iq=n
      call velo1(y(1,1),iq,1-jdown,ktseg(1),c,dc,d2c)
      p=y(1,1)*sin(y(2,1))/c     ! slowness
      if(jdebug.gt.1) write(13,*) 'p=',p
      rayvel(1)=c
      if(ktseg(1).gt.1) rayvel(1)=-c
      eps=geteps(y(1,1))       ! ellipticity at r
      z=e2g(3,1)               ! cartesian z coordinate of source
      rr=y(1,1)
      deltar=eps*rr*(third-z*z)
      dphi=y(3,2)-y(3,1)
      sum=(rr/c)**3*(deltar/rr)*dc*dphi    ! postpone factor 0.5/p
      sumd=deltar*abs(cos(y(2,1)))/c

      if(jdebug.gt.1) then
      write(13,98)
98    format(5x,'i',3x,'jsg',4x,'kd',5x,'phi',9x,'r',8x,'dr',3x,'dtell',
     &       3x,'sumd',4x,'lat',6x,'y2',7x,'c   clast')
      write(13,99) 1,jsg(1),kdwn(1),y(3,1)*r2d,y(1,1),deltar,
     &       0.5*sum/p,sumd,90.-acos(z)*r2d
99    format(3i6,f8.2,2f10.2,2f7.3,f8.1,3f8.2)
      endif

      ! integrate with phi as independent variable, unless a
      ! discontinuity is crossed (when we apply boundary conditions)
      ! we use trapezoidal rule (but omit factor 0.5 while summing)
      do i=2,nray-1                       ! step from i-1 to i
        jdownold=jdown
        jdown=1
        if(jsg(i).lt.0) jdown=0
        dphi=y(3,i+1)-y(3,i-1)
        rr=y(1,i)
        clast=c
        iseg=abs(jsg(i))
        if(iseg.gt.2) jghost=0                  ! no ghost leg any more
        ! choose correct velocity when crossing discontinuity
        kupr=jdown
        dr=y(1,i)-y(1,i-1)
        if(abs(dr).lt.0.0001) kupr=1-jdown
        call velo1(rr,iq,kupr,ktseg(iseg),c,dc,d2c)
        rayvel(i)=c
        if(ktseg(iseg).gt.1) rayvel(i)=-c
        eps=geteps(y(1,i))
        xx=cos(y(3,i))                     ! ray at equator
        yy=sin(y(3,i))
        z=e2g(3,1)*xx+e2g(3,2)*yy           ! z geographic
        deltar=eps*rr*(third-z*z)
        sum=sum+(rr/c)**3*(deltar/rr)*dc*dphi 
        ! test for discontinuity
        if(abs(dr).lt.0.0001) then              ! if interface
          if(jghost.eq.1.and.rr.ge.rsrc) goto 10        ! skip ghost
          if(jdown.eq.jdownold) then     ! transmission
            sumd=sumd-deltar*(cos(y(2,i))/c-cos(y(2,i-1))/clast)
          else                           ! reflection 
            sumd=sumd-2.*deltar*cos(y(2,i))/c ! works for top&bottom
          endif  
10        continue
        endif
          
        if(jdebug.gt.1) then
          write(13,99) i,jsg(i),kdwn(iseg),y(3,i)*r2d,y(1,i),deltar,
     &       0.5*sum/p,sumd,90.-acos(z)*r2d,y(2,i)*r2d,c,clast
        endif
      enddo  

      ! treat last ray node
      jdown=0
      iq=n
      call velo1(y(1,nray),iq,jdown,ktseg(iseg),c,dc,d2c)
      rayvel(nray)=c
      if(ktseg(iseg).gt.1) rayvel(nray)=-c
      eps=geteps(y(1,nray))
      xx=cos(y(3,nray))                     ! ray at equator
      yy=sin(y(3,nray))
      z=e2g(3,1)*xx+e2g(3,2)*yy           ! z geographic
      rr=y(1,nray)
      deltar=eps*rr*(third-z*z)
      dphi=y(3,nray)-y(3,nray-1)
      sum=sum+(rr/c)**3*(deltar/rr)*dc*dphi    ! postpone factor 0.5/p
      sumd=sumd+deltar*abs(cos(y(2,nray)))/c   
      if(jdebug.gt.1) then
        write(13,99) i,jsg(i),kdwn(iseg),y(3,i)*r2d,y(1,i),deltar,
     &     0.5*sum/p,sumd,90.-acos(z)*r2d
      endif

      ecorr=0.5*sum/p+sumd
      if(jdebug.gt.0) write(13,*) 'ecorr=',ecorr

      return
      end

      subroutine euler1(ts,ps,tr,pr,e2g)

! calculates the Euler matrix e2g to go from
! equatorial coordinates (cartesian) along source-receiver path to 
! geocentric (again Cartesian). Use euler_rot if you also need the
! inverse, but e2g here is the same as in that routine.

! input: spherical coordinates (ts,ps): theta and phi for vector s and
!        (tr,pr) for vector r, both in radians
! output: matrices g2e and e2g

! (ts,ps) and (tr,pr) define a great circle n the sphere that acts
! as the equator in a new coordinate system, where (ts,ps) is the
! origin and the longitude increases in the direction of (tr,pr)

! Usage of g2e: compute (x,y,z) geocentric, g2e*(x,y,z) gives
! cartesian coordinates equatorial.  Eg, (ts,pr) -> (0,0) and
! (tr,ps) -> (0,Delta) when expressed back in angles.

! e2g maps new (cartesian) coordinate back to the old system.

! validated with tsteuler.f

      dimension e2g(3,3)

      call cart(ts,ps,xi,yi,zi)
      call cart(tr,pr,xr,yr,zr)

      xk=yi*zr-yr*zi
      yk=xr*zi-xi*zr
      zk=xi*yr-xr*yi

      dk=sqrt(xk*xk+yk*yk+zk*zk)
      xk=xk/dk
      yk=yk/dk
      zk=zk/dk

      xj=yk*zi-yi*zk
      yj=xi*zk-xk*zi
      zj=xk*yi-xi*yk

      dj=sqrt(xj*xj+yj*yj+zj*zj)
      xj=xj/dj
      yj=yj/dj
      zj=zj/dj

c rotation matrix cos(x',x)

      e2g(1,1)=xi
      e2g(2,1)=yi
      e2g(3,1)=zi
      e2g(1,2)=xj
      e2g(2,2)=yj
      e2g(3,2)=zj
      e2g(1,3)=xk
      e2g(2,3)=yk
      e2g(3,3)=zk

      return
      end

      subroutine cart(theta,phi,x,y,z)
      ! used in euler_rot and euler1
      ! colatatitude theta and longitude phi -> Cartesian x,y,z
      s=sin(theta)
      x=s*cos(phi)
      y=s*sin(phi)
      z=cos(theta)
      return
      end

      function geteps(x)

      ! returns the ellipticity of the Earth

      data ricb/1229.48/, rcmb/3484.3/

      if(x.gt.5971.) then
        geteps=0.003224+(x-5971.0)*3.2e-7
        return
      else if(x.gt.5701.) then
        geteps=0.003200+(x-5701.)*8.8889e-8
        return
      else if(x.gt.rcmb) then
        geteps=0.002666+(x-rcmb)*2.40899e-7
        return
      else if(x.gt.ricb) then
        geteps=0.002454+(x-ricb)*9.4021e-8
        return
      else
        geteps=0.002454
      endif

      return
      end

      subroutine mkrtable2(rseg,ktseg,kdwn,nseg,tableflag,
     &  hmf,hmb,ypq,rayvel,rayq)

      ! traces rays for a suite of slownesses and prints ray table
      ! alternative to mkrtable for the case that tableflag>0.

      ! input:
      ! ray geometry in rseg,ktseg,kdwn,nseg (see subroutine rdray)
      ! tableflag: 1 for standard table computation. >1 divides the
      !    standard initial angle spacing by the value of tableflag.
      !    <0 T-X curve computation with spacing scaled by abs(tableflag)
      ! model in common /modl/ (see subroutine model)

      ! output:
      ! r, i, delta, time, h11, h22, R
      ! where r=radius, i=angle with vertical (down is 0), delta=epicentral
      ! distance, h11,h22 are secon derivatives of time w.r.t. orthogonal
      ! ray coordinates q1 and q2 (in- and out of the ray plane, resp.).
      ! All angles in radians, r in km, T in sec.

      ! output to file raydyntrace.xxxx - see write(2,.) statements

      integer tableflag
      dimension rseg(20),ktseg(20),kdwn(20)
      parameter (LRAY=20000)
      dimension rayvel(LRAY),rayq(LRAY),plist(3000)
      dimension hmf(2,LRAY),hmb(2,LRAY),ypq(10,LRAY)
      dimension ystart(4),y(4,LRAY),jsg(LRAY)
      complex zz,zztot,zmas(0:3),rps(3,3),tps(3,3)              ! refl/transm coeff
      ! background model
      common /modl/r(999),vp(999),vs(999),rho(999),
     &      qvp(3,999),qvs(3,999),qrh(999),
     &      Qs(999),qqs(3,999),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/, eps/0.0001/, epsdr/0.1/
      data zmas/(1.,0.),(0.,-1.),(-1.,0.),(0.,1.)/

      ! set jdebug =1 for short debug, =2 for endless debug
      jdebug=0

      if(jdebug.gt.0) then
        write(13,*) 'mkrtable2 called with nseg:',nseg
        write(13,*) 'rseg=',(rseg(i),i=1,nseg)
        write(13,*) 'ktseg=',(ktseg(i),i=1,nseg)
        write(13,*) ' kdwn=',(kdwn(i),i=1,nseg)
        write(13,*) 'flag=',tableflag
      endif  

      ! careful... step size h must be same in mkrtable2 and main program
      h=10.    ! integration step size (keep equal to h in main program)


      ! initialize
      rsrc=rseg(1)
      ksrc=ktseg(1)
      iq=nmoh
      call findiq(iq,0,rsrc)
      iqmax=min(nmoh-1,iq)
      call velo1(rseg(1),iq,1-kdwn(1),ktseg(1),vsrc,dc,d2c)
      pmax=0.99999*rsrc/vsrc
      if(jdebug.gt.0) write(13,*) 'vsrc, pmax=',vsrc,pmax
      iqmin=iqmax
      pmin=0.
      rturn=0.
      rrefl=r(n)
      krefl=1

      ! find min and max slowness and the layers where they turn
      r1=rseg(1)
      k1=ktseg(1)
      kturn=k1
      if(jdebug>0) write(13,*) 'Finding min/max p:'
      do i=2,nseg
        if(kdwn(i).eq.2) then           ! if end of segment is a turning point
          krefl=0
          call findiq(iq,1,rseg(i))
          v=vp(iq)
          if(ktseg(i).eq.2) v=vs(iq)
          if(pmin.lt.rseg(i)/v) then
            pmin=rseg(i)/v
            iqmin=iq
            kturn=ktseg(i)
          endif  
        else if(kdwn(i).eq.3) then      ! if end of segment is reflection
          rrefl=min(rseg(i),rrefl)
          kup=1
          if(rseg(i).gt.r1) kup=0
          call findiq(iq,kup,rseg(i))
          v=vp(iq)
          if(k1.eq.2.or.ktseg(i).eq.2) v=vs(iq)
          if(pmax.gt.rseg(i)/v) then
            pmax=rseg(i)/v
            iqmax=iq
          endif  
        endif
        k1=ktseg(i)
        r1=rseg(i)
        if(jdebug>0) write(13,*) 'seg ',i,rseg(i),kdwn(i),pmin,pmax
      enddo  

      ! compose list of slownesses, depending on turning/reflection/up
      if(krefl.eq.1) then               ! if pmin reflects at deepest point
        np=200                          ! for PcP etc, 200 nodes is enough
        if(rrefl.ge.r(n)) np=30         ! case up (p or s) only
        dp=pmax/(np-1.0)
        if(jdebug.gt.0) write(13,*) 'krefl=1, np=',np,' dp=',dp
        do i=1,np
          plist(i)=(i-1)*dp
        enddo
        plist(np)=plist(np)-0.1*dp      ! avoids cmb problems I hope
        if(jdebug>0) write(13,*) 'list:',(plist(i),i=1,np)
      else                              ! turning rays
        np=0
        drold=99999.
        plold=99999.
        if(jdebug.gt.0) write(13,*) 'iqmax,iqmin=',iqmax,iqmin
        do iq=iqmax-1,iqmin,-1          ! loop over layers downwards
          dr=(r(iq+1)-r(iq))/(abs(tableflag)+3.0)    !at least 4 turning points/layer
          if(jdebug>0) write(13,*) 'loop iq=',iq,' dr=',dr
          if(dr.lt.0.01) then
            if(jdebug>0) write(13,*) 'disc:',r(iq+1),r(iq)
            ii=iq
            call velo1(r(iq),ii,1,kturn,c1,dc,d2c)      ! velocity top of interface
            call velo1(r(iq),ii,0,kturn,c0,dc,d2c)      ! and at bottom
            dc=(c1-c0)/(abs(tableflag)+3.0)
          endif
          do j=0,abs(tableflag)+2            ! subdivide layer at least 4*
            rr=max(r(iq),r(iq+1)-j*dr)  ! work your way down in r, start at iq+1
            ii=iq
            if(abs(dr).gt.0.01) then
              kup=1
              if(j.eq.0.and.drold.eq.0.) kup=0    ! if we start from discont
              call velo1(rr,ii,kup,kturn,c,dc,d2c)
              pl=rr/c
            else
              c=c1-j*dc                 
              pl=rr/c
            endif
            if(dc>0. .or. abs(pl-plold).lt.0.001) cycle
            if(pl.gt.pmax) cycle
            plold=pl
            np=np+1
            if(np.gt.3000) stop 'Increase dimension plist in mkrtable2'
            plist(np)=pl
            if(jdebug.gt.1) write(13,*) 'iq,rr,c,pl=',iq,rr,c,pl
          enddo
          drold=dr
        enddo
        if(jdebug.gt.0) write(13,*) 'krefl=0, np=',np
      endif  
            
      plast=1.0e10
      ntbl=0
      nskip=0

      ! Now loop over the slownesses in the list      

      if(jdebug.gt.0) write(13,*) 'ip   p'
      do ip=1,np

        p=plist(ip)
        if(jdebug.gt.0) write(13,*) 'ip=',ip,plist(ip)

        if(p.gt.plast.and.krefl.eq.0) then        ! shadow zone, new branch
          if(nskip.eq.0) then
            if(tableflag.eq.0) then
              write(8,fmt='(2i6,12f10.2,i5,a)') ntbl,0,p,
     &          0.,0.,0.,0.,0.,0.,(0.,0.),0.,0.,0.,0,' shadow zone'
            else
              write(8,fmt='(2i6,6f10.2,a)') ntbl,0,p,0.,0.,0.,0.,0.,
     &          '  shadow zone'         
            endif
          endif
          nskip=nskip+1
          cycle
        endif
        plast=p
        nskip=0

        ystart(1)=rseg(1)         ! source radius
        ystart(2)=asin(plist(ip)*vsrc/rseg(1))
        if(kdwn(1).eq.1) ystart(2)=pi-ystart(2)
        ystart(3)=0.              ! source delta (phi in Dahlen 2000)
        ystart(4)=0.              ! source travel time

        ! loop over decreasing turning depths (decreasing slowness p)
  
        call tracer(ystart,rseg,ktseg,kdwn,nseg,y,jsg,rmin,tstar,h,nray)
  
        if(nray.le.0) then
          if(jdebug.gt.0) then
            write(13,*) 'tracer returns nray=',nray,' i=',ystart(2)*r2d
            write(13,*) 'ystart=',ystart,' vsrc=',vsrc,' kdwn=',kdwn(1)
          endif
          cycle
        endif  

        if(tableflag<0) then
          write(12,*) y(3,nray)*r2d,y(4,nray)
          write(8,fmt='(2i6,5f10.2,f10.3,g12.3,5f8.4,i5)') ntbl,nray,p,
     &      rmin,y(2,1)*r2d,y(3,nray)*r2d,y(4,nray),tstar
          cycle
        endif  
  
        call hessian(y,ypq,nray,jsg,rseg,ktseg,kdwn,hmf,hmb,p,
     &        rayvel,rayq,Rxs,Mxs)
  
        ntbl=ntbl+1
        Rrs=sqrt(abs(ypq(3,nray)*ypq(4,nray)))/abs(rayvel(1))      ! (A55)
        Trs=y(4,nray)
        cs=abs(rayvel(1))
        Vr=abs(rayvel(nray))
        if(tableflag.gt.0)
     &    write(2,fmt='(2i6,2f10.4,e14.6,2f10.6,f10.3,a)')ntbl,nray,p,
     &    ystart(2)*r2d,Rrs,cs,Vr,Trs,
     &    '     ray nr, # of nodes, p, starting angle,Rrs,Vs,Vr,Trs'
        raylen=0.
        zz=cmplx(1.,0.)
        zztot=cmplx(1.,0.)
        if(tableflag.gt.0) write(2,40)
40      format(9x,'r',9x,'i     Delta',9x,'T',9x,'c',8x,'q0',11x,'h11',
     &    11x,'h22',11x,'Rxs',11x,'Rxr',9x,
     &    's     ReRT     ImRT    phase')

        if(tableflag.gt.0)
     &    write(2,20) (y(j,1),j=1,4),rayvel(1),rayq(1),hmf(1,1),
     &          hmf(2,1),0.,abs(cs*Rrs/rayvel(nray)),0.,zz,0.
        iii=2
        if(tableflag.lt.0) iii=nray
        y(1,nray+1)=0.          ! avoids error in if statement at the end
        y(3,nray+1)=0.          ! avoids error in if statement at the end
        jsg(nray+1)=jsg(nray)
        hmfold=hmf(1,1)
        do i=iii,nray
          c=abs(rayvel(i))
          q0=rayq(i)
          pkm=p/y(1,i)          ! slowness in s/km
          zz=(1.,0.)            ! refl/transm coefficient
          k=jsg(i)              ! k<0 if ray travels upward
          k1=abs(jsg(i+1))
          ktrans=1
          if(sign(1,jsg(i)).ne.sign(1,jsg(i+1))) ktrans=0  ! reflection
          kk=abs(k)
          ktype=ktseg(kk)
          ktype1=ktseg(k1)
          if(abs(y(1,i)-y(1,i+1)) < epsdr.and.abs(y(3,i)-y(3,i+1)) 
     &             < eps) then
            if(ktrans.eq.0) then   ! reflection
              rr=y(1,i)  
              if(k.gt.0) then              ! down->up
                call velo2(rr,iq,1,a1,b1,d1)
                call velo2(rr,iq,0,a2,b2,d2)
                call refl(rps,a1,a2,b1,b2,d1,d2,pkm)
                zz=rps(ktype,ktype1)
              else                         ! up > down  
                call velo2(rr,iq,0,a1,b1,d1)
                if(rr.lt.r(n)-0.001) then
                  call velo2(rr,iq,1,a2,b2,d2)
                else              ! surface
                  a2=0.
                  b2=0.
                  d2=0.
                endif  
                call refl(rps,a1,a2,b1,b2,d1,d2,pkm)
                zz=rps(ktype,ktype1)
              endif
            else                          ! transmission
              rr=y(1,i)  
              if(k.gt.0) then              ! down->down
                call velo2(rr,iq,1,a1,b1,d1)
                call velo2(rr,iq,0,a2,b2,d2)
                call trans(tps,a1,a2,b1,b2,d1,d2,pkm)
                zz=tps(ktype,ktype1)
              else                         ! up -> up
                call velo2(rr,iq,0,a1,b1,d1)
                call velo2(rr,iq,1,a2,b2,d2)
                call trans(tps,a1,a2,b1,b2,d1,d2,pkm)
                zz=tps(ktype,ktype1)
              endif
            endif  
          endif
          if(zz.eq.(0.,0.)) then           ! debug
            write(13,*) 'BUG? zz=0, at i=',i,' r=',y(1,i),' parameters:'
            write(13,*) a1,a2,b1,b2,d1,d2
            write(13,*) 'ktrans, k,ktype,ktype1=',ktrans,k,ktype,ktype1
          endif  
          zztot=zztot*zz
          if(i.eq.nray) then       ! at surface zz is surface Ux,Uz factor
            call surface(pkm,ktype,ux,uz,phase)
            zz=cmplx(ux,-uz)        ! only a trick to get it printed on unit2
          endif
  
          ! sum forward and backward hessian:
          h11=hmf(1,i)+hmb(1,i)
          h22=hmf(2,i)+hmb(2,i)
          if(jdebug>2) write(13,*) 'i,hmf=',i,hmf(1,i),hmf(2,i)
          hmfold=hmf(1,i)
          dr=y(1,i)-y(1,i-1)
          rddel=y(1,i)*(y(3,i)-y(3,i-1))
          ds=sqrt(dr*dr+rddel*rddel)    ! ray step length
          raylen=raylen+ds
          Rxs=sqrt(abs(ypq(3,i)*ypq(4,i)))/abs(rayvel(1))      ! (A56)
          Rxr=sqrt(abs(ypq(7,i)*ypq(8,i)))/abs(rayvel(nray))      ! (A57)
          ! write 
          ! y(1-4,i), i=1,...,nray: r,angle,delta,time of ray
          if(tableflag.gt.0)
     &     write(2,20) (y(j,i),j=1,4),rayvel(i),q0,h11,h22,Rxs,Rxr,
     &        raylen,zz,phase
  20      format(f10.3,2f10.5,f10.3,f10.5,f10.6,4e14.6,f10.1,2f9.3,f9.5)
        enddo
        Mxs=mod(Mxs,4)
        zztot=zztot*zmas(Mxs)
        write(8,fmt='(2i6,5f10.2,f10.3,g12.3,5f8.4,i5)') ntbl,nray,p,
     &    rmin,y(2,1)*r2d,y(3,nray)*r2d,y(4,nray),tstar,Rxs,zztot,ux,
     &    -uz,phase,Mxs
        if(jdebug.gt.0)
     &  write(13,fmt='(2i6,5f10.2,4g13.4)') ntbl,nray,p,rmin,y(2,1)*r2d,
     &        y(3,nray)*r2d,y(4,nray),tstar,Rxs,Rxr
  
        plast=plist(ip)

      enddo

      return

      end
