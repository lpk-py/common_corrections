subroutine refl(rps,a1,a2,b1,b2,d1,d2,p)

! computes reflection coefficients at an interface with input:
! P,S velocities a1,b1 in medium 1 and a2,b2 in medium 2 (km/s)
! for P and S waves traveling/reflecting in medium 1.
! Densities are respectively d1 and d2 (g/ccm), 
! the slowness is p s/km for both P and SV
! output: rps(1,1)= PP, rps(1,2)=PSv, rps(2,1)=SvP, rps(2,2)=SvSv
! and rps(3,3)=ShSh [others 0 in the absence of P or Sv/Sh conversion]

! if free surface set a2=b2=0, if fluid solid set the S velocity 0
! we compute fluid coefficients by a limiting process Vs->0 (Kennett)

! Using Aki and Richards displacement coefficients, chapter 5

complex rps(3,3),sini1,sini2,cosi1,cosi2,sinj1,sinj2,cosj1,cosj2
complex ccp,denom,E,F,G,H

jdebug=0

sini1=a1*p
sinj1=b1*p
cosi1=sqrt(1.0-sini1*sini1)
cosj1=sqrt(1.0-sinj1*sinj1)
do j=1,3
  do i=1,3
    rps(i,j)=0.
  enddo
enddo  
! avoid NaN
if(abs(a1-a2)<0.001.and.abs(b1-b2)<0.001.and.abs(d1-d2)<0.001) return
pp=p*p

if(a2.eq.0.) then               ! free surface
  bb=b1*b1
  b2p=(1.0/bb-2*pp)
  b2p2=b2p*b2p
  pab=4*pp/(a1*b1)
  ccp=cosi1*cosj1*pab
  denom=b2p2+ccp
  rps(1,1)=(-b2p2+ccp)/denom
  rps(1,2)=(4*p*cosi1/b1)*b2p/denom
  rps(2,1)=(4*p*cosj1/a1)*b2p/denom
  rps(2,2)=(b2p2-ccp)/denom
  rps(3,3)=1.0
else                            ! solid-solid or
  if(b2.eq.0.) b2=0.00001       ! solid-fluid or
  if(b1.eq.0.) b1=0.00001       ! fluid-solid
  a=d2*(1.0-2*b2*b2*pp)-d1*(1.0-2*b1*b1*pp)
  b=d2*(1.0-2*b2*b2*pp)+2*d1*b1*b1*pp
  c=d1*(1.0-2*b1*b1*pp)+2*d2*b2*b2*pp
  d=2*(d2*b2*b2-d1*b1*b1)
  sini2=a2*p
  sinj2=b2*p
  cosi2=sqrt(1.0-sini2*sini2)
  cosj2=sqrt(1.0-sinj2*sinj2)
  E=b*cosi1/a1+c*cosi2/a2
  F=b*cosj1/b1+c*cosj2/b2
  G=a-d*(cosi1/a1)*(cosj2/b2)
  H=a-d*(cosi2/a2)*(cosj1/b1)
  denom=E*F+G*H*pp
  rps(1,1)=((b*cosi1/a1-c*cosi2/a2)*F-(a+d*(cosi1/a1)*(cosj2/b2))*H*pp)/denom
  rps(1,2)=-2*(cosi1/a1)*(a*b+c*d*(cosi2/a2)*(cosj2/b2))*p*a1/(b1*denom)
  rps(2,1)=-2*(cosj1/b1)*(a*b+c*d*(cosi2/a2)*(cosj2/b2))*p*b1/(a1*denom)
  rps(2,2)=-((b*cosj1/b1-c*cosj2/b2)*E-(a+d*(cosi2/a2)*(cosj1/b1))*G*pp)/denom
  rps(3,3)=(d1*b1*cosj1-d2*b2*cosj2)/(d1*b1*cosj1+d2*b2*cosj2)
endif

!debug
if(jdebug.gt.0.) then
  write(13,*) 'Debug - input subr refl:',a1,a2,b1,b2,d1,d2,p
  write(13,99) rps(1,1),rps(1,2),rps(2,1),rps(2,2)
  99 format("R=",2f10.5,3x,2f10.5,/,2x,2f10.5,3x,2f10.5)
endif  

return
end

subroutine trans(tps,a1,a2,b1,b2,d1,d2,p)

! computes transmission coefficients at an interface with input:
! P,S velocities a1,b1 in medium 1 and a2,b2 in medium 2 (km/s)
! for P and S waves traveling/reflecting in medium 1.
! Densities are respectively d1 and d2 (g/ccm), 
! the slowness is p s/km for both P and SV
! output: tps(1,1)= PP, tps(1,2)=PSv, tps(2,1)=SvP, tps(2,2)=SvSv
! and tps(3,3)=ShSh [others 0 in the absence of P or Sv/Sh conversion]

! if fluid-solid set the S velocity 0
! we compute fluid coefficients by a limiting process Vs->0 (Kennett)

! Using Aki and Richards displacement coefficients, chapter 5

complex tps(3,3),sini1,sini2,cosi1,cosi2,sinj1,sinj2,cosj1,cosj2
complex ccp,denom,E,F,G,H

jdebug=0

sini1=a1*p
sinj1=b1*p
cosi1=sqrt(1.0-sini1*sini1)
cosj1=sqrt(1.0-sinj1*sinj1)
do j=1,3
  do i=1,3
    tps(i,j)=0.
  enddo
  tps(j,j)=1.0
enddo  
if(abs(a1-a2)<0.001.and.abs(b1-b2)<0.001.and.abs(d1-d2)<0.001) return
pp=p*p

if(b2.eq.0.) b2=0.00001       ! solid-fluid or
if(b1.eq.0.) b1=0.00001       ! fluid-solid
a=d2*(1.0-2*b2*b2*pp)-d1*(1.0-2*b1*b1*pp)
b=d2*(1.0-2*b2*b2*pp)+2*d1*b1*b1*pp
c=d1*(1.0-2*b1*b1*pp)+2*d2*b2*b2*pp
d=2*(d2*b2*b2-d1*b1*b1)
sini2=a2*p
sinj2=b2*p
cosi2=sqrt(1.0-sini2*sini2)
cosj2=sqrt(1.0-sinj2*sinj2)
E=b*cosi1/a1+c*cosi2/a2
F=b*cosj1/b1+c*cosj2/b2
G=a-d*(cosi1/a1)*(cosj2/b2)
H=a-d*(cosi2/a2)*(cosj1/b1)
denom=E*F+G*H*pp
tps(1,1)=2*d1*(cosi1/a1)*F*a1/(a2*denom)
tps(1,2)=2*d1*(cosi1/a1)*H*p*a1/(b2*denom)
tps(2,1)=-2*d1*(cosj1/b1)*G*p*b1/(a2*denom)
tps(2,2)=2*d1*(cosj1/b1)*E*b1/(b2*denom)
tps(3,3)=2*d1*b1*cosj1/(d1*b1*cosj1+d2*b2*cosj2)

!debug
if(jdebug.gt.0.or.abs(tps(1,1)).eq.0.) then
  write(13,*) 'BUG? input subr trans:',a1,a2,b1,b2,d1,d2,p
  write(13,99) tps(1,1),tps(1,2),tps(2,1),tps(2,2)
  99 format("T=",2f10.5,3x,2f10.5,/,2x,2f10.5,3x,2f10.5)
endif  

return
end
