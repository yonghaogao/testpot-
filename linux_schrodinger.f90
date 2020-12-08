!  Parameters :
!  xmin(xmax):the lower(upper) limit of integral
!  tol:  error range
   module para
      implicit none
      real*8::xmin,xmax,Tol
      real*8::hbarc,rm
      integer::n
      end module
!-----------------------------------------
!  xm: match point
!  e:  Experimental energy value
!  em: Energy eigenvalue
!  rm=mn*mp/(mn+mp)
!    = 469.5662127985447 Mev
      program schrodinger
      use para
      implicit none
      integer::i
      real*8,dimension(1:2000)::V,k2
      real*8,dimension(1:2000)::x,y1,y2,y
      real*8::xm,h,c
      real*8::templ,tempr
      real*8::E
      n=500
      xmin=0.
      xmax=10.
      tol=1.E-5
      e=-2.224d0
      hbarc=19.7327d0
      rm=469.566213d0
      xm=5.d0
      delta=1.
      call left(e,xm,delta,xl,yl)
      call right(e,xm,delta,xr,yr)
      do i=1,n
!         y(i)=yl(i)
      end do
      do i=n+1,2*n
         y(i)=yr(i-n)
      end do
      
      do while()
        call equel(e,xm,x,y1)
          ddelta=
          delta=delta+ddelta
          call equel(e,xm,x,y2)
          templ=(y1(n+1)-y1(n))/(x(n+1)-x(n))-(y1(n+3)-y1(n+2))/(x(n+3)-x(n+2))
          tempr=(y2(n+1)-y2(n))/(x(n+1)-x(n))-(y2(n+3)-y2(n+2))/(x(n+3)-x(n+2))
       ! the normalization
       c=0
       do i=1,2*n+1
          c=c+(x(i+1)-x(i))/2*(y2(i)**2+y2(i+1)**2)
       end do
       do i=1,2*n+1
          y(i)=y2(i)/sqrt(c)
          write(*,*) y(i)
       end do
      end program
   !------------------------------------------------------------------------
   !
   !
      subroutine equel(e,xm,x,yyl)
      use para
      implicit none
      integer::i
      real*8,dimension(1:2000)::V,k2
      real*8,dimension(1:2000)::x,yl,yr,yyl
      real*8::h
      real*8::delta,ddelta,xm
      real*8::E
      call left(e,xm,delta,x,yl)
      call right(e,xm,delta,x,yr)
      do while(abs(ddelta)>tol)
         delta=delta-ddelta
         call left(e,xm,delta,x,yyl)
         if((yyl(n+1)-yr(1))*(yl(n+1)-yr(1))>0) then
         else
         delta=delta+ddelta
         end if
      end do
      end subroutine
   
   
   !------------------------------------------------------------------------
   !
   !   Numerov algorithm
   !   LEFT
   !
   !   INPUT:
   !   e: Test eigenvalues
   !   xm: match point
   !   del:adjustable parameter
   !   x(i):ridus
   !   k2:  k2=k**2
   
   !   OUTPUT:
   !   y(i): wave function
   
       subroutine left(e,xm,delta,x,y)
       use para
       implicit none
       integer::i
       real*8::delta,h,xm
       real*8,dimension(1:2000)::x,y
       real*8,dimension(1:2000)::V,k2
       real*8::gausspot
       real*8::E
       h=(xm-xmin)/n
       y(1)=0
       y(2)=h
       do i=1,n
           x(i)=xmin+(i-1)*h
           V(i)=gausspot(x(i))
           k2(i)=2*rm*(e-v(i))/hbarc**2
       end do
       do i=2,n
           y(i+1)=2*(1-5*h**2*k2(i)/12.*(e-v(i)))*y(i)-(1+h**2*k2(i-1)/12.*(e-v(i-1)))*y(i-1)
           y(i+1)=y(i+1)/(1+h**2*k2(i+1)/12.*(e-v(i+1)))
       end do
       end subroutine
       
   !--------------------------------------------------------------------------
   !   RIGHT	
   !
   !   INPUT:
   !   e: Test eigenvalues
   !   xm: match point
   !   delta: adjustable parameter
   !   mp: mass of proton
   !   mn: mass of neutron
   !   rm: reduced mass
   !   HETA : Sommerfeld parameter
   !   hbar: h*c/(2*pi)
   !   ec: e=-1.60218E-19
   !   OUTPUT:
   !   x(i):ridus
   !   y(i): wave function
       subroutine right(e,xm,delta,x,y)
       use para
       implicit none
       integer::i
       integer::zero
       real*8::delta,h,xm
       real*8::dy1,dy2
       real*8::HETA
       real*8,dimension(1:2000)::x,y,dyl,dyr
       real*8,dimension(1:2000)::V,k2
       real*8::gausspot
       real*8::temp,kn,knn
       real*8::E
       heta=0.
       zero=0
       h=(xmax-xm)/n
       temp=xmax-h
       !heta=z1*z2*e2*niu*/(hbar*sqrt(w*niu*Energy))
       do i=1,n
           x(i)=xm+(i-1)*h
           V(i)=gausspot(x(i))
           k2(i)=2*rm*(e-v(i))/hbarc**2
       end do
       kn=sqrt(abs(k2(n)))
       knn=sqrt(abs(k2(n-1)))
       call WHIT(HETA,xmax,kn,e,zero,y(n+1),dyl(n+1),zero)
       call WHIT(HETA,temp,knn,e,zero,y(n),dyr(n),zero)
       do i=n,2,-1
           y(i-1)=2*(1-5*h**2*k2(i)/12.*(e-V(i)))*y(i)-(1+h**2*k2(i+1)/12.*(e-v(i+1)))*y(i+1)
           y(i-1)=y(i-1)/(1+h**2*k2(i-1)/12.*(e-v(i-1)))
       end do
       write(*,*) y(1)
       end subroutine
   !-------------------------------------------------------------------------
   !     c *** Gaussian Potential
       function gausspot(r)
         implicit none
          real*8 r,v0,r0,gausspot,a
          v0=72.15d0
          r0=0.
          a=1.484d0
            if (a.gt.1e-6) then
              gausspot=-V0*exp(-(r-r0)**2/a**2)
                else
                  write(*,*)'a too small in gausspot!'
                  stop
            endif
            return
       end function
         
   !-------------------------------------------------------------------------- 
     
         SUBROUTINE WHIT(HETA,R,XK,E,LL,F,FD,IE)
   !
   !     CALCULATES  WHITTAKER  FUNCTION  WL(K,R)  WITH
   !     ASYMPTOTIC  FORM  EXP(-(KR + ETA(LOG(2KR)))
   !     E  IS  NEGATIVE
   !     If IE = 0, allowed to return result e**IE larger than Whittaker,
   !                for the IE value returned.
   !     If IE > 0, must scale results by that amount.
   !
   !   input : 
   !           HETA : Sommerfeld parameter
   !           R : radius 
   !           XK: module of wavenumber in fm^{-1}
   !           E :  C.M. energy in MeV 
   !           LL :  partial wave 
   !           IE :  normally set to 0 
   !   output
   !           F(LL+1) : WHITTAKER  FUNCTION
   !           FD(LL+1) : derivative WHITTAKER  FUNCTION
   
   !	use drier !  AMoro
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION F(LL+1),FD(LL+1) ,T(12),S(7)
   !! AMoro: to replace drier module
         REAL*8 FPMAX
   !	acc8 = epsilon(acc8);
         fpmax = huge(acc8)**0.8d0 
   !! -------------------------
   ! added by yhgao
         IE=0
         L = LL+1
   !              NOW L = NO. OF VALUES TO FIND
         EE=-1.0
         AK=XK
         ETA=HETA
         LP1=L+1
         RHO=AK*R
      S(:) = 0
         IF(L-50)1,1,2
       1 LM=60
         GO TO 3
       2 LM=L+10
       3 LMP1=LM+1
         IS=7
         PJE=30.0*RHO+1.0
         H=max(INT(PJE),4)
         H=RHO/H
         RHOA=10.0*(ETA+1.0)
         IF(RHOA-RHO)13,13,14
      13 IFEQL=1
         RHOA=RHO
         GO TO 15
      14 IFEQL=0
      15 PJE=RHOA/H+0.5
         RHOA=H*INT(PJE)
         IF(IFEQL)16,16,18
      16 IF(RHOA-RHO-1.5*H)17,18,18
      17 RHOA=RHO+2.0*H
      18 IF(EE)55,55,19
      19 STOP 'WHIT'
      27 A=2.0-10.0/12.0*H*H*EE
         B=1.0/6.0*H*ETA
         C=1.0+1.0/12.0*H*H*EE
         M1=INT(RHOA/H-0.5)
         M2=INT(RHO/H-1.5)
         T(2)=B/FLOAT(M1+1)
         T(3)=B/FLOAT(M1)
         JS=M1
         DO 29 IS=M2,M1
         DO 28 I=1,6
         S(I)=S(I+1)
      28 CONTINUE
         T(1)=T(2)
         T(2)=T(3)
         T(3)=B/FLOAT(JS-1)
         S(7)=((A+10.0*T(2))*S(6)-(C-T(1))*S(5))/(C-T(3))
         JS=JS-1
         IF(ABS(S(7)).LE.FPMAX) GO TO 29
          DO 285 I=2,7
     285   S(I) = S(I) / FPMAX
      29 CONTINUE
         T(1)=S(4)
         T(2)=(1.0/60.0*(S(1)-S(7))+0.15*(S(6)-S(2))+0.75*(S(3)-S(5)))/H
         GO TO 60
      55 C=1.0/RHOA
         A=1.0
         B=1.0-C*ETA
         F(1)=A
         FD(1)=B
         DO 56 M=1,26
         D=0.5*(ETA+FLOAT(M-1))*(ETA+FLOAT(M))*C/FLOAT(M)
         A=-A*D
         B=-B*D-A*C
         F(1)=F(1)+A
         FD(1)=FD(1)+B
      56 CONTINUE
         A=-ETA*LOG(2.0*RHOA)-RHOA
         FPMINL = -LOG(FPMAX)
         if(IE.eq.0.and.A.LT.FPMINL) IE = INT(FPMINL-A)
         A=EXP(A+IE)
         F(1)=A*F(1)
   !      FD(1)=A*FD(1)
         FD(1)=A*FD(1) * (-1d0 - 2*ETA/(RHOA))
         IF(IFEQL)57,57,61
      57 S(IS)=F(1)
         IF(IS-7)27,58,27
      58 IS=6
         RHOA=RHOA+H
         GO TO 55
      60 F(1)=T(1)
         FD(1)=T(2)
      61 C=1.0/RHO
         DO 63 M=1,L-1
         A=ETA/FLOAT(M)
         B=A+C*FLOAT(M)
         F(M+1)=(B*F(M)-FD(M))/(A+1.0)
         FD(M+1)=(A-1.0)*F(M)-B*F(M+1)
      63 CONTINUE
         DO 65 M=1,L
         FD(M)=AK*FD(M)
      65 CONTINUE
         RETURN
         END SUBROUTINE
         
   