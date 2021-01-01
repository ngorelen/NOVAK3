c Here we evaluate the dispersion relation for stability analysis of 
c M=1 internal kink mode see F.Cheng's paper in Phys.Reports.92, Eq.(6.17).
C
C     In this program we compute w2_t and w2_s - squars of w_t and w_s
C     as defined in that paper. The frequencies are normalized, as all
C     frequencies in NOVA, by (v_A/ q1 R0), where q1 is edge q value, R0
C     is the major radius of the geometric center of last surface -
C     Cheng's choice.
C
C     On input:
C     omreal = sqrt(|omega_mhd^2|)(sign(omega_mhd^2) where omega_mhd^2
C     is ideal MHD frequency (negative for unstable ideal modes)
C     omkom  is deltaWk / deltaK in the dispersion relation. 
C     omkomflr is the same but with the flr effects.
C     
C     On output it computes w2_t,w2_s and spits out on screen (pipe)
C     values of growth rates with fast ions as a function of fast ion
C     betas as follows from Cheng dispersion relation. For a more
C     complete dispersion relation see file NOVAK_m1dispersion.pdf in
C     the source directory. 
C
C     This program also gives values of deltaWk/omegaA
C     (NOVA normalization) for comparison with the values from
C     Mcclements and Porcelli papers (see below). This is usefull for
C     self benchmark for a given equilibrium.
C     
      subroutine w_tw_s(b0,deltk,ihsps,omreal,omkom,omkomflr,w2_t,w2_s
     &     ,rdw_k,rdw_kflr)
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      include 'orb2d.par'
      common/functn/omreal1,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      dimension ind(nn),xrt(nn),frt(nn)
      complex gamma_k,omega(2)
      q_m1=1.
      gam_s=5./3.
      pcrat=0.
      w2_t=0.
      do itheta=1,ict
         pcrat=pcrat+pc(1,itheta)
      enddo
      do itheta=1,iht
         w2_t=w2_t+ph(ihsps,itheta)
      enddo
      pcrat=pcrat/w2_t
      pcrat=pcrat/(1.+pcrat)
c this fix is for the case of nonmonotonic q-profiles and n=m=1 modes
      if((q1-q_m1)*(qoo(1)-q_m1).gt.0.)then
         q_m1=100.
         do i=1,nosurf
            q_m1=amin1(q_m1,qoo(i))
            print *,i,q_m1,qoo(i)
         enddo
         q_m1=q_m1*1.01
      endif
      call roots(qoo(1),ppsi(1),ind(1),nosurf,q_m1,1.e-6,0.,xrt(1),nroot
     &     )
      w2_t=0.
      w2_s=0.
      w2_a=0.
      if((q1-q_m1)*(qoo(1)-q_m1).gt.0.)then
      write(*,*) '------------------------'
         write(*,*) 
     &   'Warning: your surface is out of the safety factor limits. ',
     &   'Will use q_m1',q_m1,'and r_root',xrt(nroot),nroot
c         return
      endif
      pcrat=gam_s*ynn(ppsi(ind(nroot)),pprf(ind(nroot)),xrt(nroot))*
     &     (0.5*(xjacob2d(1,1)+xjacob2d(1,mth/2))*q1/qoo(1))**2
     &     /ynn(ppsi(ind(nroot)),rhoprf(ind(nroot)),xrt(nroot))*pcrat
      do itheta=1,mth
         hhh=ynn(ppsi(ind(nroot)),b2d(ind(nroot),itheta),xrt(nroot))**2
         w2_t=w2_t+hhh*pcrat*ynn(ppsi(ind(nroot)),bnom(ind(nroot),itheta
     &        ),xrt(nroot))**2/ynn(ppsi(ind(nroot)),grpssq2d(ind(nroot)
     &        ,itheta),xrt(nroot))
         w2_s=w2_s+pcrat/hhh/ynn(ppsi(ind(nroot)),xjacob2d(ind(nroot)
     &        ,itheta),xrt(nroot))**2
         w2_a=w2_a+1./ynn(ppsi(ind(nroot)),rhoprf(ind(nroot)),xrt(nroot)
     &        )/ynn(ppsi(ind(nroot)),xjacob2d(ind(nroot),itheta)
     &        ,xrt(nroot))**2
      enddo
      w2_t=w2_t/real(mth)
      w2_s=w2_s/real(mth)
c
      w2_a=q1*rax*sqrt(w2_a/real(mth))
      call deriv4(ppsi(ind(nroot)-1),qoo(ind(nroot)-1),xrt(nroot),dyn)
      call fun4(ppsi(ind(nroot)-1),rr(ind(nroot)-1,1),xrt(nroot),yn)
      call fun4(ppsi(ind(nroot)-1),rr(ind(nroot)-1,mth2/2),xrt(nroot)
     &     ,yn1)
      write(*,*) ' R_q=1=',yn
      r0_q1=(yn+yn1)*0.5
      r_q1=(yn-yn1)*0.5
c      write(*,*) ' r_q=1=',r_q1,' xi^2=',(ynn(ppsi(ind(nroot))
c     &     ,eigfun(ind(nroot),2-minm,1),xrt(nroot))**2
c     &     /ynn(ppsi(ind(nroot)),grpssq(ind(nroot)),xrt(nroot)))
c     &     ,' grpssq=',ynn(ppsi(ind(nroot)),grpssq(ind(nroot)),xrt(nroot
c     &     )),' xi_psi=',ynn(ppsi(ind(nroot)),eigfun(ind(nroot),2-minm,1
c     &     ),xrt(nroot))
      yn=((yn/2.-yn1/2.)/rax)**2
c      rdw_k=omkom*abs(omreal)
      rdw_k=omkomflr*abs(omreal)
ckg that was introduced for testing purposes
c*4./q1/2./xrt(nroot)/dyn*yn
c      rdw_kflr=omkomflr*abs(omreal)*4./q1/2./xrt(nroot)/dyn*yn
c      rdw_kflr=deltk*1.6726/rax*20./1.e6/b0**2/3.1415926**2/(2.
c     &     *xrt(nroot)*dyn)/(ynn(ppsi(ind(nroot)),eigfun(ind(nroot),2
c     &     -minm,1),xrt(nroot))**2/ynn(ppsi(ind(nroot)),grpssq(ind(nroot
c     &     )),xrt(nroot)))*w2_a
c      rdw_k=omkom*rdw_kflr
      write(*,*) ' q`r=',2.*sqrt(xrt(nroot))*dyn,'bet_norm',betah0(ihsps
     &     )
      write(*,*) ' r*q`r=',2.*xrt(nroot)*dyn
      err=1.e-5
      hb=.001
      beta=0.
      gamma_k=-cmplx(rdw_k,0.)/betah0(ihsps)
      write(*,*)'Cheng m=1 dispersion predicts stabilization dependence'
 1    call frequency(abs(omreal),gamma_k,sqrt(w2_s),sqrt(w2_t),beta
     &     ,omega,nomega)
      write(*,*)'beta,omega_im',beta,aimag(omega(2))
      beta=beta+hb
      if(beta.le..02+err) goto 1
c
ckg Evaluate Lenin criterium on q=1 surface 
ckg   local value of grad psi * R at q=1 surface
      transf=sqrt(0.5*(abs(ynn(ppsi(ind(nroot)),grpssq2d(ind(nroot)
     &     ,1),xrt(nroot)))+abs(ynn(ppsi(ind(nroot)),grpssq2d(ind(nroot)
     &     ,mth/2.),xrt(nroot)))))/psitot*rax
      r1qp1=sqrt(yn)*transf*dyn
ckg   Lenin critical value is
ckg    calculate first effective mass of the plasma 
      r1qp1cr=rho0/(denc0(2)+denc0(3)+denc0(4))
c      write(*,*) 'effective mass is',r1qp1cr,ynn(ppsi(ind(nroot))
c     &     ,zeff(ind(nroot)),xrt(nroot))
      call funder4(ppsi(ind(nroot)),pprf(ind(nroot)),xrt(nroot),pp1,p1)
      call funder4(ppsi(ind(nroot)),denc(ind(nroot),1),xrt(nroot),rhop1
     &     ,rho1)
      write(*,*) 'pp1,p1,rhop1,rho1= ',pp1,p1,rhop1,rho1
      r1qp1cr=1.4*(r1qp1cr*0.5/ynn(ppsi(ind(nroot)),zeff(ind(nroot))
     &     ,xrt(nroot)))**(1./6.)*(beta0*p1*abs(rhop1)*transf/(rho1
     &     *pprf(1)))**(2./3.)*(abs(pp1)*transf/p1)**(1./3.)
      write(*,*) 'd q/ d psi_norm=',dyn
      write(*,*) 'd q/ d r_norm_polflux=',2.*sqrt(xrt(nroot))*dyn
      write(*,*) 'r d q/ d r=',r1qp1
      write(*,*)'According to Lenin criterion m=1 is stable if crit'
     &     ,r1qp1cr,' > local r1qp1',r1qp1
c----------- compare with analytical stuff -------------
      write(*,*) 'oms/om to comp. with anal. model (+FLR)=gam/om_ANOVA='
c     &     ,omkom*abs(omreal)/q1,omkomflr*abs(omreal)/q1
     &     ,omkom*abs(omreal),omkomflr*abs(omreal)
c
      do ir=1,ind(nroot)
         xrt(ir)=(abs(rr(ir,mth2/2)-rr(ir,1))*0.5/r_q1)**2
      enddo
      call ASIMP(xrt(1),pprf(1),ind(nroot),res)
      beta_p1=(res/(xrt(ind(nroot))-xrt(1))-pprf(ind(nroot
     &     )))/pprf(1)*beta0/(r_q1/r0_q1)**2
      beta_pc=0.3*(1.-1.6666666*r_q1/(abs(rr(nn
     &     ,mth2/2)-rr(nn,1))*0.5))
      write(*,*) 'beta_p1=',beta_p1,'beta_pc=',beta_pc,'dW_buss=',-3.*3
     &     .1415926*(r_q1/r0_q1)**2*(beta_p1**2-beta_pc**2)
c-----------employ Mcclements formula from NF,v35,p.1761(1995).-----
ckg note that php is array of the fast ion pressure derivative over 
ckg the poloidal flux, which goes from 0 to psitot at the edge and
ckg ppsi is the array normalized and goes from 0 to 1
      sa=abs(rr(nn,mth2/2)-rr(nn,1))*0.5
      ir=1
      xrt(ir)=(abs(rr(ir,mth2/2)-rr(ir,1))*0.5/r_q1)
      do ir=1,ind(nroot)
         xrt(ir+1)=(abs(rr(ir+1,mth2/2)-rr(ir+1,1))*0.5/r_q1)
c
         frt(ir)=xrt(ir)**1.5*php(1,ihsps)*(0.6+3.2*(1-qoo(ir)-0.5
     &        *(xrt(ir+1)+xrt(ir))*(qoo(ir+1)-qoo(ir))/(qoo(ir+1)+qoo(ir
     &        ))/ (xrt(ir+1)-xrt(ir)) ))
      enddo
      call ASIMP(ppsi(1),frt(1),ind(nroot),res)
      res=res*sqrt(3.)*3.1415926/8.*sqrt(r0_q1/r_q1)
     &     *betah0(ihsps)*psitot/(ph(1,ihsps)*r1qp1)
      call elongs(elong,ind(nroot))
      write(*,*) 'Elongation at q=1 is ',elong
ckg elongation is added as per Fu,PoP,v.13,052517-1(2006)
      write(*,*) 'Mcclements dW_f=gam/om_ANOVA=',
     &     -res*q1*2*sqrt(3.)/sqrt((1+elong**2)*(1+5.*elong**2))
c-----------employ Porcelli model PPCF,v38,p.2163(1996).---------------
      do ir=1,ind(nroot)
         frt(ir)=xrt(ir)**1.5*php(1,ihsps)
      enddo
      call ASIMP(ppsi(1),frt(1),ind(nroot),res)
      res=res/sqrt(3.)*sqrt(r0_q1/r_q1)
     &     *betah0(ihsps)*psitot/(ph(1,ihsps)*r1qp1)
      write(*,*) 'Porcelli   dW_f=gam/om_ANOVA=',
     &     -res*q1*2*sqrt(3.)/sqrt((1+elong**2)*(5.+elong**2))
      return
      end
c**************************************************************************
