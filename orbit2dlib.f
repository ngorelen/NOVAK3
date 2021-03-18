c     nontest function for bounce averaging
c  Here:
c        ff(1)  is    Re(ksi_psi) x A + Re(dp) x C
c        ff(2)  is    Re(ksi_psi) x A + Re(dp) x C with FLR
c        ff(3)  is    Re(ksi_srf) x A 
c        ff(4)  is    Re(ksi_srf) x A with FLR             
c     for asymmetric case in addition we have imaginary contribution to
c     the previous
c        ff(5)  for   Re(ksi_psi) x A + Re(dp) x C         
c        ff(6)  for   Re(ksi_psi) x A + Re(dp) x C with FLR
c        ff(7)  for   Re(ksi_srf) x A                      
c        ff(8)  for   Re(ksi_srf) x A with FLR
      subroutine fgm(ff,iav,r,symm,vres,z)
      common/cgm/b0ax,brz,pmu,qpf,z0,imtae,imt
      dimension ff(iav)
      call gmoo(b0ax,brz,imtae,pmu,qpf,r,symm,vres,z,z0,ff(1),ff(2),ff(3
     &     ),ff(4),ff(5))
      return
      end
c*********************************************************
      subroutine gmoo(b0ax,brz,imtae,pmu,qpf,r1oo,symm,vres,z1oo,z0,gmo
     &     ,gmoflr,gmoi,gmoflri,ff)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      common/asymm/asym /fast/xav,xavint,x2r,x2i,psioo,kfast
      common/functn/omreal,eigenf(nn,mt,3),eigenf_s(nn,mt,3)
      common/funct/xkpop(nn,mt,nts,3)
      common/wdr/alpha,b0axxx,psedge,vioo,sigvpr
      dimension cl(2),cu(2),ff(4)
      common/saveij/i,j,cl,cu
      gmo=0.
      gmoflr=0.
      gmoi=0.
      gmoflri=0.
      call zeroo(ff(1),4)
      if(abs(asym).lt..1) asym=1.
      if(imtae.le.0) imtae=1
      roo=r1oo
      zoo=z0+abs(z1oo-z0)*asym
      if(kfast.gt.0.and.i.gt.0) goto 10
      call ij(i,j,roo,zoo,rax,z0,kfast)
      if(i.le.0) return
      if(i.le.1) then
         zoo=zz(i+1,j)
         roo=rr(i+1,j)
      endif
      z1=zoo
      r1=roo
      if(i.eq.1) then
         if(sqrt((r1-rax)**2+(z1-z0)**2).lt.1.05*abs(rr(1,1)-rr(1,mth/2)
     &        )/2)then
            if(cl(1).lt.0..or.cl(2).lt.0..or.cu(1).lt.0..or.cu(2).lt.0.
     &        )then
c               write(*,*)'1cl,cu,r1,z1,ru,zu,rl,zl=',cl,cu,r1,z1,rr(i,j)
c     &              ,rr(i+1,j),zz(i,j),zz(i+1,j),rr(i,j+1),rr(i+1,j+1)
c     &              ,zz(i,j+1),zz(i+1,j+1)
               cl(1)=0.5
               cu(1)=0.5
               cl(2)=0.
               cu(2)=0.
            endif
         else
            call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      else
         jc1=(j+jd0)/jdj
         if(jc1.eq.1.or.jc1.eq.3) then
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),zz(i,j),zz(i,j+1),rr(i,j),rr(i,j+1),z1,r1,cl,cu)
         else
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      endif
call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),roo,zoo,cl,cu)
 10   hhh=pmu*brz/rax/b0ax
c      gmo=f2dperp(divxsi(i,imtae),divxsi(i,imtae),cl,cu)*
c     $     f2dperp(dnom(i,j),dnom(i,j+1),cl,cu)
c Here       eigenf(i,im,1) is perturbed Ksi_n,
c            eigenf(i,im,2) is perturbed pressure,
c            eigenf(i,im,3) is perturbed Ksi_s
      besarg0=alpha/vioo*0.2*sqrt(hhh)*vres*b0ax/brz
c      write(*,*) '---------------------------',besarg0
c     &     ,alpha,vioo,vres,b0ax,brz,hhh
      besarg=besarg0*f2dperp(xkpop(i,imtae,j,2),xkpop(i,imtae,j+1
     &     ,2),cl,cu)
      bes2=FBES(besarg,bes0,2)
      gmo=f2dperp(eigenf(i,imtae,2),eigenf(i,imtae,2),cl,cu)/brz**2
      if(symm.eq.0.)then
         ff(1)=f2dperp(eigenf_s(i,imtae,2),eigenf_s(i,imtae,2
     &        ),cl,cu)/brz**2
         ff(1)=ff(1)*hhh
         ff(2)=ff(1)*(bes2+bes0)
      endif
      gmo=gmo*hhh
      gmoflr=gmo*(bes2+bes0)
      hhh=(1.-0.5*hhh)
c
c      write(*,*) 'besarg,bes0,bes2,gmo=',besarg,bes0,bes2,gmo,gmoflr
      besarg=besarg0*f2dperp(xkpop(i,imtae,j,1),xkpop(i,imtae,j+1
     &     ,1),cl,cu)
      bes2=FBES(besarg,bes0,0)
cnng18 recovering Larmor radius from kperp*rho to be used in radial FLR correction
c besarg0 m.b. in [m] but then becomes in psi 
      flrprsi=besarg0*sqrt(grpssq2d(i,j))   !i.e. this is the Larmor radius in psi
      psiflr=f2dperp(ppsi(i),ppsi(i),cl,cu) !this is psi at the center of the circle
      iflr=max(indx(ppsi,nosurf,psiflr)-1,1)
      call fun4(ppsi(iflr),eigenf(iflr,imtae,1),psiflr,eignfc0)
      Nflr=8
      aint=0.
      do icc=1,Nflr               !integrating only in the upper part of the Larmor circle шаг по потоку не сглажен
         psiflrc=psiflr+flrprsi/abs(psimin-psilim)* !psi on the integration path for the Larmor circle
     &        cos(3.1415926*(icc-1)/Nflr)
         iflr=max(indx(ppsi,nosurf,psiflrc)-1,1)
         iflr=min(iflr,nosurf-3)
         call fun4(ppsi(iflr),eigenf(iflr,imtae,1),psiflrc,eignfc)
         aint=aint+eignfc
      enddo
      if(eignfc0.ge.0.)then
         aint=aint/(Nflr*amax1(eignfc0,1.e-4))
      else
         aint=-aint/(Nflr*amax1(-eignfc0,1.e-4))
      endif
c      write(*,*) '2nd besarg,bes0=',besarg,bes0
      besargflr=f2dperp(eigenf(i,imtae,1),eigenf(i,imtae,1),cl,cu)
c      besargflr=aint
      besarg=hhh*besargflr
     &     *f2dperp(anom(i,j),anom(i,j+1),cl,cu)*aint
cnng18
      gmo=gmo+besarg
      gmoflr=gmoflr+bes0*besarg
      if(symm.eq.0.)then
         besarg=hhh*f2dperp(eigenf_s(i,imtae,1),eigenf_s(i,imtae,1),cl
     &        ,cu)*f2dperp(anom(i,j),anom(i,j+1),cl,cu)
         ff(1)=ff(1)+besarg
         ff(2)=ff(2)+bes0*besarg
      endif
c
c      write(*,*) '2nd gmo,gmoflr',gmo,gmoflr
c
c      write(*,*) 'besarg=',besarg
      besarg=besarg0*f2dperp(xkpop(i,imtae,j,3),xkpop(i,imtae,j+1
     &     ,3),cl,cu)
      bes2=FBES(besarg,bes0,0)
c
c      write(*,*) '3rd besarg,bes0=',besarg,bes0
c 
      if(symm.eq.0.) ff(3)=hhh*f2dperp(eigenf_s(i,imtae,3),eigenf_s(i
     &     ,imtae,3),cl,cu)
      hhh=hhh*f2dperp(eigenf(i,imtae,3),eigenf(i,imtae,3
     &     ),cl,cu)
      bes2=f2dperp(bnom(i,j),bnom(i,j+1),cl,cu)
      if(symm.eq.0.)then
         ff(3)=ff(3)*bes2
         ff(4)=ff(3)*bes0
      endif
      bes2=bes2*hhh
      gmoi=bes2
      gmoflri=bes0*bes2
      return
      end
c***********************************************************************
c wvpr is output for the coefficient 
      subroutine wdif(alpha1,b01,eps,epsp,iter,ntgr,r0,rmi,rma,pmu,
     $     nmu,psiav,psmi,psma,qpf,npf,rf,rvpar0,symm,sigm,t0,vi,v,nv
     &     ,wvpr,wqvpr,wd,wb,xhp,z0)
      parameter (iav=4)
      common/wdr/alpha,b0,psedge,vio,sigvpr
      dimension fo(iav),pmu(nmu),qpf(npf),v(nv),wvpr(nv,nmu,npf),
     $     wqvpr(nv,nmu,npf),wd(nv,nmu,npf),wb(nv,nmu,npf)
     $     ,rf(2,6,nmu,nv),rfcur(2,6),psiav(nv,nmu,npf),rvpar0(nmu)
      external fwd
      alpha=alpha1
      vio=vi
      b0=b01
      psedge=psma
      do ipf=1,npf
      if(sigm.gt.0..or.qpf(ipf).gt.0.) then
ctem
         do imu=1,nmu  
c         do imu=48,48
         if(sigm.gt.0..or.pmu(ipf).lt.r0)then
ctem
            do iv=1,nv 
c            do iv=19,20
               do i2=1,2
                  do i6=1,6
                     rfcur(i2,i6)=rf(i2,i6,imu,iv)
                  enddo
               enddo
               if(pmu(ipf).lt.r0.and.pmu(ipf).gt.rmi)then
                  ntgr1=ntgr*2
               else
                  ntgr1=ntgr
               endif
               call average(iav,iter,eps,epsp,fwd,fo,ntgr1,psma,psmi
     $              ,qpf(ipf),pmu(imu),r0,rma,rmi,rfcur(1,1)
     $              ,rvpar1,sigm,symm,t0,tb,v(iv),z0,xhp)
               if(tb.gt.1.e-10)then
                  wb(iv,imu,ipf)=2.*3.1415926/tb
               else
                  wb(iv,imu,ipf)=0.
               endif
c               if(iv.eq.nv) then
c                  write(*,*)fo,tb,imu,iv
c               endif
c               write(*,*)'---',iv,imu,tb,fo
c     $              ,qpf(ipf),pmu(imu),tb
c               write(*,*)'--------------------------------'
c               write(*,*)wd(iv,imu,ipf),wqvpr(iv,imu,ipf),wvpr(iv,imu
c     &              ,ipf),tb
               wd(iv,imu,ipf)=fo(3)
               wqvpr(iv,imu,ipf)=fo(2)
               wvpr(iv,imu,ipf)=fo(1)
               psiav(iv,imu,ipf)=fo(4)
c               write(*,*)wd(iv,imu,ipf),wqvpr(iv,imu,ipf),wvpr(iv,imu
c     &              ,ipf),tb
            enddo
            rvpar0(imu)=rvpar1
         endif
         enddo
      endif
      enddo
ctem       stop
      return
      end
c*********************************************************
      function f2d(ro1,ro2,ro3,ro4,r12,r13,r14,r32,r42,r34,f12,f34)
      dimension f12(2),f34(2)
      f2d=ro3*ro4
c      f2d=(f12(1)*ro2*f2d+f12(2)*ro1*f2d+f34(1)*ro4*ro1*ro2+
c     +     f34(2)*ro3*ro1*ro2)/(r12*r13*r14+r12*r32*r42+r13*r32*r34
c     +     +r14*r42*r34)
      f2d=(f12(1)*ro2*f2d+f12(2)*ro1*f2d+
     +     f34(1)*ro4*ro1*ro2+f34(2)*ro3*ro1*ro2)/(ro2*f2d+ro1*f2d
     +     +ro4*ro1*ro2+ro3*ro1*ro2)
      return
      end
c*********************************************************************
c     iav is number of averaged functions
c     ntgr is number of integration grid ~10
c     iter is iteration for orbit evaluation ~2-3
c
c
c
c
c
c
      subroutine average(iav,iter1,eps1,epsp,f,fo,ntgr,psma,psmi,qpf1
     $     ,pmu11,r01,rma,rmi,rf,rvpar1,sig1,symm,t0,tb,v1,z01,xhp1)
      external f
      dimension fo(iav),fo1(20),rf(2,6)
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      common/asymm/asym /fast/xav,xavint,x2r,x2i,psioo,kfast
      common/wdr/alpha,b0,psedge,vio,sigvpr
      data pi/3.1415926/
      pmu1=pmu11
 13   tb=0.
      tb1=0.
      do i=1,iav
         fo(i)=0.
         fo1(i)=0.
      enddo
      xav1=0.
      xav=0.
      xavint=0.
      x2r1=0.
      x2r=0.
      x2i=0.
c      write(*,*)' pmu11',pmu11,pmu,rma
      if(pmu1.ge.rma) return
      asym=1.
      iter=iter1
      eps=eps1
      qpf=qpf1
      pmu=pmu1
      r0=r01
      z0=z01
      call rvpar0(eps,qpf,pmu,rma,rmi,r0,rvpar,z0,zvpar0,iter)
c      write(*,*)'0 ', eps,qpf,pmu,rma,rmi,r0,rvpar,z0,zvpar0,iter
      sig=sig1
      v=v1
      xhp=xhp1
      rvpar1=rvpar
      if(rf(1,5).gt.0.) then
         rf(1,3)=amax1(rf(1,6)-rf(1,5)+rvpar,rvpar)
      else
         rf(1,3)=amax1(rf(1,6),rvpar)
      endif
      call ppcorn(a,b,iconf,psma,psmi,rma,rmi,rf,sig1)
      if(rvpar.gt.a.or.rvpar.gt.b.or.a.eq.b)return
c      write(*,*) 'a,b,rf,qpf,pmu,v1',a,b,rf,qpf,pmu,v1
      sig=sig1
c      if(sig1.lt.0..and.qpf.gt..5.and.abs(pmu-2.7).lt.1.) then
c         write(*,*)'a,b ',a,b,rmi,qpf,pmu,
c     $     sig,v,iconf,psma,psmi,iter,eps,r0,xhp,rma
c      endif
 10   if((a.lt.rmi.or.b.lt.rmi)
c.or.abs(a-r0).lt.r0*.01.or.abs(b-r0).lt.r0*.01
     &     ) then
         return
cccc         if(iconf.lt.1.or.sig.gt.0.) return
cccc         sig=1.
cccc         call ppcorn(a,b,iconf,psma,psmi,rma,rmi,rf,1.)
cccc         write(*,*)'0a,b,rmi,pmu,rvpar,sig,sig1',a,b,rmi,pmu,rvpar,
cccc     $        sig,sig1
cccc         goto 10
      endif
cccc      if(a.lt.rmi) return 
      sigvpr=sig1
      if(sig1.lt.0..and.sig.lt.0.)then
         call summation(a,b,ntgr,f,pi,eps,epsp,iter,qpf,pmu,r0,rvpar
     $        ,sig1,v,z0,xhp,iav,t0,tb,fo)
         xav1=xav
         x2r1=x2r
         if(tb.eq.0.)return
         if(symm.eq.0.) then
            asym=-1.
            call summation(a,b,ntgr,f,pi,eps,epsp,iter,qpf,pmu,r0,rvpar
     $           ,sig1,v,z0,xhp,iav,t0,tb1,fo1(1))
            if(tb1.eq.0.)return
            tb1=abs(tb1)
            do i=1,iav
               fo(i)=(fo(i)*tb+fo1(i)*tb1)/(tb+tb1)
            enddo
            xav1=(xav1*tb+xav*tb1)/(tb+tb1)
            x2r1=(x2r1*tb+x2r*tb1)/(tb+tb1)
            tb=tb+tb1
            asym=1.
         else
            tb=2.*tb
         endif
         return
      endif
      amh=epsp*rvpar
      fpsi=psio(dpsi,ddpsi,dzpsi,rvpar,r0,z0,z0,am)
      sigvpr=sig1
      if(fpsi.gt.qpf.or.abs(a-rvpar).le.amh.or.rvpar.le.rmi) then
c         write(*,*) '----',a,b,amh
         if(abs(a-rvpar).le.amh) then
            if(fpsi.gt.qpf)then
               sigvpr=sig1
            else
               sigvpr=1.
            endif
            am=a+amh
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,am
     $           ,r0,rvpar,sig1,v,zvpar,z0,xhp)
            tb1=abs(zvpar-z0)*2.*bm/b0*am*t0/xhp
            kfast=-1
            call f(fo1(1),iav,a+amh/2.,amh/2.)
            if(symm.eq.0.) then
               asym=-1.
               call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,
     $              psior,am,r0,rvpar,sig1,v,zvpar,z0,xhp)
               tb=abs(zvpar-z0)*2.*bm/b0*am*t0/xhp
               kfast=-1
               call f(fo(1),iav,a+amh/2.,amh/2.)
               if(tb.eq.0.) return
               do i=1,iav
                  fo1(i)=(fo(i)*tb+fo1(i)*tb1)/(tb+tb1)
               enddo
               tb1=tb+tb1
               asym=1.
            endif
            if(fpsi.lt.qpf.and.pmu.gt.r0) am=rvpar
         else
            am=a
         endif
corb         am=a
c         write(*,*) am,b,amh
         sigvpr=sig1
         if(am.lt.b-amh*2.) then
            call summation(am,b,ntgr,f,pi,eps,epsp,iter,qpf,pmu,r0
     $           ,rvpar,1.,v,z0,xhp,iav,t0,tb,fo)
c
            if(tb.eq.0.)return
            if(symm.eq.0.)then
               tb1=abs(tb1)
               do i=1,iav
                  fo1(i)=(fo(i)*tb+fo1(i)*tb1)/(tb+tb1)
               enddo
               xav1=(xav*tb+xav1*tb1)/(tb+tb1)
               x2r1=(x2r*tb+x2r1*tb1)/(tb+tb1)
               tb1=tb+tb1
               asym=-1.
               call summation(am,b,ntgr,f,pi,eps,epsp,iter,qpf,pmu,r0
     $              ,rvpar,1.,v,z0,xhp,iav,t0,tb,fo)
               if(tb.eq.0.)return
               tb=abs(tb)
               asym=1.
            endif
         else
c     may be improved for the bounce time of well trapped particles
            if(pmu1.eq.pmu11.and.pmu11.gt.r0+amh*2.)
     $           then
               pmu1=amax1(pmu11-amh*5./r0,0.)
               write(*,*)' more step',pmu1,pmu11
               goto 13
            endif
            r=(a+b)/2.
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior
     $           ,r,r0,rvpar,sig1,v,zvpar,z0,xhp)
            if(zvpar-z0.ge.0.) then
               zvpar=z0+max(a*epsp,zvpar-z0)
            else
               zvpar=z0+min(-a*epsp,zvpar-z0)
            endif
            tb=psio(dpsi,ddpsi,dzpsi,r,r0,zvpar,z0,am)
c            if(rvpar.ge.r)write(*,*) 'Negative',rvpar,r,dzpsi,a,b
            tb=sqrt(1.-rvpar/r)*abs(dzpsi)
            tb=max(.001,tb)
            tb=bm*r/b0/v/tb*abs(a-b)*t0*2.
            kfast=-1
            call f(fo(1),iav,r,amh/2.)
            return
         endif
      else
c trapped particle contribution
         sigvpr=-1.
         call summation(rvpar,a,ntgr,f,pi,eps,epsp,iter,qpf,pmu,r0
     $        ,rvpar,-1.,v,z0,xhp,iav,t0,tb1,fo1(1))
         xav1=xav
         x2r1=x2r
         if(tb1.eq.0.)return
c     write(*,*) (fo1(i),i=1,iav),tb1
         sigvpr=1.
         call summation(rvpar,b,ntgr,f,pi,eps,epsp,iter,qpf,pmu,r0
     $        ,rvpar,1.,v,z0,xhp,iav,t0,tb,fo)
         if(tb1.eq.0..or.tb.eq.0.)then
            tb=0.
            return
         endif
         if(symm.eq.0.)then
            do i=1,iav
               fo(i)=(fo(i)*tb+fo1(i)*tb1)/(tb+tb1)
            enddo
            xav1=(xav*tb+xav1*tb1)/(tb+tb1)
            x2r1=(x2r*tb+x2r1*tb1)/(tb+tb1)
            tb=tb+tb1
            asym=-1.
            sigvpr=-1.
            call summation(rvpar,a,ntgr,f,pi,eps,epsp,iter,qpf,pmu,r0
     $           ,rvpar,-1.,v,z0,xhp,iav,t0,tb1,fo1(1))
            if(tb1.eq.0.)return
            tb1=abs(tb1)
            do i=1,iav
               fo(i)=(fo(i)*tb+fo1(i)*tb1)/(tb+tb1)
            enddo
            xav1=(xav1*tb+xav*tb1)/(tb+tb1)
            x2r1=(x2r1*tb+x2r*tb1)/(tb+tb1)
            tb=tb+tb1
            sigvpr=1.
            call summation(rvpar,b,ntgr,f,pi,eps,epsp,iter,qpf,pmu,r0
     $           ,rvpar,1.,v,z0,xhp,iav,t0,tb1,fo1(1))
            if(tb1.eq.0.)return
            hhh=xav1
            xav1=xav
            xav=hhh
            hhh=x2r1
            x2r1=x2r
            x2r=hhh
            tb1=abs(tb1)
            asym=1.
         endif
c         write(*,*) fo,tb
      endif
      do i=1,iav
         fo(i)=(fo(i)*tb+fo1(i)*tb1)/(tb+tb1)
      enddo
      xav=(xav*tb+xav1*tb1)/(tb+tb1)
      x2r=(x2r*tb+x2r1*tb1)/(tb+tb1)
      tb=tb+tb1
      if(symm.eq.1.) tb=tb*2.
      return
      end
c******************************************************************
        SUBROUTINE deriv4(x,y,xn,dyn)
c Four points Lagrange' formula for derivatives
        dimension x(4),y(4)
        p=(xn-x(2))/(x(3)-x(2))
        p2=p*p
        dyn=(-(0.5*p2-p+1./3.)*y(1)+(1.5*p2-2.*p-0.5)*y(2)-
     -  (1.5*p2-p-1.)*y(3)+(0.5*p2-1./6.)*y(4))/(x(3)-x(2))
        return
        end
C***********************************************************************
        SUBROUTINE funder4(x,f,xn,dyn,yn)
c Four points Lagrange' formula for derivatives on nonequal grid
c and for the function itself
        dimension x(4),f(4)
        x1=xn-x(1)
        x2=xn-x(2)
        x3=xn-x(3)
        x4=xn-x(4)
        yn=f(1)/(x(1)-x(3))/(x(1)-x(4))
        dyn=yn*(x3*x4+x2*x3+x2*x4)
        hhh=f(2)/(x(2)-x(3))/(x(2)-x(4))
        dyn=(dyn-hhh*(x3*x4+x1*x3+x1*x4))/(x(1)-x(2))
        yn=(yn*x2-hhh*x1)*x3*x4/(x(1)-x(2))
        hhh=f(3)/(x(3)-x(2))/(x(3)-x(4))/(x(3)-x(1))
        dyn=dyn+hhh*(x2*x4+x1*x2+x1*x4)
        yn=yn+hhh*x2*x1*x4
        hhh=f(4)/(x(4)-x(2))/(x(4)-x(3))/(x(4)-x(1))
        dyn=dyn+hhh*(x2*x3+x1*x2+x1*x3)
        yn=yn+hhh*x2*x3*x1
        return
        end
c************************************************************************
      FUNCTION FINTR(F,XF,NPF,NINTR,X)
c if npf eq. nintr, then their is'nt dependence from nintr
c and if abscissas is equally spaced then interpolation is Lagrange
c interpolation
      DIMENSION F(NPF) ,XF(NPF)
      IF(X.GT.XF(1)) GOTO 10
      fINTR=F(1)
      GOTO 100
 10   IF(X.LT.XF(NPF)) GOTO 20
      FINTR=F(NPF)
      GOTO 100
 20   IF(NPF.GT.NINTR) GOTO 25
         I1=1
         I2=NPF
         NN=NPF
         GOTO 60
 25         DO 30 I=2,NPF
            IF(.NOT.(XF(I-1).LE.X.AND.XF(I).GE.X)) GOTO 30
            I1=I-1
            I2=I 
            GOTO 35 
 30         CONTINUE
 35   ID=int((real(NINTR)+0.5)/2.)-1
      I1=I1-ID
      I2=I1+NINTR-1
      IF(I1.LT.1)    I1=1 
      IF(I2.GT.NPF)  I2=NPF 
      NN=I2-I1+1
 60   CALL POLINT(F(I1),XF(I1),NN,X,R)
      FINTR=R
  100 RETURN
      END
C***********************************************************************
	SUBROUTINE POLINT(F,A,N,X,R)
c if abscissas is equally spaced then interpolation is Lagrange
c N-points interpolation
	DIMENSION F(N), A(N)
	R=0.0
	DO 1 J=1,N
	   AL=1.0
	   DO 2 I=1,N
	      IF(I-J) 3,2,3
 3	      AL=AL*(X-A(I))/(A(J)-A(I))
 2	   CONTINUE
 1	R=R+AL*F(J)
	RETURN
	END
c************************************************************************
c  
        subroutine rev(yy1,n,m)
        dimension yy1(n,m)
        do j=1,m
           do i=1,n/2
              hhh=yy1(i,j)
              yy1(i,j)=yy1(n-i+1,j)
              yy1(n-i+1,j)=hhh
           enddo
        enddo
        return
        end
c*********************************************************
c test function for bounce averaging
      subroutine f(ff,iav,r,z)
      common/asymm/asym
      dimension ff(iav)
      ff(1)=(r-3.52)/sqrt((r-3.52)**2+z**2)
      if(iav.gt.1)ff(2)=r
csqrt((r-3.)**2+z**2)
      if(iav.gt.2)ff(3)=z
      return
      end
c*********************************************************
      SUBROUTINE ROOT(xmin,xmax,xeps,feta,xroot,f)
      external f
      common/root2/kfast
      kfast=0
      call c05ade(xmin,xmax,xeps,feta,f,xr,ifail)
      if(ifail.ne.-13) then
         xroot=(xr)
      else
         xroot=0.
      endif
      return
      end
c*********************************************************
      SUBROUTINE ROOT1(xmin,xmax,xeps,feta,xroot,f)
      external f
      common/root2/kfast
      kfast=0
      call c06ade(xmin,xmax,xeps,feta,f,xr,ifail)
      if(ifail.ne.-13) then
         xroot=(xr)
      else
         xroot=0.
      endif
      return
      end
c*********************************************************
      SUBROUTINE C05ADE(A, B, EPS, ETA, F, X, IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     .. SCALAR ARGUMENTS ..
      REAL A, B, EPS, ETA, X
      INTEGER IFAIL
C     .. FUNCTION ARGUMENTS ..
      external F
      REAL F
C     ..
C     .. LOCAL SCALARS ..
      character*8 SRNAME
c      DOUBLE PRECISION SRNAME
      REAL FX, Y
      INTEGER IFAIL1, IND, IR
C     .. LOCAL ARRAYS ..
      REAL C(17)
C     .. FUNCTION REFERENCES ..
      INTEGER P01AAE
C     .. SUBROUTINE REFERENCES ..
C     C05AZE
C     ..
      DATA SRNAME /'  C05ADE'/
C     DRIVER FOR C05AZE
      IFAIL1 = 1
C     INPUT ERROR
      IF (A.EQ.B .OR. EPS.LE.0.) GO TO 120
      X = A
      FX = F(X)
      IFAIL1 = 0
C     ZERO AT INITIAL POINT
      IF (ABS(FX).LT.ETA .OR. FX.EQ.0.) GO TO 120
      Y = X
      C(1) = FX
      X = B
      FX = F(X)
C     ZERO AT INITIAL POINT
      IF (ABS(FX).LT.ETA .OR. FX.EQ.0.) GO TO 120
      IFAIL1 = 1
C     NO ROOT IN RANGE
      IF (SIGN(1.,FX).EQ.SIGN(1.,C(1))) then
         ifail=-13
         return
      endif
cfastGO TO 120
      IR = 0
      IND = -1
   20 CALL C05AZE(X, Y, FX, EPS, IR, C, IND, IFAIL1)
      IF (IND.EQ.0) GO TO 40
      FX = F(X)
      IF (ABS(FX).GE.ETA .AND. FX.NE.0.) GO TO 20
C     ZERO HIT EXACTLY
      IFAIL1 = 0
      GO TO 120
   40 IF (IFAIL1.EQ.0) GO TO 120
      GO TO (60, 60, 60, 100, 80), IFAIL1
C     IMPOSSIBLE EXIT
   60 IFAIL1 = 4
      GO TO 120
C     TOO MUCH ACCURACY REQUESTED
   80 IFAIL1 = 2
      GO TO 120
C     PROBABLY A POLE
  100 IFAIL1 = 3
  120 IFAIL = P01AAE(IFAIL,IFAIL1,SRNAME)
      RETURN
      END
c*********************************************************
      SUBROUTINE C06ADE(A, B, EPS, ETA, F, X, IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     .. SCALAR ARGUMENTS ..
      REAL A, B, EPS, ETA, X
      INTEGER IFAIL
C     .. FUNCTION ARGUMENTS ..
      external F
      REAL F
C     ..
C     .. LOCAL SCALARS ..
      character*8 SRNAME
c      DOUBLE PRECISION SRNAME
      REAL FX, Y
      INTEGER IFAIL1, IND, IR
C     .. LOCAL ARRAYS ..
      REAL C(17)
C     .. FUNCTION REFERENCES ..
      INTEGER P01AAE
C     .. SUBROUTINE REFERENCES ..
C     C05AZE
C     ..
      DATA SRNAME /'  C05ADE'/
C     DRIVER FOR C05AZE
      IFAIL1 = 1
C     INPUT ERROR
      IF (A.EQ.B .OR. EPS.LE.0.) GO TO 120
      X = A
      FX = F(X)
      IFAIL1 = 0
C     ZERO AT INITIAL POINT
      IF (ABS(FX).LT.ETA .OR. FX.EQ.0.) GO TO 120
      Y = X
      C(1) = FX
      X = B
      FX = F(X)
C     ZERO AT INITIAL POINT
      IF (ABS(FX).LT.ETA .OR. FX.EQ.0.) GO TO 120
      IFAIL1 = 1
C     NO ROOT IN RANGE
      IF (SIGN(1.,FX).EQ.SIGN(1.,C(1)))  then
         ifail=-13
         return
      endif
cfastGO TO 120
      IR = 0
      IND = -1
   20 CALL C05AZE(X, Y, FX, EPS, IR, C, IND, IFAIL1)
      IF (IND.EQ.0) GO TO 40
      FX = F(X)
      IF (ABS(FX).GE.ETA .AND. FX.NE.0.) GO TO 20
C     ZERO HIT EXACTLY
      IFAIL1 = 0
      GO TO 120
   40 IF (IFAIL1.EQ.0) GO TO 120
      GO TO (60, 60, 60, 100, 80), IFAIL1
C     IMPOSSIBLE EXIT
   60 IFAIL1 = 4
      GO TO 120
C     TOO MUCH ACCURACY REQUESTED
   80 IFAIL1 = 2
      GO TO 120
C     PROBABLY A POLE
  100 IFAIL1 = 3
  120 IFAIL = P01AAE(IFAIL,IFAIL1,SRNAME)
      RETURN
      END
c********************************************************************
      SUBROUTINE C05AZE(X, Y, FX, TOLX, IR, C, IND, IFAIL)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1979.
C     .. SCALAR ARGUMENTS ..
      REAL FX, TOLX, X, Y
      INTEGER IFAIL, IND, IR
C     .. ARRAY ARGUMENTS ..
      REAL C(17)
C     ..
C     .. LOCAL SCALARS ..
      character*8 SRNAME
c      DOUBLE PRECISION SRNAME
      REAL AB, DIFF1, DIFF2, DIFF, MAX, REL, TOL
      INTEGER I
      LOGICAL T
C     .. FUNCTION REFERENCES ..
      REAL SQRT, X02AAE, X02ABE
      INTEGER P01AAE
C     ..
      DATA SRNAME /'  C05AZE'/
      I = 0
      IF ((IND.GT.0 .AND. IND.LE.4) .OR. IND.EQ.-1) GO TO 20
C     USER NOT CHECKED IND OR CHANGED IT
      I = 2
      IND = 0
      GO TO 640
   20 IF (TOLX.GT.0. .AND. (IR.EQ.0 .OR. IR.EQ.1 .OR. IR.EQ.2)) GO
     * TO 40
      I = 3
      IND = 0
      GO TO 640
   40 REL = 1.
      AB = 1.
      IF (IR.EQ.1) REL = 0.
      IF (IR.EQ.2) AB = 0.
      IF (IND.EQ.-1) GO TO 80
      GO TO (60, 100, 180, 480), IND
   60 C(3) = X
      IND = 2
      RETURN
   80 C(3) = X
  100 IF (FX.NE.0.) GO TO 140
  120 Y = X
      IND = 0
      I = 0
      GO TO 640
  140 C(4) = FX
      C(15) = ABS(FX)
      C(16) = 0.
      X = Y
      Y = C(3)
      C(2) = C(4)
      C(5) = X
      IF (IND.EQ.-1) GO TO 160
      IND = 3
      RETURN
  160 FX = C(1)
      IND = 3
  180 IF (FX.EQ.0.) GO TO 120
      IF (SIGN(1.,FX).NE.SIGN(1.,C(2))) GO TO 200
      IND = 0
      I = 1
      GO TO 640
  200 C(6) = FX
      C(13) = SQRT(X02AAE(0.0))
      C(15) = AMAX1(C(15),ABS(FX))
      C(14) = X02ABE(0.0)
      C(16) = 0.0
  220 C(1) = C(5)
      C(2) = C(6)
      C(17) = 0.
  240 IF (ABS(C(2)).GE.ABS(C(4))) GO TO 280
      IF (C(1).EQ.C(5)) GO TO 260
      C(7) = C(5)
      C(8) = C(6)
  260 C(5) = C(3)
      C(6) = C(4)
      X = C(1)
      C(3) = X
      C(4) = C(2)
      C(1) = C(5)
      C(2) = C(6)
  280 TOL = 0.5*TOLX*AMAX1(AB,REL*ABS(C(3)))
      DIFF2 = 0.5*(C(1)-C(3))
      C(12) = DIFF2
      DIFF2 = DIFF2 + C(3)
      IF (C(12).EQ.0.) GO TO 340
      IF (ABS(C(12)).LE.TOL) GO TO 580
      IF (C(17).LT.2.5) GO TO 300
      C(11) = C(12)
      GO TO 460
  300 TOL = TOL*SIGN(1.,C(12))
      DIFF1 = (C(3)-C(5))*C(4)
      IF (C(17).GT.1.5) GO TO 320
      DIFF = C(6) - C(4)
      GO TO 380
  320 IF (C(7).NE.C(3) .AND. C(7).NE.C(5)) GO TO 360
  340 IND = 0
      I = 5
      GO TO 640
  360 C(9) = (C(8)-C(4))/(C(7)-C(3))
      C(10) = (C(8)-C(6))/(C(7)-C(5))
      DIFF1 = C(10)*DIFF1
      DIFF = C(9)*C(6) - C(10)*C(4)
  380 IF (DIFF1.GE.0.) GO TO 400
      DIFF1 = -DIFF1
      DIFF = -DIFF
  400 IF (ABS(DIFF1).GT.C(14) .AND. DIFF1.GT.DIFF*TOL) GO TO 420
      C(11) = TOL
      GO TO 460
  420 IF (DIFF1.GE.C(12)*DIFF) GO TO 440
      C(11) = DIFF1/DIFF
      GO TO 460
  440 C(11) = C(12)
  460 C(7) = C(5)
      C(8) = C(6)
      C(5) = C(3)
      C(6) = C(4)
      C(3) = C(3) + C(11)
      X = C(3)
      Y = C(1)
      IND = 4
      RETURN
  480 IF (FX.EQ.0.) GO TO 120
      C(4) = FX
      MAX = ABS(FX)
      IF (C(13)*MAX.LE.C(15)) GO TO 500
      IF (C(16).EQ.1.) C(16) = -1.
      IF (C(16).EQ.0.) C(16) = 1.
      GO TO 520
  500 C(16) = 0.
  520 IF (C(2).GE.0.) GO TO 540
      T = C(4).LE.0.
      GO TO 560
  540 T = C(4).GE.0.
  560 IF (T) GO TO 220
      I = IFIX(C(17)+0.1)
      I = I + 1
      IF (C(11).EQ.C(12)) I = 0
      C(17) = FLOAT(I)
      GO TO 240
  580 IF (C(16).GE.0.) GO TO 600
      I = 4
      GO TO 620
  600 Y = C(1)
      I = 0
  620 IND = 0
  640 IFAIL = P01AAE(IFAIL,I,SRNAME)
      RETURN
      END
c********************************************************************
      INTEGER FUNCTION P01AAE(IFAIL, ERROR, SRNAME)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 3 REVISED
C     MARK 4A REVISED, IER-45
C     MARK 4.5 REVISED
C     MARK 7 REVISED (DEC 1978)
C     RETURNS THE VALUE OF ERROR OR TERMINATES THE PROGRAM.
      INTEGER ERROR, IFAIL, NOUT
c+SELF,IF=IBM.
      character*8 SRNAME
c      DOUBLE PRECISION SRNAME
c+SELF.
C     TEST IF NO ERROR DETECTED
      IF (ERROR.EQ.0) GO TO 20
C     DETERMINE OUTPUT UNIT FOR MESSAGE
      CALL X04AAE (0,NOUT)
C     TEST FOR SOFT FAILURE
      IF (MOD(IFAIL,10).EQ.1) GO TO 10
C     HARD FAILURE
c      WRITE (NOUT,99999) SRNAME, ERROR
C       STOPPING MECHANISM MAY ALSO DIFFER
	P01AAE = ERROR
	return
c      STOP
C     SOFT FAIL
C     TEST IF ERROR MESSAGES SUPPRESSED
   10 IF (MOD(IFAIL/10,10).EQ.0) GO TO 20
      WRITE (NOUT,99999) SRNAME, ERROR
   20 P01AAE = ERROR
      RETURN
99999 FORMAT (1H0, 38HERROR DETECTED BY NAG LIBRARY ROUTINE , A8,
     * 11H - IFAIL = , I5//)
      END
c********************************************************************
      SUBROUTINE X04AAE(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     *** NOTE ***
C     THIS ROUTINE ASSUMES THAT THE VALUE OF NERR1 IS SAVED
C     BETWEEN CALLS.  IN SOME IMPLEMENTATIONS IT MAY BE
C     NECESSARY TO STORE NERR1 IN A LABELLED COMMON
C     BLOCK /AX04AA/ TO ACHIEVE THIS.
C
C     .. SCALAR ARGUMENTS ..
      INTEGER I, NERR
C     ..
C     .. LOCAL SCALARS ..
      INTEGER NERR1
C     ..
      DATA NERR1 /5/
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
c********************************************************************
      REAL FUNCTION X02AAE(X)
C     NAG COPYRIGHT 1975
C     MARK 4.5 RELEASE
      REAL X, Z
c      DATA Z/Z3C100000/
      DATA Z/2.e-38/
c      DATA Z/1641 4000 0000 0000 0000B/
C     * EPS *
C     RETURNS THE VALUE EPS WHERE EPS IS THE SMALLEST
C     POSITIVE
C     NUMBER SUCH THAT 1.0 + EPS > 1.0
C     THE X PARAMETER IS NOT USED
      X02AAE = x
      X02AAE = Z
      RETURN
      END
c+DECK,X02ABE.
      REAL FUNCTION X02ABE(X)
C     NAG COPYRIGHT 1975
C     MARK 4.5 RELEASE
      REAL X, Z
c+SELF,IF=R01.
c       DATA Z/Z00100000/
       DATA Z/2.e-38/
c+SELF,IF=R04.
c      DATA Z/0001 4000 0000 0000 0000B/
c+SELF.
C     * RMIN *
C     RETURNS THE VALUE OF THE SMALLEST POSITIVE REAL FLOATING-
C     POINT NUMBER EXACTLY REPRESENTABLE ON THE COMPUTER
C     THE X PARAMETER IS NOT USED
	X02ABE = x
      X02ABE = Z
      RETURN
      END
c************************************************************************
      subroutine t0xhp(alpha,t0,xhp,ai,b0,b0oo,vi,zi,psitot)
      xhp=.104563*vi*ai/zi/psitot/b0oo
      t0=0.1*b0/vi/psitot
      alpha=.52192*vi**2*ai/zi/b0/b0oo
      return
      end
C***********************************************************************
        function f2dperp(fu,fl,cu,cl)
        dimension fu(2),fl(2),cu(2),cl(2)
c        call weight(ru,rl,zu,zl,r,z,cu,cl)
        f2dperp=fl(1)*cl(1)+fl(2)*cl(2)+fu(1)*cu(1)+fu(2)*cu(2)
        return
        end
C***********************************************************************
	subroutine ini(y,n,yl,yh)
	dimension y(n)
	y(n)=(yh-yl)/real(n-1)
	y(1)=yl
	do i=2,n
	   y(i)=y(i-1)+y(n)
	enddo
	return
	end
c*********************************************************
	subroutine inistr(y,n,yl,yh,pk)
	dimension y(n)
c make linear in index first
	y(n)=(yh**(1./pk)-yl**(1./pk))/real(n-1)  ! **pk
	do i=1,n
	   y(i)=yl**(1./pk)+real(i-1)*y(n)
	enddo
c and then take a sqrt of it
	do i=1,n
	   y(i)=(y(i))**pk
	enddo
	return
	end
c*********************************************************
c This subr. makes interface between TAE (C.Z.Cheng) and
c New part ORBIT2D, which is to calculate perturbatively 
c the growth rates of TAE by direct 2D averaging along fast particle 
c     drift orbit
c
c     BE SURE TO INCLUDE XA ZA etc. INPUT bypassed BEFORE
      subroutine interfa
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      data dz0/.00001/
      save b0ax
      if(nsrf.eq.1) then
         rax=(xa(1)+xa(mth/2+1))/2.
         b0ax=(bf(1)+bf(mth/2+1))/2.
         z0=(za(1)+za(mth/2+1))/2.
      endif
      do i=1,mth+1
         rr(nsrf,i)=xa(i)
         zz(nsrf,i)=za(i)
         grthsq2d(nsrf,i)=grthsq(i)
         grpssq2d(nsrf,i)=grpssq(i)
         xjacob2d(nsrf,i)=xjacob(i)
         b2d(nsrf,i)=bf(i)
         shear2d(nsrf,i)=shear(i)
         if(abs(xtheta(i)).lt.dz0) then
            xthhh=dz0
         else
            xthhh=xtheta(i)
         endif
         dzpsi(nsrf,i)=abs((zpsi(i)-ztheta(i)*xpsi(i)/xthhh)
     $        /rr(nsrf,i))*psitot
         dpsio(nsrf,i)=rr(nsrf,i)*
     $        ztheta(i)/(xpsi(i)*ztheta(i)
     $        -xtheta(i)*zpsi(i))/psitot
c drift frequency and transit-bounce freq. staff
         wdo(nsrf,i,3)=rg/xjacob(i)
         wdo(nsrf,i,1)=-curvs(i)*bf(i)*gptgpp(i)
         wdo(nsrf,i,2)=wdo(nsrf,i,1)*xjacob(i)*rg/xsq(i)
     $        +curvn(i)*bf(i)-curvs(i)*bf(i)*qdelp(i)
         wdo(nsrf,i,1)=wdo(nsrf,i,1)+wdo(nsrf,i,3)*curvn(i)/bf(i)
         wdo(nsrf,i,4)=bsq(i)
         anom(nsrf,i)=b0ax/bf(i)
         wdo(nsrf,i,1)=wdo(nsrf,i,1)*anom(nsrf,i)
         wdo(nsrf,i,2)=wdo(nsrf,i,2)*anom(nsrf,i)
         anom(nsrf,i)=pp/rax/bf(i)**3
         wdo(nsrf,i,3)=wdo(nsrf,i,3)*anom(nsrf,i)
         wdo(nsrf,i,4)=wdo(nsrf,i,4)*anom(nsrf,i)
c     bounce average over the orbit of G_m coeff.
         anom(nsrf,i)=curvn(i)
         bnom(nsrf,i)=curvs(i)
         cnom(nsrf,i)=1./bsq(i)
         dnom(nsrf,i)=1.66666667/bsq(i)*p
      enddo
c This part defines the drift and bounce frequencies of alphas
      do i=1,lam
         wdchng(i,nsrf,1)=wd(i)
         wdchng(i,nsrf,2)=tkbou(i)
         alchng(i,nsrf,1)=gal(i)
      enddo
      do i=1,lamc
         wdcchng(i,nsrf,1)=wdc(lamc-i+1)
         wdcchng(i,nsrf,2)=tkcir(lamc-i+1)
         alchng(i,nsrf,2)=galc(lamc-i+1)
      enddo
      rr(nsrf,mth+2)=rr(nsrf,2)
      zz(nsrf,mth+2)=zz(nsrf,2)
      grthsq2d(nsrf,mth+2)=grthsq(2)
      grpssq2d(nsrf,mth+2)=grpssq(2)
      xjacob2d(nsrf,mth+2)=xjacob(2)
      b2d(nsrf,mth+2)=bf(2)
      shear2d(nsrf,mth+2)=shear(2)
      dzpsi(nsrf,mth+2)=dzpsi(nsrf,2)
      dpsio(nsrf,mth+2)=dpsio(nsrf,2)
      rg2d(nsrf)=rg
      qoo(nsrf)=q
      if(nsrf.eq.1) then
         ppsi(nsrf)=0.
      else
         ppsi(nsrf)=rrg**2
      endif
      r0=(xinf(1)+xinf(mth/2+1))/2.
      return
      end
c******************************************************************
c assume no more then two roots
      SUBROUTINE roots(y,x,ind,n,y0,xeps,feta,xroot,nroot)
C     subroutines
c        fun4
c        root
c        c05ade
c        c05aze
c        p01aae
c        x04aae
c        x02aae
c        x02abe
c     and function
c        froot
c     are used here.
      external froot
      parameter (nrot=600,nr2=2)
      dimension y(n),x(n),xroot(nr2),wy(nrot),wx(nrot),ind(n),xwr(nrot)
      common /rots/wy,wx,iym,nnn
      nnn=n
      call zero(xroot,nr2)
      call zero(xwr,nrot)
      do i=1,n
         ind(i)=0
         wx(i)=x(i)
         wy(i)=y(i)-y0
c         write(*,*)i,wx(i),wy(i)
      enddo
      nroot=0
      i=0
 1    i=i+1
 2    if(i.ge.n) then
         if(i.eq.n) then
            if(wy(i).eq.0.) then
               nroot=nroot+1
               xwr(nroot)=x(i)
               ind(nroot)=i
            endif
         endif
         call equiv(xroot(1),xwr(1),min0(nroot,nr2))
         return
      endif
      if(wy(i).eq.0.) then
         nroot=nroot+1
         xwr(nroot)=x(i)
         ind(nroot)=i
      goto 1
      endif
      if(wy(i+1).eq.0.) then
         nroot=nroot+1
         xwr(nroot)=x(i+1)
         ind(nroot)=i+1
         i=i+2
         goto 2
      endif
      if(sign(1.,wy(i)).eq.sign(1.,wy(i+1))) then
         goto 1
      else
         nroot=nroot+1
         iym=max0(i-1,1)
         iym=min0(i-1,n-3)
         call root(x(i),x(i+1),xeps,feta,xwr(nroot),froot)
         ind(nroot)=i
         goto 1
      endif
      call equiv(xroot(1),xwr(1),min0(nroot,nr2))
      return
      end
c**************************************************************************
      FUNCTION froot(xr)
      parameter (nrot=600)
      dimension wy(nrot),wx(nrot)
      common /rots/wy,wx,iym,nnn
c      write(*,*) iym,(wx(i),i=iym,4),(wy(i),i=iym,4)
      if(iym.le.0) then
         iym=1
         wx(nnn+2)=wx(nnn)-wx(nnn-1)
         wx(nnn+1)=wx(nnn)+wx(nnn+2)
         wy(nnn+1)=wy(nnn)
         wx(nnn+2)=wx(nnn+1)+wx(nnn+2)
         wy(nnn+2)=wy(nnn)
      endif
      call fun4(wx(iym),wy(iym),xr,f)
      froot=f
      end
c********************************************************************
        SUBROUTINE fun4(x,f,xn,yn)
c Four points Lagrange' formula for function on nonequal grid
        dimension x(4),f(4)
        x1=xn-x(1)
        x2=xn-x(2)
        x3=xn-x(3)
        x4=xn-x(4)
        yn=f(1)/(x(1)-x(3))/(x(1)-x(4))
        hhh=f(2)/(x(2)-x(3))/(x(2)-x(4))
        yn=(yn*x2-hhh*x1)*x3*x4/(x(1)-x(2))
        hhh=f(3)/(x(3)-x(2))/(x(3)-x(4))/(x(3)-x(1))
        yn=yn+hhh*x2*x1*x4
        hhh=f(4)/(x(4)-x(2))/(x(4)-x(3))/(x(4)-x(1))
        yn=yn+hhh*x2*x3*x1
        return
        end
c************************************************************************
	integer function indx(f,n,x)
	dimension f(n)
        integer jd, jsmesch
	indx=(x-f(1))*real(n)/(f(n)-f(1))
	indx=min(indx,n-1)
	indx=max(1,indx)
        jd=sign(1.,f(n)-f(1))
        jsmesch=(jd-1)/2
 1	continue
	if(x.lt.f(indx-jsmesch).and.indx-jsmesch.gt.1.and.indx.le.n-1
     &       +jsmesch) then
	   indx=indx-jd
	   go to 1
	endif
 2	continue
	if(x.gt.f(indx+1+jsmesch).and.indx+jsmesch.ge.1.and.indx.lt.n-1
     &       -jsmesch) then
	   indx=indx+jd
	   go to 2
	endif
	return
	end
C************************************************************************
        SUBROUTINE der4(x,f,xn,dyn)
c Four points Lagrange' formula for derivatives on nonequal grid
        dimension x(4),f(4)
        x1=xn-x(1)
        x2=xn-x(2)
        x3=xn-x(3)
        x4=xn-x(4)
        yn=f(1)/(x(1)-x(3))/(x(1)-x(4))
        dyn=yn*(x3*x4+x2*x3+x2*x4)
        hhh=f(2)/(x(2)-x(3))/(x(2)-x(4))
        dyn=(dyn-hhh*(x3*x4+x1*x3+x1*x4))/(x(1)-x(2))
        hhh=f(3)/(x(3)-x(2))/(x(3)-x(4))/(x(3)-x(1))
        dyn=dyn+hhh*(x2*x4+x1*x2+x1*x4)
        hhh=f(4)/(x(4)-x(2))/(x(4)-x(3))/(x(4)-x(1))
        dyn=dyn+hhh*(x2*x3+x1*x2+x1*x3)
        return
        end
c************************************************************************
      function wdcchn(i,j,l)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      wdcchn=wdcchng(i,j,l)
      return
      end
c*********************************************************
      function ppchn(i)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      ppchn=ppsi(i)
      return
      end
c**********************************************************************        
      function wdchn(i,j,l)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      wdchn=wdchng(i,j,l)
      return
      end
c*********************************************************
	SUBROUTINE ASIMP(x,y,nns,res)
c*** this subroutine for calculation of integral by the parabola with
c*** nonequal step.
c*** input parameters: x(nns)-massive of variables:
c***                   y(nns)-massive of function;
c***                   nns-numbers  of variables.
c*** output parameter: res - result.
	dimension x(nns),y(nns)
	res=0.
	if (nns.eq.2) then
	   res=res+(y(nns-1)+y(nns))*(x(nns)-x(nns-1))/2.
	else
	   do i=1,nns-2,2
c	      res=res+(y(i)+y(i+1))*(x(i+1)-x(i))/2.+
c	1	   (y(i+2)+y(i+1))*(x(i+2)-x(i+1))/2.
	    res1=(x(i+2)-x(i))*(y(i)+y(i+1)+y(i+2)+(y(i+1)-y(i+2))*
	1	   (x(i+1)-x(i))/(2.*(x(i+2)-x(i+1)))+(y(i+1)-y(i)
	2	   )*(x(i+2)-x(i+1))/(2.*(x(i+1)-x(i))))/3.
	    res=res+res1
	   enddo
	   rm=mod(real(nns-1),2.)
	   if(rm.lt.0.1) then
	      return
	   else
	      res=res+(y(nns-1)+y(nns))*(x(nns)-x(nns-1))/2.
	   endif
	endif
 	return
	end
C******************************************************
      subroutine abij(rii1,rj1,zii1,zj1,aint1,ajnt,bjnt,sqrajnt
     $     ,asqrajnt)
      dimension rii1(2),zii1(2)
      ajnt=(zii1(1)-zii1(2))/(rii1(1)-rii1(2))
      aint1=(rii1(1)-rj1)/(zii1(1)-zj1)
      bjnt=zii1(1)-rii1(1)*ajnt
      sqrajnt=sqrt(ajnt**2+1)
      asqrajnt=sqrajnt*ajnt
      return
      end
c***********************************************************************
      subroutine derivative(y,x,n,x0,dy,wtae)
      dimension y(n),x(n)
      data epsv/1.e-3/
      dy=0.
      if(n.gt.3) then
c         write(*,*) x,n,x0
         i0=max(indx(x,n,x0)-1,1)
         i0=min(i0,n-3)
         call der4(x(i0),y(i0),x0,dy)
c         dy=abs(dy)
         if(abs(dy).lt.abs(wtae/x0*epsv)) dy=abs(wtae/x0*epsv)
         return
      endif
      if(n.gt.2) then
         dy=(y(2)-y(1))/(x(2)-x(1))
         dy2=(y(3)-y(2))/(x(3)-x(2))
         dy=dy+(x0-x(1))*(dy2-dy)/(x(3)-x(1))
c         dy=abs(dy)
         if(abs(dy).lt.abs(wtae/x0*epsv)) dy=abs(wtae/x0*epsv)
         return
      endif
      if(n.eq.2) then
         dy=(y(2)-y(1))/(x(2)-x(1))
c         dy=abs(dy)
         if(abs(dy).lt.abs(wtae/x0*epsv)) dy=abs(wtae/x0*epsv)
         return
      endif
      return
      end
c***********************************************************************

      subroutine weightk(aii1,aj,aj1,bj,bj1,sqrajnt,sqraj1nt,asqrajnt
     $     ,asqraj1nt,rii1,rjj1,zii1,zjj1,r1,z1,ci,cj)
      dimension aii1(2),rii1(2),rjj1(2),zjj1(2),zii1(2),ci(2),cj(2)
      cj(1)=abs(z1-r1*aj-bj)
      cj(2)=abs(z1-r1*aj1-bj1)
      ax=(asqrajnt*cj(2)+asqraj1nt*cj(1))/(cj(1)*sqraj1nt+cj(2)*sqrajnt)
      haxi=ax*aii1(1)
      haxi1=ax*aii1(2)
      h1axi=1.-haxi
      h1axi1=1.-haxi1
      cj(1)=z1-zii1(1)
      cj(2)=z1-zii1(2)
      h3=rii1(2)-r1
      h4=rii1(1)-r1
c      ci(2)=(cj(2)*aii1(2)+rii1(2)-haxi1*r1)/h1axi1-(cj(1)*aii1(1)
c     $     +rii1(1)-haxi*r1)/h1axi
      ci(1)=(cj(2)*aii1(2)+h3)/h1axi1
c     $     /ci(2)
      ci(2)=-(cj(1)*aii1(1)+h4)/h1axi
c     $     /ci(2)
      h1=zjj1(1)-zii1(1)
      h2=zjj1(2)-zii1(2)
      cj(1)=abs((cj(1)+ax*h4)/h1axi/h1*ci(1))
      cj(2)=abs((cj(2)+ax*h3)/h1axi1/h2*ci(2))
      ci(1)=abs(ci(1)*(zjj1(1)-z1-ax*(rjj1(1)-r1))/h1axi/h1)
      ci(2)=abs(ci(2)*(zjj1(2)-z1-ax*(rjj1(2)-r1))/h1axi1/h2)
      ax=ci(1)+ci(2)+cj(1)+cj(2)
c      return
      ci(1)=ci(1)/ax
      ci(2)=ci(2)/ax
      cj(1)=cj(1)/ax
      cj(2)=cj(2)/ax
      return
      end
c***********************************************************************
      subroutine ij(i,j,r1,z1,r0o,z0,kfast)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      data epsi/.01/
      if(abs(r0o-r1)+abs(z1-z0).lt.sqrpsi0) then
         i=1
c         j=1
         hh2=(r1-r0o)
         hh22=hh2**2
         hhz2=(z1-z0)**2
         hh=max((hhz2+hh22),1.e-6)
         if(hh2.lt.0.) then
            hh=amax1(hh2/sqrt(hh),-1.)
         else
            hh=amin1(hh2/sqrt(hh),1.)
         endif
         if(z1-z0.lt.0) then
            j=1+acos(hh)*mth/6.2831853
         else
            j=1+(1-acos(hh)/6.2831853)*mth
         endif
         j=max(1,j)
         j=min(j,mth)
         return
      endif
      if(kfast.eq.0.and.i.gt.5) goto 10
      i=0
      if(z1.gt.zinf(jzmax).or.z1.lt.zinf(jzmin).or.r1.lt.xinf(jrmin)
     $     .or.r1.gt.xinf(jrmax)) return
      hh2=(r1-r0o)
      hh22=hh2**2
      hhz2=(z1-z0)**2
      hh=(hhz2+hh22)
      if(hh2.lt.0.) then
         hh=amax1(hh2/sqrt(hh),-1.)
      else
         hh=amin1(hh2/sqrt(hh),1.)
      endif
      if(z1-z0.lt.0) then
         j=1+acos(hh)*mth/6.2831853
      else
         j=1+(1-acos(hh)/6.2831853)*mth
      endif
      j=max(1,j)
      j=min(j,mth)
      hh2=(zmz02*hh22+rmr02*hhz2)/zmrm2
      hh=sqrt(hh2)
      i=nosurf*hh
      i=max(1,i)
      i=min(i,nosurf-1)
c     Here i,j are initial guess for the cell where r1,z1 are supposed to lay
c      write(*,*) '---------------',i,j
 10   continue
      if(j.lt.jzmin.and.z1.gt.z0) j=mth
      do iii=1,2
         jd=1
 11      continue
         if((r1-rr(i,j))*(zz(i+1,j)-zz(i,j)).gt.(z1-zz(i,j))*
     $        (rr(i+1,j)-rr(i,j)).and.j.lt.mth-1)then
            j=j+jd
            goto11
         endif
 1       continue
         if((r1-rr(i,j))*(zz(i,j+1)-zz(i,j)).lt.(z1-zz(i,j))*
     $        (rr(i,j+1)-rr(i,j)).and.i.lt.nosurf-1)then
            i=i+jd
            goto1
         endif
         jd=-1
 2       continue
         if((r1-rr(i,j))*(zz(i+1,j)-zz(i,j)).lt.(z1-zz(i,j))*
     $        (rr(i+1,j)-rr(i,j)).and.j.gt.1)then
            j=j+jd
            goto2
         endif
 3       continue
         if((r1-rr(i,j))*(zz(i,j+1)-zz(i,j)).gt.(z1-zz(i,j))*
     $        (rr(i,j+1)-rr(i,j)).and.i.gt.1)then
            i=i+jd
            goto3
         endif
      enddo
      i=min(i,nosurf-1)
      j=max(1,j)
      j=min(j,mth-1)
      if(i.ge.nosurf-1) then
         if((r1-rr(i+1,j))*(zz(i+1,j+1)-zz(i+1,j)).lt.(z1-zz(i+1,j))*
     $        (rr(i+1,j+1)-rr(i+1,j))) then
            i=0
            j=0
         endif
      endif
c      write(*,*)'__________',i,j
c     write(*,*) '_______________',i,j
      return
       end
c********************************************************************
      subroutine weight(ru,rl,zu,zl,r,z,cu,cl)
      dimension ru(2),rl(2),zu(2),zl(2),cu(2),cl(2)
 13   err=1.e-9
      do i=1,2
         cl(i)=0.
         cu(i)=0.
      enddo
      if((ru(1)-rl(1))**2+(zu(1)-zl(1))**2.ge.(rl(2)-rl(1))**2+(zl(2)
     &     -zl(1))**2) then
         klrotat=0
         if(abs(ru(1)-rl(1)).lt.abs(zu(1)-zl(1))) then
            z1=zl(1)
            r1=rl(1)
            z2=zl(2)
            r2=rl(2)
            z3=zu(1)
            r3=ru(1)
            z4=zu(2)
            r4=ru(2)
            rx=r
            zx=z
         else
            r1=zl(1)
            z1=rl(1)
            r2=zl(2)
            z2=rl(2)
            r3=zu(1)
            z3=ru(1)
            r4=zu(2)
            z4=ru(2)
            rx=z
            zx=r
         endif
      else
         klrotat=1
         if(abs(rl(2)-rl(1)).lt.abs(zl(2)-zl(1))) then
            z1=zl(2)
            r1=rl(2)
            z2=zu(2)
            r2=ru(2)
            z3=zl(1)
            r3=rl(1)
            z4=zu(1)
            r4=ru(1)
            rx=r
            zx=z
         else
            r1=zl(2)
            z1=rl(2)
            r2=zu(2)
            z2=ru(2)
            r3=zl(1)
            z3=rl(1)
            r4=zu(1)
            z4=ru(1)
            zx=r
            rx=z
         endif
      endif
      xl=sqrt((z1-zx)**2+(r1-rx)**2)
      xu=sqrt((r3-rx)**2+(z3-zx)**2)
      h=sqrt((r1-r3)**2+(z1-z3)**2)
      xlu=xl+xu-h
      if(xlu.lt.err) xlu=0.
      xlh=xl-xu+h
      if(xlh.lt.err) xlh=0.
      xuh=xu-xl+h
      if(xuh.lt.err) xuh=0.
      h13=sqrt(xlh*xlu*xuh*(xl+xu+h))/2./h
      xl=sqrt((z2-zx)**2+(r2-rx)**2)
      xu=sqrt((r4-rx)**2+(z4-zx)**2)
      h=sqrt((r2-r4)**2+(z2-z4)**2)
      xlu=xl+xu-h
      if(xlu.lt.err) xlu=0.
      xlh=xl-xu+h
      if(xlh.lt.err) xlh=0.
      xuh=xu-xl+h
      if(xuh.lt.err) xuh=0.
      h24=sqrt(xlh*xlu*xuh*(xl+xu+h))/2./h
      if(h13+h24.le.0.)then
         if(err.gt.1.e-10)then
            err=1.e-15
            goto 13
         endif
         write(*,*)'h13+h24=',h13+h24
     &        ,'Warn. Check the tolerance in programm weight'
         write(*,*) 'ru,zu,rl,zl,r,z=',ru,zu,rl,zl,r,z
         return
         cu(1)=1.
c         stop
      endif
      zll=h13*(z1-z3)*(r1*(z4-z2)+r2*(z1-z4)+r4*(z2-z1))
      zuu=h24*(z2-z4)*(r1*(z3-z2)+r2*(z1-z3)+r3*(z2-z1))
      zll=(zll*z2+zuu*z1)/(zll+zuu)
      zlll=h13*(z1-z3)*(r2*(z3-z4)+r3*(z4-z2)+r4*(z2-z3))
      zuu=h24*(z2-z4)*(r1*(z3-z4)+r3*(z4-z1)+r4*(z1-z3))
      zuu=(zlll*z4+zuu*z3)/(zlll+zuu)
      p=(zx-zll)/(zuu-zll)
      ax=((r1-r3)/(z1-z3)*h24+(r2-r4)/(z2-z4)*h13)/(h13+h24)
      pl=(ax*(zll-zx)+rx-r1)/(r2-r1)
      pu=(ax*(zuu-zx)+rx-r3)/(r4-r3)
      if(klrotat.eq.0) then
         cl(1)=(1.-pl)*(1.-p)
         cl(2)=pl*(1.-p)
         cu(1)=p*(1.-pu)
         cu(2)=p*pu
      else
         cl(2)=(1.-pl)*(1.-p) 
         cl(1)=p*(1.-pu)
         cu(1)=p*pu
         cu(2)=pl*(1.-p)
      endif
      if(k_pest(1.).ne.-1)return
      if(cl(1).lt.0.)then
         cl(2)=cl(2)+cl(1)/3.
         cu(1)=cu(1)+cl(1)/3.
         cu(2)=cu(2)+cl(1)/3.
         cl(1)=0.
         if(cl(2).lt.0.)then
            cu(1)=cu(1)+cl(2)/2.
            cu(2)=cu(2)+cl(2)/2.
            cl(2)=0.
            if(cu(2).lt.0.) then
               cu(1)=1.
               cu(2)=0.
            else if(cu(1).lt.0.) then
               cu(2)=1.
               cu(1)=0.
            endif
         endif
         return
      else if(cl(2).lt.0.) then
         cl(1)=cl(1)+cl(2)/3.
         cu(1)=cu(1)+cl(2)/3.
         cu(2)=cu(2)+cl(2)/3.
         cl(2)=0.
         if(cl(1).lt.0.)then
            cu(1)=cu(1)+cl(1)/2.
            cu(2)=cu(2)+cl(1)/2.
            cl(1)=0.
            if(cu(1).lt.0.) then
               cu(2)=1.
               cu(1)=0.
            else if(cu(2).lt.0.) then
               cu(1)=1.
               cu(2)=0.
            endif
         endif
         return
      else if(cu(2).lt.0.) then
         cu(1)=cu(1)+cu(2)/3.
         cl(1)=cl(1)+cu(2)/3.
         cl(2)=cl(2)+cu(2)/3.
         cu(2)=0.
         if(cu(1).lt.0.) then
            cl(1)=cl(1)+cu(1)/2.
            cl(2)=cl(2)+cu(1)/2.
            cu(1)=0.
            if(cl(1).lt.0.)then
               cl(2)=1.
               cl(1)=0.
            else if(cl(2).lt.0.) then
               cl(1)=1.
               cl(2)=0.
            endif
         endif
         return
      else if(cu(1).lt.0.) then
         cu(2)=cu(2)+cu(1)/3.
         cl(1)=cl(1)+cu(1)/3.
         cl(2)=cl(2)+cu(1)/3.
         cu(1)=0.
         if(cu(2).lt.0.) then
            cl(1)=cl(1)+cu(2)/2.
            cl(2)=cl(2)+cu(2)/2.
            cu(2)=0.
            if(cl(1).lt.0.)then
               cl(2)=1.
               cl(1)=0.
            else if(cl(2).lt.0.) then
               cl(1)=1.
               cl(2)=0.
            endif
         endif
      endif
      return
      end
c****************************************************************
c function for bounce averaging
      subroutine fwd(ff,iav,r,z)
      common/wdr/alpha,b0ax,psedge,vioo,sigvpr
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,voo,z0,xhp
      common/fast/xav,xavint,x2r,x2i,psioo,kfast
      dimension ff(iav)
      call wdoo(alpha,b0ax,pmu,psedge,qpf,r,
     $       rvpar,sigvpr,z,z0,vioo,voo,ff(1),ff(2),ff(3))
      if(iav.gt.3) ff(4)=psioo
      return
      end
c*********************************************************
c function for bounce averaging
      subroutine fwdf(ff,iav,r,z)
      common/wdr/alpha,b0ax,psedge,vioo,sigvpr
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,voo,z0,xhp
      common/fast/xav,xavint,x2r,x2i,psioo,kfast
      dimension ff(iav)
      ff(iav)=psioo
      call wdoof(alpha,b0ax,pmu,psedge,qpf,r,
     $       rvpar,sigvpr,z,z0,vioo,voo,ff(1),ff(2),ff(3))
      return
      end
c*********************************************************
      subroutine minorrad(htb)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      common/saveij/i,j,cl(2),cu(2)
      htb=f2dperp(grthsq2d(i,j),grthsq2d(i,j+1),cl,cu)
      if(htb.le.0.)write(*,*)'Warning in minorrad',i,j,cl,cu,grthsq2d(i
     &     ,j)
      htb=1./sqrt(abs(htb))
      return
      end
c*********************************************************
      subroutine summation(a,b,ntgr1,f,pi,eps,epsp,iter,qpf,pmu,r0
     $     ,rvpar,sig,v,z0,xhp,iav,t0,tb,fo)
      parameter(ntgr2=30)
      common/fowint/wdintfow(ntgr2,2),rintfow(ntgr2,2),tbintfow(ntgr2,2)
     &     ,klyu
      common /fast/xav,xavint,x2r,x2i,psioo,kfast
      common/saveij/ik,jk,cl(2),cu(2)
      external f
      save ntgr0,rf
      dimension ff(20),fo(iav),rf(150)
      ntgr=max(ntgr1,8)
      if(ntgr0.ne.ntgr) then
         ntgr0=ntgr
         do i=1,ntgr
            rf(i)=cos((real(2*i)-1)*pi/real(2*ntgr))
         enddo
      endif
c      ntgr=abs(a-b)/abs(a+b)*4*ntgr1
      tb=0.
c      write(*,*) 'a,b,rvpar= ',a,b,rvpar
      do i=1,ntgr
c integration is from b to a
         r=.5*(a+b+(b-a)*rf(i))
c         write(*,*)'============================',i
         if(i.le.2)then
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,r,r0
     $           ,rvpar,sig,v,z,z0,xhp)
         else
            call psorbf(bm,b0,eps,iter,qpf,pmu,psior,r,r0
     $           ,rvpar,sig,v,z,z0,xhp,zh)
         endif
         if(k_pest(psior).eq.1.or.(ik.eq.1.and.k_pest(psior).eq.-1)
     &        .or.psior.ge.1.)then
            tb=0.
            return
         endif
         zh=z
         kfast=1
         psioo=psior
cnng others averaged integrands are computed inside f, which is fwd here
         call f(ff(1),iav,r,z)
c Note that here dzpsi is R*dPsi/dZ and dpsi= R*dPsi/dR
         htb=psiof(dpsi,ddpsi,dzpsi,r,r0,z,z0,am,1)
         if(r0*b0.lt..01) goto 11
         pss=1.-pmu*bm/r0/b0
c         write(*,*) pss,pmu,bm,r0,b0
         if(pss.gt.5.*r*epsp) then
            htb1=min(abs(rvpar-a),abs(rvpar-b))
            htb=min(abs(r-a),abs(b-r))
            if(htb.lt.abs(a-b)*.063.and.(htb1.gt.1.e-6*abs(a-b).or.max(a
     &           ,b)-r.lt.r-min(a,b)))then
               call minorrad(htb)
               htb=sqrt(pss*2./htb*(dzpsi**2+dpsi**2))
     &              /sqrt(max(abs(r-a),abs(b-r)))
c
c               write(*,*) '1htb,pss,dzpsi,r,a,b,rvpar,htb1  ',htb,pss
c     &              ,dzpsi,r,a,b,rvpar,htb1
               htb=max(abs(htb),.003)
               htb=r**2/htb*bm
            else
               htb=dzpsi*sqrt(pss)/sqrt((r-a)*(b-r))
c
c               write(*,*) '2htb,pss,dzpsi,r,a,b,rvpar,htb1  ',htb,pss
c     &              ,dzpsi,r,a,b,rvpar,htb1
               htb=max(abs(htb),.003)
               htb=r**2/htb*bm
            endif
         else
            htb=min(abs(rvpar-a),abs(rvpar-b))
            if(htb.lt.1.e-6*abs(a-b))then
               if(abs(rvpar-a).lt.(abs(rvpar-b)))then
                  htb=abs(r-b)
               else
                  htb=abs(r-a)
               endif
               htb=dzpsi/sqrt(htb*r)
c
c               write(*,*) '3htb,pss,dzpsi,r,a,b  ',htb,pss,dzpsi,r,a
c     &              ,b
               htb=max(abs(htb),0.003)
               htb=bm/htb*r**2
            else
               htb=min(abs(r-a),abs(b-r))
               if(htb.lt.abs(a-b)*.063)then
                  call minorrad(htb)
                  htb=sqrt((1-rvpar/r)*2./htb*(dzpsi**2+dpsi**2))
     &                 /sqrt(max(abs(r-a),abs(b-r)))
c
c               write(*,*) '4htb,pss,dzpsi,r,a,b  ',htb,pss,dzpsi,r,a
c     &              ,b
                  htb=max(abs(htb),.003)
                  htb=r**2/htb*bm
               else
                  htb=dzpsi*sqrt(abs((r-rvpar)/(b-r)/(r-a)/r))
c
c                  write(*,*) '5htb,pss,dzpsi,r,a,b  ',htb,pss,dzpsi
c     &                 ,r,a,b
                  htb=max(abs(htb),0.003)
                  htb=r**2/htb*bm
               endif
            endif
         endif            
         do ii=1,iav
            fo(ii)=fo(ii)+ff(ii)*htb
         enddo
         xav=xav+xavint*htb
         x2r=x2r+x2i*htb
c         write(*,*) 'fo(3),ff(3),r=',fo(3),ff(3),r
c         write(*,*) 'fo(1),ff(1),r=',fo(1),ff(1),r
         tb=tb+htb
         if(klyu.eq.1)then
            if(sig.lt.0.)then
c               write(*,*)a,b,rvpar,r,htb
               wdintfow(i,1)=ff(3)*htb
               rintfow(i,1)=r
               tbintfow(i,1)=htb
            else
               wdintfow(i,2)=ff(3)*htb
               rintfow(i,2)=r
               tbintfow(i,2)=htb
            endif
         endif
 11      continue
      enddo
c
c      write(*,*) '----------------------------',tb
      if(b0.lt..01) b0=1.
      htb=pi/real(ntgr)/b0/v
c
c      write(*,*) htb,tb,fo
      if(tb.eq.0.) then
         do ii=1,iav
            fo(ii)=0.
         enddo
         return
      endif
      do ii=1,iav
         fo(ii)=fo(ii)/tb
      enddo
      xav=xav/tb
      x2r=x2r/tb
      tb=tb*htb*t0
c      write(*,*) 'wb= ',2.*pi/tb,fo(2)
      return
      end
c**********************************************************************
c This subr. gives the instant frequency to be averaged for the 
c tabulation of the resonance condition frequencies
c WKPRVPR is the frequency which stands in front of poloidal mode number M
c WQKPR   is the -- // -- N toroidal number and is WKPRVPR*Q
c WD1     is the only drift in front of N
        subroutine wdoo(alpha,b0ax,pmu,psedge,qpf,r1oo,
     $     rvpar,sigvpr,z1oo,z0,vioo,voo,wkprvpr
     &     ,wqkpr,wd1)
        include 'clich1'
        include 'clich2'
        include 'orb2d.par'
        common/asymm/asym /fast/xav,xavint,x2r,x2i,psioo,kfast
        dimension cl(2),cu(2)
        common/saveij/i,j,cl,cu
        if(abs(asym).lt..1) asym=1.
        roo=r1oo
        zoo=z0+abs(z1oo-z0)*asym
c        zoo=z1oo
        if(kfast.gt.0.and.i.gt.0) goto 10
        call ij(i,j,roo,zoo,rax,z0,kfast)
c        write(*,*) '1 ',zoo,z0,z1oo,i,j,roo
        if(i.le.0) then
           wd1=0.
           wkprvpr=0.
           wqkpr=0.
           xavint=0.
           x2i=0.
           return
        endif
        if(i.le.1) then
           zoo=zz(i+1,j)
           roo=rr(i+1,j)
        endif
      z1=zoo
      r1=roo
      if(i.eq.1) then
         if(sqrt((r1-rax)**2+(z1-z0)**2).lt.1.05*abs(rr(1,1)-rr(1,mth/2)
     &        )/2)then
            if(cl(1).lt.0..or.cl(2).lt.0..or.cu(1).lt.0..or.cu(2).lt.0.
     &        )then
c               write(*,*)'2cl,cu,r1,z1,ru,zu,rl,zl=',cl,cu,r1,z1,rr(i,j)
c     &              ,rr(i+1,j),zz(i,j),zz(i+1,j),rr(i,j+1),rr(i+1,j+1)
c     &              ,zz(i,j+1),zz(i+1,j+1)
               cl(1)=0.5
               cu(1)=0.5
               cl(2)=0.
               cu(2)=0.
            endif
         else
            call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      else
         jc1=(j+jd0)/jdj
         if(jc1.eq.1.or.jc1.eq.3) then
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),zz(i,j),zz(i,j+1),rr(i,j),rr(i,j+1),z1,r1,cl,cu)
         else
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      endif
call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),roo,zoo,cl,cu)
 10   psioo=f2dperp(ppsi(i),ppsi(i),cl,cu)
      qqq=f2dperp(qoo(i),qoo(i),cl,cu)
      bcur=f2dperp(b2d(i,j),b2d(i,j+1),cl,cu)
      xjcb=f2dperp(xjacob2d(i,j),xjacob2d(i,j+1),cl,cu)
      wkprvpr=pmu*bcur/rax/b0ax
      wqkpr=(1.-0.5*wkprvpr)
      wd1=f2dperp(wdo(i,j,2),wdo(i,j+1,2),cl,cu)*wqkpr
      wqkpr=f2dperp(wdo(i,j,1),wdo(i,j+1,1),cl,cu)
     $     *wqkpr
      wd1=wd1-(f2dperp(wdo(i,j,4),wdo(i,j+1,4),cl,cu))*pmu
      wqkpr=wqkpr-f2dperp(wdo(i,j,3),wdo(i,j+1,3),
     $     cl,cu)*pmu
c Note that for well trapped particles WD1 and WQKPR must have the same sign
      wd1=wd1*voo**2*alpha
      wqkpr=wqkpr*voo**2*alpha
      wkprvpr=(1.-wkprvpr)
      if(wkprvpr.lt.0.) wkprvpr=abs(roo-rvpar)/roo
      wkprvpr=sqrt(wkprvpr)*vioo*voo/xjcb/bcur*10.
c Note that sigvpr is different from usual definition
c Here it is positive if particle is copassing to the plasma current
      if(sigvpr.gt.0.) wkprvpr=-wkprvpr
c      if(psioo.lt.qpf) wkprvpr=-wkprvpr
c For well trapped VPR < 0 and WQKPR < 0
      wkprvpr=wkprvpr-wqkpr
      wqkpr=wkprvpr*qqq
      xavint=-100.
      x2i=-100.
c
      return
      end
c***********************fast version************************************
c This subr. gives the instant frequency to be averaged for the 
c tabulation of the resonance condition frequencies
c It similar to the WDOO except it is faster and requires certain more values
c WKPRVPR is the frequency which stands in front of poloidal mode number M
c WQKPR   is the -- // -- N toroidal number and is WKPRVPR*Q
c WD1     is the only drift in front of N
      subroutine wdoof(alpha,b0ax,pmu,psedge,qpf
     &     ,r1oo,rvpar,sigvpr,z1oo,z0,vioo,voo
     &     ,wkprvpr,wqkpr,wd1)
        include 'clich1'
        include 'clich2'
        include 'orb2d.par'
        common/asymm/asym /fast/xav,xavint,x2r,x2i,psioo,kfast
        dimension cl(2),cu(2)
        common/saveij/i,j,cl,cu
        if(abs(asym).lt..1) asym=1.
        roo=r1oo
        zoo=z0+abs(z1oo-z0)*asym
c        zoo=z1oo
        if(kfast.gt.0.and.i.gt.0) goto 10
        call ij(i,j,roo,zoo,rax,z0,kfast)
c        write(*,*) '1 ',zoo,z0,z1oo,i,j,roo
        if(i.le.0) then
           wd1=0.
           wkprvpr=0.
           wqkpr=0.
           xavint=0.
           x2i=0.
           return
        endif
        if(i.le.1) then
           zoo=zz(i+1,j)
           roo=rr(i+1,j)
        endif
      z1=zoo
      r1=roo
      if(i.eq.1) then
         if(sqrt((r1-rax)**2+(z1-z0)**2).lt.1.05*abs(rr(1,1)-rr(1,mth/2)
     &        )/2)then
            if(cl(1).lt.0..or.cl(2).lt.0..or.cu(1).lt.0..or.cu(2).lt.0.
     &        )then
c               write(*,*)'3cl,cu,r1,z1,ru,zu,rl,zl=',cl,cu,r1,z1,rr(i,j)
c     &              ,rr(i+1,j),zz(i,j),zz(i+1,j),rr(i,j+1),rr(i+1,j+1)
c     &              ,zz(i,j+1),zz(i+1,j+1)
               cl(1)=0.5
               cu(1)=0.5
               cl(2)=0.
               cu(2)=0.
            endif
         else
            call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      else
         jc1=(j+jd0)/jdj
         if(jc1.eq.1.or.jc1.eq.3) then
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),zz(i,j),zz(i,j+1),rr(i,j),rr(i,j+1),z1,r1,cl,cu)
         else
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      endif
call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),roo,zoo,cl,cu)
 10   continue
cnng
      psioo=f2dperp(ppsi(i),ppsi(i),cl,cu)
      qqq=f2dperp(qoo(i),qoo(i),cl,cu)
      bcur=f2dperp(b2d(i,j),b2d(i,j+1),cl,cu)
      xjcb=f2dperp(xjacob2d(i,j),xjacob2d(i,j+1),cl,cu)
      wkprvpr=pmu*bcur/rax/b0ax
      wqkpr=(1.-0.5*wkprvpr)
      wd1=f2dperp(wdo(i,j,2),wdo(i,j+1,2),cl,cu)*wqkpr
      wqkpr=f2dperp(wdo(i,j,1),wdo(i,j+1,1),cl,cu)
     $     *wqkpr
      wd1=wd1-(f2dperp(wdo(i,j,4),wdo(i,j+1,4),cl,cu))*pmu
      wqkpr=wqkpr-f2dperp(wdo(i,j,3),wdo(i,j+1,3),
     $     cl,cu)*pmu
      wd1=wd1*voo**2*alpha
      wqkpr=wqkpr*voo**2*alpha
c this is v||/v **2
      wkprvpr=amax1(1.-wkprvpr,0.)
      x2i=sqrt(wkprvpr)
      bfof=f2dperp(rg2d(i),rg2d(i),cl,cu)
c this is v_varphi*R/v **2
      xavint=wkprvpr*(bfof/bcur)**2
      xavint=roo**2-xavint
      x2i=x2i*bfof/bcur
c
      if(wkprvpr.lt.0.) wkprvpr=abs(roo-rvpar)/roo
      wkprvpr=sqrt(wkprvpr)*vioo*voo/xjcb/bcur*10.
c Note that sigvpr is different from usual definition
c Here it is positive if particle is copassing to the plasma current
      if(sigvpr.gt.0.) wkprvpr=-wkprvpr
c      if(psioo.lt.qpf) wkprvpr=-wkprvpr
      wkprvpr=wkprvpr-wqkpr
      wqkpr=wkprvpr*qqq
      return
      end
c***********************************************************************
      function psio(dpsi1,ddpsi1,dzpsi1,r1,r0o,z11,z0,a1)
c     psio is the magnetic flux in the equatorial plane
c     dpsi is output parameter (logariphm derivative of psi)
c     ddpsi is output parameter (logariphm second dirivative of psi)
c     dzpsi is output parameter (r*dirivative of psi on z)
c     r,z are input parameters
      psio=psiof(dpsi1,ddpsi1,dzpsi1,r1,r0o,z11,z0,a1,-1)
      return
      end
c*********************************************************
      function psiof(dpsi1,ddpsi1,dzpsi1,r1,r0o,z11,z0,a1,kfast)
c     dpsi is output parameter  R*d Psi / dR
c     ddpsi is output parameter (logariphm second derivative of psi)
c     dzpsi is output parameter R*d Psi / dZ
c     r,z are input parameters
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      common/asymm/asym
      dimension xx1(4),yy1(4,2),cl(2),cu(2)
      common/saveij/i,j,cl,cu
      if(abs(asym).lt..1) asym=1.
      a1=(xinf(1)-xinf(mth/2+1))/2
      z1=z0+abs(z11-z0)*asym
      if(kfast.gt.0.and.i.gt.0)goto 1
      call ij(i,j,r1,z1,r0o,z0,kfast)
      if(i.le.0) then
         psiof=1.
         dpsi1=0.
         ddpsi1=0.
         dzpsi1=0.
         return
      endif
      if(i.eq.1) then
         if(sqrt((r1-rax)**2+(z1-z0)**2).lt.1.05*abs(rr(1,1)-rr(1,mth/2)
     &        )/2)then
            if(cl(1).lt.0..or.cl(2).lt.0..or.cu(1).lt.0..or.cu(2).lt.0.
     &        )then
c               write(*,*)'4cl,cu,r1,z1,z11,ru,zu,rl,zl=',i,j,kfast,cl,cu
c     &              ,r1,z1,z11,rr(i,j),rr(i+1,j),zz(i,j),zz(i+1,j),rr(i
c     &              ,j+1),rr(i+1,j+1),zz(i,j+1),zz(i+1,j+1)
               cl(1)=0.5
               cu(1)=0.5
               cl(2)=0.
               cu(2)=0.
            endif
         else
            call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      else
         jc1=(j+jd0)/jdj
         if(jc1.eq.1.or.jc1.eq.3) then
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),zz(i,j),zz(i,j+1),rr(i,j),rr(i,j+1),z1,r1,cl,cu)
         else
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      endif
 1    continue
      psiof=
     c     f2dperp(ppsi(i),ppsi(i),cl,cu)
      dpsi1=
     c     f2dperp(dpsio(i,j),dpsio(i,j+1),cl,cu)
      if(i.ne.1) then
         if(j.eq.1) then
            dzpsi1=1./dzpsi(i,j+1)
         else
            dzpsi1=1./
     c        f2dperp(dzpsi(i,j),dzpsi(i,j+1),cl,cu)
         endif
      else
care .001 mozhet dat' oshibku okolo z=0
         if(abs(zz(2,j)-z0).gt.1.e-4) then
            dzpsi1=1./dzpsi(2,j)*abs(z1-z0)/abs(zz(2,j)-z0)+.001
         else if(abs(zz(2,j+1)-z0).gt.1.e-4) then
            dzpsi1=1./dzpsi(2,j+1)*abs(z1-z0)/abs(zz(2,j+1)-z0)+.001
        endif
      endif
      if(abs(z0-z1).le.1.e-4) then
         if(r1.lt.r0o) then
            j=mth/2+1
            if(r1.gt.rr(i,j).and.i.gt.1) then
               i=i-1
               if(r1.gt.rr(i,j).and.i.gt.1) i=i-1
            endif
            if(r1.lt.rr(i+1,j).and.i.lt.nosurf-1) then
               i=i+1
               if(r1.lt.rr(i+1,j).and.i.lt.nosurf-1) i=i+1
            endif
         else
            j=1
            if(r1.lt.rr(i,j).and.i.gt.1) then
               i=i-1
               if(r1.lt.rr(i,j).and.i.gt.1) i=i-1
            endif
             if(r1.gt.rr(i+1,j).and.i.lt.nosurf-1) then
               i=i+1
               if(r1.gt.rr(i+1,j).and.i.lt.nosurf-1) i=i+1
            endif
        endif
         if(i.ge.nosurf-1) then
c            psiof=fintr(ppsi(nosurf-3),rr(nosurf-3,j),4,4,r1)
            call funder4(rr(nosurf-3,j),ppsi(nosurf-3),r1,dpsi1
     $           ,psiof)
            call funder4(rr(nosurf-3,j),dpsio(nosurf-3,j),r1,ddpsi1
     $           ,dpsi1)
c            write(*,*)r1,z1,psiof,dpsi1,ddpsi1,j,ppsi(nosurf),
c     $           psi2,dpsi2*r1
c            dpsi1=fintr(dpsio(nosurf-3,j),ppsi(nosurf-3),4,4,psiof)
c            call deriv4(ppsi(nosurf-3),dpsio(nosurf-3,j),psiof,ddpsi1)
         else
            if(i.lt.2) then
               if(j.lt.3) then
                  j1 = mth/2+1
               else
                  j1=1
               endif
               xx1(1)=rr(2,j1)
               yy1(1,1)=-ppsi(2)
               yy1(1,2)=dpsio(2,j1)
               xx1(2)=r0o
               yy1(2,2)=(dpsio(1,j1)+dpsio(1,j))/2.
               yy1(2,1)=0.
               xx1(3)=rr(2,j)
               yy1(3,2)=dpsio(2,j)
               yy1(3,1)=ppsi(2)
               xx1(4)=rr(3,j)
               yy1(4,2)=dpsio(3,j)
               yy1(4,1)=ppsi(3)
c     write(*,*) xx1,yy1
            else
               do iii=1,4
                  xx1(iii)=rr(i+iii-2,j)
                  yy1(iii,1)=ppsi(i+iii-2)
                  yy1(iii,2)=dpsio(i+iii-2,j)
               enddo
            endif
            if(xx1(1).gt.xx1(4))then
               call rev(xx1,4,1)
               call rev(yy1,4,2)
            endif
            psiof=fintr(yy1(1,1),xx1,4,4,r1)
            call funder4(xx1,yy1(1,2),r1,ddpsi1,dpsi1)
         endif
c     dpsi1=fintr(yy1(1,2),yy1(1,1),4,4,psiof)
c     call deriv4(yy1(1,1),yy1(1,2),psiof,ddpsi1)
         ddpsi1=r1*ddpsi1-dpsi1
c         if(i.eq.1) write(*,*) 'tri',psiof,ddpsi1,dpsi1
      else
         ddpsi1=0.
      endif
      return
      end
c******************************************************************
      function psiofpp(r1,r0o,z11,z0,a1,kfast)
c     dpsi is output parameter (logarighmical derivative of psi)
c     ddpsi is output parameter (logarighmical second dirivative of psi)
c     dzpsi is output parameter (r*dirivative of psi on z)
c     r,z are input parameters
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      common/asymm/asym
      dimension xx1(4),yy1(4,2),cl(2),cu(2)
      common/saveij/i,j,cl,cu
      if(abs(asym).lt..1) asym=1.
      a1=(xinf(1)-xinf(mth/2+1))/2
      z1=z0+abs(z11-z0)*asym
      if(kfast.gt.0.and.i.gt.0)goto 1
      call ij(i,j,r1,z1,r0o,z0,kfast)
      if(i.le.0) then
         psiofpp=1.
         return
      endif
      if(i.eq.1) then
         if(sqrt((r1-rax)**2+(z1-z0)**2).lt.1.05*abs(rr(1,1)-rr(1,mth/2)
     &        )/2)then
            if(cl(1).lt.0..or.cl(2).lt.0..or.cu(1).lt.0..or.cu(2).lt.0.
     &        )then
c               write(*,*)'5cl,cu,r1,z1,ru,zu,rl,zl=',cl,cu,r1,z1,rr(i,j)
c     &              ,rr(i+1,j),zz(i,j),zz(i+1,j),rr(i,j+1),rr(i+1,j+1)
c     &              ,zz(i,j+1),zz(i+1,j+1)
               cl(1)=0.5
               cu(1)=0.5
               cl(2)=0.
               cu(2)=0.
            endif
         else
            call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      else
         jc1=(j+jd0)/jdj
         if(jc1.eq.1.or.jc1.eq.3) then
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),zz(i,j),zz(i,j+1),rr(i,j),rr(i,j+1),z1,r1,cl,cu)
         else
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      endif
 1    continue
      psiofpp=f2dperp(ppsi(i),ppsi(i),cl,cu)
      return
      end
c******************************************************************
      function psiofpp2(r1,r0o,z0,a1)
c     dpsi is output parameter (logarighmical derivative of psi)
c     ddpsi is output parameter (logarighmical second dirivative of psi)
c     dzpsi is output parameter (r*dirivative of psi on z)
c     r,z are input parameters
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      dimension xx1(4),yy1(4,2),cl(2),cu(2)
      a1=(xinf(1)-xinf(mth/2+1))/2
      z1=z0
      if(r1.lt.r0o) then
         j=mth/2+1
      else
         j=1
      endif
      i=indx(rr(1,j),nosurf,r1)
      if(i.ge.nosurf-1) then
         call funder4(rr(nosurf-3,j),ppsi(nosurf-3),r1,dpsi1
     $        ,psiofpp2)
      else
         if(i.lt.2) then
            if(j.lt.3) then
               j1 = mth/2+1
            else
               j1=1
            endif
            xx1(1)=rr(2,j1)
            yy1(1,1)=-ppsi(2)
            xx1(2)=r0o
            yy1(2,1)=0.
            xx1(3)=rr(2,j)
            yy1(3,1)=ppsi(2)
            xx1(4)=rr(3,j)
            yy1(4,1)=ppsi(3)
c     write(*,*) xx1,yy1
         else
            do iii=1,4
               xx1(iii)=rr(i+iii-2,j)
               yy1(iii,1)=ppsi(i+iii-2)
            enddo
         endif
         if(xx1(1).gt.xx1(4))then
            call rev(xx1,4,1)
            call rev(yy1(1,1),4,1)
         endif
         psiofpp2=fintr(yy1(1,1),xx1,4,4,r1)
      endif
      return
      end
c******************************************************************
      subroutine intini
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      jdj=mth/4+1
      jd0=mth/8
      do i=2,nosurf
         jc0=1
         do j=1,mth+1
            jc1=(j+jd0)/jdj
            if(jc1.eq.jc0)then
               jc0=jc1+1
               if(jc1.eq.2.or.jc1.eq.4)then
                  call abij(zz(i,j),zz(i,j+1),rr(i,j),rr(i,j+1)
     $                 ,aint1(i,j+jc1-1),ajnt(i,j+jc1-1),bjnt(i,j+jc1
     $                 -1),sqrajnt(i,j+jc1-1),asqrajnt(i,j+jc1-1))
               else
                  call abij(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),
     $                 aint1(i,j+jc1-1),ajnt(i,j+jc1-1),bjnt(i,j+jc1-1)
     $                 ,sqrajnt(i,j+jc1-1),asqrajnt(i,j+jc1-1))
               endif
            endif
            if(jc1.eq.1.or.jc1.eq.3)then
               call abij(zz(i,j),zz(i,j+1),rr(i,j),rr(i,j+1)
     $              ,aint1(i,j+jc1),ajnt(i,j+jc1),bjnt(i,j+jc1),
     $              sqrajnt(i,j+jc1),asqrajnt(i,j+jc1))
            else
               call abij(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),
     $              aint1(i,j+jc1),ajnt(i,j+jc1),bjnt(i,j+jc1),
     $              sqrajnt(i,j+jc1),asqrajnt(i,j+jc1))
            endif
         enddo
      enddo
      zzmin=0.
      zzmax=0.
      rrmax=0.
      rrmin=1000.
      psi0=ppsi(2)
      do j=1,mth
         if(zzmin.gt.zz(nosurf,j)) then
            zzmin=zz(nosurf,j)
            jzmin=j
         endif
         if(zzmax.lt.zz(nosurf,j)) then
            zzmax=zinf(j)
            jzmax=j
         endif
         if(rrmin.gt.rr(nosurf,j)) then
            rrmin=rr(nosurf,j)
            jrmin=j
         endif
        if(rrmax.lt.rr(nosurf,j)) then
            rrmax=xinf(j)
            jrmax=j
         endif
      enddo
      sqrpsi0=sqrt(psi0)*(rrmax-rrmin)/5.
      zmz02=(zzmax-zzmin)**2/4.
      rmr02=(rrmax-rrmin)**2/4.
      zmrm2=rmr02*zmz02
c      write(*,*) 'zmax=z(',jzmax,')=',zzmax,' jzmax=',jzmax
c      write(*,*) 'zmin=z(',jzmin,')=',zzmin,' jzmin=',jzmin
c      write(*,*) 'rmax=r(',jrmax,')=',rrmax,' jrmax=',jrmax
c      write(*,*) 'rmin=r(',jrmin,')=',rrmin,' jrmin=',jrmin
      return
      end
c***********************************************************************
      function bfo(b,b0,r1,r0o,z11,z0)
c bfo,b,b0(T),rax,z0 are output parameter, b0 is vacuum field
c r1,z1 are input parameters 
      bfo=bfof(b,b0,r1,r0o,z11,z0,-1)
      return
      end
c*********************************************************************
      function bfof(b,b0,r1,r0o,z11,z0,kfast)
c bfof,b,b0(T),rax,z0 are output parameter, b0 is vacuum field
c r1,z1 are input parameters 
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      common/asymm/asym
      dimension cl(2),cu(2)
      common/saveij/i,j,cl,cu
      z0=zz(1,1)
c      write(*,*) 'z1,z11=',z1,z11,z0,asym
      z1=z0+asym*abs(z11-z0)
c      write(*,*) 'z1,z11=',z1,z11,z0,asym
      bfof=0.
      b=0.
      b0=0.
      r0o=rax
      if(kfast.gt.0.and.i.gt.0) goto 10
c     write(*,*) 'z1,z11=',z1,z11
      call ij(i,j,r1,z1,r0o,z0,kfast)
c      write(*,*) 'z1,z11=',z1,z11
c      write(*,*)'i,iz,ir,j,r1,z1,r0o,z0',i,iz,ir,j,r1,z1,r0o,z0
      if(i.le.0) return
      if(i.eq.1) then
         if(sqrt((r1-rax)**2+(z1-z0)**2).lt.1.05*abs(rr(1,1)-rr(1,mth/2)
     &        )/2)then
            if(cl(1).lt.0..or.cl(2).lt.0..or.cu(1).lt.0..or.cu(2).lt.0.
     &        )then
c               write(*,*)'6cl,cu,r1,z1,ru,zu,rl,zl=',cl,cu,r1,z1,rr(i,j)
c     &              ,rr(i+1,j),zz(i,j),zz(i+1,j),rr(i,j+1),rr(i+1,j+1)
c     &              ,zz(i,j+1),zz(i+1,j+1)
               cl(1)=0.5
               cu(1)=0.5
               cl(2)=0.
               cu(2)=0.
            endif
         else
            call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      else
         jc1=(j+jd0)/jdj
         if(jc1.eq.1.or.jc1.eq.3) then
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),zz(i,j),zz(i,j+1),rr(i,j),rr(i,j+1),z1,r1,cl,cu)
         else
            call weightk(aint1(i,j+jc1),ajnt(i,j+jc1),ajnt(i,j+jc1+1)
     $           ,bjnt(i,j+jc1),bjnt(i,j+jc1+1),sqrajnt(i,j+jc1),
     $           sqrajnt(i,j+jc1+1),asqrajnt(i,j+jc1),asqrajnt(i,j+jc1
     $           +1),rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
         endif
      endif
call weight(rr(i,j),rr(i,j+1),zz(i,j),zz(i,j+1),r1,z1,cl,cu)
 10   bfof=f2dperp(rg2d(i),rg2d(i),cl,cu)/r1
      b=f2dperp(b2d(i,j),b2d(i,j+1),cl,cu)
      b0=(b2d(1,1)+b2d(1,mth/2+1))/2.
      return
      end
c*********************************************************
c This subr. does not need common blocks
      function zatpsir(eps,psii,r,r00,z00)
      external fzatpsir
      common/zatpsi/psi0,z0,r0,r1,izat
      common/asymm/asym
      r1=r
      psi0=psii
      r0=r00
      z0=z00
      zatpsir=psiofpp2(r,r0,z0,a)
cpsio(dpsi,ddpsi,dzpsi,r,r0,z00,z0,a)
      if(psii.le.zatpsir) then
         zatpsir=z0
      else
care  This *2.5 may be changed 
         izat=0
         call ROOT(z0,a*2.5*asym,eps,0.,zatpsir,fzatpsir)
         if(abs(zatpsir-z0).le.eps) zatpsir=z0
c         write(*,*)'++++',z0,a*2.5*asym,eps,zatpsir
      endif
      return
      end
c*********************************************************
      function fzatpsir(z)
      common/zatpsi/psii,z0,r0,r1,izat
      izat=izat+1
      if(izat.lt.40.or.mod(izat,4).eq.0)then
         fzatpsir=psiof(dpsi1,ddpsi1,dzpsi1,r1,r0,z,z0,a,-1)-psii
cpp(r1,r0,z,z0,a,-1)-psii
      else
         fzatpsir=psiofpp(r1,r0,z,z0,a,0)-psii
      endif
c      write(*,*) 'fzat=',z,fzatpsir,izat
      return
      end
c*********************************************************
      function d1eq0(r)
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      call psorb(b,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,r,r0
     $     ,rvpar,sig,v,z,z0,xhp)
      d1eq0=psio(dpsi,ddpsior,dzpsi,r,r0,z0,z0,am)
      d1eq0=dpsior-dpsi
      return
      end
c*********************************************************
      function d2eq0(r)
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      call psorb(b,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,r,r0
     $     ,rvpar,-1.,v,z,z0,xhp)
      d2eq0=psio(dpsior,ddpsi,dzpsi,r,r0,z0,z0,am)
      d2eq0=ddpsior-ddpsi
      return
      end
c*********************************************************
      function psieq0(r)
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      common/root2/kfast
      save zh
      kfast=kfast+1
c      call psorb(b,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,r,r0
c     $     ,rvpar,sig,v,z,z0,xhp)
      if(kfast.lt.3) then
         call psorb(b,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,r,r0
     $        ,rvpar,sig,v,z,z0,xhp)
         psieq0=psiofpp2(r,r0,z0,am)-psior
      else
         call psorbf(b,b0,eps,iter,qpf,pmu,psior,r,r0
     $        ,rvpar,sig,v,z,z0,xhp,zh)
         psieq0=psiofpp2(r,r0,z0,am)-psior
      endif
      zh=z
cpsiofpp(r,r0,z0,z0,am,0)
      return
      end
c*********************************************************
      FUNCTION FBES(X,fbes0,N)
      fbes0=BESJ0(X)
      N1=N-1
      IF(N1) 1,2,3
 1             FBES=fbes0
      RETURN
 2             FBES=BESJ1(X)
      RETURN
 3             X1=X-1.E-5
      IF(X1) 4,4,5
 4             FBES=0.0
      RETURN
 5             FBES=(2.*BESJ1(X))/X-fbes0
      RETURN
      END
c***********************************************************      
      FUNCTION BESJY(X)

      LOGICAL L

      ENTRY BESJ0(X)

      L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 4
    8 F=0.0625*X**2-2.0
      A =           - 0.00000 00000 000008
      B = F * A     + 0.00000 00000 000413
      A = F * B - A - 0.00000 00000 019438
      B = F * A - B + 0.00000 00000 784870
      A = F * B - A - 0.00000 00026 792535
      B = F * A - B + 0.00000 00760 816359
      A = F * B - A - 0.00000 17619 469078
      B = F * A - B + 0.00003 24603 288210
      A = F * B - A - 0.00046 06261 662063
      B = F * A - B + 0.00481 91800 694676
      A = F * B - A - 0.03489 37694 114089
      B = F * A - B + 0.15806 71023 320973
      A = F * B - A - 0.37009 49938 726498
      B = F * A - B + 0.26517 86132 033368
      A = F * B - A - 0.00872 34423 528522
      A = F * A - B + 0.31545 59429 497802
      BESJY=0.5*(A-B)
      IF(L) RETURN

      A =           + 0.00000 00000 000016
      B = F * A     - 0.00000 00000 000875
      A = F * B - A + 0.00000 00000 040263
      B = F * A - B - 0.00000 00001 583755
      A = F * B - A + 0.00000 00052 487948
      B = F * A - B - 0.00000 01440 723327
      A = F * B - A + 0.00000 32065 325377
      B = F * A - B - 0.00005 63207 914106
      A = F * B - A + 0.00075 31135 932578
      B = F * A - B - 0.00728 79624 795521
      A = F * B - A + 0.04719 66895 957634
      B = F * A - B - 0.17730 20127 811436
      A = F * B - A + 0.26156 73462 550466
      B = F * A - B + 0.17903 43140 771827
      A = F * B - A - 0.27447 43055 297453
      A = F * A - B - 0.06629 22264 065699
      BESJY=0.636619772367581*log(X)*BESJY+0.5*(A-B)
      RETURN

    4 F=256.0/X**2-2.0
      B =           + 0.00000 00000 000007
      A = F * B     - 0.00000 00000 000051
      B = F * A - B + 0.00000 00000 000433
      A = F * B - A - 0.00000 00000 004305
      B = F * A - B + 0.00000 00000 051683
      A = F * B - A - 0.00000 00000 786409
      B = F * A - B + 0.00000 00016 306465
      A = F * B - A - 0.00000 00517 059454
      B = F * A - B + 0.00000 30751 847875
      A = F * B - A - 0.00053 65220 468132
      A = F * A - B + 1.99892 06986 950373
      P=A-B
      B =           - 0.00000 00000 000006
      A = F * B     + 0.00000 00000 000043
      B = F * A - B - 0.00000 00000 000334
      A = F * B - A + 0.00000 00000 003006
      B = F * A - B - 0.00000 00000 032067
      A = F * B - A + 0.00000 00000 422012
      B = F * A - B - 0.00000 00007 271916
      A = F * B - A + 0.00000 00179 724572
      B = F * A - B - 0.00000 07414 498411
      A = F * B - A + 0.00006 83851 994261
      A = F * A - B - 0.03111 17092 106740
      Q=8.0*(A-B)/V
      F=V-0.785398163397448
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432/SQRT(V)
      IF(L) GO TO 6
      BESJY=F*(Q*A+P*B)
      RETURN
    6 BESJY=F*(P*A-Q*B)
      RETURN

      ENTRY BESJ1(X)

      L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 5
    3 F=0.0625*X**2-2.0
      B =           + 0.00000 00000 000114
      A = F * B     - 0.00000 00000 005777
      B = F * A - B + 0.00000 00000 252812
      A = F * B - A - 0.00000 00009 424213
      B = F * A - B + 0.00000 00294 970701
      A = F * B - A - 0.00000 07617 587805
      B = F * A - B + 0.00001 58870 192399
      A = F * B - A - 0.00026 04443 893486
      B = F * A - B + 0.00324 02701 826839
      A = F * B - A - 0.02917 55248 061542
      B = F * A - B + 0.17770 91172 397283
      A = F * B - A - 0.66144 39341 345433
      B = F * A - B + 1.28799 40988 576776
      A = F * B - A - 1.19180 11605 412169
      A = F * A - B + 1.29671 75412 105298
      BESJY=0.0625*(A-B)*X
      IF(L) RETURN

      B =           - 0.00000 00000 000244
      A = F * B     + 0.00000 00000 012114
      B = F * A - B - 0.00000 00000 517212
      A = F * B - A + 0.00000 00018 754703
      B = F * A - B - 0.00000 00568 844004
      A = F * B - A + 0.00000 14166 243645
      B = F * A - B - 0.00002 83046 401495
      A = F * B - A + 0.00044 04786 298671
      B = F * A - B - 0.00513 16411 610611
      A = F * B - A + 0.04231 91803 533369
      B = F * A - B - 0.22662 49915 567549
      A = F * B - A + 0.67561 57807 721877
      B = F * A - B - 0.76729 63628 866459
      A = F * B - A - 0.12869 73843 813500
      A = F * A - B + 0.04060 82117 718685
      BESJY=0.636619772367581*log(X)*BESJY-0.636619772367581/X
     1     +0.0625*(A-B)*X
      RETURN

    5 F=256.0/X**2-2.0
      B =           - 0.00000 00000 000007
      A = F * B     + 0.00000 00000 000055
      B = F * A - B - 0.00000 00000 000468
      A = F * B - A + 0.00000 00000 004699
      B = F * A - B - 0.00000 00000 057049
      A = F * B - A + 0.00000 00000 881690
      B = F * A - B - 0.00000 00018 718907
      A = F * B - A + 0.00000 00617 763396
      B = F * A - B - 0.00000 39872 843005
      A = F * B - A + 0.00089 89898 330859
      A = F * A - B + 2.00180 60817 200274
      P=A-B
      B =           + 0.00000 00000 000007
      A = F * B     - 0.00000 00000 000046
      B = F * A - B + 0.00000 00000 000360
      A = F * B - A - 0.00000 00000 003264
      B = F * A - B + 0.00000 00000 035152
      A = F * B - A - 0.00000 00000 468636
      B = F * A - B + 0.00000 00008 229193
      A = F * B - A - 0.00000 00209 597814
      B = F * A - B + 0.00000 09138 615258
      A = F * B - A - 0.00009 62772 354916
      A = F * A - B + 0.09355 55741 390707
      Q=8.0*(A-B)/V
      F=V-2.356194490192345
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432/SQRT(V)
      IF(L) GO TO 7
      BESJY=F*(Q*A+P*B)
      RETURN
    7 BESJY=F*(P*A-Q*B)
      IF(X .LT. 0.0) BESJY=-BESJY
      RETURN

      ENTRY BESY0(X)

      IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 4
      GO TO 8

      ENTRY BESY1(X)

      IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 5
      GO TO 3

    9 BESJY=0.
      write(*,100)X
      RETURN
  100 FORMAT(1X,36HBESJY ... NON-POSITIVE ARGUMENT X = ,E15.4)

      END
c*************************
      subroutine equiv(fb,fb1,n)
      dimension fb(n),fb1(n)
      do i=1,n
         fb(i)=fb1(i)
      enddo
      return
      end
c******************************************************************
      subroutine zero(f,n)
      dimension f(n)
      do i=1,n
         f(i)=0.
      enddo
      return
      end
c******************************************************************
      subroutine zeroo(f,n)
      dimension f(n)
      do i=1,n
         f(i)=0.
      enddo
      return
      end
c*********************fast version************************
c This subr. does not need common blocks
      function zatpsirf(eps,psii,r,r00,z00,zgues,h2)
      external fzatpsir
      common/zatpsi/psi0,z0,r0,r1,izat
      common/asymm/asym
      r1=r
      psi0=psii
      r0=r00
      z0=z00
      zatpsirf=psiofpp2(r,r0,z0,a)
      if(psii.le.zatpsirf) then
         zatpsirf=z0
      else
         h=amax1(abs(z0-zgues)-h2*a,0.)*asym
care  This *2.5 may be changed 
         izat=0
         call ROOT(z0+h,zgues+h2*a*asym,eps,0.,zatpsirf,fzatpsir)
c         write(*,*)'------',z0,a*2.5*asym,zatpsirf,h,zgues
cczgues+.1*asym
         if(abs(zatpsirf-z0).le.eps) zatpsirf=z0
      endif
      return
      end
c***********************************************************************
c     Output are
c
c     dpsior= dPsi/dlnR, ddpsior= RRd^2Psi/dR^2,
c     psior= Psi, z
c
c     Inputs are eps,iter,qpf,pmu,r,r0,rvpar,sigm,v,z0,xhp
c
      subroutine psorbf(b,b0,eps,iter,qpf,pmu,psior,r,r0
     $     ,rvpar,sigm,v,z,z0,xhp,zgues)
      common /kfow/xfow
      psior=0.
      z=0.
c
      if(r.lt.rvpar) return
      if(rvpar.eq.r)then
         pss=0.
      else
         pss=sqrt(1.-rvpar/r)
      endif
c      write(*,*)'1',pss
      dpsior=sigm*xfow*xhp*v*r
      if(iter.gt.1)then
         psior=max(qpf,0.)+dpsior*pss/2.
      else
         psior=qpf+dpsior*pss
      endif
      z=zatpsirf(eps,psior,r,r0,z0,zgues,.5)
c      write(*,*) '1--- psorb:: ',psior,z,pss,b,b0
      hp=bfof(b,b0,r,r00,z,z00,1)
      ddpsior=dpsior
c      write(*,*)'1 psorb:: ',b,b0,r,r00,z,z00,ddpsior
cvac change this line if vacuum included
      if(b0.lt..1) then
         psior=1.
         return
      endif
ckg      if(abs(rvpar-r)/r.lt.((b*r/b0/r00)-1.)**2) return
c      write(*,*)'2',pss
      i2=1
      do i=i2,iter
         if(i.ne.i2) then
            psior=bfof(b,b0,r,r00,z,z00,1)
         else
            psior=hp
         endif
c         write(*,*) '2 psorb:: ',b0,r,z,b,psior
         if(b0.lt..1) then
            psior=1.
            return
         endif
         pss=(1.-pmu*b/b0/r00)
         if(rvpar.eq.r)then
            pss=0.
         else
            pss=(1.-pmu*b/b0/r00)
            if(pss.lt.0.) pss=max((1.-rvpar/r),1.e-10)
         endif
         pss=sqrt(pss)
         dpsior=ddpsior*psior/b
         psior=qpf+dpsior*pss
         hp=z
         z=zatpsirf(eps,psior,r,r0,z0,hp,1.)
c         write(*,*) '3--- psorb:: ',psior,z,pss,b,b0
      enddo
      return
      end
c******************************************************************
        function ynn(x,f,xn)
c Two points Lagrange' formula for a function
        dimension x(2),f(2)
        ynn=f(1)+(f(2)-f(1))*(xn-x(1))/(x(2)-x(1))
        return
        end
c******************************************************************
	SUBROUTINE Atra(x,y,nns,res)
c*** this subroutine for calculation of integral by the trapetsii with
c*** nonequal step.
c*** input parameters: x(nns)-massive of variables:
c***                   y(nns)-massive of function;
c***                   nns-numbers  of variables.
c*** output parameter: res - result.
	dimension x(nns),y(nns)
	res=0.
	if (nns.eq.2) then
	   res=res+(y(nns-1)+y(nns))*(x(nns)-x(nns-1))/2.
	else
	   do i=1,nns-2,2
	      res=res+(y(i)+y(i+1))*(x(i+1)-x(i))/2.+
	1	   (y(i+2)+y(i+1))*(x(i+2)-x(i+1))/2.
	   enddo
	   rm=mod(real(nns-1),2.)
	   if(rm.lt.0.1) then
	      return
	   else
	      res=res+(y(nns-1)+y(nns))*(x(nns)-x(nns-1))/2.
	   endif
	endif
 	return
	end
C******************************************************
      subroutine omstari_calc(ai,omstari,vi,zi)
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      include 'orb2d.par'
      dimension ind(nn),xrt(nn)
      q_m1=1.
c assume no more then two roots
      call roots(qoo(1),ppsi(1),ind(1),nosurf,q_m1,1.e-6,0.,xrt(1),nroot
     &     )
      if((q1-q_m1)*(qoo(1)-q_m1).gt.0.)then
      write(*,*) '------------------------'
         write(*,*) 
     &   'Warning: your surface is out of the safety factor limits'
ckg
         return
      endif
c     
c      nroot=1
c      xrt(nroot)=ppsi(1)+0.32**2*(ppsi(nosurf)-ppsi(1))
c      ind(nroot)=indx(ppsi,nosurf,xrt(nroot))
      call funder4(ppsi(ind(nroot)-1),pc(ind(nroot)-1,1),xrt(nroot)
     &       ,dyn,yn)
c This one is for D ions
      omstari=omstari*zi/ai*(vc0(2)/vi/1.e+9)**2*dyn/yn
c This one is for electrons
c      omstari=omstari*zi/ai*(vc0(1)/vi/1.e+9)**2*dyn/yn*rmcp(1)/2.
      write(*,*) ' D plasma rot. frequency =',omstari
c      write(*,*) ppsi,pc,ind(nroot),xrt(nroot),dyn,yn
      return
      end
c************************************************************************
      subroutine fun4_eq(f,p,yn)
c Four points Lagrange' formula for function on equal step grid
c p must be: p=(x-x0)/h, where x0 is the second point
      dimension f(4)
      yn=p-1.
      yn=(f(1)*(2.-p)+f(4)*(p+1.))/6.*yn*p+(f(2)*yn-p*f(3))/2.
     &     *(p-2.)*(p+1.)
      return
      end
c**************************************************************
cHere the error function is calculated with the accuracy $\le 2.5\ 10^{-5}$
cN. Gorelenkov, 24 Sept. 2002
      function erf1(x)
      if(x.gt.5.) then
         erf1=1.
      elseif(x.lt.-5.)then
         erf1=-1.
      else
         erf1=1./(1.+0.47047*abs(x))
         erf1=sign(1.,x)*(1.-(0.3480242-(0.0958798-0.7478556*erf1)*erf1)
     &        *erf1*exp(-x**2))
      endif
      return
      end
c**************************************************************
      SUBROUTINE E04CCF(N,X,FMIN,EPS,N1,PDSTAR,PSTAR,PBAR,STEP,Y,P,
     *                  FUNCT,MONIT,MAXIT,IFAIL)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 3 REVISED.
C     MARK 4.5 REISSUE. LER-F7
C     MARK 8 REVISED. IER-220 (MAR 1980)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-716 (DEC 1989).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04CCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, FMIN
      INTEGER           IFAIL, MAXIT, N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  P(N1,N), PBAR(N), PDSTAR(N), PSTAR(N), STEP(N),
     *                  X(N), Y(N1)
C     .. Subroutine Arguments ..
      EXTERNAL          FUNCT, MONIT
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, C, CENT, COEFF, DERIV, DERIV2, F1, F2, F3,
     *                  FMAX, R, SERROR, X1, X2, X3, XMIN, YDSTAR,
     *                  YMEAN, YSTAR, F1FIX
      INTEGER           H, I, IV, J, K, L, LASTMX, MCOUNT, NCALL, NP1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT, DBLE
C     .. Executable Statements ..
      FMAX = 1000.D0
      FMIN = 0.D0
      NCALL = 0
      R = 0
      IF (N.LT.1 .OR. EPS.LT.X02AJF()
     *    .OR. N1.NE.N+1 .OR. MAXIT.LT.1 .OR.
     *    (IFAIL.NE.1 .AND. IFAIL.NE.0)) GO TO 1220
      CALL FUNCT(N,X,F1FIX)
      NCALL = NCALL + 1
   20 DO 300 I = 1, N
         F1 = F1FIX
         X1 = 0.D0
         COEFF = 1.D0
         DO 80 J = 1, N
            IF (I-J) 60, 40, 60
   40       PSTAR(J) = X(J) + COEFF
            GO TO 80
   60       PSTAR(J) = X(J)
   80    CONTINUE
         CALL FUNCT(N,PSTAR,F2)
         NCALL = NCALL + 1
         X2 = 1.D0
         PSTAR(I) = PSTAR(I) + COEFF
         CALL FUNCT(N,PSTAR,F3)
         NCALL = NCALL + 1
         X3 = 2.D0
  100    IF (NCALL.GT.MAXIT) GO TO 1120
         DERIV = (X2-X3)*F1 + (X3-X1)*F2 + (X1-X2)*F3
         IF (ABS(DERIV)-X02AJF()) 140, 120, 120
  120    DERIV2 = DERIV/(X1-X2)/(X2-X3)/(X3-X1)
         IF (DERIV2) 240, 140, 140
  140    IF (F1-F3) 160, 200, 200
  160    IF (X1.LE.-5.0D0) GO TO 180
         F3 = F2
         X3 = X2
         F2 = F1
         X2 = X1
         X1 = X1 - COEFF
         PSTAR(I) = X(I) + X1
         CALL FUNCT(N,PSTAR,F1)
         NCALL = NCALL + 1
         GO TO 100
  180    XMIN = -5.0D0
         GO TO 280
  200    IF (X3.GE.5.0D0) GO TO 220
         F1 = F2
         X1 = X2
         F2 = F3
         X2 = X3
         X3 = X3 + COEFF
         PSTAR(I) = X(I) + X3
         CALL FUNCT(N,PSTAR,F3)
         NCALL = NCALL + 1
         GO TO 100
  220    XMIN = 5.0D0
         GO TO 280
  240    XMIN = .5D0*((X2**2-X3**2)*F1+(X3**2-X1**2)*F2+(X1**2-X2**2)
     *          *F3)/DERIV
         IF (XMIN.NE.0.0D0) GO TO 260
         XMIN = 0.1D0
         GO TO 280
  260    IF (ABS(XMIN).LT.0.1D0) XMIN = SIGN(0.1D0,XMIN)
         IF (ABS(XMIN).GT.5.0D0) XMIN = SIGN(5.0D0,XMIN)
  280    STEP(I) = XMIN
  300 CONTINUE
      NP1 = N + 1
      DO 380 I = 1, NP1
         DO 360 J = 1, N
            IF (I-J-1) 320, 340, 320
  320       PSTAR(J) = X(J)
            P(I,J) = X(J)
            GO TO 360
  340       PSTAR(J) = X(J) + STEP(J)
            P(I,J) = X(J) + STEP(J)
  360    CONTINUE
         NCALL = NCALL + 1
         CALL FUNCT(N,PSTAR,Y(I))
  380 CONTINUE
      A = 1.D0
      B = .5D0
      C = 2.D0
      LASTMX = 0
      MCOUNT = 0
      K = 0
  400 K = K + 1
      FMAX = Y(1)
      FMIN = Y(1)
      H = 1
      L = 1
      DO 480 I = 2, N1
         IF (Y(I)-FMAX) 440, 440, 420
  420    FMAX = Y(I)
         H = I
         GO TO 480
  440    IF (Y(I)-FMIN) 460, 480, 480
  460    FMIN = Y(I)
         L = I
  480 CONTINUE
      IF (LASTMX-H) 580, 500, 580
  500 MCOUNT = MCOUNT + 1
      IF (MCOUNT-5) 600, 520, 600
  520 IF (H-1) 560, 540, 560
  540 H = 2
      FMAX = Y(H)
      GO TO 600
  560 H = 1
      FMAX = Y(H)
      GO TO 600
  580 LASTMX = H
      MCOUNT = 0
  600 CONTINUE
      CALL MONIT(FMIN,FMAX,P,N,N1,NCALL)
      IF (K.EQ.1) GO TO 620
      IF (SERROR-EPS) 1140, 620, 620
  620 IF (NCALL-MAXIT) 640, 1160, 1160
  640 DO 700 J = 1, N
         CENT = 0.D0
         DO 680 I = 1, N1
            IF (I-H) 660, 680, 660
  660       CENT = CENT + P(I,J)
  680    CONTINUE
         PBAR(J) = CENT/DBLE(N)
  700 CONTINUE
C     REFLECTION
      DO 720 I = 1, N
         PSTAR(I) = (1.D0+A)*PBAR(I) - A*P(H,I)
  720 CONTINUE
      CALL FUNCT(N,PSTAR,YSTAR)
      NCALL = NCALL + 1
      IF (YSTAR-FMIN) 740, 780, 780
C     EXPANSION
  740 DO 760 I = 1, N
         PDSTAR(I) = (1.D0+C)*PSTAR(I) - C*PBAR(I)
  760 CONTINUE
      CALL FUNCT(N,PDSTAR,YDSTAR)
      NCALL = NCALL + 1
      IF (YDSTAR-YSTAR) 980, 1020, 1020
C     CONTRACTION
  780 DO 820 I = 1, N
         IF (I-H) 800, 820, 800
  800    IF (YSTAR-Y(I)) 1020, 820, 820
  820 CONTINUE
      IF (FMAX-YSTAR) 880, 840, 840
  840 DO 860 I = 1, N
         P(H,I) = PSTAR(I)
  860 CONTINUE
  880 CONTINUE
      DO 900 I = 1, N
         PDSTAR(I) = B*P(H,I) + (1.D0-B)*PBAR(I)
  900 CONTINUE
      CALL FUNCT(N,PDSTAR,YDSTAR)
      NCALL = NCALL + 1
      IF (YDSTAR-FMAX) 980, 980, 920
  920 DO 960 I = 1, N1
         DO 940 J = 1, N
            PBAR(J) = (P(I,J)+P(L,J))*0.5D0
            P(I,J) = PBAR(J)
  940    CONTINUE
         CALL FUNCT(N,PBAR,Y(I))
         NCALL = NCALL + 1
  960 CONTINUE
      GO TO 1060
  980 DO 1000 J = 1, N
         P(H,J) = PDSTAR(J)
 1000 CONTINUE
      Y(H) = YDSTAR
      GO TO 1060
 1020 DO 1040 J = 1, N
         P(H,J) = PSTAR(J)
 1040 CONTINUE
      Y(H) = YSTAR
 1060 YMEAN = 0.D0
      SERROR = 0.D0
      DO 1080 I = 1, N1
         YMEAN = YMEAN + Y(I)
 1080 CONTINUE
      YMEAN = YMEAN/DBLE(N+1)
      DO 1100 I = 1, N1
         SERROR = SERROR + (Y(I)-YMEAN)**2
 1100 CONTINUE
      SERROR = SQRT(SERROR/DBLE(N+1))
      GO TO 400
 1120 IV = 2
      GO TO 1240
 1140 IV = 0
      GO TO 1180
 1160 IV = 2
 1180 DO 1200 I = 1, N
         X(I) = P(L,I)
 1200 CONTINUE
      GO TO 1240
 1220 IV = 1
 1240 IFAIL = P01ABF(IFAIL,IV,SRNAME,0,P01REC)
      RETURN
      END
c**************************************************************
      subroutine elongs(elong,isrf)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      zmin=10000.
      zmax=-zmin
      rmin=zmin
      rmax=zmax
      do its=1,nts0
         zmin=min(zmin,zz(isrf,its))
         zmax=max(zmax,zz(isrf,its))
         rmin=min(rmin,rr(isrf,its))
         rmax=max(rmax,rr(isrf,its))
      enddo
c      write (*,*)'zmin,zmax,rmin,rmax',zmin,zmax,rmin,rmax
      elong=(zmax-zmin)/(rmax-rmin)
      return
      end
c**************************************************************
        SUBROUTINE boundtoc(xb,y,xc,nx)
        dimension work(200),xc(nx),xb(nx),y(nx)
c change from the boundary grid  TRANSP to the central
         do i=1,nx
            i0=max(1,i-2)
            i0=min(i0,nx-3)
            call fun4(xb(i0),y(i0),xc(i),work(i))
         enddo
         do i=1,nx
            y(i)=work(i)
         enddo
         return
         end
c**************************************************************
c this is not used actually
      subroutine datin
c
c      parameter(ndat=20)
      include 'gridparam'
c      parameter(ndat=40)
      common/dat/rsdat(ndat),rdat(ndat),rhodat(ndat)
     >,ptotdat(ndat),qdat(ndat),betadat(ndat)
     >,dedat(ndat),dddat(ndat),dtdat(ndat),dzdat(ndat),dhdat(ndat)
     >,dhmindat(ndat),dbddat(ndat),dbtdat(ndat),dadat(ndat)
     >,tedat(ndat),tddat(ndat),ttdat(ndat),tzdat(ndat),thdat(ndat)
     >,thmindat(ndat),pbddat(ndat),pbtdat(ndat),padat(ndat)
     >,fluxt(ndat),fluxp(ndat),rtsdat(ndat),shift(ndat)
     >,diondat(ndat),dbeamdat(ndat)
     >,ue(ndat),ui(ndat),ubdpar(ndat),ubdper(ndat),ubtpar(ndat)
     >,ubtper(ndat),ufipa(ndat),ufipp(ndat),uphi(ndat)
     >,ufastpa(ndat),ufastpp(ndat),utot(ndat)
     >,ftotdt(ndat),ttntx(ndat)
     >,betath(ndat),betaadat(ndat),pthdat(ndat),condic(ndat)
     >,condec(ndat)
c2d     >,dvol(ndat),darea(ndat),surfa(ndat)
      common / chnel / iodat,iomode,iomap1,ioequ1,mp0,mp1,mp2,iotty

      character*80 title,sc
c
c
      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c      write(*,*) title,t
c Toroidal flux (Webers) at zone bdry
      read(10,99) (fluxt(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Poloidal Flux (Webers) at zone bdry
      read(10,99) (fluxp(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c x=zone centr grid = sq_root(normlz toroid flux) (~r/a)
      read(10,99) (rtsdat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Shafranov shift (cm) of zone bdry
      read(10,99) (shift(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Electron density at x (/cm3)
      read(10,99) (dedat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Hydrogen density at x (/cm3)
      read(10,99) (dhdat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Thermal deuterium density at x (/cm3)
      read(10,99) (dddat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Thermal tritium density at x (/cm3)
      read(10,99) (dtdat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Impurity density at x (/cm3)
      read(10,99) (dzdat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Total ion density (NH+ND+NT+NIMP) at x (/cm3)
      read(10,99) (diondat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Fast alpha density at x (/cm3)
      read(10,99) (dadat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Total beam ion density at x (/cm3)
      read(10,99) (dbeamdat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c  D-Beam ion density at x (/cm3)
      read(10,99) (dbddat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c T-Beam ion density at x (/cm3)
      read(10,99) (dbtdat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Impurity Temperature at x (keV)
      read(10,99) (tzdat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Hydrogenic (ions) Temperature at x (keV)
      read(10,99) (tddat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Electron Temperature at x (keV)
      read(10,99) (tedat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c q_psi at zone bdry
      read(10,99) (qdat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Electron energy density 3*NE*TE/2 (Joules/cm3)
      read(10,99) (ue(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Ion energy density 3*NI*TI/2 (Joules/cm3)
      read(10,99) (ui(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c D-beam parallel energy density (Joules/cm3)
      read(10,99) (ubdpar(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c D-beam perpendicular energy density (Joules/cm3)
      read(10,99) (ubdper(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c T-beam parallel energy density (Joules/cm3)
      read(10,99) (ubtpar(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c T-beam perpendicular energy density (Joules/cm3)
      read(10,99) (ubtper(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Fusion ion parallel energy density (Joules/cm3)
      read(10,99) (ufipa(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Fusion ion perpendicular energy density (Joules/cm3)
      read(10,99) (ufipp(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c  Rotation energy density (Joules/cm3)
      read(10,99) (uphi(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Total energy den (UE+UI+UPHI+UFASTPA+UFASTPP)
      read(10,99) (utot(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c DT neutron emission rate (#/cm3/sec)
      read(10,99) (ftotdt(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Total neutron emission rate (#/cm3/sec)
      read(10,99) (ttntx(k),k=1,ndat)
 
      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Total Toroidal Beta (=BTPL+BTBM+BTMIN+BTFI) at x
      read(10,99) (betadat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Thermal Toroidal Beta (=BTE+BTI) at x
      read(10,99) (betath(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Alpha Toroidal Beta at x
      read(10,99) (betaadat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c Thermal Pressure 2(UI+UE)/3 (Pascals)
      read(10,99) (pthdat(k),k=1,ndat)

      read(10,8) sc
c      write(*,*) sc
      read(10,9) title,t
c      write(*,9) title,t
c Total Pressure (thermal+3*UPHI/2+UFASTPA+UFASTPP/2)
      read(10,99) (ptotdat(k),k=1,ndat)

c      read(10,8) sc
c      write(*,*) sc
c      read(10,9) title,t
c Z_eff profile
c      read(10,99) (zeffdat(k),k=1,ndat)

      rewind 10
      do 21 k=1,ndat
cc  thermal Tritium and Hydrogen have same temperature as thermal Deuterium 
      ttdat(k)=tddat(k)
      thdat(k)=tddat(k)
cc  define beam-D and beam-T pressure
      pbddat(k)=(ubdpar(k)+ubdper(k)*0.5)*1.e6
      pbtdat(k)=(ubtpar(k)+ubtper(k)*0.5)*1.e6
cc  alpha pressure = beta_alpha (arbitrary magnitude)
      padat(k)=(ufipa(k)+ufipp(k)*0.5)*1.e6
   21 continue

    8 format(a80)
    9 format(a59,f8.4)
   99 format(5(2x,1pe11.4))
      call smooth(ptotdat,ndat)
      call smooth(betadat,ndat)

      call smooth(dedat,ndat)
      call smooth(dddat,ndat)
      call smooth(dtdat,ndat)
      call smooth(dhdat,ndat)
      call smooth(dzdat,ndat)
cc smooth beam density data twice 
      call smooth(dbddat,ndat)
      call smooth(dbtdat,ndat)
      call smooth(dadat,ndat)
      call smooth(tedat,ndat)
      call smooth(tddat,ndat)
      call smooth(ttdat,ndat)
      call smooth(thdat,ndat)
      call smooth(tzdat,ndat)
cc smooth beam pressure data twice 
      call smooth(pbddat,ndat)
      call smooth(pbtdat,ndat)

      call smooth(padat,ndat)
      call smooth(betaadat,ndat)
c
ccc assume impurity (Carbon and other heavy metal) mass/proton ratio = 13
c
      do 13 i=1,ndat
      rhodat(i)=2.*dddat(i)+3.*dtdat(i)
     &         +2.*dbddat(i)+3.*dbtdat(i)
     &         +4.*dadat(i)+13.0*dzdat(i)
   13 continue
c
c
cc extrapolate the quantities to normalized sqrt(tor_flux)=0 and 1
      call extrap1(fluxp,rtsdat,ndat)
      call extrap1(dedat,rtsdat,ndat)
      call extrap1(rhodat,rtsdat,ndat)
      call extrap1(ptotdat,rtsdat,ndat)
      call extrap1(betadat,rtsdat,ndat)
      call extrap1(qdat,rtsdat,ndat)

      call extrap1(dddat,rtsdat,ndat)
      call extrap1(dtdat,rtsdat,ndat)
      call extrap1(dzdat,rtsdat,ndat)
      call extrap1(dhdat,rtsdat,ndat)

      call extrap1(dbddat,rtsdat,ndat)
      call extrap1(dbtdat,rtsdat,ndat)
      call extrap1(dadat,rtsdat,ndat)
c      call extrap1(drdat,rtsdat,ndat)

      call extrap1(tedat,rtsdat,ndat)
      call extrap1(tddat,rtsdat,ndat)
      call extrap1(ttdat,rtsdat,ndat)
      call extrap1(thdat,rtsdat,ndat)
      call extrap1(tzdat,rtsdat,ndat)

      call extrap1(pbddat,rtsdat,ndat)
      call extrap1(pbtdat,rtsdat,ndat)
      call extrap1(padat,rtsdat,ndat)
      call extrap1(betaadat,rtsdat,ndat)
c
      do 20 k=1,ndat
cc  rsdat is the normalized poloidal flux
      rsdat(k)=(fluxp(k)-fluxp(1))/(fluxp(ndat)-fluxp(1))
      rdat(k)=sqrt(rsdat(k))
   20 continue
c
c   rescale rhodat profile so that rhodat(1)=1.0
      dedat0=dedat(1)
      rhodat0=rhodat(1)
c2d      do 977 i=1,ndat
ccc choose electron density to be proportional to mass density profile
c      rhodat(i)=dedat(i)/dedat0
ccc use mass density profile
c2d      rhodat(i)=rhodat(i)/rhodat0
c2d  977 continue
c   prepare for spline fitting
      call splprp(rsdat,ndat)
c2d      call depose(rhodat,spcd)
      return
      end
c**************************************************************
c dummy subroutine to silent the compiler
        SUBROUTINE ij2rz(x8,y8,xw8,yw8)
        return
        end
c************************************************************************
        SUBROUTINE fun2(x,f,xn,yn)
c Two points Lagrange' formula for a function
        dimension x(2),f(2)
        yn=f(1)+(f(2)-f(1))*(xn-x(1))/(x(2)-x(1))
        return
        end
c************************************************************************
c     The real part of the integral (int) involved in the steady
c     solutions existence criterion was calculated for several values of
c     ratio=nu_diff/nu_drag using Mathematica. The results were
c     interpolated in order to be used in NOVA by V.Duarte
c
      subroutine chirpcritvd(ratio, xint)
      real*8 xint, ratio
      xint=0.
      if (0. <= ratio .and. ratio <= 0.2) then
         xint=-0.5
      else if (0.2 < ratio .and. ratio< 1.4) then
c
         xint=1.252 +(-25.430 +(144.393+(-415.267+(648.395+
     &        (-552.896+(242.236-42.715*ratio)*ratio)*ratio)
     &        *ratio)*ratio)*ratio)*ratio

      else if (1.4 <= ratio .and. ratio <=  7.28) then
         xint=0.3686 +(-0.3383 +(0.12596+(-0.023475+(0.0021759-
     &        0.00008*ratio)*ratio)*ratio)*ratio)*ratio
      else 
         xint=0.
      end if
      
      return
      end
c************************************************************************
c     dump the data into the u4mat 
      subroutine dump24mat(x, y, z, nx, ny, runid, time, egnid
     &     , coms)
      dimension x(nx),y(ny),z(nx,ny)
      character runid*(*), time*(*), egnid*(*), coms*(*)
      open(unit=19,file='Out/contmap.u4m',status='unknown')
      write(19,*) 
     &     'arrays to plot contour map; runid ',runid,
     &     '; time ',time,'; mode ',egnid,'; ',coms,';'
      write(19,*) '3'
      write(19,*) '0 2; Sclrs; Nsclrs;'
      write(19,*) nx,ny,'; ndrag nscatt;'
      write(19,*) '1 2; nudrag; nuscatt;'
      write(19,*) nx,'; nnudrag;'
      write(19,*) x
      write(19,*) ny,'; nnuscatt;'
      write(19,*) y
      write(19,*) '2 1; nuscatt over nudrag;'
      write(19,*) nx*ny,'; nuratios;'
      write(19,*) z
c   fprintf(fpn,"2 2; %s; %s;\n",l1,l2);
c   fprintf(fpn,"%d; %s;\n", (*nx)*(*nlc), l2);
c   while(ic<(*nx)*(*nlc)){
c     fprintf(fpn, "%11.5e ",xin[ic]);
c     ic++;
c   }
c   ic=0; 
c   fprintf(fpn,"\n %d; v/v0; \n", (*nx)*(*nlc));
c   while(ic<(*nx)*(*nlc)){
c     fprintf(fpn, "%11.5e ",yin[ic]);
      close(19)
      return
      end
c************************************************************************
C Subroutine to be incorporated in NOVA-K to calculate Coulomb logarithms accurately. This subroutine is adapted from the one used in TRANSP (developed by Greg Hammett) - see  /p/transpusers/mgorelen/transp/trunk/codesys/source/comput/r8_coulog.for. Here the units (eV,cm**-3,T) are different from the ones in TRANSP (keV,m**-3,T).
C
C The index t refers to test particle while b refers to background particles. The index i accounts for the background particles
C
      subroutine coulog(et,at,zt,denb,tempb,ab,zb,cln,nb,Bmag)

C This calculates the coulomb logarithm for a test particle of energy et
C (ev), atomic mass at and atomic charge zt, colliding with NB
C maxwellian background species with densities denb (/cm**3),
C temperatures tempb (eV), atomic masses ab, atomic charges zb.  The
C electrons should be included as one of the background species, with
C Z=-1 and A=(1./1836.1).  The list of coulomb logarithms is returned
C in cln.  The magnetic field should be input as bmag (Tesla).

      INTEGER nb,i

      REAL*8 sum,omega2,vrel2,rmax,rmincl,rminqu,rmin
      REAL*8 et,at,zt,Bmag
      REAL*8 denb(nb),tempb(nb),ab(nb),zb(nb),cln(nb)
c
c first calculate the maximum impact parameter, rmax:
c
c rmax=( sum_over_species omega2/vrel2 )**(-1/2)
c
c where,
c omega2 = omega_p**2 + omega_c**2
c vrel2=T/m+2.*E_t/m_t
c here rmax is being calculated in cm
      sum=0.
      do 100 i=1,nb
         omega2=1.74*(10.**6.)*(zb(i)**2.)/ab(i)*denb(i)
     1	+9.18*(10.**15.)*zb(i)**2./ab(i)**2.*Bmag**2.
         vrel2=9.58*(10.**11.)*(tempb(i)/ab(i) + 2.*et/at)
         sum=sum+omega2/vrel2
 100  continue
      rmax=sqrt(1./sum)
 
c next calculate rmin, including quantum corrections.  The classical
C rmin is:
c
C rmincl = e_alpha e_beta / (m_ab vrel**2)
C
C where m_ab = m_a m_b / (m_a+m_b) is the reduced mass.
c vrel**2 = 3 T_b/m_b + 2 E_a / m_a
c (Note:  the two different definitions of vrel2 used in this code
c are each correct for their application.)
c
C The quantum rmin is:
C
C rminqu = hbar/( 2 exp(0.5) m_ab vrel)
C
C and the proper rmin is the larger of rmincl and rminqu
C
      do 200 i=1,nb
         vrel2=9.58*(10.**11.)*(3.*tempb(i)/ab(i)+2.*et/at)
         rmincl=137964.*abs(zb(i)*zt)*(ab(i)+at)/ab(i)/at/vrel2
         rminqu=1.9121*(10**(-4))*(ab(i)+at)/ab(i)/at/sqrt(vrel2)
         rmin=max(rmincl,rminqu)
         cln(i) =log(rmax/rmin)
200   continue
      return
      end
c************************************************************************
c dump the Xir component into the gxm.. file in u4m format
      subroutine savexi(iop)
      include 'clich1'
      include 'clich2'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      write(iop,*) nn*mt,'; eigfun;'
      write(iop,*) ((eigfun(in,im,1),in=1,nn),im=1,mt)
      return
      end
