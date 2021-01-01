c 
      subroutine orbit2d(ai,b0oo,deltk,dpga,gamkom,gamkomflr,im1,iter
     &     ,immax,isrfmax,ntgr,omkom,omkomflr,p0ga,pk,scalvavf,symm,
     &     tip2,vi,xfow,zi,istart)
      include 'clich1'
      include 'clich1b'
      include 'clich2'
      include 'orb2d.par'
c no sense to make nres>2 since only first 2 are used
      parameter(iav=6,mtl=mt+1,nres=2)
      common/asymm/asym
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      common/funct/xkpop(nn,mt,nts,3)
      common /kfow/xfow1
      character tip2*1
c
      dimension rrplot(nosurf),yyplot(nosurf,10),zzplot(nts),dintgr(
     &     -mtl:mtl,nres,4,2),fo(iav,mt,-mtl:mtl,nres),vres(nres,
     &     -mtl:mtl),q_res(nres,-mtl:mtl),nroot(-mtl:mtl),dfalp(nres,
     &     -mtl:mtl),dfalpdv(nres,-mtl:mtl),aeqv(nres,-mtl:mtl)
     &     ,psiavres(nres,-mtl:mtl),rfres(2,6,nres,-mtl:mtl),hwb1(nres,
     &     -mtl:mtl),hnub(nres,-mtl:mtl),hnub_se(nres,-mtl:mtl)
      tip=tip2
      if(im1.ne.1.and.im1.ne.0.and.im1.ne.10.and.im1.ne.2.
     &     and.im1.ne.3)then
         write(*,*)
     &        ' Define the value of variable IM1 to be 0 or 1 or 10'
         stop
      endif
      rb=1.6
      asym=1.
      rf=2.5
      zl=-1.
      zu=1.
      r1=rf
      z1=0.
      call intini
      r1=bfo(b,b0,r1,r0o,z1,z0)
      r1=bfo(b,b0,r,r0o,z1,z0)
c      write(*,*) 'b,b0,r,bfo,r0o,z1,z0',b,b0,r,r1,r0o,z1,z0
c      write(*,*) bsqav
c      write(*,*) r1,psitot,rax,r0o,z0,r,b0,rr(nosurf,1)-rr(1,1)
c      do i=1,nosurf
      do i=30,30
         aaamx=0.
         aabmx=0.
         aacmx=0.
         do imt=1,mt
            aaamx=max(abs(eigfun(i,imt,1)),aaamx)
            aabmx=max(abs(eigfun(i,imt,2)),aabmx)
            aacmx=max(abs(eigfun(i,imt,3)),aacmx)
         enddo
c         write(*,*) (rr(1,j),j=1,66)
c         write(*,*) ppsi(i),aaamx,aabmx,aacmx
         r1=rb+real(i-1)*(Rf-rb)/real(nosurf-1)
         rrplot(i)=r1
c         do j=1,mth2
         do j=90,90
            z1=zl+real(j-1)*(zu-zl)/real(mth2)
            if(z1-z0.lt.0.) then
               asym=-1.
            else
               asym=1.
            endif
            zzplot(j)=z1
c            z1=(zz(i,j)+zz(i+1,j))/2.
c            z1=(zz(i,j))
c            if(j.eq.1.or.j.eq.61) then
            ddpsio(i,j)=psio(dpsi1,ddpsi1,dzpsi1,r1,r0o,z1,z0,a1)
c            else
c               ddpsio(i,j)=psiof(dpsi1,ddpsi1,dzpsi1,r1,r0o,z1,z0,a1,0)
c            endif
c            ddpsio(i,j)=bfo(b,b0,r1,r0o,z1,z0)
            ddpsio(i,j)=dzpsi1
c            ddpsio(i,j)=dpsi1
c            write(*,*)i,j,r1,z1,ddpsio(i,j)
            if(j.eq.(60)) then
               yyplot(i,1)=psio(dpsi1,ddpsi1,dzpsi1,r1,r0o,z0,z0,a1)
               yyplot(i,2)=abs(dpsi1)
               yyplot(i,3)=(ddpsi1)
               rrplot(i)=r1
            endif
         enddo
      enddo
c      open(19,file='pldens.dat',status='unknown')
      do i=1,nosurf
c         write(19,*) rgrid(i),rhoprf(i)*rho0
         do j=1,mt
            if(abs(eigfun(i,j,2)).gt.1.)then
               write(*,*) ' i=',i,' j=',j,' phi=',(eigfun(i,j,k),k=1,3)
               eigfun(i,j,2)=0.
            endif
            if(abs(eigfun(i,j,3)).gt.100.)then
               write(*,*) ' i=',i,' j=',j,' phi=',(eigfun(i,j,k),k=3,3)
c               eigfun(i,j,3)=0.
            endif
         enddo
c         rrplot(i)=(-rr(1,1)+rr(i,1))/(-rr(1,1)+rr(nosurf,1))
      enddo
c      close(19)
c      eigfun(100,6,1)=abs(eigfun(100,6,1))
c      eigfun(100,6,3)=abs(eigfun(100,6,3))
c      write(*,*) ((eigfun(i,j,3),i=1,nosurf),j=1,6)
      moffset=-2
      call grafinit(1)
ckg
      call twodgraf(rgrid,eigfun(1,immax+moffset,1),nn,7,'white','solid'
     &     ,'red',5,'rXi_vsr','\gq\d\gh\u\u1/2','r\gc\dn\u(a.u.)','ts',6
     &     ,6,minm+immax-1+moffset,yyplot,0)
      call twodgraf(rr(1,1),eigfun(1,immax+moffset,1),nn,7,'white'
     &     ,'solid','red',5,'rXi_vsR','R(m)','r\gc\dn\u(a.u.)','ts',6,6
     &     ,minm+immax-1+moffset,yyplot,0)
ckg      call frame(0)
ckg
      call twodgraf(rgrid,eigfun(1,immax+moffset,3),nn,7,'white','solid'
     &     ,'red',5,'Xis_vsr','\gq\d\gh\u\u1/2','\gc\ds\u(a.u.)','ts',6
     &     ,6,minm+immax-1+moffset,yyplot,0)
c fast ion beta, species 1-3
      call twodgraf(rr(1,1),ph(1,1),nn,3,'white','solid','red',5
     &     ,'gbfivsR','R(m)','p\dfi','ts',6,6,1,yyplot,0)
c
      call twodgraf(rgrid,ph(1,1),nn,3,'white','solid','red',5
     &     ,'gbfivsr','\gq\d\gh\u\u1/2','p\dfi','ts',6,6,1
     &     ,yyplot,0)
c fast ion temperatures
      call twodgraf(rgrid,th(1,1),nn,3,'white','solid','red',5
     &     ,'sTfivsr','\gq\d\gh\u\u1/2','T\dfi','ts',6,6,1
     &     ,yyplot,0)
ckg      call frame(0)
ckg      call twodgraf(rgrid,eigfun(1,immax-3,2),nn,7,'white','solid',
ckg     $     'red',5,'delta_p','r/a','delta p.','ts',6,6,minm+immax-1
ckg     &     ,yyplot,0)
ckg      call frame(0)
c      a0=2.5
c      a1=.3
c      do i=1,nosurf
c         rrplot(i)=a0-ppsi(i)**1.2*(a0-a1)
c      enddo
c      call twodgraf(rr(1,1),rrplot(1),nosurf,1,'white','solid',
c     $     'red',5,'ts','R,m','Rotation_w','ts',6,6,1,yyplot,0)
c      call frame(0)
c      call twodgraf(rr(1,1),qoo(1),nosurf,1,'white','solid',
c     $     'red',5,'ts','R,m','q profile','ts',6,6,1,yyplot,0)
c      call frame(0)
c      if(symm.eq.0)then
c      call twodgraf(rgrid,eigfun_s(1,2,1),nn,4,'white','solid',
c     $     'red',5,'ts','r/a','Displ,norm','ts',6,6,1,yyplot,0)
c      call frame(0)
c      call twodgraf(rgrid,eigfun_s(1,2,3),nn,4,'white','solid',
c     $     'red',5,'ts','r/a','Displ,surf','ts',6,6,1,yyplot,0)
c      call frame(0)
c      call twodgraf(rgrid,eigfun_s(1,2,2),nn,4,'white','solid',
c     $     'red',5,'ts','r/a','Per.Pres.','ts',6,6,1,yyplot,0)
c      call frame(0)
c      endif
c      call contgraf(ddpsio,nosurf,nts,rb,rf,
c     $     zl,zu,1,0,xx,0,rr,zz,2,0,
c     $     'test','R','Z','Mol',6,6)
c      call frame(0)
c      call contgraf(ddpsio,nosurf,nts,rb,rf,
c     $     zl,zu,1,0,xx,0,rr,zz,2,0,
c     $'test','R','Z','Mol',6,6)
c c     call frame(0)
c c     call twodgraf(rrplot,yyplot(1,1),nosurf,3,'white','solid',
c     $     'red',5,'ts','R','R2d2Ps/d2R','ts',6,6,1,yyplot,0)
c      call frame(0)
c      call twodgraf(rrplot,yyplot(1,1),nosurf,1,'white','solid',
c     $     'red',5,'ts','R','R2d2Ps/d2R','ts',6,6,1,yyplot,0)
c     call frame(0)
cc      call twodgraf(rrplot,yyplot(1,3),nosurf,1,'white','solid',
cc     $     'red',5,'ts','R','R2d2Ps/d2R','ts',6,6,1,yyplot,0)
c      call frame(0)
c this is recommended values for TFTR
      eps=1.e-5
      epsp=1.e-4
c 
c if plasma is up down symmetric then symm=1. otherwise symm=0.
c      symm=1.
c these parameters are critical
      xfow1=amax1(0.001,xfow)
c this IF is for JET test case
      if(tip.eq.'j'.or.tip.eq.'e') then
         do i=1,nosurf
            rhoprf(i)=1.
            tc(i,1)=1.
         enddo
         tc0(1)=5.e3
         dene0=5.0e13
         aeff=2.00
         rho0=dene0*aeff
         xfow1=1.
         ai=1.
         zi=25.
         b0oo=3.0
         talpha=893.e3 
         vi=1.3845e-3*sqrt(talpha/ai)/2.*1.
         ntor=2
         im1=0
ctem         rax=3.0
c 8  $inp1 te0=5.e3 ti0=5.0e3 dene0=5.0e13 b0=3.0e4 zalpha=100.
c  9  iphd=0 iph=8 ihdi=2 aeff=2.00 rmaj=300. rmip=2.0 rmap=1.0
c 10  amin=90. talpha=893.e3 ipbmax=7 zeff=1.00 pawds=1.0 $
      endif
      rmi=rr(nosurf,mth2/2)*1.0001
      rma=rr(nosurf,1)*0.9999
      b0oo1=b0oo
c      b0oo1=b0oo*r/r0o
      write(*,*)' b0oo, r, r0, Iter = ',b0oo,r,r0o, iter
c 1 tests orbits of particles, works as of 12/2004 --NNG
c      call orbittest(rr(nosurf,mth2/2),rr(nosurf,1),1.,1.)
c 2 shows the instant drift frequencies
c      call testwd(rr(nosurf,mth2/2),rr(nosurf,1),1.,1.,ai,
c     $     b0oo1,vi,zi,psitot)
c 3 shows averaged drift frequencies vs v
c      call testorbitv(rr(nosurf,mth2/2),rr(nosurf,1),1.,1.,tip,ai,
c     $     b0oo1,vi,zi,psitot,vivalf)
c 3.a shows averaged drift frequencies vs Pfi
      psmi=1.
      psma=1.
c      call testorbitpf(rr(nosurf,mth2/2),rr(nosurf,1),psmi,psma,tip,ai,
c     $     b0oo1,vi,zi,psitot,vivalf)
c 4 compares with Frank's frequencies      
c      call testcomp(rr(nosurf,mth2/2),rr(nosurf,1),psmi,psma,ai,
c     $     b0oo,b2d,nn,nts,mth,vi,rr,r,zi,zz,psitot)
c 5 test profiles
c       call testprof(b0oo1,ai,b0,ntor,tip,vi,zi)
cc 6 The growth rate calculator thyself
      if(im1.eq.0.or.im1.eq.2.or.im1.eq.1.or.im1.eq.3)then
         call growth_prpl(ai,aeqv,b0oo1,deltk,dfalp,dfalpdv,
     &        dpga,dintgr,fo,hnub,hnub_se,hwb1,mt,mtl,nres,eps
     &        ,epsp,gamkom,gamkomflr,iav,im1,iter,immax,isrfmax,ntgr,
     &        nosurf,nroot,omkom,omkomflr,p0ga,pk,psiavres,psitot,
     &        rfres,rma,rmi,symm,vi,vres,scalvavf,tip,xfow1,zi,istart)
      else if(im1.eq.10)then
c 7 The particle redistribution due to TAE activity
         call nlin_mix(ai,aeqv,b0oo1,deltk,dfalp,dfalpdv,dpga,dintgr,fo
     &        ,hnub,hwb1,mt,mtl,nres,eps,epsp,gamkom,gamkomflr
     &        ,im1,iter,ntgr,nroot,omkom,omkomflr,p0ga,psiavres,psitot
     &        ,q_res,rfres,rma,rmi,symm,vi,vres,scalvavf,tip,xfow1,zi)
      else
         call grafinit(2)
         stop 'something is wrong'
      endif
      call grafinit(2)
      return
      end
c**********************************************************************        
c This subroutine is to calculate the distribution function of trapped
c      particles in the presence of near threshold TAE activity
c      January 1998 (Gorelenkov)
      subroutine nlin_mix(ai,aeqv,b0oo,deltk,dfalp,dfalpdv,dpga
     &     ,dintgr,fo,hnub,hwb1,mt,mtl,nnres,eps,epsp,gamkom,gamkomflr
     &     ,im1,iter,ntgr,nroot,omkom,omkomflr,p0ga,psiavres,psitot
     &     ,q_res,rfres,rma,rmi,symm,vi,vres,scalvavf,tip,xfow,zi)
      parameter (nv=1,nmu=25,npf=20)
      parameter (iav=6,iberk=30,mt_m1=20)
      dimension qpf(npf),pmu(nmu),v(nv),wvpr(nv,nmu,npf,2),fo(iav,mt,
     &     -mtl:mtl,nnres),wqvpr(nv,nmu,npf,2),wd(nv,nmu,npf,2),wb(nv
     &     ,nmu,npf,2),wres(npf),vres(nnres,-mtl:mtl),q_res(nnres,
     &     -mtl:mtl),gamintpfi(npf,6,2),gamintmu(nmu,8),rf(2,6,nmu,nv)
     &     ,vplot(nv,nmu,4),ind(npf),rfres(2,6,nnres,-mtl:mtl)
     &     ,gamkomis(6,2),psiav(nv,nmu,npf,2),gtrap(npf,6),dintgr(
     &     -mtl:mtl,nnres,4,2),gtrap1(6),nroot(-mtl:mtl),dfalp(nnres,
     &     -mtl:mtl),dfalpdv(nnres,-mtl:mtl),aeqv(nnres,-mtl:mtl)
     &     ,psiavres(nnres,-mtl:mtl),dbob(iberk)
     &     ,gberkintmu(nmu,iberk,2),gberkintpfi(npf,iberk,4)
     &     ,gberkis(iberk,4),hwb1(nnres,-mtl:mtl),hnub(nnres,-mtl:mtl)
     &     ,yy(npf,nmu),q_mix(nmu),i_mix(nmu),psi_aver(nmu)
     &     ,grps2(1000),rgrid(1000)
c stuff for m1 study
      dimension rf_m1(2,6,nv),fo_m1(iav,nv,mt_m1),dwv_m1(iav,
     &     -mt_m1:mt_m1),dfalpdv_m1(nv),omintmu(nmu,2),omintpfi(npf
     &     ,2),omkomis(2)
      common/wdr/alpha,b0ax,psedge,vioo,sigvpr
      common/cgm/b0,brz,pmui,qpfi,z0,imtae0,imt
      save rf
c_________graph\/staff
      character col(npf),pat(npf),tip*1,fow*11,filename*11
      complex wtae
      gam=1.
      db_b=1.e-5
      call zero(yy(1,1),nmu*npf)
      write(fow(1:4),'(f4.2)') scalvavf
      if(ai.eq.4.)then
         filename='gam'//fow(1:1)//fow(3:4)//'.'//fow(5:5)//fow(7:8)/
     &        /'a'
         ihsps=3
      else if(ai.eq.1.)then
         filename='gam'//fow(1:1)//fow(3:4)//'.'//fow(5:5)//fow(7:8)/
     &        /'h'
         ihsps=3
      else if(ai.eq.2.)then
         filename='gam'//fow(1:1)//fow(3:4)//'.'//fow(5:5)//fow(7:8)/
     &        /'d'
         ihsps=1
      else if(ai.eq.3.)then
         filename='gam'//fow(1:1)//fow(3:4)//'.'//fow(5:5)//fow(7:8)/
     &        /'t'
         ihsps=2
      else 
         filename='gam'//fow(1:1)//fow(3:4)//'.'//fow(5:5)//fow(7:8)/
     &        /'?'
         ihsps=3
         write(*,*)
     &        '?????????warning:  Your particle weight is uncommon'
      endif
      write(fow(5:8),'(f4.2)') xfow
      do i=1,npf
         col(i)='white'
         pat(i)='solid'
      enddo
c---------graph/\temporary
c========== definitions for integrations
      psma=bfo(b,b0,(rmi+rma)/2.,r0,0.,z0)
      psma=psio(dpsi,ddpsi,dzpsi,rma,r0,z0,z0,am)
      psmi=psio(dpsi,ddpsi,dzpsi,rmi,r0,z0,z0,am)
      call eigenmod(b0oo,b0,dbobmax,dbobmaxs,dnonmax,grps2(1),rgrid(1),
     &     min,min_f,max,max_f,ntae,scalvavf,valfv,wtae,deltk,ximax)
      call rescondif(ai,alpha,b0,b0oo,eps,epsp,im1,iter,ntgr
     &     ,psiav,psitot,psma,psmi,pmu,nmu,qpf,npf,r0,rf,rma
     &     ,rmi,symm,t0,tip,valfv,v,nv,vi,z0,zi,wvpr,wqvpr,wd,wb
     &     ,xfow,xhp)
      call eigenmod(b0oo,b0,dbobmax,dbobmaxs,dnonmax,grps2(1),rgrid(1),
     &     min,min_f,max,max_f,ntae,scalvavf,valfv,wtae,deltk,ximax)
      do ipf=1,npf
         do imu=1,nmu
            yy(ipf,imu)=-real(ntae)*(wqvpr(1,imu,ipf,1)+wd(1,imu,ipf,1))
     &           +real(wtae)
         enddo
      enddo
      call twodgraf(qpf(1),yy(1,1),npf,nmu,col(1),pat(1),'white'
     $     ,2,'test','P_fi','w_d','orbit',6,6,1,wd,0)
      call frame(0)
      chi=-.06
      do iv=1,nv
         do imu=1,nmu
            q_mix(imu)=rf(1,5,imu,iv)*(1.-chi**2)
            yy(imu,4)=q_mix(imu)
            q_mix(imu)=psio(dpsi1,ddpsi1,dzpsi1,q_mix(imu),r0o,z0,z0,a1)
     &           -xfow*chi*xhp*v(iv)*q_mix(imu)*bfof(b,b0,q_mix(imu),r0
     &           ,z0,z0,1)
            i_mix(imu)=indx(qpf(1),npf,q_mix(imu))
            yy(1,3)=psiav(iv,imu,i_mix(imu),1)
            yy(2,3)=psiav(iv,imu,i_mix(imu)+1,1)
            yy(imu,1)=ynn(qpf(i_mix(imu)),yy(1,3),q_mix(imu))
            call distrfun(ai,b0*b0oo,yy(imu,2),dfalpdv(1,1),dfalpdpsi
     &       ,dpga,eps,yy(imu,1),gam,hnub(1,1),hh2,ihsps,iter,q_mix(imu)
     &           ,p0ga,pmu(imu),rf(1,5,imu,iv),1.,tb,tip,v(1),vi,v(nv),0
     &           .,z0,xhp,zi,dwbdpphi,wb_res,istart)
            yy(imu,2)=yy(imu,2)/1.e+9
            write(*,*) yy(imu,4),q_mix(imu),yy(imu,2)
         enddo
      enddo
      call twodgraf(yy(1,4),yy(1,2),nmu,1,col(1),pat(1),'white'
     $     ,2,'test','R_mjr','F_0','Initial_distribution',6,6,1,wd,0)
      call frame(0)
c      stop 'Stopped grid definition'
c
      write(*,*) ' max of Xi[m],dB/B,dB/B_sin and dn/n are',ximax
     &     ,dbobmax,dbobmaxs,dnonmax
      write(*,*) 'm range ',min,max,'ntae ',ntae,'Signif.m`s are',min_f
     &     ,max_f
      write(*,*) 'om= ',wtae,'f= ',real(wtae)/2./3
     &     .1415926
      write(*,*) ' Delta K= ',deltk,' Valfven=',valfv
      write(*,*) ' Va/Valfv=',vi/valfv,' xFOW=',xfow,' Type=',tip
      write(*,*) ' Alpha`s A,Z,V_0', ai,zi,vi
      write(*,*) 'B_0 (T)',b0,' R_ax',r0,' Z_ax',z0
      write(*,*) 'Psi_edge=',psitot
      b0ax=b0
      psedge=psma
      vioo=vi
      gamkom=16.18282*ai**2*vi**4*sign(abs(wtae),real(wtae))/xhp/zi
     &     /b0oo/b0/4./deltk/r0
      gamkomhlp=gamkom
      omstari=real(ntae)*xhp*vi*10.
      omstar=omstari/sign(abs(wtae),real(wtae))
      is1=1
      ipf1=1
      open(unit=18,err=10,file=filename,status='unknown')
      read(18,*,err=10,end=10) is1,ipf1
      read(18,*,err=10,end=10) gamkomis,omkomis,gtrap1
      read(18,*,err=10,end=10) gamintpfi
      read(18,*,err=10,end=10) omintpfi
      read(18,*,err=10,end=10) gtrap
      read(18,*,err=10,end=10) gberkintpfi
      read(18,*,err=10,end=10) gberkis
 10   continue
      close(18)
      itest=0
      call zero(dintgr(-mtl,1,1,1),nnres*(2*mtl+1)*8)
ctem
      if(itest.eq.1)then
         is1t=1
         is2t=1
         ipf1t=13
c        indx(qpf(1),npf,0.16)
         ipf2t=ipf1t
         imu1t=1
         imu2t=nmu
c To test for each bounce harmonics:
c make it =10 to see all the resonance condition or =1 for simpler
c output or it = 20 for m1 stability
         itestil=0
c To test for each mode harmonics:
         itestm=0
c allows to work only with harmonics at such numbers, which are the
c poloidal numbers, which is differs from Max M value printed!!!!
         im1t=1
         im2t=2
      else
         im1t=1
         im2t=2
         is1t=is1
         is2t=2
         ipf1t=ipf1
         ipf2t=npf
         imu1t=1
         imu2t=nmu
         itestil=0
         itestm=0
      endif
      is2t=1
      do is=is1t,is2t
         sigm=1.
         if(is.eq.2) sigm=-1.
ctem
         do iv=1,nv
            imlmax=0
            imtamax=0
            vii=v(iv)
            do imu=1,nmu
               pmui=pmu(imu)
               call zero(aeqv,nnres*(mtl*2+1))
               min1=min_f
               max1=max_f
               if(itestm.eq.1)then
                  min1=im1t
                  max1=min0(im2t,max_f)
               endif
               call zero(q_res(1,1),nnres*(mtl*2+1))
               do imtae=min1,max1
                  imtae0=imtae+1-min
c     hereafter the summ over the bounce harmonics
c
                  il=0
                  if(imtae.eq.min1) then
                     nres=0
                     ires0=0
                     if(itestil.eq.10)write(*,*)
     &                    ' V, <V_||>, W_d, W_b, W_res',imtae,il
                     do ipf=1,npf
                        if(ires0.le.i_mix(imu).and.abs(wb(iv,imu,ipf,is)
     &                       ).gt.1.e-10) then
                           ires0=ipf
                        endif
                        if(ires0.ne.0.and.abs(wb(iv,imu,ipf,is))
     $                       .gt.1.e-10) nres=nres+1
                        wres(ipf)=-ntae*(wqvpr(iv,imu,ipf,is)+wd(iv
     &                       ,imu,ipf,is))+real(wtae)+il*wb(iv
     &                       ,imu,ipf,is)
                        if(itestil.eq.10)then
                           write(*,*) qpf(ipf),wqvpr(iv,imu,ipf,is)
     &                             ,wd(iv,imu,ipf,is),wb(iv,imu,ipf,is
     &                          ),wres(ipf)
                        endif
                     enddo
                     nroot(il)=0
c     defining the resonance points
                     if(nres.gt.1) then
                        call roots(wres(ires0),qpf(ires0),ind(1),nres
     $                       ,0.,eps,0.,q_res(1,il),nroot(il))
                     endif
c     summation over the resonance points
                  endif
c     /\ end of definition of the resonance velocity for a give IL
c     \/ Start of resonance calculation for nroot(il) resonances
                  do ires=1,1
                     qpfi=q_mix(imu)
c     define critpoints for a given integrals
                     if(imtae.eq.min1) then
                        do i1=1,2
                           do i2=1,4
                              rfres(i1,i2,ires,il)=rf(i1,i2,imu,iv)
                           enddo
                        enddo
                        wvpr_res=wvpr(iv,imu,ind(ires),is)
                        psiavres(ires,il)=yy(imu,1)
                     endif
c     calculation of Gm -> Re Gm = fo(2); Im Gm = fo(3);
c     write(*,*) imu,imtae,il
c     aeqv is an indicator, which shows that average refers to trapped
c     particles if it is -100.
                     call xaverage(hhh1,beqv,iav,iter,eps,epsp,fo(1
     &                    ,imtae0,il,ires),imtae,ntae,ntgr,psma,psmi
     &                    ,qpfi,pmui,r0,rma,rmi,rvpar,rfres(1,1,ires,il)
     &                    ,sigm,symm,t0,tb,vii,z0,xhp
     &                    ,real(wtae),0.,wvpr_res)
                     if(tb.gt.1.e-10.and.imtae.eq.min1) then
                        if(itestil.eq.1.or.itestil.eq.10)then
                           write(*,*) 'a,b=',hhh1,beqv
                           write(*,*) 'at m,l=',imtae,il
     &                          ,'ind,q_res=',ind(ires)+ires0-1
     &                          ,q_res(ires,il),ires
c     
                        endif
                        aeqv(ires,il)=hhh1
c     if calculation is succesfull de rivative is next -> dwv, which is
c     eventually dwv == |dwv|
                        if(nroot(il).eq.0)then
                           call derivative(wres(ires0),qpf(ires0),nres
     &                          ,qpfi,dwv,real(wtae))
                           call sderiv(qpf(i_mix(imu)+1),wres(i_mix(imu)
     &                          +1),ddwv)
                        else
                           call derivative(wres(ires0),qpf(ires0),nres
     &                          ,q_res(ires,il),dwv,real(wtae))
                           call sderiv(qpf(ind(ires)+1),wres(ind(ires)+1
     &                          ),ddwv)
                        endif
                        dwv=dwv*omstar*vii
                        ddwv=ddwv*(omstar*vii)**2
c     following call is for equilibr. distr. function of alphas
c     derivatives dfalpdv,dfalpdpsi are normalized to V_al0,
c     Psi_edge
c     and are in cm^-3
                        call distrfun(ai,b0*b0oo,dfalp(ires,il)
     &                       ,dfalpdv(ires,il),dfalpdpsi,dpga,eps
     &                    ,psiavres(ires,il),gam,hnub(ires,il),hh2,ihsps
     &                       ,iter,qpfi,p0ga,pmui,rvpar,sigm,tb,tip
     &                       ,vii,vi,v(nv),0.,z0,xhp,zi,dwbdpphi,wb_res
     &                       ,istart)
                        dintgr(il,ires,4,2)=-1.01
                        if(itestil.eq.1.or.itestil.eq.10)then
                           write(*,*)'tb,qpfi,dwv,dfalpdv,dfalpdpsi'
     &                          ,',fo,psiaver',tb,qpfi,dwv,dfalpdv(ires
     &                          ,il),dfalpdpsi,fo(2,imtae0,il,ires),fo(3
     &                          ,imtae0,il,ires),psiavres(ires,il)
                        endif
                     endif
                  enddo
c        /\ end of 'as if' resonance integration
               enddo
c     /\ end of poloidal harmonic definition - imtae
c     we assemble here all quantities
c......................................................................
               xmlmax=0.
               hhh=0.
               imlmax=0.
               imtamax=0
               xmtamax=0.
               il=0.
               hhh11=0.
               do ires=1,1
c     do ires=1,min0(1,nroot(il)),1
                  if(abs(dintgr(il,ires,4,2)+1.01).lt.0.01)then
                     hhh1=0.
                     hhh2=0.
                     do imtae=min1+1-min,max1+1-min
                        hhh1=hhh1+fo(2,imtae,il,ires)
                        hhh2=hhh2+fo(3,imtae,il,ires)
                        hhh111=sqrt(fo(2,imtae,il,ires)**2+fo(3,imtae
     &                       ,il,ires)**2)
                        if(xmtamax.lt.hhh111) then
                           xmtamax=hhh111
                           imtamax=imtae-1+min
                        endif
                     enddo
                     hhh1=hhh1**2+hhh2**2
                     hhh11=0.
                     hhh2=0.
                     do imtae=min1+1-min,max1+1-min
                        hhh11=hhh11+fo(5,imtae,il,ires)
                        hhh2=hhh2+fo(4,imtae,il,ires)
                     enddo
                     hhh2=hhh2**2+hhh11**2
                     disp1=sqrt(db_b/dbobmax*hhh2*vii**3/abs(dwv
     &                    *real(wtae)))*real(ntae)*xhp*vi*10.
                     disp2=(db_b/dbobmax*3.*hhh2*vii**4/abs(ddwv)
     &                    /real(wtae)**2)**(1./3.)*real(ntae)*xhp*vi*10.
                     write(*,*)' ',yy(imu,4),q_mix(imu),sqrt(hhh2)
                     write(*,*)'Some = ',dwv,ddwv,q_res(1,0),disp1,disp2
c     berk
                     hhb=dfalpdv(ires,il)*hhh2
                  endif
               enddo
c /\ End of bounce harmonic and multiple resonances assembly
c......................................................................
c
            enddo
c /\ end of MU - pitch angle integration loop
         enddo
c /\ end of integration in velocity
      enddo
c /\ end of integration for co(+trapped) and contragoing particles. 
      call twodgraf(yy(1,4),yy(1,2),nmu,1,col(1),pat(1),'white'
     $     ,2,'test','R_mjr','F_0','Initial_distribution',6,6,1,wd,0)
      call frame(0)
      stop
c passing particles only without FLR and cutoff
      call twodgraf(qpf,gamintpfi(1,1,1),npf,1,col(1),pat(1),'white'
     $     ,2,'Passing_w/o_cutoff','sqrt(P_phi)'
     &     ,'d(gam/om)/dsqrt(P_phi)','FOW',6,6,1,wd,0)
      call frame(0)
c trapped particles only without FLR and cutoff
      call twodgraf(qpf,gtrap(1,1),npf,1,col(1),pat(1),'white'
     $     ,2,'Trapped_w/o_cutoff','sqrt(P_phi)'
     &     ,'d(gam/om)/dsqrt(P_phi)','FOW',6,6,1,wd,0)
      call frame(0)
      do i=1,npf
         gamintpfi(i,1,1)=abs(gamintpfi(i,3,1)-gtrap(i,3)+
     &        gamintpfi(i,3,2))*2.*abs(qpf(i)*gamkomhlp)
         gtrap(i,1)=abs(gtrap(i,3))*2.*abs(qpf(i)*gamkomhlp)
      enddo
      if(tip.eq.'i'.or.tip.eq.'h'.or.tip.eq.'r'.or.tip.eq.'x'.or.
     &     tip.eq.'m'.or.tip.eq.'e')return
c      stop ' stop in growth '
c passing particles only without FLR at cutoff
      call twodgraf(qpf,gamintpfi(1,1,1),npf,1,col(1),pat(1),'white'
     $     ,2,'Passing_cutoff','sqrt(P_phi)','d(gam/om)/dsqrt(P_phi)'
     &     ,'FOW',6,6,1,wd,0)
      call frame(0)
c trapped particles only without FLR at cutoff
      call twodgraf(qpf,gtrap(1,1),npf,1,col(1),pat(1),'white'
     $     ,2,'Trapped_cutoff','sqrt(P_phi)','d(gam/om)/dsqrt(P_phi)'
     &     ,'FOW',6,6,1,wd,0)
      call frame(0)
c      return
      do i=1,npf
         gamintpfi(i,1,1)=abs(gamintpfi(i,5,1)-gtrap(i,5)+
     &        gamintpfi(i,5,2))*2.*abs(qpf(i))
         gtrap(i,1)=abs(gtrap(i,5))*2.*abs(qpf(i))
      enddo
c      stop ' stop in growth '
c passing particles only without FLR after cutoff
      call twodgraf(qpf,gamintpfi(1,1,1),npf,1,col(1),pat(1),'white'
     $     ,2,'Passing_za_cutoff','sqrt(P_phi)','d(gam/om)/dsqrt(P_phi)'
     &     ,'FOW',6,6,1,wd,0)
      call frame(0)
c trapped particles only without FLR after cutoff
      call twodgraf(qpf,gtrap(1,1),npf,1,col(1),pat(1),'white'
     $     ,2,'Trapped_za_cutoff','sqrt(P_phi)','d(gam/om)/dsqrt(P_phi)'
     &     ,'FOW',6,6,1,wd,0)
      bmax=0.
      return
      end
c*********************************************************
c simple second derivative program
      subroutine sderiv(x,y,d2)
      dimension x(3),y(3)
      d2=-2.*((y(2)-y(1))/(x(2)-x(1))-(y(3)-y(2))/(x(3)-x(2)))/(x(3)-x(1
     &     ))
      return
      end
c*********************************************************
