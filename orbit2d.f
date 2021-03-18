c this subroutine gives the averaged drift frequencies as a function of P_phi
c You need to provide it with correct RM and R1 and with the 
c file 'wd****' with the grids and change the parameters here
      subroutine testorbitpf(rmi,rma,psmi,psma,tip,ai,b0oo,vi,zi,psitot
     &     ,vivalf)
      include 'orbit.grid'
c      parameter (nv=30,nmu=70,npf=70)
      parameter (nr=npf,nz=3,iav=3)
      dimension vpsi(nr,nz),vb(nr,nz),vd(nr,nz),vbv(nv,nz),vdv(nv,nz)
     &     ,fo(iav),rr(nr),rf(2,6,nmu,nv),pmu(nmu),qpf(npf),v0(nv)
     &     ,wvpr(nv,nmu,npf,2),wqvpr(nv,nmu,npf,2),wd(nv,nmu,npf,2)
     &     ,wb(nv,nmu,npf,2),rvpar0(nmu),psiav(nv,nmu,npf,2)
      character col(nz),pat(nz),tip
      external fwd
      common /kfow/xfow
      do i=1,nz
         col(i)='white'
         pat(i)='solid'
      enddo
      ntgr=60
      rm=rma
      z1=0.01
      zm=.85
      eps=1.e-5
      iter=2
      hz=(zm-z1)/(nr-1)
      v=1.
      vio=vi
      moo=1
      noo=1
      r1=(rmi+rma)/2.
      hr=bfo(b,b0,r1,r0,z1,z0)
      r1=2.
      hr=(rm-r1)/(nr-1)
c      rmi=rmi
c      rma=rma
      psmi=psio(dpsi,ddpsi,dzpsi,rmi,r0,z0,z0,am)
      psma=psio(dpsi,ddpsi,dzpsi,rma,r0,z0,z0,am)
      psedge=psma
      sigm=1.
      symm=1.
      iv=25
      qpf1=.15
      pmu1=r0
      epsp=1.e-4
      valfv=1.
c
      call rescondif(ai,alpha,b0,b0oo,eps,epsp,0,iter,ntgr,psiav,psitot
     $     ,psma,psmi,pmu,nmu,qpf,npf,r0,rf,rma,rmi,symm,t0,tip,valfv
     &     ,v0,nv,vi,z0,zi,wvpr,wqvpr,wd,wb,xfow,xhp)
      write(*,*) ' Psmi,ma=',psmi,psma,xfow
      write(*,*) ' Type ntgr,sigm,symm',ntgr,sigm,symm
      read(*,*) ntgr,sigm,symm
call t0xhp(alpha,t0,xhp,ai,b0,b0oo,vi,zi,psitot)
c
      xfow=.3
      nv1=23
      iv1=1
      nr1=1
      j1=4
      qpf1=qpf(34)
c
      do j=1,3
         if(j.eq.1)then
            nv1=29
         else if(j.eq.2)then
            nv1=20
         else
            nv1=13
         endif
         xv0=v0(nv1)
c         xv0=xv0*1.298/vi
         write(*,*) ' Mu=',pmu(j1),' r0=',r0,' v=',xv0,nv1
         write(*,*) ' Hereafter W_d & W_b'
         do i=1,nr
            rr(i)=z1+(i-1)*hz
            vb(i,3)=0.
            vd(i,3)=0.
            call critpnts(eps,iter,1,1,qpf(1),pmu(j1),psmi,psma,r0,rf(1
     &           ,1,j1,nv1),rmi,rma,xv0,xhp,z0)
            do ivv=1,1
               call wdif(alpha,b0,eps,epsp,iter,ntgr,r0,rmi,rma,pmu(j1)
     &              ,nr1,vd(i,3),psmi,psma,rr(i),1,rf(1,1,j1,nv1)
     &              ,rvpar0(1),symm,sigm,t0,vi,xv0,1,vd(i,1),vb(i,1)
     &              ,vb(i,2),vd(i,2),xhp,z0)
c                      wvpr____wqvpr___wd_____wb
c               vd(i,1)=4.*abs(vb(i,1)+vb(i,2))
            enddo
            write(*,*)'rf===',rf(1,1,j1,nv1),vb(i,1),vb(i,2),vd(i,2),i
         enddo
      enddo
c     write(*,*)((vb(i,2),vd(i,2),i),i=j1,j1-1+nr1)
      call twodgraf(rr,vd(1,1),nr,2,col(1),pat(1),'white'
     $     ,2,'test','R,m','W_t,W_b,*10^6_1/sec','orbit',6,6,1
     &     ,vd,0)
      call frame(0)
      call twodgraf(rr,vb(1,1),nr,2,col(1),pat(1),'white'
     $     ,2,'test','R,m','qW_t,W_d,*10^6_1/sec','orbit',6,6,1,vd,0)
      call frame(0)
      call grafinit(2)
      stop
      end
c**********************************************************************        
c     this subroutine gives the averaged drift frequencies as a function
c     of velocity
c You need to provide it with correct RM and R1 and with the 
c file 'wd****' with the grids and change the parameters here
      subroutine testorbitv(rmi,rma,psmi,psma,tip,ai,
     $     b0oo,vi,zi,psitot,vivalf)
c      parameter (nv=15,nmu=40,npf=70)
      parameter (nr=40,nv=15,nz=3,iav=3,nmu=nr,npf=70)
      dimension vpsi(nr,nz),vb(nr,nz),vd(nr,nz),vbv(nv,nz),vdv(nv,nz)
     &     ,fo(iav),rr(nr),rf(2,6,nmu,nv),pmu(nmu),qpf(npf),v0(nv)
     &     ,wvpr(nv,nmu,npf,2),wqvpr(nv,nmu,npf,2),wd(nv,nmu,npf,2)
     &     ,wb(nv,nmu,npf,2),rvpar0(nmu),psiav(nv,nmu,npf,2)
      character col(nz),pat(nz),tip
      external fwd
      common /kfow/xfow
      do i=1,nz
         col(i)='white'
         pat(i)='solid'
      enddo
      ntgr=20
      rm=.7
      z1=-0.1
      zm=.5
      eps=1.e-5
      iter=0
      hz=(zm-z1)/(nz-1)
      v=3.
      vio=vi
      moo=1
      noo=1
      r1=(rmi+rma)/2.
      hr=bfo(b,b0,r1,r0,z1,z0)
      r1=2.
      hr=(rm-r1)/(nr-1)
c      rmi=rmi
c      rma=rma
      psmi=psio(dpsi,ddpsi,dzpsi,rmi,r0,z0,z0,am)
      psma=psio(dpsi,ddpsi,dzpsi,rma,r0,z0,z0,am)
      psedge=psma
      sigm=1.
      symm1=1.
      symm=symm1
      iv=25
      write(*,*) ' Psmi,ma=',psmi,psma,xfow
      write(*,*) ' Type ntgr,sigm,symm',ntgr,sigm,symm
      read(*,*) ntgr,sigm,symm
      qpf1=.15
      pmu1=r0
      epsp=1.e-4
      valfv=1.
c
      call rescondif(ai,alpha,b0,b0oo,eps,epsp,0,iter,ntgr,psiav,psitot
     $     ,psma,psmi,pmu,nmu,qpf,npf,r0,rf,rma,rmi,symm1,t0,tip,valfv
     &     ,v0,nv,vi,z0,zi,wvpr,wqvpr,wd,wb,xfow,xhp)
call t0xhp(alpha,t0,xhp,ai,b0,b0oo,vi,zi,psitot)
c
      nv1=nv
      iv1=1
      nr1=1
c mu grid number
      j1=15
      qpf1=qpf(25)
c
      do j=1,1
c     do j=1,nz
         z=z1+(j-1)*hz
c     do i=3,5
         do i=1,nr
            rr(i)=pmu(i)
            vb(i,3)=0.
            vd(i,3)=0.
         enddo
         do ivv=1,nv
c         do ivv=12,14
         call wdif(alpha,b0,eps,epsp,iter,ntgr,r0,rmi,rma,rr(j1),nr1
     $           ,vdv(ivv,3),psmi,psma,qpf1,1,rf(1,1,j1,ivv),rvpar0,symm
     &           ,sigm,t0,vi,v0(ivv),1,vdv(ivv,1),vbv(ivv,1),vbv(ivv,2)
     &           ,vdv(ivv,2),xhp,z0)
         enddo
c                      wvpr____wqvpr___wd_____wb
c
         write(*,*)(vbv(i,2),vdv(i,2),i,i=iv1,iv1-1+nv1)
      enddo
c      do i=1,nr
c         vb(i,1)=abs(vb(i,1))
c         vb(i,2)=abs(vb(i,2))
c         vd(i,1)=abs(vd(i,1))
c         vd(i,2)=abs(vd(i,2))
c      enddo
c      write(*,*)vb
c      write(*,*)vd
c      return
      call twodgraf(v0(1),vbv(1,1),nv,2,col(1),pat(1),'white'
     $     ,2,'test','v/v_alp0','<qW_t>,W_d,*10^6_1/sec','orbit',6,6,1
     &     ,vd,0)
      call frame(0)
c      call chng(50)
      call twodgraf(v0(1),vdv(1,1),nv,nz,col(1),pat(1),'white'
     $     ,2,'test','v/v_alp0','W_b,<W_t>,<Psi>*10^6_1/sec','orbit',6,6
     &     ,1,vd,0)
      call frame(0)
      return
      end
c**********************************************************************
c shows the instant drift frequencies without bounce average
      subroutine testwd(rmi,rma,psmi,psma,ai,
     $     b0oo,vi,zi,psitot)
      parameter (nr=70,nz=70,iav=3,ntgr=20)
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      common/cgm/b0ax,brz,pmui,qpfi,z0i,imtae,imt
      common/asymm/asym /fast/xav,xavint,x2r,x2i,psioo,kfast
      common/wdr/alpha11,b011,psedge11,vio11,sigvpr1
      dimension vpsi(nr,nz),vb(nr,nz),vd(nr,nz),fo(iav),rr(nr)
      character col(nz),pat(nz)
      external fwd
      do i=1,nz
         col(i)='white'
         pat(i)='solid'
      enddo
      rm=rma
      z1=-1.
      zm=1.
      eps=1.e-5
      iter=0
      hz=(zm-z1)/(nz-1)
      v=1.
      vio=vi
      vio11=vi
      moo=0
      noo=1
      sigm=-1.
      sigvpr=sigm
      r1=(rmi+rma)/2.
      hr=bfo(b,b0,r1,r0,(z1+zm)/2.,z0)
      b011=b0
      r1=rmi
      hr=(rm-r1)/(nr-1)
c      rmi=rmi
c      rma=rma
      psmi=psio(dpsi,ddpsi,dzpsi,rmi,r0,z0,z0,am)
      psma=psio(dpsi,ddpsi,dzpsi,rma,r0,z0,z0,am)
      psedge=psma
      psedge11=psedge
      write(*,*) ' Psmi,ma=',psmi,psma
      qpf1=.1
      pmu1=r0
      epsp=1.e-4
      call t0xhp(alpha,t0,xhp,ai,b0,b0oo,vi,zi,psitot)
      write(*,*) 'alpha,t0,xhp',alpha,t0,xhp
      alpha11=alpha
      symm=1.
      imtae=1
c      istat=ishell('ps')
      do i=1,nr
         r=r1+(i-1)*hr
         rr(i)=r
c         do j=nr/2+10,nr/2+10
         do j=1,nz
            z=z1+(j-1)*hz
            if(z-z0.lt.0.) then
               asym=-1
            else
               asym=1.
            endif
            pmu=r1
c            qpf=z1
            pmui=pmu
            qpfi=qpf
            z0i=z0
            call  rvpar0(eps,qpf,pmu,rma,rmi,r0,rvpar,z0,zvpar0,iter)
c            write(*,*) 'rvpar0 ',eps,qpf,pmu,rma,rmi,r0,rvpar,z0,
c     $           zvpar0,iter
            kfast=-1
            call wdoo(alpha,b0,pmu,psma,qpf,r,rvpar,sigvpr,z,z0,vi,v
     &           ,wkprvpr,wdkpr,wd1)
cb0ax,brz,imtae,pmu,qpf,z0
            vb(i,j)=bfof(brz,b0ax,r,r0,z,z0,1)
c            call gmoo(b0,brz,1,pmu,qpf,r,vres,z,z0,vb(i,j)
c     &     ,gmoflr)
c            vb(i,j)=wd1*1.e+6
c            call fwd(fo,iav,r,z)
            vb(i,j)=wkprvpr
c            vb(i,j)=fo(1)
            vpsi(i,j)=wdkpr
c            vd(i,j)=fo(2)
            vd(i,j)=wd1
c            vpsi(i,j)=fo(3)
cfo(2)
c     
         enddo
      enddo
c      istat=ishell('ps')
c      write(*,*)vpsi
c      call twodgraf(rr,vpsi(1,1),nr,4,col(1),pat(1),'white'
c     $     ,2,'test','R,m','W_b,*10^6_1/sec','orbit',6,6,1,vd,0)
c      call frame(0)
      call contgraf(vb,nr,nz,r1,rm,z1,zm,2,0,6.,0,x,y,3,
     *     0,'W_vpr','R','z','BR',6,6)
      call frame(0)
      call contgraf(vpsi,nr,nz,r1,rm,z1,zm,1,0,6.,0,x,y,3,
     *     0,'W_qvpr','R','z','BR',6,6)
      call frame(0)
      call contgraf(vd,nr,nz,r1,rm,z1,zm,1,0,6.,0,x,y,3,
     *     0,'W_d','R','z','BR',6,6)
      call frame(0)
      return
      end
c*********************************************************
c if you plot, chech the numbers below such as r1
      subroutine orbittest(rmi,rma,psmi,psma)
      parameter (nr=201,nz=61)
      dimension rf(2,6)
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      dimension vpsi(nr,nz),vb(nr,nz),vd(nr,nz)
      character col(nz),pat(nz)
      do i=1,nz
         col(i)='white'
         pat(i)='solid'
      enddo
      rmi=rmi
      rma=rma
      rm=rma
      z1=.05
      zm=(rma-rmi)/4.
      eps=1.e-5
      iter=2
      hz=(zm-z1)/(nr-1)
c      xhp=1.98
      xhp=.1
      v=1.
      sigm=1.
      r1=(rmi+rma)/2
c compute the bounce point
      pmu=r1*0.1
      r1=bfo(b,b0,r1,r0,z1,z0)
      call rvpar0(eps,z1,pmu,rma,rmi,r0,rvpar1,z0,zvpar0,iter)
      qpf=z1
      rvpar=rvpar1
      call turnpnts(psma,psmi,rma,rmi,rf)
      pat(5)='dashed'
      pat(6)='dotted'
      write(*,*) 'rf,rvpar',rf,rvpar,'rmi,rma',rmi,rma
c      r1=rmi+.01
c      r1=7.3
      r1=pmu
      hr=(rm-r1)/(nr-1)
      nj=3
      do j=1,nj
      qpf=z1
      call rvpar0(eps,qpf,pmu,rma,rmi,r0,rvpar,z0,zvpar0,iter)
      call turnpnts(psma,psmi,rma,rmi,rf)
         do i=1,nr
c         do i=98,nr
c compute the point to investigate on the orbit
            r=r1+(i-1)*hr
            z=z1+(j-1)*(zm-z1)/real(max((nj-1),1))
            qpf=z
            vb(i,j)=bfo(b,b0,r,r0,z,z0)*r
            rvpar=rvpar1
            call rvpar0(eps,qpf,pmu,rma,rmi,r0,rvpar,z0,zvpar0,iter)
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,r
     $           ,r0,rvpar,sigm,v,z,z0,xhp)
c            if(j.eq.1) 
            vd(i,j)=r
            vb(i,j)=psior
            vb(i,j+5)=dpsior
            vb(i,j+1)=z
            vb(i,j+2)=psio(dpsi,ddpsi,dzpsi,r,r0,z0,z0,a)
            vb(i,j+2)=psiofpp2(r,r0,z0,a)
            if(abs(ddpsior).ge.30.) ddpsior=0.
            vb(i,j+10)=ddpsior
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,r
     $           ,r0,rvpar,-sigm,v,z,z0,xhp)
ckg            write(*,*) i,r,vb(i,j),psior,vb(i,j+1),z
            Vb(i,j+3)=psior
            vb(i,j+4)=z
            vb(i,j+6)=dpsior
            vb(i,j+7)=dpsi
            vb(i,j+8)=ddpsi
            if(abs(ddpsior).ge.30.) ddpsior=0.
            vb(i,j+9)=ddpsior
            qpf=z1+(i-1)*hz
            call rvpar0(eps,qpf,pmu,rma,rmi,r0,rvpar,z0,zvpar0,iter)
            call ppcorn(a,b,iconf,psma,psmi,rma,rmi,rf,sigm)
c            if(j.eq.1) 
c            write(*,*) 'i=',i,a,b,rvpar,qpf,iconf
            vd(i,j+1)=qpf
            vb(i,j+11)=a
            vb(i,j+12)=b
            vb(i,j+13)=iconf
         enddo
         write(*,*) 'rf,rvpar',rf,rvpar,'nr',nr
         call twodgraf(vd(1,1),vb(1,j),nr,5,col(1),pat(1),'white'
     $        ,2,'test','R,m','Psi,Z','orbit',6,6,1,vd,0)
         call frame(0)
         call twodgraf(vd(1,1),vb(1,j+5),nr,3,col(1),pat(1),'white'
     $        ,2,'test','R,m','dPsi(R,Z=0)/dR,dPsi_orb/dR','orbit',6,6,1
     &        ,vd,0)
         call frame(0)
         call twodgraf(vd(1,1),vb(1,9),nr,3,col(1),pat(1),'white'
     $        ,2,'test','R,m','d2Psi(R,Z=0)/dR2,d2Psi_orb/dR2','orbit',6
     &        ,6,1,vd,0)
         call frame(0)
      enddo
      qpf=z1
      call ppcorn(a,b,iconf,psma,psmi,rm,2.,rf,sigm)
      write(*,*) a,b,rvpar,qpf
c to run this part you need to have 2 DO loops upstairs c call
contgraf(vb,nr,nz,r1,rm,z1,zm,1,0,6.,0,x,y,3, 
c *0,'BR','R','z','BR',6,5)
c      return
      call twodgraf(vd(1,nj+1),vb(1,nj+11),nr,3,col(1),pat(1),'white'
     $     ,2,'test','P_phi','R,m','orbit',6,6,1,vd,0)
      call frame(0)
      call grafinit(2)
      stop
      end
c*********************************************************
c define arrays for NOVA-K code, qpf, pmu, and v
      subroutine rescondif(ai,alpha,b0,b0oo,eps,epsp,im1,iter,ntgr,psiav
     &     ,psitot,psma,psmi,pmu,nmu,qpf,npf,r0,rf,rma,rmi,symm,t0,tip1
     &     ,valfv,v,nv,vi,z0,zi,wvpr,wqvpr,wd,wb,xfow,xhp)
      dimension qpf(npf),pmu(nmu),v(nv),wvpr(nv,nmu,npf,2),
     $     wqvpr(nv,nmu,npf,2),wd(nv,nmu,npf,2),wb(nv,nmu,npf,2)
     $     ,rf(2,6,nmu,nv),rvpar0(150),psiav(nv,nmu,npf,2)
c      external fwd
common /kfow/xfow1
      character fow*8,tip,tip1,filename*10
      tip=tip1
      if(tip1.eq.'d')tip='s'
      is1=0
      ipf1=0
      xfow1=0.
      krf=0
      write(fow(1:4),'(f4.2)')xfow
      write(fow(6:8),'(i3)')npf
c      write(*,*) fow,xfow,npf
      call t0xhp(alpha,t0,xhp,ai,b0,b0oo,vi,zi,psitot)
      if(ai.eq.4.)then
         filename='wd'//fow(7:8)//'.'//fow(1:1)//fow(3:4)//tip//'a'
      else if(ai.eq.1.)then
         filename='wd'//fow(7:8)//'.'//fow(1:1)//fow(3:4)//tip//'h'
      else if(ai.eq.2.)then
         filename='wd'//fow(7:8)//'.'//fow(1:1)//fow(3:4)//tip//'d'
      else if(ai.eq.3.)then
         filename='wd'//fow(7:8)//'.'//fow(1:1)//fow(3:4)//tip//'t'
      else 
         filename='wd'//fow(7:8)//'.'//fow(1:1)//fow(3:4)//tip//'?'
         write(*,*)
     &        '?????????warning:  Your particle weight is uncommon'
      endif
      write(*,*) '2open'
      open(unit=11,err=10 ,file=filename,form='unformatted',status
     &     ='unknown')
c      write(*,*) 'opened'
      read(11,err=10,end=10) nv1,npf1,nmu1,alpha1,t01,xhp1,xfow1,symm1
c      write(*,*) 'read1',nv1,npf1,nmu1,nv,npf,nmu
      read(11,err=10,end=10) krf,rf,psiav
c      write(*,*) 'read2'
      read(11,err=10,end=10) v,pmu,qpf
c      write(*,*) 'read3'
      read(11,err=10,end=10) wvpr,wqvpr,wd,wb
c      write(*,*) 'read4'
      read(11,err=10,end=10) is1,ipf1
c      write(*,*) 'read5'
ctem
ctem      is1=2
ctem      ipf1=1
      if(is1*ipf1.lt.2*npf) goto 10
      if(nv1.ne.nv.or.npf1.ne.npf.or.nmu1.ne.nmu.or.
     $     abs(alpha-alpha1)+abs(t0-t01).gt..1.or.xfow1.ne.xfow.or
     &     .symm1.ne.symm.or.abs(xhp-xhp1).ge.1.e-3*abs(xhp+xhp1)) then
         write(*,*) is1*ipf1,nv1,npf1,nmu1,alpha1,t01,xfow1,symm1,xhp1
         write(*,*) 2*npf,nv,npf,nmu,alpha,t0,xfow,symm,xhp
         is1=0
         ipf1=0
         goto 10
      endif
      close(11)
      write(*,*) '1 alpha,t0,xhp',alpha,t0,xhp,' FOW=',xfow
ctem
      return
 10   continue
      krf=0
      write(*,*) 'is1,ipf1,npf= ',is1,ipf1,npf
      write(*,*) '1 alpha,t0,xhp',alpha,t0,xhp,' FOW=',xfow
ctem      return
      xfow1=xfow
      symm1=symm
      is1=max(is1,1)
      ipf1=max(ipf1,1)
      close(11)
      write(*,*) 'didn`t read satisfact.'
      imudk=2
      pmudk=0.6
      if(imudk.eq.1.or.imudk.eq.2)pmudk=1.
c pitch angle grid definition
      pmu(nmu)=rma-.02
      pmu(1)=.02
cnng14      if(tip.eq.'h'.or.tip.eq.'i')then
c make it more non ICRH like sl. down DF 
cnng16      if(tip.eq.'i')then
      if(tip.eq.'h')then
         pmu(1)=rmi*1.03
         pmu(1)=rmi*0.9
ckg         pmu(nmu)=rma-0.2
         pmu(nmu)=rma*0.98
      endif
      if(im1.eq.1)then
cben this is used to include only trapped particle effects
c         pmu(1)=rmi
      else if(im1.eq.10)then
         pmu(1)=r0
         pmu(nmu)=(rma+r0)/2.
      endif
      if(imudk.eq.1)then
         pmu(1)=0.005*r0
         pmu(nmu)=1.5*r0
      endif
      do imu=2,nmu-1
         pmu(imu)=pmu(1)+(pmu(nmu)-pmu(1))*
     $        (real(imu-1)/real(nmu-1))**pmudk
      enddo
      write(*,*) ' Grid in Mu is',pmu
c velocity grid definition
      if(tip.eq.'s'.or.tip.eq.'r'.or.tip.eq.'t'.or.tip.eq.'j'.or.tip.eq.
     ^     'c'.or.tip.eq.'g'.or.tip.eq.'l'.or.tip.eq.'b'.or.tip.eq.'.'
     ^     )then
         vmax=1.07
      else
         vmax=max(valfv/vi,1.)*4.
      endif
      pk=0.5
      if(im1.eq.0.or.im1.eq.2.or.im1.eq.3)then
         call inistr(v,nv,.1,vmax,pk)
      else if(im1.eq.1)then
cben         call ini(v,nv,.3,vmax)
         call inistr(v,nv,.1,vmax,pk)
      else if(im1.eq.10)then
         v(1)=vi
      endif
      write(*,*) ' Grid in V is',v
c Pphi grid definition
      if(im1.ne.10)then
         qpf(1)=-sqrt(r0*abs(-pmu(1)+r0))*xfow*xhp*v(nv)
         qpf(npf)=sqrt(rmi*abs(-pmu(1)+rmi))*xfow*xhp*v(nv)
     $        +psmi
      endif
c This stuff is to speed up the process
      if(imudk.eq.1)then
         qpf(1)=-1.
         qpf(npf)=2.
      endif
      if(im1.eq.0.or.im1.eq.3)then
c         if(tip.eq.'r')then
c            call ini(qpf,npf,qpf(1)+(qpf(npf)*.7-qpf(1))*3./7.,qpf(1)
c     &           +(qpf(npf)*.7-qpf(1))*4.5/7.)
c         else
            call ini(qpf,npf,qpf(1),qpf(npf)*pmudk)
c         endif
      else if(im1.eq.1.or.im1.eq.2)then
         call ini(qpf,npf,qpf(1),qpf(npf))
      else if(im1.eq.10)then
         qpf(1)=psio(dpsi1,ddpsi1,dzpsi1,pmu(1),r0,z0,z0,a1)
         qpf(npf)=.5
         call ini(qpf,npf,qpf(1),qpf(npf))
      endif
ctem       krf=0
      if(krf.eq.0) then
         call critpnts(eps,iter,nmu,nv,qpf(1),pmu,psmi
     $     ,psma,r0,rf,rmi,rma,v,xhp,z0)
         krf=1
         open(unit=11,file=filename,form='unformatted',status='unknown')
         write(11) nv,npf,nmu,alpha,t0,xhp,xfow,symm
         write(11) krf,rf,psiav
         write(11) v,pmu,qpf
         write(11) wvpr,wqvpr,wd,wb
         write(11) is1,ipf1
         close(11)
      endif
ctem
      do is=is1,2
c      do is=1,1
ctem
         do ipf=ipf1,npf
c         inpf=indx(qpf(1),npf,0.14)
c         do ipf=inpf,inpf
c         do ipf=31,31
         sigm=1.
         if(is.eq.2) sigm=-1.
         if(ipf.eq.1.and.is.eq.1) then
            do iii=1,nmu
               do iiii=1,nv
                  rf(2,6,iii,iiii)=rf(2,3,iii,iiii)
                  rf(1,5,iii,iiii)=-1.
                  rf(1,6,iii,iiii)=rf(1,3,iii,iiii)
               enddo
            enddo
         endif
c         write(*,*) 'ipf=',ipf
         call wdif(alpha,b0,eps,epsp,iter,ntgr,r0,rmi,rma,pmu,
     $        nmu,psiav(1,1,ipf,is),psmi,psma,qpf(ipf),1,rf,rvpar0(1)
     &        ,symm,sigm,t0,vi,v,nv,wvpr(1,1,ipf,is),wqvpr(1,1,ipf,is)
     &        ,wd(1,1,ipf,is),wb(1,1,ipf,is),xhp,z0)
         if(ipf.eq.1.and.is.eq.1) then
            do iii=1,nmu
               do iiii=1,nv
                  rf(1,5,iii,iiii)=rvpar0(iii)
               enddo
            enddo
         endif
         if(ipf.eq.npf.or.mod(ipf,30).eq.0) then
            open(unit=11,file=filename,form='unformatted',status
     &           ='unknown')
         write(11) nv,npf,nmu,alpha,t0,xhp,xfow,symm
         write(11) krf,rf,psiav
         write(11) v,pmu,qpf
         write(11) wvpr,wqvpr,wd,wb
         write(11) is,ipf+1
ctem         write(11) is,npf+1
         close(11)
         write(*,*) 'I just calculated Sigm=',sigm,'case ipf=',ipf
         endif
      enddo
         ipf1=1
      enddo
c write precessional frequencies for Joe Snipes
      open(unit=11,file='phidot'//filename,status
     &           ='unknown')
      write(11,*) 'tables of <phi> dot in 10^6 sec'
      iv=max(indx(v,nv,1.),1)
      do ipf=ipf1,npf,3
         write(11,*)'P_phi =',qpf(ipf)
         write(11,*)'      lambda,      R,      <psi>,     <phi> dot'
         do imu=1,nmu
            write(11,'(4e14.5)')pmu(imu)/r0,pmu(imu),psiav(iv,imu,ipf,1
     &           ),wd(iv,imu,ipf,1)
         enddo
      enddo
      close(11)
      return
ctem      stop
      end
c**********************************************************************        
c     Output are
c
c     dpsior= dPsi/dlnR, ddpsior= RRd^2Psi/dR^2,
c     psior= Psi, z
c
c     Inputs are eps,iter,qpf,pmu,r,r0,rvpar,sigm,v,z0,xhp
c
      subroutine psorb(b,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,r,r0
     $     ,rvpar,sigm,v,z,z0,xhp)
c      write(*,*)'psorb',b,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,r,r0
c     $     ,rvpar,sigm,v,z,z0,xhp
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
      if(iter.gt.1.and.r.gt.r0*.85)then
         psior=max(qpf,0.)+dpsior*pss/2.
      else
         psior=qpf+dpsior*pss
      endif
      z=zatpsir(eps,psior,r,r0,z0)
      hp=bfof(b,b0,r,r00,z,z00,-1)
      ddpsior=dpsior
c      write(*,*)'1 psorb:: ',b,b0,r,r00,z,z00,ddpsior
cvac change this line if vacuum included
      if(b0.lt..1) then
         psior=1.
         return
      endif
c Ne znayu kakogo khrena ya eto vstavil
ckg      if(abs(rvpar-r)/r.lt.((b*r/b0/r00)-1.)**2) goto 10
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
c         z=zatpsir(eps,psior,r,r0,z0)
         if(z.eq.z0)then
            z=zatpsirf(eps,psior,r,r0,z0,hp,1.)
         else
            z=zatpsirf(eps,psior,r,r0,z0,hp,.5)
         endif
c         write(*,*) '3 psorb:: ',psior,z,pss,b,b0
      enddo
 10   if(pss.lt.1.e-3) pss=1.e-3
c      write(*,*)'4',pss,b0,r0
      ddpsior=-dpsior*(b*pmu/b0/r0)**2/4./pss**3
      dpsior=dpsior*(pss+b*pmu/b0/r0/2./pss)
c      write(*,*)'5',pss
      return
      end
c*********************************************************
c     iav is number of averaged functions
c     ntgr is number of integration grid ~10
c     iter is iteration for orbit evaluation ~2-3
c
c
c
c
c     rvpar is output
c
      subroutine xaverage(a,b,iav,iter1,eps1,epsp,fo,mtae,ntae
     $     ,ntgr,psma,psmi,qpf1,pmu11,r01,rma,rmi,rvpar1,rf,sig1
     $     ,symm,t0,tb,v1,z01,xhp1,wtae,wb_res,wres,wvpr_res)
      dimension fo(iav),fo1(8),fo2(6),rf(2,6)
      common/cgm/b0ax,bm,pmui,qpfi,z0i,imtae,imt
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      common/asymm/asym /fast/xav,xavint,x2r,x2i,psioo,kfast
      common/wdr/alpha,b0axx,psedge,vio,sigvpr
      data pi/3.1415926/
      save zvpar,pmusav,qpfsav
      pmu1=pmu11
 13   tb=0.
      tb1=0.
      tb2=0.
      wdw=0.
      wdw6=0.
      wdwass=0.
      const=1.
      b0ax=b0axx
      ntgr2=ntgr
      call zero(fo(1),iav)
      call zero(fo1(1),8)
      call zero(fo2(1),6)
      xav1=0.
      xav=0.
      xavint=0.
      x2r1=0.
      x2r=0.
      x2i=0.
      if(pmu1.ge.rma) return
      asym=1.
      iter=iter1
      eps=eps1
      qpf=qpf1
      pmu=pmu1
      r0=r01
      z0=z01
      if(pmu.ne.pmusav.or.qpf.ne.qpfsav)then
         pmusav=pmu
         qpfsav=qpf
         call rvpar0(eps,qpf,pmu,rma,rmi,r0,rvpar,z0,zvpar0,iter)
         rf(1,3)=amax1(rf(1,6)-rf(1,5)+rvpar,rvpar)
      endif
      rvpar1=rvpar
      sig=sig1
      v=v1
      xhp=xhp1
      call ppcorn(a,b,iconf,psma,psmi,rma,rmi,rf,sig1)
      if(rvpar.gt.a.or.rvpar.gt.b.or.a.eq.b)return
c      write(*,*) 'xxx',a,b,qpf,pmu,v
      sig=sig1
      sigvpr=sig1
      if(a.lt.rmi.or.b.lt.rmi.or.kenter_pest(a,b).eq.1) return
      wvpr=0.
      if((a-r0)*(b-r0).ge.0.) then
         wvpr=wvpr_res
      else
         wvpr=wvpr_res-sign(wb_res,wvpr_res)
      endif
c Negative passing particles are averaged over the time here
      if(sig1.lt.0.)then
         const=(-1)**mtae
         ntgr2=ntgr
 111     call xsummation(b,a,const,mtae,ntae,ntgr2,pi,eps,epsp
     &        ,iter,qpf,pmu,r0,rvpar,sig1,symm,v,z0,xhp,iav,t0,tb
     &        ,fo,wtae,warn,wres,wvpr)
         xav1=xav
         x2r1=x2r
         if(warn.eq.1..and.ntgr2.lt.ntgr*6)then
            ntgr2=ntgr2*3
            tb=0.
            do i=1,6
               fo(i)=0.
            enddo
            xav=0.
            x2r=0.
            goto 111
         else if(warn.eq.1)then
            tb=0.
            do i=1,6
               fo(i)=0.
            enddo
            xav=0.
            x2r=0.
            return
         endif
c     write(*,*) 'Up tb,fo',tb,fo
c to be changed for antisymmetric case
         if(symm.eq.0.) then
            asym=-1.
            tb1=tb
            do i=1,iav
               fo1(i)=fo(i)*tb
            enddo
            xav=xav*tb
            x2r=x2r*tb
            call xsummation(a,b,const,mtae,ntae,ntgr2,pi,eps,epsp
     &           ,iter,qpf,pmu,r0,rvpar,sig1,symm,v,z0,xhp,iav,t0
     &           ,tb1,fo1(1),wtae,warn,wres,wvpr)
c            write(*,*) 'Down tb,fo',tb,fo
            asym=1.
            if(tb1.eq.0..or.tb.eq.0.)then
               tb=0.
               return
            else
               tb=tb1
            endif
            call equiv(fo(1),fo1(1),6)
         else
            tb=2.*tb
         endif
         return
      endif
      amh=epsp*rvpar
      fpsi=psio(dpsi,ddpsi,dzpsi,rvpar,r0,z0,z0,am)
      if(fpsi.gt.qpf.or.abs(a-rvpar).le.amh.or.rvpar.le.rmi) then
c
c this part for copassing particles
c
c     The main part of copassing particle contribution
corb 
c         write(*,*) ' am,b,amh,symm=',am,b,amh,symm
         if(abs(a-rvpar).le.amh) then
            if(a.lt.b)then
               am=a+amh
            else
               am=a-amh
            endif
         else
            am=a
         endif
         sigvpr=sig1
         if(am.lt.b-amh*2.) then
            ntgr2=ntgr
            const=1.
 112        continue
            call xsummation(am,b,const,mtae,ntae,ntgr2,pi,eps,epsp
     &           ,iter,qpf,pmu,r0,rvpar,1.,symm,v,z0,xhp,iav,t0
     &           ,tb,fo(1),wtae,warn,wres,wvpr)
            if(warn.eq.1..and.ntgr2.lt.ntgr*6)then
               ntgr2=ntgr2*3
               tb=0.
               call zero(fo(1),iav)
               xav=0.
               x2r=0.
               goto 112
            else if(warn.eq.1)then
               tb=0.
               do i=1,6
                  fo(i)=0.
               enddo
               xav=0.
               x2r=0.
               return
            endif
c            write(*,*) '2am,b,epsp,amh=',am,b,epsp,amh
            if(abs(a-rvpar).le.amh.and.tb.ne.0.) then
               if(fpsi.gt.qpf)then
                  sigvpr=sig1
               else
                  sigvpr=1.
               endif
               call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,am
     $              ,r0,rvpar,sig1,v,zvpar,z0,xhp)
               tb1=abs(zvpar-z0)*2.*bm/b0*am*t0/xhp
               kfast=-1
               tb2=tb1+tb
               call fwd(fo1(1),3,a+amh/2.,amh/2.)
               ss=mtae*(fo1(1)-wvpr)
               cs=ntae*(fo1(2)+fo1(3))-wtae+wres
               fo1(1)=ss-cs
               fo1(6)=-ss-cs
               fo(1)=(fo(1)*tb+fo1(1)*tb1)/tb2
               fo(6)=(fo(6)*tb+fo1(6)*tb1)/tb2
               kfast=1
               call fgm(fo1(1),8,a+amh/2.,symm,v1,amh/2.)
               ss=const*tb1
               if(symm.ne.0.) then
c                  cs=ss*(cos(fo(1))+cos(fo(6)))
c                  ss=ss*(sin(fo(1))-sin(fo(6)))
                  cs=ss*cos(fo(1))
                  ss=ss*sin(fo(1))
c     FLR included
                  fo(4)=(fo(4)*tb+fo1(2)*cs+fo1(4)*ss)/tb2
                  fo(5)=0.
c     FLR not included
                  fo(2)=(fo(2)*tb+fo1(1)*cs+fo1(3)*ss)/tb2
                  fo(3)=0.
                  tb=tb2
               else
c     Change it for nonsymetric case
                  cs=ss*cos(fo(1))
c                  cs1=ss*cos(fo(6))
c                  ss1=ss*sin(fo(6))
                  cs1=0.
                  ss1=0.
                  ss=ss*sin(fo(1))
c     FLR included
                  fo(4)=(fo(4)*tb+fo1(2)*(cs+cs1)+fo1(4)*(ss-ss1)
     &                 +fo1(6)*(cs-cs1)-fo1(8)*(ss+ss1))/tb2
                  fo(5)=(fo(5)*tb+fo1(2)*(ss+ss1)+fo1(4)*(cs1-cs)
     &                 +fo1(6)*(ss-ss1)+fo1(8)*(cs1+cs))/tb2
c     FLR not included
                  fo(2)=(fo(2)*tb+fo1(1)*(cs+cs1)+fo1(3)*(ss-ss1)
     &                 +fo1(5)*(cs-cs1)-fo1(7)*(ss+ss1))/tb2
                  fo(3)=(fo(3)*tb+fo1(1)*(ss+ss1)+fo1(3)*(cs1-cs)
     &                 +fo1(5)*(ss-ss1)+fo1(7)*(cs1+cs))/tb2
                  tb=tb2
c     
                  asym=-1.
                  call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,
     $                 psior,am,r0,rvpar,sig1,v,zvpar,z0,xhp)
                  tb1=abs(zvpar-z0)*2.*bm/b0*am*t0/xhp
                  kfast=-1
                  call fwd(fo1(1),3,a+amh/2.,amh/2.)
                  tb2=tb1+tb
                  ss=mtae*(fo1(1)-wvpr)
                  cs=ntae*(fo1(2)+fo1(3))-wtae+wres
                  fo1(1)=ss-cs
                  fo1(6)=-ss-cs
                  fo(1)=(fo(1)*tb+fo1(1)*tb1)/tb2
                  fo(6)=(fo(6)*tb+fo1(6)*tb1)/tb2
                  kfast=1
                  call fgm(fo1(1),8,a+amh/2.,symm,v1,amh/2.)
                  ss=const*tb1
                  cs=ss*cos(fo(1))
c                  cs1=ss*cos(fo(6))
c                  ss1=ss*sin(fo(6))
                  cs1=0.
                  ss1=0.
                  ss=ss*sin(fo(1))
c     FLR included
                  fo(4)=(fo(4)*tb+fo1(2)*(cs+cs1)+fo1(4)*(ss-ss1)
     &                 +fo1(6)*(cs-cs1)-fo1(8)*(ss+ss1))/tb2
                  fo(5)=(fo(5)*tb+fo1(2)*(ss+ss1)+fo1(4)*(cs1-cs)
     &                 +fo1(6)*(ss-ss1)+fo1(8)*(cs1+cs))/tb2
c     FLR not included
                  fo(2)=(fo(2)*tb+fo1(1)*(cs+cs1)+fo1(3)*(ss-ss1)
     &                 +fo1(5)*(cs-cs1)-fo1(7)*(ss+ss1))/tb2
                  fo(3)=(fo(3)*tb+fo1(1)*(ss+ss1)+fo1(3)*(cs1-cs)
     &                 +fo1(5)*(ss-ss1)+fo1(7)*(cs1+cs))/tb2
                  asym=1.
               endif
            endif
c     For nonsymetric case
            if(symm.eq.0.)then
               do i=1,iav
                  fo(i)=fo(i)*tb
               enddo
               xav=xav*tb
               x2r=x2r*tb
               asym=-1.
c               write(*,*) '3am,b,epsp,amh=',am,b,epsp,amh
               call xsummation(b,am,const,mtae,ntae,ntgr2,pi,eps
     &              ,epsp,iter,qpf,pmu,r0,rvpar,1.,symm,v,z0,xhp
     &              ,iav,t0,tb,fo,wtae,warn,wres,wvpr)
c               write(*,*)' Here=',tb,fo,wdwass
               asym=1.
               if(tb.eq.0.)return
               tb=abs(tb)
            else
               tb=2.*tb
            endif
c            write(*,*) ' after tb,fo1,wdw=',tb,fo1,wdw
            return
         else if(abs(a-b).gt.amh*.1)then
c     may be improved for the bounce time of well trapped particles
            if(pmu1.eq.pmu11.and.pmu11.gt.r0+amh*2.)
     $           then
               pmu1=amax1(pmu11-amh*5./r0,0.)
               write(*,*)' more step',pmu1,pmu11
               goto 13
            endif
            write(*,*)' more step',pmu1,pmu11,a,b,amh
            r=(a+b)/2.
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior
     $           ,r,r0,rvpar,sig1,v,zvpar,z0,xhp)
            if(zvpar-z0.ge.0.) then
               zvpar=z0+max(a*epsp,zvpar-z0)
            else
               zvpar=z0+min(-a*epsp,zvpar-z0)
            endif
            tb=psio(dpsi,ddpsi,dzpsi,r,r0,zvpar,z0,am)
            tb=sqrt(1.-rvpar/r)*abs(dzpsi)
            tb=max(.001,tb)
            tb=bm*r/b0/v/tb*abs(a-b)*t0*2.
            kfast=-1
            call fgm(fo1(1),8,r,symm,v1,amh/2.)
            xav=-1000.
            x2r=-1000.
cpmu
cpsioo
            return
         else
            return
         endif
      else
c this part for trapped particles
         sigvpr=1.
         const=1.
         call xsummation(rvpar,b,const,mtae,ntae,ntgr,pi,eps,epsp
     &        ,iter,qpf,pmu,r0,rvpar,1.,symm,v,z0,xhp,iav,t0,tb,fo1
     &        ,wtae,warn,wres,wvpr)
         xav=xav*tb
         x2r=x2r*tb
         do i=1,iav
            fo(i)=fo1(i)*tb
         enddo
c      write(*,*) 'bef.fo1,tb ',(fo1(i),i=1,iav),tb
         sigvpr=-1.
         ntgr2=ntgr
         tb2=tb
         if(tb.eq.0.)return
 113     call xsummation(a,rvpar,const,mtae,ntae,ntgr2,pi,eps,epsp
     &        ,iter,qpf,pmu,r0,rvpar,-1.,symm,v,z0,xhp,iav,t0,tb2
     &        ,fo(1),wtae,warn,wres,wvpr)
         if(tb2.eq.0.)then
            tb=0.
            return
         endif
         if(warn.eq.1..and.ntgr2.lt.ntgr*6)then
            ntgr2=ntgr2*3
            tb2=tb
            do i=1,iav
               fo(i)=fo1(i)*tb
            enddo
            xav=xav*tb
            x2r=x2r*tb
            goto 113
         else if(warn.eq.1)then
            tb=0.
            do i=1,6
               fo(i)=0.
            enddo
            xav=0.
            x2r=0.
            return
         endif
         tb=tb2
c      write(*,*) 'bef.fo,tb ',(fo(i),i=1,iav),tb
         if(symm.eq.0.)then
c asymmetric case
            do i=1,6
               fo1(i)=fo(i)*tb
            enddo
            xav=xav*tb
            x2r=x2r*tb
            tb1=tb
            asym=-1.
            sigvpr=-1.
            call xsummation(rvpar,a,const,mtae,ntae,ntgr,pi,eps
     &           ,epsp,iter,qpf,pmu,r0,rvpar,-1.,symm,v,z0,xhp,iav
     &           ,t0,tb1,fo1(1),wtae,warn,wres,wvpr)
            tb=abs(tb1)
            sigvpr=1.
            do i=1,6
               fo(i)=fo1(i)*tb
            enddo
            xav=xav*tb
            x2r=x2r*tb
c 
            call xsummation(b,rvpar,const,mtae,ntae,ntgr,pi,eps
     &           ,epsp,iter,qpf,pmu,r0,rvpar,1.,symm,v,z0,xhp,iav
     &           ,t0,tb,fo(1),wtae,warn,wres,wvpr)
c 
            asym=1.
            if(tb1.eq.0..or.tb.eq.0.)then
               tb=0.
               return
            endif
         else
            tb=tb*2.
         endif
         a=-100.
c         write(*,*) '-----tb,fo= ',tb,fo
c         write(*,*) fo,tb
      endif
      return
      end
c*********************************************************
      subroutine xsummation(a,b,const,mtae,ntae,ntgr1,pi,eps,epsp,
     $     iter,qpf,pmu,r0,rvpar,sig,symm,v,z0,xhp,iav,t0,tb,fo
     &     ,wtae,warn,wres,wvpr)
      common/cgm/b0,bm,pmui,qpfi,z0i,imtae,imt
      dimension ff(8),fo(iav),rf(270)
      common /fast/xav,xavint,x2r,x2i,psioo,kfast
      common/saveij/ik,jk,cl(2),cu(2)
      save ntgr0,rf
      warn=0.
c      ntgr=abs(a-b)/abs(a+b)*4*ntgr1
      ntgr=max(ntgr1,8)
      if(ntgr0.ne.ntgr) then
         ntgr0=ntgr
         do i=1,ntgr
            rf(i)=cos((real(2*i)-1)*pi/real(2*ntgr))
         enddo
      endif
c      write(*,*) 'a,b,rvpar= ',a,b,rvpar
      deli=pi/real(ntgr)
      delr=.5*(b-a)
      elr=.5*(a+b)
c
      htb0=t0*pi/real(ntgr)/b0/v
      do 11 i=1,ntgr
         r=elr+delr*rf(i)
copt for optimization delete calculation of derivatives
         if(i.le.2)then
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior
     $           ,r,r0,rvpar,sig,v,z,z0,xhp)
         else
            call psorbf(bm,b0,eps,iter,qpf,pmu,psior,r,r0
     $           ,rvpar,sig,v,z,z0,xhp,zh)
         endif
         if(k_pest(psior).eq.1.or.(ik.eq.1.and.k_pest(psior).eq.-1)
     &        .or.psior.ge.1.)then
            tb=0.
            warn=0.
            return
         endif
         zh=z
         kfast=1
         psioo=psior
         call fwdf(ff(1),3,r,z)
         htb=psiof(dpsi,ddpsi,dzpsi,r,r0,z,z0,am,1)
         if(r0.lt..01.or.b0.lt..01) goto 11
         pss=1.-pmu*bm/r0/b0
c
         if(pss.gt.5.*r*epsp) then
            htb1=min(abs(rvpar-a),abs(rvpar-b))
            htb=min(abs(r-a),abs(b-r))
            if(htb.lt.abs(a-b)*.063.and.(htb1.gt.1.e-6*abs(a-b).or.max(a
     &           ,b)-r.lt.r-min(a,b)))then
               call minorrad(htb)
               htb=sqrt(pss*2./htb*(dzpsi**2+dpsi**2))
     &              /sqrt(max(abs(r-a),abs(b-r)))
c
c               write(*,*) '1htb,pss,dzpsi,r,a,b,m  ',htb,pss,dzpsi,r,a
c     &              ,b,imtae
               htb=max(abs(htb),.003)
               htb=r**2/htb*bm
            else
               htb=dzpsi*sqrt(pss)/sqrt((r-a)*(b-r))
c
c               write(*,*) '2htb,pss,dzpsi,r,a,b,m  ',htb,pss,dzpsi,r,a
c     &              ,b,imtae
               htb=max(abs(htb),.003)
               htb=r**2/htb*bm
            endif
         else if(pss.lt.-5.*r*epsp)then
            goto 11
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
c
c       write(*,*) '2a,rvpar,r,z  ',a,rvpar,r,z
c            htb=dzpsi*sqrt((r-rvpar)/(r-a)/r/(b-r))
c         write(*,*) '2dzp,htb,r,z  ',dzpsi,htb,r,z
c            htb=max(abs(htb),.2)
c            htb=1./htb*bm*r**2
         endif
         htb=htb*htb0
         ss=real(mtae)*(ff(1)
     &        -wvpr)
         cs=real(ntae)*(ff(2)+ff(3))-wtae
     &        +wres
         hwarn=(ss-cs)*htb
         if(abs(hwarn).gt.1.5)warn=1.
         fo(1)=fo(1)+hwarn
         hwarn=(-ss-cs)*htb
         fo(6)=fo(6)+hwarn
         if(abs(hwarn).gt.1.5)warn=1.
         xav=xav+xavint*htb
         x2r=x2r+x2i*htb
         tb=tb+htb
         kfast=1
c         if(mtae.eq.3) then
c            write(*,*) 'cs,3xfs,r,wtae,wres='
c            write(*,*) cs,ff(1),ff(2),ff(3),r,wtae,wres
c            write(*,*) '--- i=',i,z,fo(1),fo(6),htb,dzpsi,wvpr
c         endif
         call fgm(ff(1),8,r,symm,v,z)
         if(symm.eq.0.)then
c     asymmetric case
            ss=const*htb
            cs=ss*cos(fo(1))
            cs1=0.
css*cos(fo(6))
            ss1=0.
css*sin(fo(6))
            ss=ss*sin(fo(1))
c     FLR included
            fo(4)=fo(4)+ff(2)*(cs+cs1)+ff(4)*(ss-ss1)
            fo(5)=fo(5)+ff(2)*(ss+ss1)+ff(4)*(cs1-cs)
c                      + imaginary contribution
            fo(4)=fo(4)+ff(6)*(cs-cs1)-ff(8)*(ss+ss1)
            fo(5)=fo(5)+ff(6)*(ss-ss1)+ff(8)*(cs1+cs)
c     FLR not included
            fo(2)=fo(2)+ff(1)*(cs+cs1)+ff(3)*(ss-ss1)
            fo(3)=fo(3)+ff(1)*(ss+ss1)+ff(3)*(cs1-cs)
c                      + imaginary contribution
            fo(2)=fo(2)+ff(5)*(cs-cs1)-ff(7)*(ss+ss1)
            fo(3)=fo(3)+ff(5)*(ss-ss1)+ff(7)*(cs1+cs)
c            if(mtae.eq.1)then
c               write(*,*) '--i',i,ff
c               write(*,*) ' Sum:i,Phase`s',i,fo(1),fo(6)
c               write(*,*) fo(5),fo(4),ss-ss1
c            endif
c     write(*,*) 'dfo3(2),dfo5(2)='
c     write(*,*) ff(1)*(ss+ss1),ff(3)*(cs1-cs),ff(2)*(ss+ss1
c     &           ),ff(4)*(cs1-cs)
         else
c     symmetric case
            ss=const*htb
c            cs=ss*(cos(fo(1))+cos(fo(6)))
c            ss=ss*(sin(fo(1))-sin(fo(6)))
            cs=ss*cos(fo(1))
            ss=ss*sin(fo(1))
c     FLR included
            fo(4)=fo(4)+ff(2)*cs+ff(4)*ss
            fo(5)=0.
c     FLR not included
            fo(2)=fo(2)+ff(1)*cs+ff(3)*ss
            fo(3)=0.
c            if(mtae.eq.6)then
c               write(*,*) '--i',i,ff(1),(cs),ff(3),ss
c               write(*,*) ' Sum:i,Phase`s',i,fo(1),fo(6)
c              write(*,*) ff(2),ff(4),fo(4),ss
c            endif
c      write(*,*) ' ---------tb,phase -,+ m fo2=',tb,fo(1),fo(6)
c     &           ,fo(2)
         endif
c     if(mtae.eq.1)write(*,*) htb,const*htb,fo(2),fo(3),mtae,tb
c     
c     Write(*,*) 'Dzpsi,r,z,htb,pss,fo(3),fo(1)= ',dzpsi,r,z,htb
c     &        ,pss,ff(1)
 11   continue
      if(b0.lt..01) b0=1.
c     
c      write(*,*) '---------------------------Sum ',htb,tb,fo(2),iav
      if(tb.eq.0.) then
         do ii=1,iav
            fo(ii)=0.
         enddo
         warn=0.
         return
      endif
      do ii=1,iav
         fo(ii)=fo(ii)/tb
      enddo
      xav=xav/tb
      x2r=x2r/tb
c      write(*,*) 'wb= ',2.*pi/tb,fo,wtae
      return
      end
c**********************************************************************
c On output it produces the distribution function in cm^-3,
c  i.e. it is the same spatial function but multiplied by v_0^3.
c Also here gam is the parameter so that the df is
c  f~/(v^3 + v*^3)^(2-gam) and gam=1 corresponds to the slowing down df.(BepHo)
c tb,wb_res,dwbdpphi are used only if tip='r'
      subroutine distrfun(ai,b0oo,dfalp,dfalpdv,dfalpdpsi,dpga,eps,psior
     &     ,gam,hnub,hchirp,ihsps,iter,qpf,p0ga,pmu,rvpar,sigm,tb,tip2
     &     ,vinstant1,vi,vmax,vstortkeal,z0,xhp,zi,dwbdpphi,wb_res,
     &     istart)
      include 'clich1'
      include 'clich1b'
      include 'clich2'
      include 'orb2d.par'
      character*1 tip1,tip2
      dimension vstar3(nn),vstar3s(nn),cvsld(nn),denshh(nn),hchirp(2),
     &     dtmzl(ispc,4),cln(ispc,nn)
      parameter(Nlor=26)
      real*8 sa(Nlor,2),sa2(Nlor,2),sa3(Nlor,2) ! sa, sa2 and sa3 are the sources at full, half and third of inj. energy
c fractions of beam depositions at the respective values of their energies
c comments from transp code: # If Jb is the beam current, the relationship between injection energy
c                            # (kvolt_nbi) and injection power (power_nbi) satisfies:
c                            #   power_nbi = Jb*kvolt_nbi*(FFULL + FHALF/2 + FTHIRD/3)
c values for NSTX-U devices are left according to M. Gorelenkova: xfull/0.461/, xhalf/0.386/
      data xfull/0.461/, xhalf/0.386/
      save tip1,vstar3,vstar3s,cvsld,denshh,cv,dnpsi0,betalpha0,cnorm,
     &     sa,sa2,sa3
      Tip=tip2
c Calculate the normalization every time since the normalization depends on the velocity, which changes
c every call.
      if(istart.eq.1) then
         tip1=''
         vstar3=0.
         vstar3s=0.
         cvsld=0
         denshh=0
         cv=0
         dnpsi0=0
         betalpha0=0
         cnorm=0
         sa=0
         sa2=0
         sa3=0
         istart=0
      endif
      if(tip.ne.tip1.or.tip.eq.'g')then
         if(tip.ne.tip1) then
            if(tip.eq.'x') then
               print *,'Compute ion Landau damping accurately, using -K'
               print *,'use D ions with central beta',betac0(2),
     &              'type of DF is',tip
            else
               write(*,*)' Ihsps=',ihsps,' betah0=',betah0(ihsps),tip
     &              ,tip1
               write(57,*) 'Fast ion beta, beta_h'
               write(57,*) betah0(ihsps)
            endif
c
            if(tip.eq.'l'.or.tip.eq.'b')then
               call sal(Nlor, sa(1,1), max(0.02,chidelt(2))*0.5,
     &               abs(chidelt(1)))
               call sal(Nlor, sa(1,2), max(0.02,chidelt(2))*0.5,
     &              -abs(chidelt(1)))
               sa2=sa
               sa3=sa
               do i=1,Nlor
                  if(tip.eq.'b')then
                     sa2(i,2)=0.5*(sa(i,2)+sa(i,1))
                     sa3(i,2)=0.5*(sa(i,2)+sa(i,1))
                  endif
                  sa(i,2)=0.5*(sa(i,2)+sa(i,1))
               enddo
               write(*,*) 'Using Lorentz distribution with Chi0'
     &              ,chidelt(1),'and width',sqrt(max(0.02,chidelt(2)))
c               write(*,*) sa
            endif
         endif
         cnorm=1.
         if(tip.eq.'c'.or.tip.eq.'g'.or.tip.eq.'h'.or.tip.eq.'i')then
            if(tip.ne.'g') 
     &           write(*,*)'Copassing at P0',
     &           p0ga,' width',dpga,'rax= ',rax
            do i=1,nosurf
               if(tip.eq.'g')then
                  dpga=vinstant1**3/chidelt(5)
                  chis=sqrt(1.-b2d(i,1)/b2d(i,nts0/2))
                  dpga=sqrt(max(0.02
     &                 ,chidelt(2)+sqrt(ppsi(i))*chidelt(3)-chidelt(4)
     &                 *alog(dpga*(1.+chidelt(6)**3)/(dpga+chidelt(6)**3
     &                 ))))
c this normalization is good for beams with different injections
                  if(abs(chidelt(1)).gt.chis)then
                     denshh(i)=erf1((chis+2-abs(chidelt(1)))/dpga)
     &                 +erf1((abs(chidelt(1))+chis)/dpga)
     &                 -0.5*(erf1((-abs(chidelt(1))+2-chis)/dpga)
     &                 -erf1((-abs(chidelt(1))+3-2.*chis)/dpga))
                  else
                     denshh(i)=erf1((1-chidelt(1))/dpga)
     &                 +erf1((chidelt(1)+1)/dpga)
                  endif
c here one of these 1/2 is because of integration from -1 to  1,
c i.e. the domain would give 2
                  denshh(i)=abs(denshh(i))*3.1415926**0.5/2.*dpga/2.
               else
                  call nrmlzn(b2d(i,1)/b2d(1,1),dpga,p0ga,denshh(i))
               endif
            enddo
c            if(psior.gt.0.3.and.psior.lt.0.32)write(*,*)'inside'
c     &           ,denshh(1),denshh(100),p0ga,dpga,pmu,qpf
         endif
         if(tip.eq.'j'.or.tip.eq.'e')then
            dnpsi0=0.09
            betalpha0=0.01
         else if(tip.eq.'t')then
            dnpsi0=0.1
ckg            dnpsi0=0.5
            betalpha0=0.002
         endif
         tip1=tip
         cp=1.
         do i=1,nosurf
            if(tip.eq.'d'.or.tip.eq.'t'.or.tip.eq.'j'.or.tip.eq.'e')then
               vstar3(i)=0.037037037
               z1=.5*(1.+denc(i,4)/denc(i,1))-denc(i,3)/denc(i,1)/6.
               vstar3s(i)=vstar3(i)*zeff(i)/z1
            else
c a place to call the Coulomb log subroutine
cold expr if log_e=log_i
               z1a=.5*(1.+denc(i,4)/denc(i,1))-denc(i,3)/denc(i,1)/6.
               do il=1,ispc               !loop for e,d,t,h,c
                  dtmzl(il,1)=denc(i,il)  !e density at il=1
                  dtmzl(il,2)=tc(i,il)    !e temp
                  dtmzl(il,3)=rmcp(il)    !e mass to proton mass
                  dtmzl(il,4)=zc(il)      !e charge
               enddo
               call coulog(vi**2*ai*882000/1.298**2,ai,zi,dtmzl(1,1),
     &              dtmzl(1,2),dtmzl(1,3),dtmzl(1,4),cln(1,i),5,b0oo)
               z1=0.
               ze=0.
               do il=2,ispc 
c mb logi/loge 
                  z1=z1+dtmzl(il,4)**2*dtmzl(il,1)/(dtmzl(il,3)*
     &                 dtmzl(1,1))*cln(il,i)/cln(1,i)
                  ze=ze+dtmzl(il,4)**2*dtmzl(il,1)/dtmzl(1,1)
     &                 *cln(il,i)/cln(1,i)
               enddo
c               print *,'Coulog',z1,ze,(cln(il,i),il=1,ispc)
               vstar3(i)=z1a*(.16784/vi*sqrt(tc(i,1)/1.e3))**3
               vstar3s(i)=zeff(i)*(.16784/vi*sqrt(tc(i,1)/1.e3))**3
c            print *,i,tc(i,1),zeff(i),vstar3(i)**0.33333,vstar3s(i)**0
c     &           .333333
            endif
         enddo
c pause the program for a few seconds
c         call sleep(10)
         if(tip.eq.'h'.or.tip.eq.'i'.or.tip.eq.'m'.or.tip.eq.'x'.or.
     &        tip.eq.'e')then
            cv=sqrt(pi)/4.*pi*4.
            do i=1,nosurf
               if(tip.eq.'e') then
                  denshh(i)=betalpha0*4.755e12*b0oo**2/vi**2/ai*
     &                 (1-ppsi(i))**2
c     exp(-ppsi(i)/dnpsi0)
               else if(tip.eq.'h'.or.tip.eq.'i')then
                  denshh(i)=betah0(ihsps)*ph(i,ihsps)*4.755e12*b0oo**2
     &                 *th(1,ihsps)/(ph(1,ihsps)*denshh(i)*vi**2*ai*th(i
     &                 ,ihsps))
               else if(tip.eq.'x')then
c /2 m vozniknut iz-za sravneniya s transpom
                  denshh(i)=denc(i,2)
c13 
               else
                  denshh(i)=betah0(ihsps)*ph(i,ihsps)*4.755e12*b0oo**2
     &                 /(ph(1,ihsps)*vi**2*ai)
               endif
            enddo
         else if(tip.eq.'r')then
ckg b0oo*3. byl vveden (ran'she ne bylo)
            cv=0.3195*zi*psitot*b0oo*3.
            cnorm=vi**3*zi*ai**2*0.0882*3.
         else
c               write(*,*)  'sqrt(\Psi),     rhoprf'
            do i=1,nosurf
               hhhb=vstar3(i)**0.333333
               if(tip.ne.'c'.and.tip.ne.'g')denshh(i)=1.
c                  write(*,*) 'df_norm0',denshh(i),vstar3(i)
c this comes from the velocity integration
               if(abs(gam-1.).lt..001) then
                  cvsld(i)=alog(1./vstar3(i)+1.)/3.
                  denshh(i)=denshh(i)*(0.5+hhhb**2/3.*(0.5*alog((1.+hhhb
     &                 )**2/(1.-hhhb+hhhb**2))-1.732051*(0.5236+atan((2.
     &                 /hhhb-1.)/1.732051))))/cvsld(i)
               else
                  gamc=2-gam
                  cvsld(i)=((1.+vstar3(i))**(1.-gamc)-
     +                 vstar3(i)**(1.-gamc))/3./(1.-gamc)
                  denshh(i)=denshh(i)*dfint(gam,vstar3(i))/cvsld(i)
               endif
c                  write(*,*) 'df_norm1',denshh(i),vstar3(i)
               cvsld(i)=cvsld(i)*pi*4.
               if(tip.eq.'t'.or.tip.eq.'j')then
                  denshh(i)=3./2.*betalpha0*4.755e12*b0oo**2/vi**2/ai
     &                 /denshh(i)*
ckg     &                    exp(-ppsi(i)/dnpsi0)
     &                 (1-ppsi(i))**2
               else
c     The following expression doesnot contain ph * 2 which comes from the
c     comparison of Frank's value and real Alphas pressure.
c This 3/2 comes from the definition of temperature as a 2/3 of the mean
c kinetic energy.
                  denshh(i)=3./2.*betah0(ihsps)*ph(i,ihsps)*4.755e12
     &                 *b0oo**2/(ph(1,ihsps)*vi**2*ai*denshh(i))
               endif
c                  write(*,*)  rgrid(i),rhoprf(i)
            enddo
         endif
      endif
c     Berk's stuff
      i0=max(indx(ppsi,nosurf,psior),1)
      i0=min(i0,nosurf-3)
      call funder4(ppsi(i0),vstar3s(i0),psior,dvstrdpsi,vstar)
      call funder4(ppsi(i0),denc(i0,1),psior,dvstrdpsi,hnub)
      call funder4(ppsi(i0),tc(i0,1),psior,dvstrdpsi,hhhb)
      hhhb=abs(hhhb)
c drag sl.down freq. in sec^-1
      hchirp(2)=(hnub*1.e-14)*zi**2/(hhhb/1000.)**1.5*1.e2/ai  ! plus plus
     &     *(1.+vstar)
c introducing brute force approach coefficient _xTebf for scaling Te, read mkorb to activate this option
#ifdef _xTebf
      hchirp(2)=hchirp(2)/(_xTebf)**1.5 ! plus plus
c      print *,'used scaling factor for Te ',_xTebf
c      stop 
#endif
c this is collisional scatt. frequency in 10^6 (yes) sec^-1, extra 2ai 
c is from 2 in the first equation on p.68 by Goldston_JCP81
c both sl.down and coll.scat include velocity as a factor
      hnub=0.5e-18*hnub*zi**2/ai**2*vstar/vinstant1**3/(hhhb/1000.)
     &     **1.5
c coll.scat. in sec^-1
      hchirp(1)=hnub*1.e6 ! plus plus; in many cases drag is faster as scatt is weaker,i.e. vcr < v=vA 
c
      if(psior.lt..005)then
         dfalp=0.
         dfalpdv=0.
         dfalpdpsi=0.
         return
      endif
      if(tip.eq.'r')then
         hhhb=cnorm/tb
         dfalp=dfn_fun(dfalpdv,dfalpdpsi,th0(ihsps)*vinstant1**2,cv*qpf
     &        ,pmu*th0(ihsps)*vinstant1**2/(b0oo*rax))
c         write(*,*) 'V: ',dfalpdv*th0(ihsps)*vinstant1*2.,dfalp*wb_res
         dfalpdv=hhhb*(dfalpdv*th0(ihsps)*vinstant1*2.+dfalp*wb_res)
         dfalp=hhhb*dfalp
c     dwbdpphi, wb_res
c         write(*,*) 'P_phi',-cv*dfalpdpsi,dfalp*dwbdpphi
         dfalpdpsi=hhhb*(-cv*dfalpdpsi+dfalp*dwbdpphi)
c         write(*,*)hhhb,cv,cnorm,th0(ihsps),ihsps,vinstant1,wb_res,dfalp
c         pause
         return
      endif
      call funder4(ppsi(i0),vstar3(i0),psior,dvstrdpsi,vstar)
cnng13 define vstar to be constant for FY14 milestone
cnng13      vstar=0.96**3
cnng14      vstar=0.48**3
cnng13
c     Here dfalpdpsi is derivative over the psi, 
c     dfalpdv is the density itself
      if(tip.eq.'t'.or.tip.eq.'j'.or.tip.eq.'e')then
         dfalpdv=denshh(1)*
ckg     &        exp(-psior/dnpsi0)
     &        (1-psior)**2
         dfalpdpsi=
ckg     &         -dfalpdv/dnpsi0
     &        -2.*denshh(1)*(1-psior)
      else
c Try to campare with frank's using derivative from d P_\alpha / d\psi
c         call funder4(ppsi(i0),denh(i0,ihsps),psior,dfalpdpsi,dfalpdv)
ccc        call funder4(ppsi(i0),ph(i0,ihsps),psior,dfalpdpsi,dfalpdv)
ccc
         call funder4(ppsi(i0),denshh(i0),psior,dfalpdpsi,dfalpdv)
         call fun4(ppsi(i0),denshh(i0),psior,denshhor)
         dfalpdpsi=denshhor*dfalpdpsi/dfalpdv
         dfalpdv=denshhor
      endif
      if(tip.eq.'h'.or.tip.eq.'i'.or.tip.eq.'x'.or.tip.eq.'m'
     &     .or.tip.eq.'e')then
          dfalp=dfalpdv/cv
          if(tip.eq.'x') then
c
ckg add maxwellian part for the model of thermal ions with eta_i 
contribution
             call funder4(ppsi(i0),tc(i0,2),psior,dtidpsi,tiloc)
             call funder4(ppsi(i0),denc(i0,2),psior,dnidpsi,dniloc)
cnng13
             if(abs(dnidpsi*tiloc).gt.1.e-20) etai=
     &            dtidpsi*dniloc/(dnidpsi*tiloc)
c recycle tiloc variable
             if(abs(tiloc).gt.1.e-20) then
                tiloc=vinstant1**2*ai*vi**2*1.e8/(191.69*tiloc)
             else
                tiloc=0
             endif
             dfalpdpsi=dfalpdpsi/cv*(1.+etai*(-1.5+tiloc))
             if(tiloc.gt.10.)then
                tiloc=0.
             else
                tiloc=exp(-tiloc)
             endif
             dfalpdpsi=dfalpdpsi*tiloc
             dfalp=dfalp*tiloc
          else
             dfalp=dfalp*exp(-vinstant1**2)
             dfalpdpsi=dfalpdpsi*exp(-vinstant1**2)/cv
          endif
          dfalpdv=-2.*vinstant1*dfalp
          if(tip.eq.'i')then
             call funder4(ppsi(i0),th(i0,ihsps),psior,dtidpsi,tiloc)
             call funder4(ppsi(i0),denshh(i0),psior,dnidpsi,dniloc)
             etai=dtidpsi*dniloc/(dnidpsi*tiloc)
c
             dfalpdpsi=dfalpdpsi*(1.+etai*(-1.5+vinstant1**2))
          endif
      else    ! slowing down, beam and "general" DF 
          call fun4(ppsi(i0),cvsld(i0),psior,cv)
          vinstant=min(vinstant1,1.)
          vinst2=vinstant*sqrt(2.)
          vinst3=vinstant*sqrt(3.)
          vstar2=vstar*2**1.5
          vsta3 =vstar*3**1.5
          p3xv=vstar
          p3xv2=vstar2
          p3xv3=vsta3
c
          hhhb=1.
          vstar =vinstant**3+vstar
          vstar2=vinst2**3+vstar2
          vsta3=vinst3**3+vsta3
c          hhdenom=vinstant**3/vstar
c collecting 3 sources into 1 value by similarity with dfalp-like variables
          dfalp=dfalpdv/cv
          dfalpdv=vstar
          dfalp2=dfalpdv/cv
          dfalpdv2=vstar
          dfalp3=dfalpdv/cv
          dfalpdv3=vstar
          gamc=2-gam
          hhhb=hhhb*vstar**gamc
c DFALP is alpha particle distribution function
         dfalp=dfalp/hhhb
         dfalpdv=dfalp/dfalpdv
c next thing is the derivative over the velocity
         if(vinstant1.le.1.) then
            dfalpdv=-3.*gamc*vinstant**2*dfalpdv
         else
            dfalpdv=dfalp/(1.-vmax)
         endif
         if(vinstant1.le.0.7071) then
            dfalpdv2=-3.*gamc*vinstant**2*dfalpdv2
         elseif(vinstant1.le.0.7071*vmax)then
            dfalpdv2=dfalp2/(1.-vmax)
         else
            dfalpdv2=0.
         endif
         if(vinstant1.le.0.57735) then
            dfalpdv3=-3.*gamc*vinstant**2*dfalpdv3
         elseif(vinstant1.le.0.57735*vmax)then
            dfalpdv3=dfalp3/(1.-vmax)
         else
            dfalpdv3=0.
         endif
            
            
c here we finally calculate the derivative $d F / d\psi $
         dfalpdpsi=dfalpdpsi/hhhb/cv
      endif
      if(tip.eq.'l'.or.tip.eq.'b'.or.tip.eq.'g')then
         call fun4(ppsi(i0),b2d(i0,1),psior,hhhb)        ! HFS value of B^2 is hhhb
         call fun4(ppsi(i0),b2d(i0,nts0/2),psior,chis)   ! the same is for LFS here: chis
         chis=sqrt(1.-hhhb/chis)
ckg central pitch angle calculation
         p0ga=sqrt(abs(1-pmu*hhhb/b2d(1,1)/rax))*sigm
      endif
c if Lorentz collisional scattering coefficient is used
      if(tip.eq.'l'.or.tip.eq.'b')then
c                  print *,vinstant,p3xv,vstar
c         if (p3xv.lt.0.) then
c            pause
c         endif
         tau=-alog(abs(vinstant**3/vstar)*(1.+p3xv))*0.33333333333
         tau2=-alog(abs(vinst2**3/vstar2)*(1.+p3xv2))*0.33333333333
         tau3=-alog(abs(vinst3**3/vsta3)*(1.+p3xv3))*0.33333333333
c
         if(abs(chidelt(1)).lt.chis)then
            call lorentzdf(chis, chidelt(1), abs(p0ga), dfl,
     &           ddflpdgc,ddflpdgt, 2.*tau,Nlor, sa(1,2))
            call lorentzdf(chis, chidelt(1), abs(p0ga), dfl2,
     &           ddflpdgc2,ddflpdgt2, 2.*tau2,Nlor, sa2(1,2))
            call lorentzdf(chis, chidelt(1), abs(p0ga), dfl3,
     &           ddflpdgc3,ddflpdgt3, 2.*tau3,Nlor, sa3(1,2))
c            dfl=gfL(Nlor, ddflpdgc, ddflpdgt, 2.*tau, p0ga, sa(1,2))
         else
            if(chidelt(1).gt.0.)then
               call lorentzdf(chis, chidelt(1), p0ga, dfl,
     &              ddflpdgc,ddflpdgt, 2.*tau,Nlor, sa(1,1))
               call lorentzdf(chis, chidelt(1), p0ga, dfl2,
     &              ddflpdgc2,ddflpdgt2, 2.*tau2,Nlor, sa2(1,1))
               call lorentzdf(chis, chidelt(1), p0ga, dfl3,
     &              ddflpdgc3,ddflpdgt3, 2.*tau3,Nlor, sa3(1,1))
            else
               call lorentzdf(chis, chidelt(1), -p0ga, dfl,
     &              ddflpdgc,ddflpdgt, 2.*tau,Nlor, sa(1,1))
               call lorentzdf(chis, chidelt(1), -p0ga, dfl2,
     &              ddflpdgc2,ddflpdgt2, 2.*tau2,Nlor, sa2(1,1))
               call lorentzdf(chis, chidelt(1), -p0ga, dfl3,
     &              ddflpdgc3,ddflpdgt3, 2.*tau3,Nlor, sa3(1,1))
            endif
            ddflpdgc=ddflpdgc*sigm
            ddflpdgc2=ddflpdgc2*sigm
            ddflpdgc3=ddflpdgc3*sigm
         endif
c summ all the DF from 3 sources depending on the ion energy
         dfl3=abs(dfl3)
         dfalp3=dfalp3*2.*dfl3
         dfalpdv3=2.*(dfalpdv3*dfl3-2.*ddflpdgt3*dfalp3*p3xv3/(vinst3*
     &        vsta3)+dfalp3*ddflpdgc3*(1.-p0ga*p0ga)/(vinst3*p0ga))
         dfalpdpsi3=dfalpdpsi*2.*dfl3
c
         dfl2=abs(dfl2)
         dfalp2=dfalp2*2.*dfl2
         dfalpdv2=2.*(dfalpdv2*dfl2-2.*ddflpdgt2*dfalp2*p3xv2/(vinst2*
     &        vstar2)+dfalp2*ddflpdgc2*(1.-p0ga*p0ga)/(vinst2*p0ga))
         dfalpdpsi2=dfalpdpsi*2.*dfl2
c
         dfl=abs(dfl)
         dfalp=dfalp*2.*dfl
         dfalpdv=2.*(dfalpdv*dfl-2.*ddflpdgt*dfalp*p3xv/(vinstant*vstar)
     &        +dfalp*ddflpdgc*(1.-p0ga*p0ga)/(vinstant1*p0ga))
         dfalpdpsi=dfalpdpsi*2.*dfl
c         write(*,*)'ls-chi,df__df,dfdv,dfdpsi',p0ga,dfalp*2.*dfl,2.
c     &        *(dfalpdv*dfl-2.*ddflpdgt*dfalp*p3xv/(vinstant*vstar)
c     &        +dfalp*ddflpdgc*(1.-p0ga*p0ga)/(vinstant1*p0ga)),dfalpdpsi
c     &        *2.*dfl,tau, -2.*2.*ddflpdgt*dfalp*p3xv/(vinstant*vstar)
c     &        +2.*dfalp*ddflpdgc*(1.-p0ga*p0ga)/(vinstant1*p0ga),2.
c     &        *dfalp*ddflpdgc*(1.-p0ga*p0ga)/(vinstant1*p0ga),chis
c     &        ,chidelt(1),dfl
c            write(*,*) sa
c         return
      endif
      if(tip.eq.'h'.or.tip.eq.'i'.or.tip.eq.'c'.or.tip.eq.'g')then
         if(tip.eq.'g')then
            dpga=vinstant1**3/chidelt(5)
            dpga=sqrt(max(0.02
     &           ,chidelt(2)+sqrt(psior)*chidelt(3)-chidelt(4)
     &           *alog(dpga*(1.+chidelt(6)**3)/(dpga+chidelt(6)**3
     &           ))))
c            p0ga=(1-chidelt(1)**2)*rr(i,1)/rax
c            dpga=vinstant1**3*vi**2*ai*0.88/(1.3**2*chidelt(5))
c            dpga=2.*chidelt(1)*sqrt(max(0.05,chidelt(2) +sqrt(ppsi(i))
c     &           *chidelt(3)-chidelt(4)*alog(dpga* (1.+chidelt(6)**3)
c     &           /(dpga+chidelt(6)**3))))
c         if(i0.gt.80.and.i0.lt.85)write(*,*)'inside',p0ga,dpga
c     &        ,pmu
            hhhb=(p0ga-chidelt(1))/dpga
            if((p0ga.gt.chis.and.chidelt(1).gt.chis).or.(p0ga.lt.-chis
     &           .and.chidelt(1).lt.-chis))then
               vstar=(exp(-hhhb**2)
     &              +exp(-((abs(p0ga+chidelt(1))-2)/dpga)**2)
     &              -0.5*exp(-((abs(p0ga+chidelt(1))-2.*chis)/
     &              dpga)**2))/cnorm
            elseif((p0ga.gt.chis.and.chidelt(1).lt.-chis).or.(p0ga.lt.
     &              -chis.and.chidelt(1).gt.chis))then
               vstar=0.5*exp(-((p0ga-chidelt(1)+sign(2.
     &              ,chidelt(1))*chis)/dpga)**2)/cnorm
            else
               vstar=0.5*(exp(-hhhb**2)+exp(-((p0ga+chidelt(1))
     &              /dpga)**2))/cnorm
            endif
         else
            hhhb=(pmu/rax-p0ga)/dpga
            vstar=exp(-hhhb**2)/cnorm            
         endif
         dfalp=dfalp*vstar
c Following expression includes anisotropy drive
         if(tip.eq.'g')then
            dfalpdv=dfalpdv*vstar
            dfhhh=dfalpdv
c            goto 15
            dchidv=-1.5*chidelt(4)*chidelt(6)**3/(vinstant1*dpga**2
     &           *(vinstant1**3/chidelt(5)+chidelt(6)**3))
            if((p0ga.gt.chis.and.chidelt(1).gt.chis).or.(p0ga.lt.-chis
     &           .and.chidelt(1).lt.-chis))then
c               dfalpdv=dfalpdv+dfalp*2.*hhhb**2*dchidv
c     &              -2.*dfalp*hhhb*(1.-p0ga*p0ga)/(dpga*vinstant1*p0ga)
               dfalpdv=dfalpdv+dfalp/(vstar*cnorm)*2.*dchidv*(hhhb**2
     &              *exp(-hhhb**2)+((abs(p0ga+chidelt(1))-2)/dpga)**2
     &              *exp(-((abs(p0ga+chidelt(1))-2)/dpga)**2)
     &              -0.5*((abs(p0ga+chidelt(1))-2.*chis)/dpga)**2
     &              *exp(-((abs(p0ga+chidelt(1))-2.*chis)/dpga)**2)
     &              )
     &              -2.*dfalp/(vstar*cnorm)*(hhhb*exp(-hhhb**2)
     &              +((p0ga+chidelt(1)-sign(2.,chidelt(1)))/dpga)
     &              *exp(-((abs(p0ga+chidelt(1))-2)/dpga)**2)
     &              -0.5*((p0ga+chidelt(1)-sign(2.,chidelt(1))*chis)
     &              /dpga)*exp(-((abs(p0ga+chidelt(1))-2.*chis)/dpga)**2
     &              ))*(1.-p0ga*p0ga)/(dpga*vinstant1*p0ga)
ckg               if(vinstant1.gt.0.5.and.abs(p0ga).gt.0.5) write(*,*) p0ga
ckg     &              ,vinstant1,dfhhh,dfalpdv-dfhhh,dfhhh/(dfalpdv-dfhhh)
ckg     &              ,dchidv,hhdenom
            elseif((p0ga.gt.chis.and.chidelt(1).lt.-chis).or.(p0ga.lt.
     &              -chis.and.chidelt(1).gt.chis))then
               dfalpdv=dfalpdv+dfalp*2.*((p0ga-chidelt(1)+sign(2.
     &              ,chidelt(1))*chis)
     &              /dpga)**2*dchidv-2.*dfalp*((p0ga-chidelt(1)+sign(2.
     &              ,chidelt(1))*chis)/dpga)*(1.-p0ga*p0ga)/(dpga
     &              *vinstant1*p0ga)
            else
               dfalpdv=dfalpdv+dfalp/(vstar*cnorm)*2.*dchidv*(hhhb**2
     &              *exp(-hhhb**2)+((p0ga+chidelt(1))/dpga)**2*exp(
     &              -((p0ga+chidelt(1))/dpga)**2))
     &              -2.*dfalp/(vstar*cnorm)*(hhhb*exp(-hhhb**2)
     &              +((p0ga+chidelt(1))/dpga)
     &              *exp(-((p0ga+chidelt(1))/dpga)**2))*(1.-p0ga*p0ga)
     &              /(dpga*vinstant1*p0ga)
            endif
            if(abs(chidelt(1)).gt.chis)then
               xint=erf1((-chis+2-abs(chidelt(1)))/dpga)
     &                 +erf1((abs(chidelt(1))+chis)/dpga)
               dfalpdv=dfalpdv+dfalp*dchidv/(dpga)*(
     &              2./(xint*3.1415926**0.5)*(exp(-(2-chis-abs(chidelt(1
     &              )))**2/dpga**2)*(2-chis-abs(chidelt(1)))+exp(-(chis
     &              +abs(chidelt(1)))**2/dpga**2)*(chis+abs(chidelt(1)))
     &              )-dpga)
            else
               xint=erf1((1-chidelt(1))/dpga)
     &                 +erf1((chidelt(1)+1)/dpga)
               dfalpdv=dfalpdv+dfalp*dchidv/(dpga)*(2./(xint*3.1415926
     &              **0.5)*(exp(-(1-chidelt(1))**2/dpga**2)*(1-chidelt(1
     &              ))+exp(-(1+chidelt(1))**2/dpga**2)*(1+chidelt(1)))
     &              -dpga)
            endif
         else
            dfalpdv=dfalpdv*vstar
     &           +4.*pmu*dfalp*hhhb/(dpga*rax*vinstant1)
         endif
c one can use this simplification where we neglected exact integration
c to improve accuracy but may be OK without it
c2*dfalp/vinstant1
 15      dfalpdpsi=dfalpdpsi*vstar
      endif
cnng14 to account for three energy sources in the beam distribution function
c      see above definition for xfull xhalf
      if(tip.eq.'b')then
cnng20         xfull=0.
cnng20         xhalf=0.
         hhdenom=1./(xfull+0.5*xhalf+0.3333333333*(1.-xfull-xhalf))
         dfalpdv=dfalpdv*xfull*hhdenom
         Dfalp=dfalp*xfull*hhdenom
         dfalpdpsi=dfalpdpsi*xfull*hhdenom
         if(vinstant1.lt.0.7071*vmax)then
            dfalpdv=dfalpdv+dfalpdv2*xhalf*0.5*hhdenom
            dfalp=dfalp+dfalp2*xhalf*hhdenom*0.5
            dfalpdpsi=dfalpdpsi+dfalpdpsi2*xhalf*0.5*hhdenom
            if(vinstant1.lt.0.57735*vmax)then
               dfalpdv=dfalpdv+dfalpdv3*(1-xfull-xhalf)*hhdenom/3.
               dfalp=dfalp+dfalp3*(1-xfull-xhalf)*hhdenom/3.
               dfalpdpsi=dfalpdpsi+dfalpdpsi3*(1-xfull-xhalf)*hhdenom/3.
            endif
         endif
      endif
c      dfalpdv=0.
c      write(*,*)'chi,df__df,dfdv,dfdpsi',p0ga,dfalp,dfalpdv
c     &     ,dfalpdpsi
      return
      end
c*********************************************************
      subroutine eigenmod(b0oo,b0,condec,condic,dbobmax,dbobmaxs,
     &     dnonmax,grps2,rgridc,min,min_f,max,max_f,ntae,scalvavf,
     &     valfv,wtae,deltk,ximax)
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      include 'orb2d.par'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      common/funct/xkpop(nn,mt,nts,3)
      dimension delkintgr(nn),grps2(1000),rgridc(1000),
     &     condec(1000),condic(1000)
      complex wtae
ccc  choose the sign of om to be same as alpha's wstar to get instability
ccc  for n>0, om=-omreal; for n<0, om=omreal, the former is usual
      write(*,*) 'q1,rho0,rax,r',q1,rho0,rax,r
c      write(*,*) 'scalvavf,omreal,b0oo,r,rho0',scalvavf,omreal,b0oo,r
c     &     ,rho0
      if(omreal.lt.0.)then
         wtae=cmplx(0.,-scalvavf*omreal*2.18e7*b0oo/r/q1/sqrt(rho0))
      else
         wtae=cmplx(-scalvavf*omreal*2.18e7*b0oo/r/q1/sqrt(rho0),0.)
      endif
      valfv=2.18e6*scalvavf*b0oo/sqrt(rho0)
c      wtae=-scalvavf*omreal*2.18e7*b0oo*rax/r**2/q1/sqrt(rho0)
c      valfv=2.18e6*scalvavf*b0oo*rax/r/sqrt(rho0)
      ntae=ntor
      min=minm
      max=maxm
      min_f=max
      max_f=min
      do i=1,nosurf
         if(i.gt.3) delkintgr(i)=delk(i)*rhoprf(i)*rgrid(i)
         grps2(i)=0.5*(grpssq2d(i,nts0/2)+grpssq2d(i,1))
         rgridc(i)=rgrid(i)
         condec(i)=conde(i)
         condic(i)=condi(i)
      enddo
      call asimp(rgrid(4),delkintgr(4),nosurf-3,deltk)
      deltk=deltk*abs(wtae)**2*1.e-2*psitot*rho0*2./scalvavf**2
      ximax=0.
      ximaxsum=0.
      dbobmax=0.
      dbobmaxs=0.
      dnonmax=0.
      write(*,*) 'psitot=',psitot
      do i=0.1*nosurf,nosurf
         if(i.gt.5.and.i.lt.nosurf-2)then
            do jts=1,mth/2
               hhhdbobs=0.
               hhhdbob=0.
               hhhdnon=0.
               hhhximax=0.
               do im=1,mt
                  i0=max0(1,i-1)
                  i0=min0(nosurf-3,i0)
                  theta1=(min+im-1)*2.*3.1415926*real(jts-1)/real(mth)
                  hhhdbobs=sqrt(grpssq2d(i,jts))/b2d(i,jts)**2
     &                 *((min+im-1-ntae*qoo(i))*eigfun(i,im,3)
     &                 /xjacob2d(i,jts)-eigfun(i,im,1)*shear2d(i,jts))
                  dbobmaxs=amax1(dbobmaxs,2.*abs(hhhdbobs))
                  hhhdbob=hhhdbob+sqrt(grpssq2d(i,jts))/b2d(i,jts)**2
     &                 *((min+im-1-ntae*qoo(i))*eigfun(i,im,3)
     &                 /xjacob2d(i,jts)-eigfun(i,im,1)*shear2d(i,jts))
                  hhhdnon=hhhdnon+eigfun(i,im,1)*cos(theta1)
                  hhhximax=hhhximax+eigfun(i,im,1)*cos(theta1)
     &                     /sqrt(grpssq2d(i,1))
               enddo
               call deriv4(ppsi(i0),rhoprf(i0),ppsi(i),dyn)
               hhhdnon=abs(hhhdnon*dyn/psitot/rhoprf(i))
               dnonmax=amax1(dnonmax,hhhdnon)
               dbobmax=amax1(dbobmax,2.*abs(hhhdbob))
               ximaxsum=amax1(hhhximax,ximaxsum)
ckg
            enddo
         endif
         do im=1,mt
            if(abs(eigfun(i,im,1)).gt.1.e-3.and.min_f.gt.(min+im-1)
     &           )min_f=min+im-1
            if(abs(eigfun(i,im,1)).gt.1.e-3.and.max_f.lt.(min+im-1)
     &           )max_f=min+im-1
            i0=max0(1,i-1)
            i0=min0(nosurf-3,i0)
            m1=min+im-1
            m2=m1**2
            if(i.gt.5)ximax=amax1(ximax,abs(eigfun(i0,im,1))
     &           /sqrt(grpssq2d(i,1)))
            do i3=1,3
               call funder4(ppsi(i0),eigfun(i0,im,i3),ppsi(i),dyn,yn)
               dyn=dyn/amax1(abs(yn),.01)
               dyn=(dyn/psitot)**2
               do jts=1,mth+2
cben include only poloidal contribution to the perpendicular wave vector
c    expression 
                  xkpop(i,im,jts,i3)=sqrt(grthsq2d(i,jts)*m2)
c this part includes both radial and poloidal: may be a bad approsm due to strong contributions of small Xi regions!! 
c                  xkpop(i,im,jts,i3)=sqrt(grthsq2d(i,jts)*m2+dyn
c     &                 *grpssq2d(i,jts))
               enddo
            enddo
         enddo
      enddo
      max_f=min0(max_f,abs(10*ntor))
      ximax=ximaxsum
      return
      end
c*********************************************************
      subroutine ppcorn(a,b,iconf,psimax,psimin,rmax,rmin
     $     ,rf,sigm)
      external d1eq0,psieq0,d2eq0
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      dimension rf(2,6)
      iconf=1
      a=0.
      b=0.
      sig=sigm
      if(sigm.lt.0.) go to 10
      call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rmax
     $     ,r0,rvpar,sigm,v,z,z0,xhp)
c      write(*,*) '1 psior,psimax,rmax ',psior,psimax,rmax
      if(psior.ge.psimax) then
         iconf=0
         return
      endif
      if(rvpar.gt.r0)then
         fpsi=psio(dpsior,ddpsior,dzpsi,rvpar,r0,z0,z0,am)
c         write(*,*) ' fpsi ',fpsi
         if(fpsi.lt.qpf) then
            rf(1,4)=rvpar
            sig=-1.
            call root1(rvpar,rmax,eps,0.,a,psieq0)
            sig=1.
         else
call root1(rvpar,rmax,eps,0.,rf(1,4),d1eq0)
c            write(*,*)'rvpar,rmax,eps,rf(1,4),d1eq0',rvpar,rmax,eps,rf(1,4)
c            if(rf(1,4).lt.rvpar) rf(1,4)=rvpar
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior
     $           ,rf(1,4),r0,rvpar,sigm,v,z,z0,xhp)
            fpsi=rf(2,4)
cpsio(dpsior,ddpsior,dzpsi,rf(1,4),r0,z0,z0,am)
c            write(*,*)'fpsi,psorb',fpsi,psior
            if(psior.le.fpsi)then
               iconf=0
               return
            endif
            call root1(rvpar,rf(1,4),eps,0.,a,psieq0)
         endif
c         write(*,*)'rf(1,4),rmax,b',rf(1,4),rmax,b
         call root1(rf(1,4),rmax,eps,0.,b,psieq0)
         return
      endif
      if(rvpar.le.rmin)then
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rmin
     $        ,r0,rvpar,sigm,v,z,z0,xhp)
         if(psior.ge.psimin)then
            iconf=0
            return
         endif
call root1(r0,rmax,eps,0.,rf(1,4),d1eq0)
c         if(rf(1,4).lt.r0) rf(1,4)=r0
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rf(1,4)
     $        ,r0,rvpar,sigm,v,z,z0,xhp)
         fpsi=rf(2,4)
cpsio(dpsior,ddpsior,dzpsi,rf(1,4),r0,z0,z0,am)
         if(psior.le.fpsi)then
            iconf=0
            return
         endif
         call root1(rmin,rf(1,4),eps,0.,a,psieq0)
         call root1(rf(1,4),rmax,eps,0.,b,psieq0)
         return
      endif
      fpsi=psio(dpsior,ddpsior,dzpsi,rvpar,r0,z0,z0,am)
c      write(*,*)'fpsi ',fpsi
      if(fpsi.gt.qpf) then
c         r58=amax1(rvpar,rmin)
c         write(*,*) rvpar,rmin
call root1(r58,rmax,eps,0.,rf(1,4),d1eq0)
c         write(*,*) rf(1,4)
c         if(rf(1,4).lt.r0) rf(1,4)=r0
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rf(1,4)
     $        ,r0,rvpar,sigm,v,z,z0,xhp)
         fpsi=rf(2,4)
cpsio(dpsior,ddpsior,dzpsi,rf(1,4),r0,z0,z0,am)
c         write(*,*) psior,fpsi
         if(psior.le.fpsi)then
            iconf=0
            return
         endif
         r58=amax1(rvpar,rmin)
         call root1(r58,rf(1,4),eps,0.,a,psieq0)
         call root1(rf(1,4),rmax,eps,0.,b,psieq0)
      else
         r58=amax1(rvpar,rmin)
call root1(r58,r0,eps,0.,rf(1,2),d2eq0)
         if(rf(1,2).lt.rvpar.or.d1eq0(rf(1,2)).gt.0.) then
            sig=-1.
            call root1(r58,rmax,eps,0.,a,psieq0)
            sig=sigm
            call root1(r0,rmax,eps,0.,b,psieq0)
         else
call root1(r58,rf(1,2),eps,0.,rf(1,4),d1eq0)
c            if(rf(1,4).lt.r0) rf(1,4)=r0
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior
     $           ,rf(1,4),r0,rvpar,sigm,v,z,z0,xhp)
            fpsi=rf(2,4)
cpsio(dpsior,ddpsior,dzpsi,rf(1,4),r0,z0,z0,am)
            if(psior.le.fpsi)then
c it resolves nothing
               sig=-1.
               call root1(r58,rf(1,4),eps,0.,a,psieq0)
               sig=sigm
               call root1(r0,rmax,eps,0.,b,psieq0)
               if(a.ne.0..or.b.ne.0.) write(*,*) 'nothing to resolve',a
     &              ,b
            else
c it resolves trapped orbits
               sig=-1.
               call root1(r0,rmax,eps,0.,a,psieq0)
               sig=sigm
               call root1(r0,rmax,eps,0.,b,psieq0)
            endif
c
         endif
      endif
      return
c--------------   sigm < 0  ---------------           ----
 10   if(rvpar.ge.r0) then
         sig=1.
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rmax
     $        ,r0,rvpar,sig,v,z,z0,xhp)
         if(psior.ge.psimax) then
            iconf=0
            return
         endif
c         r58=amax1(rvpar,rmin)
call root1(r58,rmax,eps,0.,rf(1,4),d1eq0)
c         if(rf(1,4).lt.rvpar) rf(1,4)=rvpar
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rf(1,4)
     $        ,r0,rvpar,sig,v,z,z0,xhp)
         fpsi=rf(2,4)
cpsio(dpsior,ddpsior,dzpsi,rf(1,4),r0,z0,z0,am)
         if(psior.le.fpsi) iconf=0
         return
      endif
      if(rvpar.le.rmin) then
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rmin
     $        ,r0,rvpar,sigm,v,z,z0,xhp)
c Check for left limiter intersection at R=Rmin
         if(psior.ge.psimin)then
            iconf=0
            return
         endif
      endif
c      r58=amax1(rvpar,rmin)
call root1(r58,r0,eps,0.,rf(1,2),d2eq0)
      if(rf(1,2).lt.rvpar.or.d1eq0(rf(1,2)).le.0.) then
c
         if(rf(1,3).lt.rvpar) rf(1,3)=rvpar
         hrf=max(rf(1,3),rmin)
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,hrf
     $        ,r0,rvpar,sigm,v,z,z0,xhp)
         fpsi=psiofpp2(hrf,r0,z0,am)
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior1,rf(1,1)
     $        ,r0,rvpar,sigm,v,z,z0,xhp)
c         fpsi=rf(2,1)
         if(psior.le.fpsi.and.psior1.gt.rf(2,1))then
            rf(1,3)=max(rf(1,3),rmin)
            call root1(rf(1,1),rmax,eps,0.,b,psieq0)
            call root1(rf(1,3),rf(1,1),eps,0.,a,psieq0)
            return
         endif
c
         sig=1.
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rmax
     $        ,r0,rvpar,sig,v,z,z0,xhp)
         if(psior.ge.psimax) then
            iconf=0
            return
         endif
call root1(r0,rmax,eps,0.,rf(1,4),d1eq0)
c         if(rf(1,4).lt.r0) rf(1,4)=r0
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rf(1,4)
     $        ,r0,rvpar,sig,v,z,z0,xhp)
         fpsi=rf(2,4)
cpsio(dpsior,ddpsior,dzpsi,rf(1,4),r0,z0,z0,am)
         if(psior.le.fpsi) iconf=0
         return
      else
call root1(rf(1,2),r0,eps,0.,rf(1,1),d1eq0)
c         if(rf(1,1).lt.rf(1,2)) rf(1,1)=r0
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,rf(1,1)
     $        ,r0,rvpar,sigm,v,z,z0,xhp)
         fpsi=rf(2,1)
cpsio(dpsior,ddpsior,dzpsi,rf(1,1),r0,z0,z0,am)
c         write(*,*) 'psior,fpsi,qpf,pmu,rf(1,1)',psior,fpsi,qpf,pmu,rf(1,1)
         if(psior.le.fpsi) then
            sig=1.
call root1(r0,rmax,eps,0.,rf(1,4),d1eq0)
c            if(rf(1,4).lt.r0) rf(1,4)=r0
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior
     $           ,rf(1,4),r0,rvpar,sig,v,z,z0,xhp)
            fpsi=rf(2,4)
cpsio(dpsior,ddpsior,dzpsi,rf(1,4),r0,z0,z0,am)
            if(psior.le.fpsi) iconf=0
            return
         endif
c         r58=amax1(rvpar,rmin)
call root1(r58,rf(1,2),eps,0.,rf(1,3),d1eq0)
         if(rf(1,3).lt.rvpar) rf(1,3)=rvpar
         hrf=max(rf(1,3),rmin)
         call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior,hrf
     $        ,r0,rvpar,sigm,v,z,z0,xhp)
ckg         fpsi=psiofpp2(rf(1,3),r0,z0,am)
         fpsi=psiofpp2(hrf,r0,z0,am)
c     write(*,*) 'psior,fpsi,qpf,pmu,r5i',psior,fpsi,qpf,pmu,hrf
c     ,rf(1,1)
         if(psior.gt.fpsi)then
           sig=1.
c           write(*,*) '-----------------------------------------------'
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior
     $           ,rmax,r0,rvpar,sig,v,z,z0,xhp)
c            write(*,*) '--------------',psior,rmax,qpf,pmu
            if(psior.ge.psimax) then
               iconf=0
               return
            endif
call root1(r0,rmax,eps,0.,rf(1,4),d1eq0)
c            if(rf(1,4).lt.r0) rf(1,4)=r0
            call psorb(bm,b0,dpsior,ddpsior,eps,iter,qpf,pmu,psior
     $           ,rf(1,4),r0,rvpar,sig,v,z,z0,xhp)
            fpsi=rf(2,4)
cpsio(dpsior,ddpsior,dzpsi,rf(1,4),r0,z0,z0,am)
            if(psior.le.fpsi) iconf=0
c            write(*,*) 'psior,fpsi,qpf,pmu,r5i',psior,fpsi,qpf,pmu
c     $           ,rf(1,3),rf(1,1)
c            write(*,*) '1a,b,iconf ',a,b,iconf
            return
         else
            rf(1,3)=max(rf(1,3),rmin)
            call root1(rf(1,1),rmax,eps,0.,b,psieq0)
            call root1(rf(1,3),rf(1,1),eps,0.,a,psieq0)
         endif
c         write(*,*) '2a,b ',a,b
      endif         
c      write(*,*) rf(1,2),rf(1,4),rf(1,3)
      return
      end
c**********************************************************************
      subroutine critpnts(eps1,iter1,nmu,nv,qpfi,pmui,psimin
     $     ,psimax,r01,rfi,rmin,rmax,vi,xhp1,z01)
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      dimension rfi(2,6,nmu,nv),pmui(nmu),vi(nv)
      iter=iter1
      eps=eps1
      qpf=0.01
      r0=r01
      z0=z01
      xhp=xhp1
      do iv=1,nv
         v=vi(iv)
         do imu=1,nmu
            pmu=pmui(imu)
            call rvpar0(eps,qpf,pmu,rmax,rmin,r0,rvpar,z0,zvpar0,iter)
            rfi(1,5,imu,iv)=rvpar
            call turnpnts(psimax,psimin,rmax,rmin,rfi(1,1,imu,iv))
c            write(*,*) 'mu,v,p1-4',pmu,v,rfi(1,1,imu,iv),rfi(1,2,imu,iv)
c     $           ,rfi(1,3,imu,iv),rfi(1,4,imu,iv)
         enddo
      enddo
      return
      end
c**********************************************************************
c Here we determine some critical points for the drift orbit trajectory
c i.e. arrays with RF(1,i) and the value of the poloidal flux at those
c      points, i.e. RF(2,i), total 6 points for up and down parts:
c      RF(i2,1) - right tangential point for negative branch (sigm < 0).
c      RF(i2,2) - separates left and right tangential points for 
c                 negative branch(sigm < 0).
c      RF(i2,3) - left tangential point for negative branch (sigm < 0).
c      RF(i2,4) - position when the orbit is tangential for
c                 sigm >0
c      RF(i2,6) - =RF(i2,3)
      subroutine turnpnts(psimax,psimin,rmax,rmin
     $     ,rf)
      external d1eq0,d2eq0
      common/ppcor/iter,eps,qpf,pmu,r0,rvpar,sig,v,z0,xhp
      dimension rf(2,6)
c--------------   sigm > 0  ---------------           ----
      sig=1.
      rf(1,3)=amax1(rvpar,r0)
      call root1(rf(1,3),rmax,eps,0.,rf(1,4),d1eq0)
      rf(1,4)=amax1(rf(1,4),rvpar)
      rf(2,4)=psio(dpsior,ddpsior,dzpsi,rf(1,4),r0,z0,z0,am)
c--------------   sigm < 0  ---------------           ----
      sig=-1.
      rf(1,3)=amax1(rvpar,rmin)
      call root1(rf(1,3),r0,eps,0.,rf(1,2),d2eq0)
      if(rf(1,2).le.rvpar.or.d1eq0(rf(1,2)).le.0.) then
c     write(*,*) ' I`m out'
         if(rvpar.lt.r0) then
            rf(2,1)=max((rvpar+r0)/2,rmin)
            call root1(rf(2,1),r0,eps,0.,rf(1,1),d1eq0)
            if(rf(1,1).le.rvpar) rf(1,1)=r0
            rf(2,1)=psio(dpsior,ddpsior,dzpsi,rf(1,1),r0,z0,z0,am)
         else
            rf(1,1)=0.
            rf(2,1)=0.
         endif
         rf(2,2)=1.
         rf(1,3)=0.
         rf(2,3)=1.
         return
      endif
      rf(2,2)=psio(dpsior,ddpsior,dzpsi,rf(1,2),r0,z0,z0,am)
      call root1(rf(1,2),r0,eps,0.,rf(1,1),d1eq0)
      if(rf(1,1).lt.rf(1,2)) rf(1,1)=r0
      rf(2,1)=psio(dpsior,ddpsior,dzpsi,rf(1,1),r0,z0,z0,am)
      r58=amax1(rvpar,rmin)
      call root1(r58,rf(1,2),eps,0.,rf(1,3),d1eq0)
c      write(*,*) 'rvpar,rf(1,3)=',rvpar,rf(1,3)
      rf(1,3)=amax1(rvpar,rf(1,3))
      rf(1,6)=rf(1,3)
      rf(2,3)=psio(dpsior,ddpsior,dzpsi,rf(1,3),r0,z0,z0,am)
      return
      end
c**********************************************************************
      function kenter_pest(a,b)
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      include 'orb2d.par'
      kenter_pest=0
      if((a-rr(1,1))*(a-rr(1,mth/2+1)).lt.0..or.(b-rr(1,1))*(b-rr(1,mth
     &     /2+1)).lt.0.)kenter_pest=1
      return
      end
c**********************************************************************
      function k_pest(psior)
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      include 'orb2d.par'
      k_pest=0
      if(psior.lt.ppsi(2).and.ishft.eq.0)then
         k_pest=1
      else if(ishft.eq.0)then
         k_pest=-1
      endif
      return
      end
c******************************************************************
c Here B0R0= B0*R0 are vacuum values at the plasma magnetic axis
c     eps,qpf,pmu,r0,z0 are input
c     r1 is bounce tip
c     iter is 1 for 1 or > 1 for 2 iterations to find the banana tip
      subroutine rvpar0(eps,qpf,pmu,rma,rmi,r0,r1,z0,z1,iter)
      parameter(ntgr2=30)
      common/fowint/wdintfow(ntgr2,2),rintfow(ntgr2,2),tbintfow(ntgr2,2)
     &     ,klyu
      common/ffratpsi/eps1,qpf1,pmu1,r01,z01
      external fratpsi
      if(pmu.lt.rmi) then
         r1=pmu
         return
      endif
ckg this check is introduced for particles in high beta plasmas
      r1mx=ratpsi(max(qpf,.001),rma,r1atpsi)
      zr1=bfof(b,b0,r1mx-eps,r0,z0,z0,-1)
      if(b0.lt.0.05)then
c.and.iter.gt.1
         r1=rma
         return
      elseif(r1mx.gt.rma-eps.or.pmu*b/(b0*r0).gt.1.)then
cif(iter.gt.1)then
            r1=rma
            return
c         endif
      endif
      z1=zatpsir(eps,qpf,pmu,r0,z0)
      r1=pmu
      if(iter.gt.1.and.pmu.gt.r1atpsi.and.r1atpsi.gt.rmi)then
         eps1=eps
         qpf1=qpf
         pmu1=pmu
         r01=r0
         z01=z0
         call ROOT1((r1atpsi+pmu)/2.,r1mx-eps,eps,0.,r1,fratpsi)
         if(r1.eq.0.) r1=rma
         return
      endif
      if(iter.gt.1)then
         iter0=iter+4
      else
         iter0=2
      endif
      do i=1,iter0
         zr1=bfof(b,b0,r1,r0,z1,z0,0)
         if(b0.lt..1)then
            r1=pmu
            return
         endif
         r1=pmu*b*r1/(b0*r0)
         zr1=z1
         z1=zatpsirf(eps,qpf,r1,r0,z0,zr1,0.5)
      enddo
      return
      end
c**********************************************************************
      function fratpsi(r)
      common/ffratpsi/eps,qpf,pmu,r0,z0
      z=zatpsir(eps,qpf,r,r0,z0)
      fratpsi=bfof(b,b0,r,r0,z,z0,0)
      fratpsi=pmu*b/(b0*r0)-1
      return
      end
c**********************************************************************
c gives Major radius at 'low field' side at given Psi
      function ratpsi(xxpsi,rma,r1atpsi)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      if(xxpsi.gt.1..or.xxpsi.lt.0.)then
         ratpsi=rma
         r1atpsi=0.
         return
      endif
      ii0=indx(ppsi(1),nosurf,xxpsi)
      ratpsi=ynn(ppsi(ii0),rr(ii0,1),xxpsi)
      r1atpsi=ynn(ppsi(ii0),rr(ii0,mth/2+1),xxpsi)
      return
      end
c**********************************************************************
c gives Major radius  given Psi index
      function ratpsii(ii0,psii,bii)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      ii00=max(1,ii0)
      ii00=min(nosurf,ii00)
      ratpsii=rr(ii00,1)
      psii=ppsi(ii00)
      bii=b2d(ii00,1)
      return
      end
c**********************************************************************
c gives distrfunction in m^-1 x T/eV x 1/eV for START shot 36484 at
c t= 26.ms provided e in eV, p=P_phi*a_i/c in m, mu in eV/T
c derivatives with the respect to e and p
      function dfn_fun(dfde,dfdp,e,p1,xmu)
      dimension a(0:6),dp(4),xwrk(4),ywrk(4)
      data a/0,-41.959,3.878e+4,-3.285e+6,-5.167e+9,1.572e+12,-1.
     &     307e+14/,eedge/31500./
c     
      data pi/3.1415927/
c     Define profiles for the constants
c     
      xe0(p)=max(0.,4.82e+5*(1.79*abs(p)**0.1-1.))
c     
      xde(p)=exp(-12.81*(1.-789.41*p+9.24e+4*p**2))
c     
      xlam(p)=1.78e+4*( 1.+6.52*cos( 2.83e+3*pi*abs(p-7.5e-4)**1.5 ) )
c     
      xdlam(p)=8.46e+6*(1.82e-2-p)*exp(-1.11e+5*(p-2.5e-3)**2)
c     Define the distribution itself
      dfn_me(e,xmu,p)= (-(xmu*eedge/E-xlam(P))**2
     &     /xdlam(P)**2-(E-xe0(P))**2/xde(P)**2)
      dfn_p(p)= (a(0)+p*(a(1)+p*(a(2)+p*(a(3)+p*(a(4)+p*(a(5)+p*a(6))
     &     )))))
      p=-p1
c
      if(p.lt.0.000003245)then
         dfn_fun=0.
         dfde=0.
         dfdp=0.
         return
      endif
      dfde= dfn_p(p)
c      write(*,*) 'p,e,xmu,dfde',p,e,xmu,dfn_me(e,xmu,p)+dfde
c      stop
      hh=max(dfn_me(e,xmu,p)+dfde,-15.)
      hh=min(hh,15.)
      dfn_fun= exp(hh)
      xwrk(4)=0.5e-5
      xwrk(2)=p
      xwrk(1)=p-xwrk(4)
      xwrk(3)=p+xwrk(4)
      xwrk(4)=xwrk(3)+xwrk(4)
      hh=max(dfn_me(e,xmu,xwrk(1))+dfn_p(xwrk(1)),-15.)
      hh=min(hh,15.)
      ywrk(1)=exp(hh)
      ywrk(2)=dfn_fun
      hh=max(dfn_me(e,xmu,xwrk(3))+dfn_p(xwrk(3)),-15.)
      hh=min(hh,15.)
      ywrk(3)=exp(hh)
      hh=max(dfn_me(e,xmu,xwrk(4))+dfn_p(xwrk(4)),-15.)
      hh=min(hh,15.)
      ywrk(4)=exp(hh)
      call deriv4(xwrk(1),ywrk(1),xwrk(2),dfdp)
      dfdp=1.23e+13*dfdp
      xwrk(4)=xde(P)/30.
      xwrk(2)=e
      xwrk(1)=e-xwrk(4)
      xwrk(3)=e+xwrk(4)
      xwrk(4)=xwrk(3)+xwrk(4)
      hh=max(dfn_me(xwrk(1),xmu,p)+dfde,-15.)
      hh=min(hh,15.)
      ywrk(1)=exp(hh)
      hh=max(dfn_me(xwrk(3),xmu,p)+dfde,-15.)
      hh=min(hh,15.)
      ywrk(3)=exp(hh)
      hh=max(dfn_me(xwrk(4),xmu,p)+dfde,-15.)
      hh=min(hh,15.)
      ywrk(4)=exp(hh)
      call deriv4(xwrk(1),ywrk(1),xwrk(2),dfde)
      dfde=1.23e+13*dfde
      dfn_fun= 1.23e+13*dfn_fun
      return
      end
c**********************************************************************
c this provides better normalization for the anistropic distribution
c function in case of coinjected beam or ICRH minority ions
      subroutine nrmlzn(h,dlambda,ylambda0,yres)
      parameter (nint=100)
      dimension fint(nint), xint(nint)
      do i=1,nint
         xint(i)=(float(i-1)/float(nint-1))**0.8
         fint(i)=exp(-((1-ylambda0*h-xint(i)**2)/(h*dlambda))**2)
      enddo
      call ASIMP(xint,fint,nint,yres)
      yres=yres/2.
      return
      end
c**************************************************************
      function dfint(gamma,vstar3)
      real dfint,gamma,vstar3
      parameter (nint=100)
      dimension fint(nint),xint(nint)
      do i=1,nint
         xint(i)=(float(i-1)/float(nint-1))**1
         dfint=xint(i)**3
ckg note that the slowing down case normalization is when this integral
ckg is taken at v^4/(vstar3+dfint)
         fint(i)=xint(i)*dfint/(vstar3+dfint)
ckg the following factor is added to account for the Krook collisional 
ckg operator in the form f/tau_loss. So that according to Cordey'81 
ckg        tau_se/3tau_loss=(gamma-1)
ckg gamma=2 is for the flat distribution in velocity. 
ckg gamma=3 is one of the inverse df. ...
         fint(i)=fint(i)*(vstar3+dfint)**(gamma-1)
      enddo
      call ASIMP(xint,fint,nint,dfint)
      return
      end
c**************************************************************
c
