c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c......................................................................
      subroutine errmes ( nout,mesage )
c......................................................................
      write(*,100)mesage
      write(nout,100)mesage
  100 format(" error in subroutine",1xa8/)
      stop
      end
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      subroutine defolt(im1,ihsps,xfow)
c
      include 'clich1'
      include 'clich2'
c      common/switch/iswitch,itransp
c default parameters to be used in distribution function
      chidelt(1)=0.5
      chidelt(2)=0.0
      chidelt(3)=0.4
      chidelt(4)=0.2
      chidelt(5)=0.13
      chidelt(6)=0.2
c
      zero = 0.0
      one = 1.0
      two = 2.0
      four = 4.0
      pi = acos ( -one )
      pi2 = pi * pi
      twopi = two * pi
ckg please keep this number equal -10 to be consistent with the scripts
      ntor = -10
      qscale=1.0
      gamma = 0.0
ckg this value is not used anywhere in the stability part      gamma = 1.6667
      rho = one
      lmin = minm
      lmax = maxm
      mp0 = 41
      mp1 = 42
      mp2=  19
ckg      mp0 = 6hmapdsk
ckg      mp1 = 6hmpout1
ckg      mp2=  6hequou1
c
      xfow=1.
      im1=0
      if(tip.ne.'x') then
         ihsps=1
c         if(im1.ne.3) im1=0
      endif
      ipbmax=jb
      ipcmin=jcmin
      ipcmax=jcmax
      idet=8
      if(idet.gt.lam/5)write(*,*) 'WARNING your basis function number',
     & ' is more that you can afford. Increase LAM'
      klamb=lam
      klamc=lamc
      nnsurf=nogrid
c
      iodat = 16
      iotty = 5
c
      ltot=lmax-lmin+1
      ltotsq=ltot*ltot
      ltots2=ltotsq*2
      ltots8=ltotsq*8
      ktot=jt
      kmin=jmin
      np=jb
      npc=jcc
c
cccc   read abscissas and weight factors for Gaussian integration
c      if(lam.eq.16) iread=7
ckg      if(lam.eq.48) 
      iread=8
      do 51 k=1,lam
      read(iread,52) xali(k),wali(k)
   51 continue
      do 53 k=1,lamc
      xalic(k)=xali(k)
      walic(k)=wali(k)
   53 continue
   52 format(2e25.15)
      rewind(iread)
cc  set defolt input parameters
      call inparam
      pk=1.
      call inparam_read(im1,ihsps,xfow,pk)
c
c
c      write(0, 10 )
   10 format ( 1x, " default values.....  ",/,
     .    1x, " lmin=minm, lmax=maxm, ntor = ?, iden=1,  ",/,
     .    1x, " mp0=6hmapdsk, mp1=6hmpout1, mp2=6hequou1, ",/,
     . 1x, " alphar=0.9, prho=1., arho=1.5, qscale=1",//)
c
c
c
      return
      end
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c......................................................................
      subroutine input
c......................................................................
c
      include 'clich1'
      include 'clich2'
cor
      namelist /inp/ ntor,nnsurf,lmin,lmax,qscale
     &,ipbmax,klamb,ipcmin,ipcmax,klamc,idet
c
c
c     call talk ( " INP input desired variables..     " )
c
      write(*,inp )
c2d      read (5,inp )
      write(iodat,inp )
      write(17,inp )
c
c
      return
      end
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     1.3   disk data input.
c.......................................................................
c
      subroutine dskin
c.......................................................................
      include 'clich1'
      include 'clich2'
c
      integer*8 nadres,ierr,lhhh,length
      integer t, surf
      real dat
      data iflag / 0 /
c
c
c  ..... iref is index to set reference values of f, g, p, etc. in
c  .....  subroutine varymo  each time the mapping is read.
c
      iref = 0
c
c     1.3.2   read equilibrium quantities from disk storage.
c              first read any       useful information.
c
      length = nts * (nsf+ishft+1)
c
c      read the equilibrium functions if mapping done before.
c      nadres = 1
c      call zrd(mp1, ntitle(1),41,nadres+ishft,lgivup,ierr)
c      lhhh = 41
c      call getwa('mpout1', ntitle(1),nadres,lhhh,ierr )
c      write(*,'(20a)') ntitle
c
      mth1 = mth + 1
      mth2 = mth1 + 1
      dth = twopi / mth
      do 10 i=1,mth2
      tgrid(i)=(i-1)
   10 continue
      dpsi = one / ( nosurf-1 )
      psiv = ( nsrf - 1 ) * dpsi
      if(nsrf.eq.1) psiv=0.001
      psitot=(psilim-psimin)/twopi
      psif=psiv*psitot
      psifsq=psif**2
c
c      start reading from mapdsk
c
c
c         read interface grid points...
c
      lgivup=1
      nadres = 50
c      call zrd(mp0, xinf(1),mth2,nadres+ishft,lgivup,ierr)
      lhhh=mth2
      call getwa('mapdsk',xinf(1),nadres,lhhh,ierr)
      nadres = nadres + mth2
c      call zrd(mp0, zinf(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',zinf(1),nadres,lhhh,ierr)
c
c         finally get poloidal flux values...
c
c     4.1.2.2      calculate and store the theta integrals on psi
      nadres = nadres + mth2 + nsrf -1
c      call zrd(mp0, psbg,1,nadres+ishft,lgivup,ierr)
      lhhh=1
      call getwa('mapdsk',psbg,nadres,lhhh,ierr)
c
c
      r2     = r * r
      r4     = r2 * r2
      r6     = r4 * r2
c
c
      nadres = 50 + nsrf - 1
c      call zrd(mp1, p,1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', p, nadres,lhhh,ierr )
      nadres = nadres + nosurf
c      call zrd(mp1, pp,1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', pp, nadres,lhhh,ierr )
      nadres = nadres + nosurf
c      call zrd(mp1, q,1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', q, nadres,lhhh,ierr )
      nadres = nadres + nosurf
c      call zrd(mp1, qp,1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', qp, nadres,lhhh,ierr )
      nadres = nadres + nosurf
c      call zrd(mp1, g,1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', g, nadres,lhhh,ierr )
      nadres = nadres + nosurf
c      call zrd(mp1, gp,1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', gp, nadres,lhhh,ierr )
      nadres = nadres + nosurf
c      call zrd(mp1, f,1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', f, nadres,lhhh,ierr )
      nadres = nadres + nosurf
c      call zrd(mp1, fp,1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', fp, nadres,lhhh,ierr )
      nadres = nadres + nosurf
       nadres_yj=nadres          !!!!! This line is added by yjhu
      if(ishft.eq.-1)then
c         call zrd(mp1, rrg,1,nadres+ishft,lgivup,ierr)
         call getwa('mpout1', rrg, nadres,lhhh,ierr )
         if(nsrf.eq.1)then
            nadres = nadres + mth2 + nsrf -1
c            call zrd(mp0, psbg,1,nadres+ishft,lgivup,ierr)
            call getwa('mapdsk', psbg, nadres,lhhh,ierr )
c            call zrd(mp1, rgrid(1),nosurf,nadres+ishft,lgivup,ierr)
            nadres=nadres_yj        !!!! This line is added by yjhu
            lhhh=nosurf
            call getwa('mpout1', rgrid(1), nadres,lhhh,ierr )
c            write(*,*) 'check rgrid',rgrid
         endif
      else
         pp=pp/f
         qp=qp/f
         gp=gp/f
         fp=fp/f
         if(nsrf.eq.1)then
            nadres = 50 + mth2 + mth2 + nsrf -1
c            call zrd(mp0, rgrid(1),nosurf ,nadres+ishft,lgivup,ierr)
            call getwa('mapdsk', rgrid(1), nadres,lhhh,ierr )
            do iii=1,nosurf
               rgrid(iii)=sqrt((rgrid(iii)-psimin+psilim)/(psilim
     &              -psimin))
            enddo
c            write(*,*) 'check rgrid',rgrid
         endif
         rrg= rgrid(nsrf)
      endif
cccc..  rescale q & g here
c      g2=g*g
c      q2=q*q
c      gnew2=g2+gconst
c      gnew=sqrt(gnew2)
c      gpnew=g*gp/gnew
c      qnew=r*gnew/f
c      qpnew=r*(f*gpnew-fp*gnew)/(f**2)
c      q=qnew
c      g=gnew
c      qp=qpnew
c      gp=gpnew
cccc...   end of rescalings of q & g
ccc  set phot=photp=0.0 in computing coeficients due to perturbative analysis  
cccc  define hot & cold components pressures
      phot=0.0
      photp=0.0
      pcore=p-phot
      pcorep=pp-photp
c
      gammapc=gamma*pcore
      g2=g*g
      ggp=g*gp
      rg=r*g
      r2g2=r2*g2
      qprf(nsrf)=q
      pprf(nsrf)=p
      gprf(nsrf)=g
      rgrid(nsrf)=rrg
      if(nsrf.eq.1) rgrid(nsrf)=0.0
c      write(*,*) 'check rgrid',rgrid(nsrf),psbg,lhhh,mth2,nsrf
c
ccc...   calculate and store theta quantities
c
      call theint
c	write(16,*)'betah profile i,betah(i)'
c	write(16,*)nsrf,betah(nsrf)
c
c     call zch(iomap1,if1,999)
c     call zch(iomode,if2,999)
c
  100 format ( 1x, 20a4, /, 1x,
     .       "equilibrium from disk calculated on date", a8, /,
     .       " nosurf = ", i4, " mth = ", i4, / )
c
      return
c  999 call errmes ( iodat, 5hdskin )
      end
c.......................................................................
c
      subroutine theint
      include 'clich1'
      include 'clich2'
      integer*8 nadres,ierr,lhhh,length
c.......................................................................
c
c        for each t construct an array of theta values on the magnetic surface
c
c
      length = nts * (nsf+ishft+1)
c      length = nts * (nsf+ishft+1)
      ladres = 50 + 2 * mth2 + nosurf
      nadres = ladres + ( nsrf - 1 ) * nts
      lgivup = 1
c      call zrd(mp0,grpssq(1),mth2,nadres+ishft,lgivup,ierr)
      lhhh=mth2
      call getwa('mapdsk',grpssq(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp0,xsq(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',xsq(1),nadres,lhhh,ierr)
      if(nosurf.eq.nsrf.and.xsq(1).eq.0.)then
c         call zrd(mp0,xsq(1),mth2,nadres+ishft-nts,lgivup,ierr)
         nadres=nadres-nts
         call getwa('mapdsk',xsq(1),nadres,lhhh,ierr)
         nadres=nadres+nts
      endif
      nadres = nadres + length
c      call zrd(mp0,grthsq(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',grthsq(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp0,grpsth(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',grpsth(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp0,gptdth(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',gptdth(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp0,xsqdps(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',xsqdps(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp0,gpsdth(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',gpsdth(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp0,xsqdth(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',xsqdth(1),nadres,lhhh,ierr)
      if(ishft.eq.-1)then
      nadres = nadres + length
c      call zrd(mp0,xjacob(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',xjacob(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp0,xjprym(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',xjprym(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp0,delta(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',delta(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp0,qdelp(1),mth2,nadres+ishft,lgivup,ierr)
      call getwa('mapdsk',qdelp(1),nadres,lhhh,ierr)
      endif
c
cccc
cccc  the followings will be bypassed
c      go to 199
      ladres = 50 + (9-1-ishft) * nosurf
      nadres = ladres + (nsrf-1) * nts
      lhhh=mth1
c      call zrd(mp1,xa(1),mth1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1',xa(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp1,za(1),mth1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1',za(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp1,xtheta(1),mth1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1',xtheta(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp1,ztheta(1),mth1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1',ztheta(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp1,xpsi(1),mth1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1',xpsi(1),nadres,lhhh,ierr)
      nadres = nadres + length
c      call zrd(mp1,zpsi(1),mth1,nadres+ishft,lgivup,ierr)
      call getwa('mpout1',zpsi(1),nadres,lhhh,ierr)
c  199 continue
c
c
ccc..   prepare surface arrays
      mthd2p1=mth/2+1
      do 5 i = 1, mth2
            grpssq(i)=abs(grpssq(i))
         if(ishft.eq.0)then
            grpssq(i)=grpssq(i)*f*f
            grpsth(i)=grpsth(i)*f
            gptdth(i)=gptdth(i)*f
            xsqdps(i)=xsqdps(i)/f
            gpsdth(i)=gpsdth(i)*f*f
            if(i.lt.mth2)then
               xpsi(i)=xpsi(i)/f
               zpsi(i)=zpsi(i)/f
            endif
            xjprym(i)=upsiln/r/twopi/f
            xjacob(i)=xsq(i)*xjprym(i)
            xjprym(i)=xsqdps(i)*xjprym(i)
            delta(i) =0.
            qdelp(i) =0.
         endif
      xsdpxs(i)=xsqdps(i)/xsq(i)
      xsdtxs(i)=xsqdth(i)/xsq(i)
      gptgpp(i)=grpsth(i)/grpssq(i)
      gpdtgp(i)=gpsdth(i)/grpssq(i)
      gptdtgp(i)=gptdth(i)/grpssq(i)
      xjsq(i)=xjacob(i)*xjacob(i)
      xjdtxj(i)=0.5*(mjx*xsdtxs(i)-njg*gpdtgp(i))
      xjprxj(i)=xjprym(i)/xjacob(i)
      bsq(i)=(grpssq(i)+r2g2)/xsq(i)
      bf(i)=sqrt(bsq(i))
      bsdtbs(i)=gpsdth(i)/(bsq(i)*xsq(i))-xsdtxs(i)
      curvs(i)=rg*bsdtbs(i)/(xjacob(i)*bsq(i))
      sigma(i)=-r*(g*pp/bsq(i)+gp)
      grgrpp(i)=xjprxj(i)+(xjdtxj(i)-gpdtgp(i))*gptgpp(i)+gptdtgp(i)
      shear(i)=r*gp+(grgrpp(i)-xsdpxs(i)-gptgpp(i)*xsdtxs(i))*rg
      shear(i)=-shear(i)/xsq(i)
      curvn(i)=r*gp-(xsdpxs(i)+gptgpp(i)*xsdtxs(i))*rg
      curvn(i)=rg*curvn(i)/xsq(i)+pp
      curvn(i)=(curvn(i)-grgrpp(i)*grpssq(i)/xsq(i))/bsq(i)
cben      if(nsrf.eq.130)write(*,'(i3,6f10.4)')i,q,curvn(i),curvs(i)
c     &     ,grpsth(i),f,fp
      qdelt(i)=rg*xjacob(i)/xsq(i)-q
    5 continue
cben      if(nsrf.eq.130)stop
ccc..   odd parity quantities
c     gptgpp(1)=0.
c     gpdtgp(1)=0.
c     xsdtxs(1)=0.
c     xjdtxj(1)=0.
c     bsdtbs(1)=0.
c     curvs(1)=0.
c     qdelp(1)=0.
c     gptgpp(mthd2p1)=0.
c     gpdtgp(mthd2p1)=0.
c     xsdtxs(mthd2p1)=0.
c     xjdtxj(mthd2p1)=0.
c     bsdtbs(mthd2p1)=0.
c     curvs(mthd2p1)=0.
c     qdelp(mthd2p1)=0.
c
ccc...  check qp by integrating over xjacob*shear
      rqp(nsrf)=0.0
      agrgrpp(nsrf)=0.
      agptdtgp(nsrf)=0.
      agptdth(nsrf)=0.0
      do 998 i=1,mth
      rqp(nsrf)=rqp(nsrf)+xjacob(i)*shear(i)
      agrgrpp(nsrf)=agrgrpp(nsrf)+grgrpp(i)*xjacob(i)
      agptdtgp(nsrf)=agptdtgp(nsrf)+gptdtgp(i)
      agptdth(nsrf)=agptdth(nsrf)+gptdth(i)
  998 continue
      rqp(nsrf)=-rqp(nsrf)*dth/twopi
ccc..  set up radial arrays
      rshear(nsrf)=shear(1)
      rcurvn(nsrf)=curvn(1)
      rsigma(nsrf)=sigma(1)
      rgrgrpp(nsrf)=grgrpp(1)
      rxsdpxs(nsrf)=xsdpxs(1)
      rgptdtgp(nsrf)=gptdtgp(1)
      rxjsq(nsrf)=xjsq(1)
      rxjprxj(nsrf)=xjprxj(1)
c
      rcurvs(nsrf)=curvs(mth/4+1)
      rxjdtxj(nsrf)=xjdtxj(mth/4+1)
      rgpdtgp(nsrf)=gpdtgp(mth/4+1)
      rgptgpp(nsrf)=gptgpp(mth/4+1)
c
cccc   compute bsqinv(nsrf) etc. on each surface for surface averaging
c
      volth=0.0
      beh=0.0
      bsqav(nsrf)=0.0
      do 3 i=1,mth
      volth=volth+xjacob(i)
      beh=beh+xjacob(i)/bsq(i)
      bsqav(nsrf)=bsqav(nsrf)+bsq(i)
    3 continue
      bsqav(nsrf)=bsqav(nsrf)/mth
      bsqinv(nsrf)=2.*beh
      volr(nsrf)=volth
      return
c  999 call errmes(iodat,6htheint)
      end
c.......................................................................
      subroutine matrix
c.......................................................................
      include 'clich1'
      include 'clich2'
ccc..  construct coefficients of the mhd stability eqns. and take fft
ccc..  to fourier harmonic space.
      rho=rhoprf(nsrf)
      omr=oma2*rho
      xx2=psitot
      xx2=xx2*2.*rrg
      do 5 i=1,mth
ccc.. slow-sound wave eqn
      as(i)=xjacob(i)*curvs(i)
      ad1(i)=xjacob(i)*(1.+gammapc/bsq(i))
      ad3(i)=gammapc/(omr*xjacob(i)*bsq(i))
      ad2(i)=-ad3(i)*(xjdtxj(i)+bsdtbs(i))
      ap(i)=-xjacob(i)/bsq(i)
      at(i)=-xjacob(i)*(curvn(i)-photp/bsq(i))
c
ccc.. shear alfven wave eqn
      ss1(i)=xjacob(i)*omr*grpssq(i)/bsq(i)
      ss3(i)=grpssq(i)/(xjacob(i)*bsq(i))
      ss2(i)=ss3(i)*(gpdtgp(i)-xjdtxj(i)-bsdtbs(i))
      sd(i)=gammapc*xjacob(i)*curvs(i)
      sp1(i)=xjacob(i)*ntor
      sp2(i)=-curvs(i)*xjacob(i)
      sp3(i)=rg/bsq(i)
      st1(i)=grpssq(i)*shear(i)/bsq(i)
      st2(i)=pcorep*sp2(i)
      st2i(i)=-ntor*photp*xjacob(i)
      st3(i)=st1(i)-sigma(i)-photp*rg/bsq(i)
c
ccc..  pressure radial derivative eqn
      gp0i(i)=-(qdelp(i)+gptgpp(i)*(qdelt(i)+q))*ntor*xx2
      gp0(i)=curvn(i)*xx2
      gp1(i)=-gptgpp(i)*xx2
      gt0(i)=omr/grpssq(i)*xx2
      gt1(i)=pcorep*curvn(i)
      gt1(i)=gt1(i)+(sigma(i)-shear(i)*grpssq(i)/bsq(i))*shear(i)
      gt1(i)=gt1(i)*xx2
      gt1i(i)=-photp*gp0i(i)
      gt3(i)=1./(xjsq(i)*grpssq(i))*xx2
      gt2(i)=-(xjdtxj(i)+gpdtgp(i))*gt3(i)+photp*gptgpp(i)*xx2
      gs(i)=(grpssq(i)*shear(i)/bsq(i)-sigma(i))/xjacob(i)*xx2
      gd(i)=gammapc*curvn(i)*xx2
c
ccc..  xipsi radial derivative eqn
      dti(i)=gp0i(i)
      dt(i)=-(grgrpp(i)+curvn(i)-photp/bsq(i))*xx2
      dt1(i)=-gptgpp(i)*xx2
      dp(i)=-1./bsq(i)*xx2
      ds1(i)=-ntor*xx2
      ds2(i)=rg/(xjacob(i)*bsq(i))*xx2
      dd(i)=-gammapc/bsq(i)*xx2
    5 continue
ccc...  fft and store arrays
ccc...  matrix coefficients for adiabatic pressure law
c
      call fourr(as,wi,mth)
      call mat0(cas,as,wi,mth,ktot)
      call matm2(cas,ktot)
c
      call fourr(ad1,wi,mth)
      call mat0(cad1,ad1,wi,mth,ktot)
      call fourr(ad2,wi,mth)
      call fourr(ad3,wii,mth)
      call mat4(cad2,ad2,wi,ad3,wii,q,ntor,mth,kmin,ktot)
c
      call fourr(ap,wi,mth)
      call mat0(cap,ap,wi,mth,ktot)
c
      call fourr(at,wi,mth)
      call mat0(cat,at,wi,mth,ktot)
ccc..  matrix coefficients for surface component eqn
c
      call fourr(ss1,wi,mth)
      call mat0(css1,ss1,wi,mth,ktot)
c
      call fourr(ss2,wi1,mth)
      call fourr(ss3,wi2,mth)
      call mat4(css2,ss2,wi1,ss3,wi2,q,ntor,mth,kmin,ktot)
c
      call fourr(sd,wi,mth)
      call mat0(csd,sd,wi,mth,ktot)
      call matm1(csd,ktot)
c
      call fourc(sp2,sp1,mth)
      call fourr(sp3,wii,mth)
      call mat2(csp,sp2,sp1,sp3,wii,q,ntor,mth,kmin,ktot)
      call matm1(csp,ktot)
c
      call fourr(st1,wi,mth)
      call fourc(st2,st2i,mth)
      call mat3(cst1,st1,wi,st2,st2i,q,ntor,mth,kmin,ktot)
      call fourr(st3,wi,mth)
      call mat1(cst,st3,wi,q,ntor,mth,kmin,ktot)
      call matadd(cst,cst1,ktot)
      call matm1(cst,ktot)
cccc...  matrix coefficients for surface components
c
      do 10 k2=1,ktot
      do 10 k1=1,ktot
      cmsd(k1,k2)=css2(k1,k2)
      cmsd(k1+ktot,k2)=cas(k1,k2)
      cmsd(k1,k2+ktot)=csd(k1,k2)
      cmsd(k1+ktot,k2+ktot)=cad1(k1,k2)
      cmtp(k1,k2)=cst(k1,k2)
      cmtp(k1+ktot,k2)=cat(k1,k2)
      cmtp(k1,k2+ktot)=csp(k1,k2)
      cmtp(k1+ktot,k2+ktot)=cap(k1,k2)
   10 continue
cccc...  matrix coefficients for normal components
c
      call fourc(gp0,gp0i,mth)
      call fourr(gp1,wii,mth)
      call mat2(cgp,gp0,gp0i,gp1,wii,q,ntor,mth,lmin,ltot)
c
      call fourr(gt0,wi,mth)
      call mat0(cgt0,gt0,wi,mth,ltot)
c
      call fourc(gt1,gt1i,mth)
      call mat0(cgt1,gt1,gt1i,mth,ltot)
      call fourr(gt2,wi,mth)
      call fourr(gt3,wii,mth)
      call mat4(cgt,gt2,wi,gt3,wii,q,ntor,mth,lmin,ltot)
      call matadd(cgt,cgt1,ltot)
c
      call fourr(gs,wi,mth)
      call mat1(cgs,gs,wi,q,ntor,mth,lmin,ltot)
      call matm2(cgs,ltot)
c
      call fourr(gd,wi,mth)
      call mat0(cgd,gd,wi,mth,ltot)
ccc..  matrix coefficients for div(xpsi) eqn
c
      call fourc(dt,dti,mth)
      call fourr(dt1,wii,mth)
      call mat2(cdt,dt,dti,dt1,wii,q,ntor,mth,lmin,ltot)
c
      call fourr(dp,wi,mth)
      call mat0(cdp,dp,wi,mth,ltot)
c
      call fourr(ds1,wi,mth)
      call fouri(wii,ds2,mth)
      call mat2(cds,ds1,wi,wii,ds2,q,ntor,mth,lmin,ltot)
c
      call fourr(dd,wi,mth)
      call mat0(cdd,dd,wi,mth,ltot)
c
ccc..  set real matrices
      call rmatrix(rad2,cad2,jt)
      call rmatrix(rss1,css1,jt)
      call rmatrix(rmsd,cmsd,jt2)
      call rmatrix(rmtp,cmtp,jt2)
      call rmatrix(rgt0,cgt0,mt)
      call rmatrix(rgp,cgp,mt)
      call rmatrix(rgt,cgt,mt)
      call rmatrix(rgs,cgs,mt)
      call rmatrix(rgd,cgd,mt)
      call rmatrix(rdp,cdp,mt)
      call rmatrix(rdt,cdt,mt)
      call rmatrix(rds,cds,mt)
      call rmatrix(rdd,cdd,mt)
c
      return
      end
c........................................................
      subroutine fourc(xr,xi,mth)
c........................................................
ccc...  take fft & store
      dimension xr(mth),xi(mth)
c     write(16,101) (xr(i),i=1,mth)
c     write(16,101) (xi(i),i=1,mth)
      call cpft(xr,xi,mth,1,-1)
c     write(16,101) (xr(i),i=1,mth)
c     write(16,101) (xi(i),i=1,mth)
  101 format(10e12.4)
      do 2 k=1,mth
      xr(k)=xr(k)/mth
      xi(k)=xi(k)/mth
    2 continue
      return
      end
c........................................................
      subroutine fouri(xr,xi,mth)
c........................................................
ccc...  take fft & store
      dimension xr(mth),xi(mth)
      do 1 i=1,mth
      xr(i)=0.0
    1 continue
c     write(16,101) (xr(i),i=1,mth)
c     write(16,101) (xi(i),i=1,mth)
      call cpft(xr,xi,mth,1,-1)
      do 2 k=1,mth
      xr(k)=xr(k)/mth
      xi(k)=xi(k)/mth
    2 continue
      return
      end
      subroutine fourr(xr,xi,mth)
c........................................................
ccc...  take fft & store
      dimension xr(mth),xi(mth)
      do 1 i=1,mth
      xi(i)=0.0
    1 continue
c     write(16,101) (xr(i),i=1,mth)
c     write(16,101) (xi(i),i=1,mth)
      call cpft(xr,xi,mth,1,-1)
      do 2 k=1,mth
      xr(k)=xr(k)/mth
      xi(k)=xi(k)/mth
    2 continue
      return
      end
c........................................................
      subroutine mat0(am,xr,xi,mth,ltot)
c........................................................
      complex am(ltot,ltot)
      dimension xr(mth),xi(mth)
      do 1 l2=1,ltot
      do 2 l1=1,ltot
      l=l1-l2+1
      if(l.lt.1) l=mth+l
      am(l1,l2)=cmplx(xr(l),xi(l))
    2 continue
    1 continue
      return
      end
c........................................................
      subroutine mat1(am,xr,xi,q,ntor,mth,lmin,ltot)
c........................................................
cccc...  am=(xr,xi)*(d/dt-i*n*q)
      complex am(ltot,ltot)
      dimension xr(mth),xi(mth)
      do 1 l2=1,ltot
      al2=l2+lmin-1
      aq1=al2-ntor*q
      do 2 l1=1,ltot
      l=l1-l2+1
      if(l.lt.1) l=mth+l
      am(l1,l2)=cmplx(xr(l),xi(l))*cmplx(0.,aq1)
    2 continue
    1 continue
      return
      end
c........................................................
      subroutine mat2(am,xr1,xi1,xr2,xi2,q,ntor,mth,lmin,ltot)
c........................................................
ccc..   am=(xr1,xi1)+(xr2,xi2)*(d/dt-i*n*q)
      complex am(ltot,ltot)
      dimension xr1(mth),xi1(mth),xr2(mth),xi2(mth)
  101 format(10e12.4)
      do 1 l2=1,ltot
      al2=l2+lmin-1
      aq1=al2-ntor*q
      do 2 l1=1,ltot
      l=l1-l2+1
      if(l.lt.1) l=mth+l
      am(l1,l2)=cmplx(xr2(l),xi2(l))*cmplx(0.,aq1)
     .          +cmplx(xr1(l),xi1(l))
    2 continue
    1 continue
      return
      end
c........................................................
      subroutine mat3(am,xr1,xi1,xr2,xi2,q,ntor,mth,lmin,ltot)
c........................................................
cccc..  am=(d/dt)(xr1,xi1)+(xr2,xi2)
      complex am(ltot,ltot)
      dimension xr1(mth),xi1(mth),xr2(mth),xi2(mth)
      do 1 l2=1,ltot
      do 2 l1=1,ltot
      l0=l1-l2
      l=l0+1
      if(l.lt.1) l=mth+l
      am(l1,l2)=cmplx(xr1(l),xi1(l))*cmplx(0.,1.)*l0
     .          +cmplx(xr2(l),xi2(l))
    2 continue
    1 continue
      return
      end
c........................................................
      subroutine mat4(am,xr1,xi1,xr2,xi2,q,ntor,mth,lmin,ltot)
c........................................................
ccc..  am=(xr1,xi1)*(d/dt-i*n*q)+(xr2,xi2)*(d/dt-i*n*q)**2
      complex am(ltot,ltot)
      dimension xr1(mth),xi1(mth),xr2(mth),xi2(mth)
      do 1 l2=1,ltot
      al2=l2+lmin-1
      aq1=al2-ntor*q
      aq2=aq1*aq1
      do 2 l1=1,ltot
      l=l1-l2+1
      if(l.lt.1) l=mth+l
      am(l1,l2)=cmplx(xr1(l),xi1(l))*cmplx(0.,aq1)
     .          -cmplx(xr2(l),xi2(l))*aq2
    2 continue
    1 continue
      return
      end
c........................................................
      subroutine mat5(am,xr1,xi1,xr2,xi2,q,ntor,mth,lmin,ltot)
c........................................................
cc..   am=i*(xr1,xi1)+(xr2,xi2)*(d/dt)
      complex am(ltot,ltot)
      dimension xr1(mth),xi1(mth),xr2(mth),xi2(mth)
      do 1 l2=1,ltot
      al2=l2+lmin-1
      do 2 l1=1,ltot
      l=l1-l2+1
      if(l.lt.1) l=mth+l
      am(l1,l2)=cmplx(xr1(l),xi1(l))*cmplx(0.,1.)
     .          +cmplx(xr2(l),xi2(l))*cmplx(0.,al2)
    2 continue
    1 continue
      return
      end
c........................................................
      subroutine mat6(am,xr1,xi1,xr2,xi2,q,ntor,mth,lmin,ltot)
c........................................................
cc..   am=(xr1,xi1)+(xr2,xi2)*(d/dt)
      complex am(ltot,ltot)
      dimension xr1(mth),xi1(mth),xr2(mth),xi2(mth)
      do 1 l2=1,ltot
      al2=l2+lmin-1
      do 2 l1=1,ltot
      l=l1-l2+1
      if(l.lt.1) l=mth+l
      am(l1,l2)=cmplx(xr1(l),xi1(l))
     .          +cmplx(xr2(l),xi2(l))*cmplx(0.,al2)
    2 continue
    1 continue
      return
      end
      subroutine cpft(r, i, n, incp, signp)
c  fortran transliteration of singleton's 6600 assembly-coded fft.
c  differs from singleton's original in that there is a special loop
c  for angle=pi/2. this should be faster on machines whose floating
c  point arithmetic is much slower than indexing (not true on cdc 6600).
c  see comments in other version.
c  a. bruce langdon, m division, l.l.l., 1971.
      real r(1), i(1)
      integer signp, span, rc
      real sines(15), i0, i1
      data sines(1)/0./
      if( sines(1).eq.1. ) go to 1
      sines(1)=1.
      t=atan(1.)
      do 2 is=2,15
      sines(is)=sin(t)
    2 t=t/2.
    1 continue
      if( n.eq.1 ) return
      inc=incp
      sgn=isign(1,signp)
      ninc=n*inc
      span=ninc
      it=n/2
      do 3 is=1,15
      if( it.eq.1 ) go to 12
    3 it=it/2

10    t=s+(s0*c-c0*s)
      c=c-(c0*c+s0*s)
      s=t
11    k1=k0+span
      r0=r(k0+1)
      r1=r(k1+1)
      i0=i(k0+1)
      i1=i(k1+1)
      r(k0+1)=r0+r1
      i(k0+1)=i0+i1
      r0=r0-r1
      i0=i0-i1
      r(k1+1)=c*r0-s*i0
      i(k1+1)=s*r0+c*i0
      k0=k1+span
      if( k0.lt.ninc ) go to 11
      k1=k0-ninc
      c=-c
      k0=span-k1
      if( k1.lt.k0 ) go to 11
      k0=k0+inc
      k1=span-k0
      if( k0.lt.k1 ) go to 10
12    continue
      span=span/2
      k0=0
13    k1=k0+span
      r0=r(k0+1)
      r1=r(k1+1)
      i0=i(k0+1)
      i1=i(k1+1)
      r(k0+1)=r0+r1
      i(k0+1)=i0+i1
      r(k1+1)=r0-r1
      i(k1+1)=i0-i1
      k0=k1+span
      if( k0.lt.ninc ) go to 13
      if( span.eq.inc ) go to 20
      k0=span/2
14    k1=k0+span
      r0=r(k0+1)
      r1=r(k1+1)
      i0=i(k0+1)
      i1=i(k1+1)
      r(k0+1)=r0+r1
      i(k0+1)=i0+i1
      r(k1+1)=sgn*(i1-i0)
      i(k1+1)=sgn*(r0-r1)
      k0=k1+span
      if( k0.lt.ninc ) go to 14
      k1=inc+inc
      if( span.eq.k1 ) go to 12
      c0=2.*sines(is)**2
      is=is-1
      s=sign( sines(is),sgn )
      s0=s
      c=1.-c0
      k0=inc
      go to 11

20    n1=ninc-inc
      n2=ninc/2
      rc=0
      ji=rc
      ij=ji
      if( n2.eq.inc ) return
      go to 22
21    ij=n1-ij
      ji=n1-ji
      t=r(ij+1)
      r(ij+1)=r(ji+1)
      r(ji+1)=t
      t=i(ij+1)
      i(ij+1)=i(ji+1)
      i(ji+1)=t
      if( ij.gt.n2 ) go to 21
22    ij=ij+inc
      ji=ji+n2
      t=r(ij+1)
      r(ij+1)=r(ji+1)
      r(ji+1)=t
      t=i(ij+1)
      i(ij+1)=i(ji+1)
      i(ji+1)=t
      it=n2
23    it=it/2
      rc=rc-it
      if( rc.ge.0 ) go to 23
      rc=rc+2*it
      ji=rc
      ij=ij+inc
      if( ij.lt.ji ) go to 21
      if( ij.lt.n2 ) go to 22

      return
      end
c
c........................................................
      subroutine matsub ( a, b, n, c )
c........................................................
c
      complex a, b, c
      dimension a(n,n), b(n,n), c(n,n)
c
      do 100 i = 1, n
      do 100 j = 1, n
c
      c(j,i)=a(j,i)-b(j,i)
c
  100 continue
      return
      end
c........................................................
      subroutine matadd ( a, b, n )
c........................................................
c
      complex a, b
      dimension a(n,n), b(n,n)
c
      do 100 i = 1, n
      do 100 j = 1, n
c
      a(j,i)=a(j,i)+b(j,i)
c
  100 continue
      return
      end
c........................................................
      subroutine matmul ( a, b, n, c )
c............................................................
      complex a, b, c, sum
      dimension a(n,n), b(n,n), c(n,n)
c
      do 100 i = 1, n
      do 100 j = 1, n
c
      sum = cmplx(0.0,0.0)
      do 60 k = 1, n
      sum = sum + a(i,k)*b(k,j)
   60 continue
c
      c(i,j) = sum
c
  100 continue
c
      return
      end
      subroutine matm2(a,n)
      complex a(n,n)
      do 1 j=1,n
      do 1 k=1,n
      a(k,j)=cmplx(0.,-1.0)*a(k,j)
    1 continue
      return
      end
      subroutine matm1(a,n)
      complex a(n,n)
      do 1 j=1,n
      do 1 k=1,n
      a(k,j)=cmplx(0.,1.0)*a(k,j)
    1 continue
      return
      end
      subroutine rmatrix(rm,cm,n)
      complex cm(n,n)
      dimension rm(n,n)
      do 1 i=1,n
      do 1 j=1,n
      rm(j,i)=real(cm(j,i))
    1 continue
      return
      end
c
c
c.......................................................................
c
      subroutine deltak(engk)
      include 'clich1'
      include 'clich2'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
c.......................................................................
cc  compute inertial kinetic energy
      engk=0.0
      nphi=50
      dphi=twopi/nphi
      do 10 j=1,nphi
      phi=(j-1.)*dphi
      do 11 i=1,mth
      th=(i-1.)*dth
      qndel=ntor*(phi-q*delta(i))
      xir=0.0
      xis=0.0
      xir1=0.0
      xis1=0.0
      do 12 l=1,mt
      thl=(l-1+minm)*th
      phase=thl-qndel
      xir=xir+eigfun(nsrf,l,1)*cos(phase)
      xis=xis+eigfun(nsrf,l,3)*sin(phase)
      if(ishft.eq.0) then
         xir1=xir1+eigfun_s(nsrf,l,1)*sin(phase)
         xis1=xis1+eigfun_s(nsrf,l,3)*cos(phase)
      endif
   12 continue
      engk=engk+xjacob(i)*((xir**2+xir1**2)/grpssq(i)
     &     +(xis**2+xis1**2)*grpssq(i)/bsq(i))
   11 continue
   10 continue
      engk=engk*dth*dphi
c      write(*,*)nsrf,engk,rgrid(nsrf),grpssq(1),eigfun(nsrf,8,1)
c     &     ,eigfun(nsrf,8,2),eigfun(nsrf,8,3)
c
      return
      end
c
c
cxxxxxxxxxxxx
      subroutine kinetic
c
      include 'clich1'
      include 'clich2'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      common/rbccol/abe(nn,lam),abp(nn,lam),abe_av(nn,lam),abe_avav(lam
     &     ,nn)
      dimension galx(lam),gale(lam),fn(lam),gn(lam),hn(lam,3)
      dimension tkbouc(lam),tpnt1c(lam),tpnt2c(lam)
c
cc..  mth is # of grid alomth b-field line
ccc..   lam is total # of grid in pitch angle space
cccc  define max. & min. of h=(b0/b)
cccc..  define quantities for bounce averages
ccc..  equilibrium quantities are defined on theta-grid which is clock-wise
ccc..  phi (toroidal angle) is in the counter clockwise direction from the top of tokamak
ccc.. so that grad(psi) x grad(theta) . grad(phi) is positive
c
c
cccc set phot=photp=0.0
      photp=0.0
cccccc
      
      do 1 i=1,mth2
      ss3(i)=xjacob(i)*bf(i)
ccc..  magnetic drift freq.
      dp(i)=curvn(i)-curvs(i)*(qdelp(i)+gptgpp(i)*(q+qdelt(i)))
      dd(i)=-(0.5*dp(i)+pp/bsq(i))*bf(i)
ccc..   the following terms involve exp(i*m*th)
ccc.. bounce averaged terms
      sp2(i)=-(photp/bsq(i)+0.5*curvn(i))*bf(i)
      sp3(i)=-0.5*curvs(i)*bf(i)
c
      st1(i)=(.5*curvn(i)+gptdtgp(i)-gptgpp(i)*gpdtgp(i))*bsq(i)
      st2(i)=ntor*(qdelp(i)+gptgpp(i)*(q+qdelt(i)))*bsq(i)
      st2i(i)=gptgpp(i)*bsq(i)
      st3(i)=-ntor*xjacob(i)*bsq(i)
      gt0(i)=rg*bsdtbs(i)-0.5*xjacob(i)*curvs(i)*bsq(i)
      gp0(i)=0.5*bsdtbs(i)/xjacob(i)
      gp1(i)=(xjdtxj(i)+0.5*bsdtbs(i))
    1 continue
c
ccc  assume there are only one local maxima and minimum in b field
      bmin=1.e10
      bmax=0.0
      do 11 i=1,mth
      if(bmin.gt.bf(i)) bmin=bf(i)
      if(bmax.lt.bf(i)) bmax=bf(i)
   11 continue
      hmin(nsrf)=1.0/bmax
      hmax(nsrf)=1.0/bmin
ccc
ccc   lamc is the number of pitch angle grid points for circulating
ccc   particle
ccc   by assuming the integrand is singular at al=hmin, we multiply sqrt(al-hmin)
ccc   and devide the same factor to the integrand and performe the
ccc   pitch angle integration by formula (25.4.37) in abramowitz et al
ccc   xalic(i), i=1,lamc is i'th zero of lengendre polynomial of order 2*lamc
ccc   walic(i) is the weighting factor.  salmc(k) will be multiplied to the
ccc   integrand before the pitch angle summation.
c
      aac=0.0
      bbc=hmin(nsrf)
      smc=bbc-aac
      salm0c=2.0*sqrt(smc)
      do 22 k=1,lamc
c     
c     galc(1) => hmin, galc(lamc) => 0.
c     
         galc(k)=smc*(1.-xalic(k)**2)
         salmc(k)=salm0c*walic(k)*sqrt(bbc-galc(k))
 22   continue
      tpnt1=-pi
      tpnt2=pi
      do 23 k=1,lamc
         call transit(galc(k),k,tpnt1,tpnt2)
 23   continue
c     
ccc..  lam is the number of trapped particle pitch angle grid points
ccc   by assuming the integrand is singular at al=hmin, we multiply sqrt(al-hmin)
ccc   and devide the same factor to the integrand and performe the
ccc   pitch angle integration by formula (25.4.37) in abramowitz et al
ccc   xali(i), i=1,lamc is i'th zero of lengendre polynomial of order 2*lamc
ccc   wali(i) is the weighting factor.  salm(k) will be multiplied to the
ccc   integrand before the pitch angle summation.
      aa=-hmax(nsrf)
      bb=-hmin(nsrf)
      sm=bb-aa
      salm0=2.*sqrt(sm)
      hlam=1.0/float(lam)
      do 12 k=1,lam
c
c gal(1) => hmin, gal(lam) => hmax
c
c next one is for coll el calculsns and introduces quadratic grid
c gale(1) => hmax, gale(lam) => hmin so that galx(1) => 0, galx(lam) => 1
c
      rlam=hlam*float(k-1)+0.5*hlam
      galx(k)=rlam-0.5*hlam
      gale(k)=-bb+sm*(1.0-rlam**2)
      gal(k)=-(aa+sm*(1.-xali(k)**2))
      salm(k)=salm0*wali(k)*sqrt(gal(k)-hmin(nsrf))
   12 continue
c      galx(lam+1)=1.0
c      write(*,*)'k,galx(k),gale(k),galx(k)'
      csss=abs(psitot)
c      write(*,*)'So far 1',nsrf
      do 13 k=lam,1,-1
ccc.. for trapped electron collisional damping calculation
         call tps(gale(k),tpnt1c(k),tpnt2c(k))
         kp=min(k+1,lam)
         call bouncee(csss,gale(k),galx(k),galx(kp),k,tpnt1c(k),tpnt2c(k
     &        ),tkbouc)
c         write(*,'(i2,5f10.6)') k,gal(k),gale(k),galx(k),tpnt1,tpnt2
ccc.. calculate turning point for each pitch angle
         call tps(gal(k),tpnt1,tpnt2)
ccc.. bounce averaged quantities
         call bounce(csss,gal(k),k,tpnt1,tpnt2)
   13 continue
c This computes the components of the parallel electric field
c this IF is to avoid some FOP you may fix that.
      if(ishft.ne.0.and.nsrf.gt.3)then
         call E_parallel(abe_avav,abe_av,abe,abp,abpe,csss,galx,gale,jb
     &     ,lam,mt,nn,nsrf,salmc,sm,smc,tkbouc,tkcir,tpnt1c,tpnt2c)
      endif
c      stop
c
c add  bounce averaged for collisonal electron terms
c
c      call tps(gal(1),tpnt1,tpnt2)
c      write(16,*)'gal(1),tpnt1,tpnt2'
c      write(16,*)gal(1),tpnt1,tpnt2
c      write(6,*)"input sl1,sl2,ns,idet"
c      read(5,*)sl1,sl2,nsx,idet
      sl1=0.1
      sl2=0.5
      nsx=1
      do 40 k=1,lam-1
	fn(k)=2.0*tdbou(k)/((galx(k)+galx(k+1))*sm)
	gn(k)=tkbouc(k)*(galx(k)+galx(k+1))*0.5
	hn(k,1)=abp(nsrf,k)
	hn(k,2)=abe(nsrf,k)
	hn(k,3)=abe_av(nsrf,k)
 40   continue
      k=lam
      fn(k)=2.0*tdbou(k)/((galx(k)+1.)*sm)
      gn(k)=tkbouc(k)*(galx(k)+1.)*0.5
      hn(k,1)=abp(nsrf,k)
      hn(k,2)=abe(nsrf,k)
      hn(k,3)=abe_av(nsrf,k)
ckg   Also above sm: \delta \Lambda = 1/b_min - 1/bmax
ckg    gn: t_bounce * \eta
ckg    fn: t_bounce * < \chi^2 > / \eta / (\delta\Lambda)
c?????????????? missed \Lambda in fn definition
c
c      write(16,*)"i,tdbou,abe"
c      write(16,*)(k,tdbou(k),hn(k,1),k=1,lam)
c
c solve the bounce-averaged collison equation for trapped
c electrons, the number of base function is idet.
c
      call solu(nsrf,lam,galx(1),fn(1),gn(1),hn(1,1),sl1,sl2,nsx,idet
     &     ,xmu)
ccc
c
ccc.. store matrices due to trapped particles
c      write(iodat,98) nsrf, lam
c      write(iodat,99) (gal(k),k=1,lam)
c      write(iodat,99) (wd(k),k=1,lam)
   98 format('nsrf=',i3,'lamda grid pts =',i3,'lamda, wd, wbounce are')
   99 format(10e12.3)
      return
      end
cxxxxxxxxxxxx
cxxxxxxxxxxxx
c
c     
      subroutine solu(nsrf,n,x,fn,gn,hn,sl1,sl2,ns,idet,xmu)
      parameter(nlam=48,ndet=12,nnsrf=501)
      dimension z(nlam),dx(nlam)
      dimension fn(n),gn(n),hn(n,3),x(n)
      dimension eig(ndet,nlam),a(nlam,nlam),b(nlam,nlam)
      common/coleig/xlam(ndet,nnsrf),cx(nnsrf,ndet,3)
      if(n.ne.nlam.or.ndet.lt.idet.or.nnsrf.lt.nsrf) stop
     &     'Check dimens in SOLU'
      h=1.0/float(n)
      c=-h**2
      do 199 i=1,n-1
         dx(i)=(x(i+1)-x(i))/h
 199  continue
      a=0
      b=0
      dx(n)=(1.-x(n))/h
      a(1,1)=-fn(1)
      a(1,2)=fn(1)
      b(1,1)=c*gn(1)*dx(1)**2
      b(1,2)=c*gn(1)*dx(1)**2
      do 200 i=2,n-1
         a(i,i-1)=fn(i-1)/dx(i-1)
         a(i,i)=-(fn(i-1)/dx(i-1)+fn(i)/dx(i))
         a(i,i+1)=fn(i)/dx(i)
         b(i,i-1)=c*gn(i-1)*dx(i-1)
         b(i,i)=c*(gn(i-1)*dx(i-1)+gn(i)*dx(i))
         b(i,i+1)=c*gn(i)*dx(i)
 200  continue
      a(n,n-1)=fn(n-1)/dx(n-1)
      a(n,n)=-(fn(n-1)/dx(n-1)+fn(n)/dx(n))
      b(n,n-1)=c*gn(n-1)*dx(n-1)
      b(n,n)=c*(gn(n-1)*dx(n-1)+gn(n-1)*dx(n-1))
c     
      eps=1.0e-4
c     slam=0.9
      hs=(sl2-sl1)/float(ns)
c     idet=0
      if(idet.eq.0)then
      do 600 i=1,ns
      slam=sl1+hs*float(i-1)
      call bansol(a(1,1),b(1,1),n,slam,eps,z,det,idet)
c      write(6,*)"eigenvalue,det"
c      write(6,*)slam,det
 600  continue
      stop "subr. solu stop"
      end if
      slx=sl1
      call bansol(a(1,1),b(1,1),n,sl1,eps,z,det,0)
      detold=det
      iroot=0
 620  continue
 610  slx=slx+hs
      slxx=slx
      call bansol(a(1,1),b(1,1),n,slxx,eps,z,det,0)
      detnew=det
      if(detold*detnew.lt.0.0)then
         call bansol(a(1,1),b(1,1),n,slxx,eps,z,det,1)
         iroot=iroot+1
         xlam(iroot,nsrf)=slxx
         do 700 i=1,n
            eig(iroot,i)=z(i)
 700     continue
c         write(6,*)"iroot,slxx"
c         write(6,*)iroot,slxx
      else
         go to 610
      end if
c
c      iroot=iroot+1
      detold=detnew
      if(iroot.lt.idet)go to 620
ckg By now eigenvalues and eigenfrequencies are defined
      if(nsrf.eq.1000)then
      write(*,*) 'nsrf,idet',nsrf,idet
      write(*,*) 'eigenvalues',(xlam(i,nsrf),i=1,idet)
ckg this small part to check the orthogonality. It seems working fine up to 2% error
c      do ii=1,idet
c      do i=1,n
c         do ir=1,idet
c            cx(i,ir,1)=eig(ir,i)
c            a(i,ir)=eig(ii,i)*eig(ir,i)*gn(i)
c            b(i,ir)=eig(ir,i)*eig(ir,i)*gn(i)
c         enddo
c      enddo
c      call ASIMP(x(1),b(1,ii),n,yn1)
c      do ir=1,idet
c         call ASIMP(x(1),b(1,ir),n,yn2)
c         call ASIMP(x(1),a(1,ir),n,yn)
c         yn=yn/sqrt(yn1*yn2)
c         write(*,*) ii,ir,yn
c      enddo
c      enddo
c...n=lam..................................................
      do i=1,n
         do ir=1,idet
            cx(i,ir,1)=eig(ir,i)
         enddo
      enddo
      do i=1,n
         i0=max(1,i-1)
         i0=min(i0,n-3)
         do ir=1,idet
            a(i,ir)=eig(ir,i)*xlam(ir,nsrf)*gn(i)
            call funder4(x(i0),cx(i0,ir,1),x(i),b(i,ir),yn)
            b(i,ir)=b(i,ir)*fn(i)
            call ASIMP(x(1),a(1,ir),i,yn)
            b(i,ir)=b(i,ir)+4.*yn
         enddo
      enddo
      call twodgraf(x(1),a(1,1),nlam,4,'white','solid',
     $     'red',5,'ts','lambda','Eigenfs','ts',6,6,1,a,0)
      call twodgraf(x(1),b(1,1),nlam,4,'white','solid',
     $     'red',5,'ts','lambda','Eigenfs','ts',6,6,1,a,0)
      call twodgraf(x(1),gn(1) ,nlam,1,'white','solid',
     $     'red',5,'ts','lambda','Eigenfs','ts',6,6,1,a,0)
      call frame(0)
      stop 'You just hit the test loop'
      endif
c So far the check gives that one needs 4 times bigger eigenvalues to
c  satisfy the eigenmode equation
ckg ---------------------------------------------------------
      do 800 i=1,idet
         cx(nsrf,i,1)=0.
         cx(nsrf,i,2)=0.
         cx(nsrf,i,3)=0.
         c2=0.0
         do 801 j=1,n
            c1=gn(j)*eig(i,j)
            cx(nsrf,i,1)=cx(nsrf,i,1)+c1*hn(j,1)
            cx(nsrf,i,2)=cx(nsrf,i,2)+c1*hn(j,2)
            cx(nsrf,i,3)=cx(nsrf,i,3)+c1*hn(j,3)
            c2=c2+c1*eig(i,j)
 801     continue
         c2=sqrt(h/c2)
         cx(nsrf,i,1)=cx(nsrf,i,1)*c2
         cx(nsrf,i,2)=cx(nsrf,i,2)*c2
         cx(nsrf,i,3)=cx(nsrf,i,3)*c2
 800  continue
      if(nsrf.eq.46)then
c hn(i,1) could be high at the suspicious point
      write(16,*)"i,gn(i),hn(i,1),hn(i,2),hn(i,3)"
      write(16,*)(i,gn(i),hn(i,1),hn(i,2),hn(i,3),i=1,n)
      end if
      if((nsrf/10)*10.eq.nsrf)then
      write(16,*)"i,xlam(i,nsrf),cx(nsrf,i,1)"
      write(16,*)(i,xlam(i,nsrf),cx(nsrf,i,1),i=1,idet)
      write(16,*)"i,eig(idet,i)"
      write(16,*)(i,eig(idet,i),i=1,n)
      end if
      return
      end
c
c
c
c
      subroutine bansol(a1,b1,n,slam,eps,z2,det,idet)
      parameter(nlam=48)
      dimension a1(n,n),b1(n,n),a(nlam,nlam),b(nlam,nlam),z1(nlam),z2(n)
      if(n.ne.nlam) stop 'Check dimens in BANSOL'
      do 100 i=1,n
      do 100 j=1,n
         a(i,j)=a1(i,j)-slam*b1(i,j)
         b(i,j)=b1(i,j)
 100  continue
c      write(6,*)"a(1,1),a(1,2),a(2,1)"
c      write(6,*)a(1,1),a(1,2),a(2,1)
c      write(6,*)"b(1,1),b(1,2),a(n,n)"
c      write(6,*)b(1,1),b(1,2),a(n,n)
      niter=0
      do 200 i=1,n-1
         c=a(i+1,i)/a(i,i)
         do 200 j=1,i+1
            a(i+1,j)=a(i+1,j)-a(i,j)*c
            b(i+1,j)=b(i+1,j)-b(i,j)*c
 200     continue
         do 300 i=1,n-1
            c=a(n-i,n+1-i)/a(n+1-i,n+1-i)
            do 300 j=1,n
               a(n-i,j)=a(n-i,j)-a(n+1-i,j)*c
               b(n-i,j)=b(n-i,j)-b(n+1-i,j)*c
 300        continue
            det=1.0
            do 301 i=1,n
                det=det*a(i,i)
 301        continue
            if(idet.eq.0)return
c            write(6,*)"i,a(i,i)"
c            write(6,*)(i,a(i,i),i=1,n)
            do 400 i=1,n
               z1(i)=1.0/sqrt(float(n))
 400        continue
c     
c     start to iterate
c     
 20         continue
            niter=niter+1
            if(niter.gt.10080)then
            write(*,*)"number of iteration"
            write(*,*)niter
            stop
            end if
            anor1=0.0
            anor2=0.0
            anor=0.0
            do 500 i=1,n
            z2(i)=0.0
               do 501 j=1,n
                  z2(i)=z2(i)+b(i,j)*z1(j)
 501           continue
               z2(i)=z2(i)/a(i,i)
               anor=anor+z2(i)**2
 500        continue
c            write(6,*)"anor"
c            write(6,*)anor
            do 600 i=1,n
               do 600 j=1,n
                  anor1=anor1+z2(i)*a(i,j)*z2(j)
                  anor2=anor2+z2(i)*b(i,j)*z2(j)
 600           continue
c               write(6,*)"anor1,anor2"
c               write(6,*)anor1,anor2
               slam1=1.0/sqrt(anor)
               slam2=anor1/anor2
c               write(6,*)"niter slam1,slam2"
c               write(6,*)niter,slam1,slam2
               crit=abs((slam1-slam2)/slam2)
               crit1=abs((slam1+slam2)/slam2)
               if(crit.lt.eps.or.crit1.lt.eps)then
                  slam=slam+slam2
                  return
               end if
               do 700 i=1,n
                  z1(i)=z2(i)/sqrt(anor)
 700           continue
            go to 20
            return
            end
c
c
      subroutine tps(al,tpnt1,tpnt2)
c
      include 'clich1'
      include 'clich2'
c
      dimension st(2), iis(2)
      common/ppar/par(4)
      external stf
ccccc...   calculate trapped particle turning pts for a given pitch angle al
ccc..   for mirror st can't be in end loss region
cccc..  b-field is symmetric about mid-plane; ns=odd
      do 1 is=1,mth2
      ss1(is)=1.-al*bf(is)
    1 continue
ccc.. assume at least one pair of turning points for each pitch angle
      nst=0
      mthd2p1=mth2/2
      do 2 is=mthd2p1,mth
      if((ss1(is+1)*ss1(is)).gt.0.) go to 2
      nst=nst+1
      iis(nst)=is
    2 continue
c      write(*,*) nst,(bf(is),is=1,mth2)
ccccc...  4-point lagrange interpolation
      do 3 ist=1,nst
      i=iis(ist)
      par(1)=ss1(i-1)
      par(2)=ss1(i)
      par(3)=ss1(i+1)
      par(4)=ss1(i+2)
      x=par(2)/(par(2)-par(3))
      if(x.eq.1.0.and.abs(par(2)).lt.abs(par(4))) x=0.3
      if(x.eq.0.0.and.abs(par(3)).lt.abs(par(1))) x=0.7
cc..  stf=ss1 defined by 4-point Lagrangian near the turning points
cccc..  zreal is a root finder of stf(x)=0 from the 'imsl' package
      eps1=1.e-9
      eps2=1.e-9
      eps=1.e-2
      eta=1.e-2
      nroot=1
      itmax=200
      xguess=x
      call c05ade(0.,1.,eps1,eps2,stf,x,ier)
c      call zreal1(stf,eps1,eps2,eps,eta,nroot,itmax,xguess,x,ier)
      st(ist)=twopi-(i-1.+x)*dth
    3 continue
cccc  the turning points are tpnt1 & tpnt2
      tpnt1=-st(1)
      tpnt2=st(1)
      if(nst.eq.1) return
      tpnt1=st(1)
      tpnt2=st(2)
      return
      end
cxxxxxxxxxxxx
cxxxxxxxxxxx
      subroutine bounce(csss,al,k,tpnt1,tpnt2)
c
      include 'clich1'
      include 'clich2'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      common/coefs/c1(nts),c2(nts),c3(nts),c4(nts),igm1(nts)
     &,ig(nts),igp1(nts),igp2(nts)
      common/cossin/cost1(nts),sint1(nts),cost2(nts),costb(nts)
     &,sintb(nts),cost(nts),sint(nts),cost3(nts),sint3(nts)
     &,cost1p(2*nts),sint1p(2*nts),cost2p(2*nts)
     &,cost1m(2*nts),sint1m(2*nts),cost2m(2*nts)
     &,costp(nts),costm(nts),sintp(nts),sintm(nts)
      dimension coste(nts,3)
c
cc..  xx2 is the radial scale factor
      real nq
      nq=ntor*q
      xx2=psitot*2.*rrg
c
      bf0=0.
      wd(k)=0.
      tpp=tpnt2+tpnt1
      tpm=tpnt2-tpnt1
      wj=pi/nthetaj
      wj0=2.0/twopi
ccc.. bounce average of magnetic drift freq. wd(k)
ccc   integrate from 0 to tpnt2 due to up-down symmetry over nthetaj points
cccc  add bounce resonant terms by fu july, 88
cccc  add tboun(theta), the bounce orbit time for trapped particle
      do 22 i=1,nthetaj
      dt(i)=0.5*(tpp+tpm*sin(rx(i)))
      dtb(i)=dt(i)
c     dt(i)=0.5*(tpp+tpm*sin((i-0.5)*wj*0.5))
c     dtb(i)=0.5*(tpp+tpm*sin((i-0.5)*wj*0.5))
   22 continue
c
ccc  set up coefficients for 4-point Lagrange interpolation
      do 11 i=1,nthetaj
      call fx4p(dt(i),i,mth2,dth)
   11 continue
c
      do 23 i=1,nthetaj
c
c      call fx4(ss2(i),ss1,mth2)
c      call fx4(at(i),ss3,mth2)
c      call fx4(fwb,dd,mth2)
c      call fx4(fwk,dp,mth2)
c     call fx4(sd(i),delta,mth2)
c     call fx4(gd(i),curvn,mth2)
c     call fx4(gs(i),curvs,mth2)
c     call fx4(ap(i),sp2,mth2)
c     call fx4(as(i),sp3,mth2)
c     call fx4(ds1(i),bf,mth2)
c     call fx4(ds2(i),bsq,mth2)
c     call fx4(ad1(i),st1,mth2)
c     call fx4(ad2(i),st2,mth2)
c     call fx4(ad3(i),st2i,mth2)
c     call fx4(wi1(i),st3,mth2)
c     call fx4(wi2(i),gt0,mth2)
c     call fx4(wi(i),gp0,mth2)
c     call fx4(wii(i),gp1,mth2)
      at(i)=c1(i)*ss3(igm1(i))+c2(i)*ss3(ig(i))
     &     +c3(i)*ss3(igp1(i))+c4(i)*ss3(igp2(i))
      ss2(i)=c1(i)*ss1(igm1(i))+c2(i)*ss1(ig(i))
     &     +c3(i)*ss1(igp1(i))+c4(i)*ss1(igp2(i))
      fwb(i)=c1(i)*dd(igm1(i))+c2(i)*dd(ig(i))
     &     +c3(i)*dd(igp1(i))+c4(i)*dd(igp2(i))
      fwk(i)=c1(i)*dp(igm1(i))+c2(i)*dp(ig(i))
     &     +c3(i)*dp(igp1(i))+c4(i)*dp(igp2(i))
      sd(i)=c1(i)*delta(igm1(i))+c2(i)*delta(ig(i))
     &     +c3(i)*delta(igp1(i))+c4(i)*delta(igp2(i))
      gd(i)=c1(i)*curvn(igm1(i))+c2(i)*curvn(ig(i))
     &     +c3(i)*curvn(igp1(i))+c4(i)*curvn(igp2(i))
      gs(i)=c1(i)*curvs(igm1(i))+c2(i)*curvs(ig(i))
     &     +c3(i)*curvs(igp1(i))+c4(i)*curvs(igp2(i))
      ap(i)=c1(i)*sp2(igm1(i))+c2(i)*sp2(ig(i))
     &     +c3(i)*sp2(igp1(i))+c4(i)*sp2(igp2(i))
      as(i)=c1(i)*sp3(igm1(i))+c2(i)*sp3(ig(i))
     &     +c3(i)*sp3(igp1(i))+c4(i)*sp3(igp2(i))
      ds1(i)=c1(i)*bf(igm1(i))+c2(i)*bf(ig(i))
     &     +c3(i)*bf(igp1(i))+c4(i)*bf(igp2(i))
      ds2(i)=c1(i)*bsq(igm1(i))+c2(i)*bsq(ig(i))
     &     +c3(i)*bsq(igp1(i))+c4(i)*bsq(igp2(i))
      ad1(i)=c1(i)*st1(igm1(i))+c2(i)*st1(ig(i))
     &     +c3(i)*st1(igp1(i))+c4(i)*st1(igp2(i))
      ad2(i)=c1(i)*st2(igm1(i))+c2(i)*st2(ig(i))
     &     +c3(i)*st2(igp1(i))+c4(i)*st2(igp2(i))
      ad3(i)=c1(i)*st2i(igm1(i))+c2(i)*st2i(ig(i))
     &     +c3(i)*st2i(igp1(i))+c4(i)*st2i(igp2(i))
      wi1(i)=c1(i)*st3(igm1(i))+c2(i)*st3(ig(i))
     &     +c3(i)*st3(igp1(i))+c4(i)*st3(igp2(i))
      wi2(i)=c1(i)*gt0(igm1(i))+c2(i)*gt0(ig(i))
     &     +c3(i)*gt0(igp1(i))+c4(i)*gt0(igp2(i))
      wi(i)=c1(i)*gp0(igm1(i))+c2(i)*gp0(ig(i))
     &     +c3(i)*gp0(igp1(i))+c4(i)*gp0(igp2(i))
      wii(i)=c1(i)*gp1(igm1(i))+c2(i)*gp1(ig(i))
     &     +c3(i)*gp1(igp1(i))+c4(i)*gp1(igp2(i))
c
   23 continue
c      write(*,*) nthetaj,tpm,(ss2(i),i=1,10)
      do 24 i=1,nthetaj
      if(ss2(i).lt.0.) ss2(i)=-ss2(i)
   24 continue
      do 25 i=1,nthetaj
      gt1(i)=sqrt((dt(i)-tpnt1)*(tpnt2-dt(i))/ss2(i))
      gt1i(i)=gt1(i)*at(i)
      wd(k)=wd(k)+gt1i(i)*(fwk(i)+fwb(i)*al)*drx(i)
   25 continue
      do 26 i=1,nthetaj
      bf0=bf0+gt1i(i)*drx(i)
   26 continue
      tboun(1)=0.5*gt1i(1)/bf0*drx(1)
c      deltb(1)=tboun(1)-dtb(1)*2.0/tpm
      do 16 i=2,nthetaj
      tboun(i)=(gt1i(i-1)+gt1i(i))/bf0*0.25*(drx(i)+drx(i-1))+tboun(i-1)
c      deltb(i)=tboun(i)-dtb(i)*2.0/tpm
   16 continue
      tkbou(k)=bf0*2.0
      wjb=2.0/tkbou(k)
      wd(k)=wd(k)*wjb
      wjb1=wjb
c
ccc..  bounce average of quantities involving exp(i*m*th)
c
      do 32 j=1,np
ccc..  j=1 is for p=0 bounce term
      jm1=j-1
      do 32 l=1,mt
      l2=l+minm-1
      alq=l2-ntor*q
c
      rbce(k,l,j,1)=0.0
      rbce(k,l,j,2)=0.0
      rbce(k,l,j,3)=0.0
c
      do 33 i=1,nthetaj
      cost(i)=cos(alq*dt(i))
      sint(i)=sin(alq*dt(i))
c
      costb(i)=cos(0.5*jm1*pi*tboun(i))*drx(i)
      sintb(i)=sin(0.5*jm1*pi*tboun(i))*drx(i)
      cost1(i)=cost(i)*cos(jm1*pi/2.0)*costb(i)
     &        +sint(i)*sin(jm1*pi/2.0)*sintb(i)
      sint1(i)=sint(i)*cos(jm1*pi/2.0)*costb(i)
     &        -cost(i)*sin(jm1*pi/2.0)*sintb(i)
      cost2(i)=gt1i(i)/ds1(i)*cost1(i)
      cost1(i)=gt1i(i)*(gd(i)+al*ap(i))*cost1(i)
      sint1(i)=gt1i(i)*(gs(i)+al*as(i))*sint1(i)
   33 continue
      do 34 i=1,nthetaj
      rbce(k,l,j,1)=rbce(k,l,j,1)+cost1(i)
      rbce(k,l,j,3)=rbce(k,l,j,3)+sint1(i)
      rbce(k,l,j,2)=rbce(k,l,j,2)+cost2(i)
   34 continue
      rbce(k,l,j,1)=rbce(k,l,j,1)*wjb1
      rbce(k,l,j,3)=rbce(k,l,j,3)*wjb1
      rbce(k,l,j,2)=rbce(k,l,j,2)*wjb1*al
c      
c            if(l.eq.2.and.j.eq.1.and.abs(q-1.3).lt.0.01)
c     +         write(*,*)'old rbce',
c     +         (rbce(k,l,j,jj),jj=1,1),tpnt1,tpnt2,al
c
   32 continue
c
      if(nsrf.eq.1) go to 98
      isf0=max(2,nsrf-1)
      isf0=min(nosurf-3,isf0)
      if(nsrf.gt.1) then
         sss=0.5/rgrid(nsrf)
      else
         sss=0.0
      end if
      do 40 l=1,mt
        call der4(rgrid(isf0),eigfun(isf0,l,3),rgrid(nsrf),dyn)
        do ijb=1,jb
        abpe(k,l,nsrf,ijb,1)=csss*rbce(k,l,ijb,1)*eigfun(nsrf,l,1)
        abpe(k,l,nsrf,ijb,4)=csss*rbcec(k,l,ijb,1)*eigfun(nsrf,l,1)
ckg this is parallel electric field term, may give large contribution to
ckg the damping if there is singularity
ckg tune     &       +rbcep(k,l,1)*dyn*sss
c
c     >              (eigfun(nsrf,l,3)-eigfun(nsrf-1,l,3))
c     >              (rgrid(nsrf)-rgrid(nsrf-1))
c define v.E terms for electrons
        enddo
        do ijb=1,jb
           habp=0
           habpc=0
           do j=1,3
              habp=habp+rbce(k,l,ijb,j)*eigfun(nsrf,l,j)
              habpc=habpc+rbcec(k,l,ijb,j)*eigfun(nsrf,l,j)
           enddo
           abpe(k,l,nsrf,ijb,2)=habp
           abpe(k,l,nsrf,ijb,5)=habpc
        enddo
 40   continue
   98 continue
c
      return
      end
c
cxxxxxxxxxxx
ccc add transit particle and include the transit resonance
      subroutine transit(al,k,tpnt1,tpnt2)
c
      include 'clich1'
      include 'clich2'
      common/coefs/c1(nts),c2(nts),c3(nts),c4(nts),igm1(nts)
     &,ig(nts),igp1(nts),igp2(nts)
      common/cossin/cost1(nts),sint1(nts),cost2(nts),costb(nts)
     &,sintb(nts),cost(nts),sint(nts),cost3(nts),sint3(nts)
     &,cost1p(2*nts),sint1p(2*nts),cost2p(2*nts)
     &,cost1m(2*nts),sint1m(2*nts),cost2m(2*nts)
     &,costp(nts),costm(nts),sintp(nts),sintm(nts)
      dimension coste(nts,3)
c
cc..  xx2 is the radial scale factor
      real nq
      nq=ntor*q
      xx2=psitot*2.*rrg
c
      bf0=0.
      wdc(k)=0.
      tpp=tpnt2+tpnt1
      tpm=tpnt2-tpnt1
      wj=pi/nthetaj
      wj0=2.0/twopi
ccc.. transit average of magnetic drift freq. wd(k)
ccc   integrate from 0 to tpnt2 due to up-down symmetry over nthetaj points
cccc  add tcirc(theta), the circulating orbit time for passing particle
      do 1 is=1,mth2
      ss1(is)=1.-al*bf(is)
    1 continue
      do 22 i=1,nthetaj
      dt(i)=0.5*(tpp+tpm*sin(rx(i)))
      dtc(i)=dt(i)
c     dt(i)=0.5*(tpp+tpm*sin((i-0.5)*wj*0.5))
c     dtc(i)=0.5*(tpp+tpm*sin((i-0.5)*wj*0.5))
   22 continue
c
ccc  set up coefficients for 4-point Lagrange interpolation

      do 11 i=1,nthetaj
      call fx4p(dt(i),i,mth2,dth)
   11 continue
c
      do 23 i=1,nthetaj
c
c      call fx4(ss2(i),ss1,mth2)
c      call fx4(at(i),ss3,mth2)
c      call fx4(fwb,dd,mth2)
c      call fx4(fwk,dp,mth2)
c     call fx4(sd(i),delta,mth2)
c     call fx4(gd(i),curvn,mth2)
c     call fx4(gs(i),curvs,mth2)
c     call fx4(ap(i),sp2,mth2)
c     call fx4(as(i),sp3,mth2)
c     call fx4(ds1(i),bf,mth2)
c     call fx4(ds2(i),bsq,mth2)
c     call fx4(ad1(i),st1,mth2)
c     call fx4(ad2(i),st2,mth2)
c     call fx4(ad3(i),st2i,mth2)
c     call fx4(wi1(i),st3,mth2)
c     call fx4(wi2(i),gt0,mth2)
c     call fx4(wi(i),gp0,mth2)
c     call fx4(wii(i),gp1,mth2)
      at(i)=c1(i)*ss3(igm1(i))+c2(i)*ss3(ig(i))
     &     +c3(i)*ss3(igp1(i))+c4(i)*ss3(igp2(i))
      ss2(i)=c1(i)*ss1(igm1(i))+c2(i)*ss1(ig(i))
     &     +c3(i)*ss1(igp1(i))+c4(i)*ss1(igp2(i))
      fwb(i)=c1(i)*dd(igm1(i))+c2(i)*dd(ig(i))
     &     +c3(i)*dd(igp1(i))+c4(i)*dd(igp2(i))
      fwk(i)=c1(i)*dp(igm1(i))+c2(i)*dp(ig(i))
     &     +c3(i)*dp(igp1(i))+c4(i)*dp(igp2(i))
      sd(i)=c1(i)*delta(igm1(i))+c2(i)*delta(ig(i))
     &     +c3(i)*delta(igp1(i))+c4(i)*delta(igp2(i))
      gd(i)=c1(i)*curvn(igm1(i))+c2(i)*curvn(ig(i))
     &     +c3(i)*curvn(igp1(i))+c4(i)*curvn(igp2(i))
      gs(i)=c1(i)*curvs(igm1(i))+c2(i)*curvs(ig(i))
     &     +c3(i)*curvs(igp1(i))+c4(i)*curvs(igp2(i))
      ap(i)=c1(i)*sp2(igm1(i))+c2(i)*sp2(ig(i))
     &     +c3(i)*sp2(igp1(i))+c4(i)*sp2(igp2(i))
      as(i)=c1(i)*sp3(igm1(i))+c2(i)*sp3(ig(i))
     &     +c3(i)*sp3(igp1(i))+c4(i)*sp3(igp2(i))
      ds1(i)=c1(i)*bf(igm1(i))+c2(i)*bf(ig(i))
     &     +c3(i)*bf(igp1(i))+c4(i)*bf(igp2(i))
      ds2(i)=c1(i)*bsq(igm1(i))+c2(i)*bsq(ig(i))
     &     +c3(i)*bsq(igp1(i))+c4(i)*bsq(igp2(i))
      ad1(i)=c1(i)*st1(igm1(i))+c2(i)*st1(ig(i))
     &     +c3(i)*st1(igp1(i))+c4(i)*st1(igp2(i))
      ad2(i)=c1(i)*st2(igm1(i))+c2(i)*st2(ig(i))
     &     +c3(i)*st2(igp1(i))+c4(i)*st2(igp2(i))
      ad3(i)=c1(i)*st2i(igm1(i))+c2(i)*st2i(ig(i))
     &     +c3(i)*st2i(igp1(i))+c4(i)*st2i(igp2(i))
      wi1(i)=c1(i)*st3(igm1(i))+c2(i)*st3(ig(i))
     &     +c3(i)*st3(igp1(i))+c4(i)*st3(igp2(i))
      wi2(i)=c1(i)*gt0(igm1(i))+c2(i)*gt0(ig(i))
     &     +c3(i)*gt0(igp1(i))+c4(i)*gt0(igp2(i))
      wi(i)=c1(i)*gp0(igm1(i))+c2(i)*gp0(ig(i))
     &     +c3(i)*gp0(igp1(i))+c4(i)*gp0(igp2(i))
      wii(i)=c1(i)*gp1(igm1(i))+c2(i)*gp1(ig(i))
     &     +c3(i)*gp1(igp1(i))+c4(i)*gp1(igp2(i))
c
   23 continue
      do 24 i=1,nthetaj
      if(ss2(i).lt.0.) ss2(i)=-ss2(i)
   24 continue
      do 25 i=1,nthetaj
      gt1(i)=sqrt((dt(i)-tpnt1)*(tpnt2-dt(i))/ss2(i))
      gt1i(i)=gt1(i)*at(i)
      wdc(k)=wdc(k)+gt1i(i)*(fwk(i)+fwb(i)*al)*drx(i)
   25 continue
      do 26 i=1,nthetaj
      bf0=bf0+gt1i(i)*drx(i)
   26 continue
      tcirc(1)=0.5*gt1i(1)/bf0*drx(1)
      deltc(1)=tcirc(i)-dtc(1)/pi
      do 16 i=2,nthetaj
      tcirc(i)=(gt1i(i-1)+gt1i(i))/bf0*0.5*(rx(i)-rx(i-1))+tcirc(i-1)
      deltc(i)=tcirc(i)-dtc(i)/pi
   16 continue
      tkcir(k)=bf0*2.0
      wjb=2.0/tkcir(k)
      wdc(k)=wdc(k)*wjb
      wjb1=wjb
c
ccc..  transit average of quantities involving exp(i*m*th)
c
c
      do 14 l=1,mt
      do 14 j=1,jcc
      l2=l+minm-1
      alq=l2-ntor*q
      alq1=j-1+jcmin-ntor*q
c
      rbcec(k,l,j,1)=0.0
      rbcec(k,l,j,2)=0.0
      rbcec(k,l,j,3)=0.0
c
      do 33 i=1,nthetaj
      cost1(i)=cos(alq*dt(i)-pi*alq1*tcirc(i))*drx(i)
      sint1(i)=sin(alq*dt(i)-pi*alq1*tcirc(i))*drx(i)
      cost2(i)=gt1i(i)/ds1(i)*cost1(i)
      cost1(i)=gt1i(i)*(gd(i)+al*ap(i))*cost1(i)
      sint1(i)=gt1i(i)*(gs(i)+al*as(i))*sint1(i)
   33 continue
      do 34 i=1,nthetaj
c
      rbcec(k,l,j,1)=rbcec(k,l,j,1)+cost1(i)
      rbcec(k,l,j,3)=rbcec(k,l,j,3)+sint1(i)
      rbcec(k,l,j,2)=rbcec(k,l,j,2)+cost2(i)
c
c
   34 continue
      rbcec(k,l,j,1)=rbcec(k,l,j,1)*wjb1
      rbcec(k,l,j,3)=rbcec(k,l,j,3)*wjb1
      rbcec(k,l,j,2)=rbcec(k,l,j,2)*wjb1*al
c
   14 continue
c
      return
      end
cxxxxxxxxxxx
      subroutine fx4p(x,i,ng,dth)
ccc...  ng=mth+2
cccc....   4 point lagrange interpolation
cccc....   fg is a periodic function with fg(1)=fg(ng-1)
cccc...     dth is the grid size, ng is total # of grids
      include 'clich1'
      common/coefs/c1(nts),c2(nts),c3(nts),c4(nts),igm1(nts)
     &,ig(nts),igp1(nts),igp2(nts)
      ix=x/dth
      dp=(x/dth-ix)
      dpm1=dp-1.
      dpp1=dp+1.
      dpm2=dp-2.
      ig(i)=ix+1
      igm1(i)=ig(i)-1
      if(ig(i).eq.1) igm1(i)=ng-2
      igp1(i)=ig(i)+1
      igp2(i)=ig(i)+2
      c1(i)=-dp*dpm1*dpm2/6.
      c2(i)=dpp1*dpm1*dpm2/2.
      c3(i)=-dp*dpp1*dpm2/2.
      c4(i)=dp*dpp1*dpm1/6.
      return
      end
c
      function stf(x)
      common/ppar/par(4)
ccc..  function stf(x)=vpar**2 is a function of b thru position x
ccc..   it is used in 'zreal' to find root of stf(x)=0
      dp=x
      dpm1=x-1.
      dpp1=x+1.
      dpm2=x-2.
      c1=-dp*dpm1*dpm2/6.
      c2=dpp1*dpm1*dpm2/2.
      c3=-dp*dpp1*dpm2/2.
      c4=dp*dpp1*dpm1/6.
      stf=c1*par(1)+c2*par(2)+c3*par(3)+c4*par(4)
      return
      end
c
      subroutine dskout
      include 'clich1'
      include 'clich2'
      integer*8 ierr,length,nadres
c
ccc.. store matrices of mhd stability eqns.
      lgivup=1
      length=(3*jb*mt+4)*lam+(3*jcc*mt+4)*lamc
      nadres=(nsrf-1)*length+10+1
cben      call zwr(mp2,rbce,length,nadres+ishft,lgivup,ierr)
c
      call putwa('equou1',rbce,nadres,length,ierr)
ccc.....
ccc  define pitch angle-flux surface variables for trapped and circulating particles
c
      do 1 k=1,lam
      gal2d(k,nsrf)=gal(k)
      tkbou2d(k,nsrf)=tkbou(k)
c      if(nsrf.eq.30)then
c        write(*,*)k,gal(k)tkbou(k),wd(k)
c      endif
    1 continue
c      if(nsrf.eq.30)then
c         stop
c      endif
      do 2 k=1,lamc
      galc2d(k,nsrf)=galc(k)
      tkcir2d(k,nsrf)=tkcir(k)
    2 continue
      do 3 it=1,mth
      bf2d(it,nsrf)=bf(it)
      xj2d(it,nsrf)=xjacob(it)
    3 continue
      return
      end
ccccc.........
      subroutine pitdist
c
ccc...  define all pitch angle distributions and the derivative for all species
c
      include 'clich1'
      include 'clich2'
      include 'clich1b'

      common/pitchang/fpht(lam,nn,isph),dfpht(lam,nn,isph)
     &,fphc(lamc,nn,isph),dfphc(lamc,nn,isph)

ckg      parameter(klam=48)
      parameter(klam=lam)
      common/intpit/ xal(klam),wal(klam),y(klam),aly(klam),ri(nts)
c
      do 51 k=1,klam
      read(8,52) xal(k),wal(k)
      y(k)=1.0-xal(k)**2
   51 continue
   52 format(2e25.15)
      rewind 8
ccc.....
cc  fp is a Gaussian  pitch angle distribution
cc  dfp = gal * (d/d_gal) fp
c
      do 21 is=1,iht
      do 10 isrf=1,nosurf
c      call fpprp(hmax(isrf),const,dhlamda(is),galh0(is))
      do 1 k=1,lam
cc trapped particle distribution in pitch angle space
cc for non-singular pitch angle distribution
      call fpdist(gal2d(k,isrf)
     &,fpht(k,isrf,is),dfpht(k,isrf,is),dhlamda(is),galh0(is))
ccc  for uniform pitch angle distribution
c      if(is.eq.3) fpht(k,isrf,is)=1.0
c      if(is.eq.3) dfpht(k,isrf,is)=0.0
    1 continue
      do 2 k=1,lamc
cc circulating particle distribution in pitch angle space
cc for non-singular pitch angle distribution
      call fpdist(galc2d(k,isrf)
     &,fphc(k,isrf,is),dfphc(k,isrf,is),dhlamda(is),galh0(is))
ccc  for uniform pitch angle distribution
c      if(is.eq.3) fphc(k,isrf,is)=1.0
c      if(is.eq.3) dfphc(k,isrf,is)=0.0
    2 continue
   10 continue

      if(ihdist(is).eq.1) then
cc for singular pitch angle  distribution
cc  Here, fp and dfp have to multiply an  delta function in pitch angle
      do 20 isrf=1,nosurf
      do 3 k=1,lam
      fpht(k,isrf,is)=1.0
      dfpht(k,isrf,is)=-1.0
    3 continue
      do 4 k=1,lamc
      fphc(k,isrf,is)=1.0
      dfphc(k,isrf,is)=-1.0
    4 continue
   20 continue
      else
      endif


cc  rkhc is the pitch angle integration and surface average of the pitch angle distribution
cc  surface-averaged alpha pressure ph = ch * th0 * ghprf * rkhc
c
      
      do 155 isrf=1,nosurf
cc for uniform pitch angle distribution, rkhc=1

      rkhc(isrf,is)=0.0

c      call fpprp(hmax(isrf),const,dhlamda(is),galh0(is))
      volrx=volr(isrf)*dth

cc for non-uniform(non-singular) pitch angle distribution
cc integrate over pitch angle and poloidal angle to obtain rkhc(isrf)
      do 14 it=1,mth
      hbf=1./bf2d(it,isrf)

      ri(it)=0.0 
      do 15 k=1,klam
      aly(k)=hbf*y(k)
      call fpdist(aly(k),fp,dfp,dhlamda(is),galh0(is))
      ri(it)=ri(it)+wal(k)*fp
   15 continue
      ri(it)=2.0*ri(it)
      rkhc(isrf,is)=rkhc(isrf,is)+xj2d(it,isrf)*ri(it)
   14 continue
      rkhc(isrf,is)=0.5*rkhc(isrf,is)/volr(isrf)
 
cc for singular pitch angle dist.
      if(ihdist(is).eq.1.and.galh0(is).ge.hmin(isrf)
     &   .and.galh0(is).le.hmax(isrf))then
      call galgrid(galh0(is),iphc(isrf,is),gal2d(1,isrf)
     &            ,tkbou2d(1,isrf),tkb,lam)
      rkhc(isrf,is)=tkb/(2.*volrx)
      else
      end if
      if(ihdist(is).eq.1.and.galh0(is).lt.hmin(isrf)) then
      call galcgrid(galh0(is),iphc(isrf,is),galc2d(1,isrf)
     &             ,tkcir2d(1,isrf),tkb,lamc)
      rkhc(isrf,is)=tkb/(2.*volrx)
      else
      end if
  155 continue
c
c     write(*,*) (rkhc(i,is),i=1,40,10)
c     write(*,*) (gal2d(k,4),k=1,40,10)
c     write(*,*) (galc2d(k,4),k=1,40,10)
c     write(*,*) (fphc(k,4,is),k=1,40,10)
c     write(*,*) (dfphc(k,4,is),k=1,40,10)
c     write(*,*) (fpht(k,4,is),k=1,40,10)
c     write(*,*) (dfpht(k,4,is),k=1,40,10)
   21 continue

      return
      end
c
      subroutine fpprp(hbf,const,dbl,bl0)

      include 'clich1'
      include 'clich2'
      include 'clich1b'
c
c const is the normalization constant of the pitch angle distribution
c
      imax=200
      totp=hbf/(imax-1.)
      ft=0.0
      do 1 i=1,imax
      al=(i-1.)*totp
      fpp=exp(-((al-bl0)/dbl)**2)
      ft=ft+fpp
    1 continue
      const=1./(ft*totp)
c     const=100.
c
      return
      end
c
      subroutine fpdist(al,fpp,dfpp,dbl,bl0)
      include 'clich1'
      include 'clich2'
      include 'clich1b'
c
c  define analytical pitch angle distributions fpp 
c  and its dfpp = al*d(fpp)/d(al)
c  al = (magnetic moment * b0)/energy
c
      fpp=exp(-((al-bl0)/dbl)**2)
      dfpp=-2.0*al*(al-bl0)/dbl**2*fpp
c      fpp=1.0
c      dfpp=0.0
c      fpm=exp(-((al-bl0)/dbl)**2)
c      dfpm=-2.0*al*(al-bl0)/dbl**2*fpm
c      fpm=1.0
c      dfpm=0.0
c
      return
      end
c
ccc
ccccc.........
      subroutine inparam
c
cccc.....   define all radial profiles for JT-60U ICRF xp's
c
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      common/plasmap/rmaj,amin,b0

c
ccc..  All quantities are in unit of cgs unit or the equilibrium scales.
ccc..  magnetic field in Gauss, temperature in eV, density in cm**-3, distance in cm
c
ccc   set  iswitch=1 then is=3 for alpha particle
ccc   if  iswitch.ne.1 then is=3 for ICRF tail ion with Maxwellian energy dist.
      iswitch=1
ccc   set  itransp=1 for reading TRANSP generated profiles
ccc   set  itransp=0 for analytical profiles
ckg   set  itransp=2 for analytical profiles with T ~ p, and prescribed 
ckg        density (used for NSTX q_0 = 2.8 case)
ckg NOTE: this will be overwriten in inparam_read
      itransp=1

ccc   central hot particle birth energy eh0 and density, and  beta
ccc   ihdist=1 for singular pitch angle distribution at galh0
ccc   ihdist=0 for nonsingular pitch angle distribution
ccc   ihdir=0 for both co- and counter-going particles
ccc   ihdir=1 for co-going particles
ccc   ihdir=2 for counter-going particles
c     TFTR
      rmaj=252.96
      amin=87.727
      b0=5.08045e4
ckg   ITER
      rmaj=814.0
      amin=277.8
      b0=3.e4
ckg   MAST
      rmaj=84.0
      amin=60.
      b0=0.39e4
ckg   NSTX q_0=2.8
      rmaj=84.8
      amin=68.
      b0=0.3e4
      denc0(1)=0.6E+14
c     DIID
      rmaj=168.1
      amin=63.4
      b0=1.87e4
      denc0(1)=0.45E+14
ckg some default values to be used in the test run
      rmaj=262.3
      amin=95
      b0=4.45e4
      denc0(1)=0.5E+14
ckg change   START
c      rmaj=25.1
c      amin=60.
c      b0=0.34e4
c      denc0(1)=0.33E+14
      zeff(1)=1.
ccc
      tc0(1)=6830.
ckg change 3.e3
      tc0(2)=1430.
ckg change 0.1128E+05 
      tc0(3)=1430.
      tc0(4)=1430.
      tc0(5)=tc0(2)
ckg change 0.1128E+05 
      denc0(2)=denc0(1)
      denc0(3)=0.000001*denc0(1)
      denc0(4)=0.000001*denc0(1)
      denc0(5)=0.000001*denc0(1)

      rmcp(1)=1./1836.0
      rmcp(2)=2.
cnng19 to work with Astra
      rmcp(2)=4.
cnng13      rmcp(2)=2.5
      rmcp(3)=3.
c
      rmcp(4)=1.0
      rmcp(5)=12.0

      zc(1)=-1.0
      zc(2)=1.0
cnng19 to work with Astra
      zc(2)=2.0
      zc(3)=1.0
      zc(4)=1.0
      zc(5)=6.0
ckg change include Zeff
c      denc0(5)=denc0(1)*(zeff(1)-zc(2))/(zc(5)*(zc(5)-zc(2)))
c      denc0(2)=denc0(1)*(zeff(1)-zc(5))/(zc(2)*(zc(2)-zc(5)))

ccc   hot particle parameters

      betah0(1)=8.000E-02
      betah0(2)=2.e-10
      betah0(3)=3.e-10

      th0(1)=0.100E+07
ckg change START      th0(1)=0.3e5
      th0(2)=0.8e5
      th0(3)=3.52e6

      ehc0(1)=th0(1)/9.
      ehc0(2)=th0(2)/9.
      ehc0(3)=th0(3)/9.

      ihdist(1)=0
      ihdist(2)=0
      ihdist(3)=0
      if(tip.ne.'x') tip='s'

      ihdir(1)=0
      ihdir(2)=0
      ihdir(3)=0

      galh0(1)=0.3
      galh0(2)=0.3
      galh0(3)=0.3

      dhlamda(1)=10000.0
      dhlamda(2)=10000.0
      dhlamda(3)=10000.0

      rmhp(1)=2.0
      rmhp(2)=3.0
      rmhp(3)=4.0

      zh(1)=1.0
      zh(2)=1.0
      zh(3)=2.0

      if(iswitch.ne.1) then
cccc  for ICRF proton tail
      hscale=3.75
      ascale=0.4
      th0(3)=1.0e6
      ihdist(3)=1
      galh0(3)=0.9
      dhlamda(3)=0.1
      rmhp(3)=1.0
      zh(3)=1.0
      endif
ccc
ccc   electron density profile parameters
      iden=1 
      alphar=-1.394
      alphar=0.
      prho=1.438
      arho=-1.643
      arho1=1.716
      arho2=-0.936
c MAST
c      alphar=1.
c      prho=1.62
c      arho=0.48
c NSTX
c      alphar=1.
c      prho=10.1
c      arho=0.122
c default values for fast ion profiles:
c      ph(i,1)=ph0(1)*(1.0001-1.*rsq**1.1)**10.1
c      exp(-abs((sqrt(rsq)-alpharh)/arhoh)**prhoh)
      alpharh=0.
      prhoh=2.
      arhoh=1.
c default values for bakcground temperature profile:
c      tc(i,1)=(1-alphart*rsq**prhot)**arhot
c
      alphart=0.
      prhot=2.
      arhot=1.

c START
ckg      alphar=1.
ckg      prho=0.15
ckg      arho=0.8499
c
      return
      end
ccccc.........
c the difference with the inprofa from NOVAK is in the analytic profiles for alphas ph(i,{1,2,3})
c it is left compatible with the NOVA_param file specification
      subroutine inprofa
c
cccc.....   define all radial profiles
c
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      integer*8 nadres,ierr,lhhh,length
      common/plasmap/rmaj,amin,b0

c
ccc..  read p-profile
      nadres=50
c      call zrd(mp1,pprf,nosurf,nadres+ishft,lgivup,ierr)
      lhhh=nosurf
      call getwa('mpout1', pprf, nadres,lhhh,ierr )
c
      nadres=nadres+nosurf
c      call zrd(mp1,ppprf,nosurf,nadres+ishft,lgivup,ierr)
      nadres=nosurf+50
      call getwa('mpout1', ppprf, nadres,lhhh,ierr )
      if(ishft.eq.0)then
         do is=1,nosurf
            ppprf(is)=ppprf(is)/f
         enddo
      endif
c
      psitot=(psilim-psimin)/twopi
      if(pprf(nosurf).le.0.0) pprf(nosurf)=-pprf(nosurf)
c
ccc   calculate central mass density
      rho0=0.0
      do is=2,ict
      rho0=rho0+denc0(is)*rmcp(is)
      end do
c
      beta0=2.*pprf(1)/bsqav(1)
cc..  Velocities are defined as sqrt(2.*T/m)
      do 16 is=1,ict
      vc0(is)=9.79e5*sqrt(2.*tc0(is)/rmcp(is))
      betac0(is)=4.03*tc0(is)*denc0(is)/(bsqav(1)*b0**2*1.0e11)
      pc0(is)=betac0(is)*pprf(1)/beta0
      wcc(is)=9.58e3*zc(is)*b0/rmcp(is)
   16 continue
c
      do 26 is=1,iht
      vh0(is)=9.79e5*sqrt(2.*th0(is)/rmhp(is))
      ph0(is)=betah0(is)*pprf(1)/beta0
      wch(is)=9.58e3*zh(is)*b0/rmhp(is)
   26 continue

cc..  valpha is the alpha birth velocity or ICRF tail ion thermal speed, etc.
c
      valpha=vh0(3)
c
c
      pie=acos(-1.)
      twopi=2.*pie
c
  678 format("** read itype")
      do 98 i=1, nosurf
      rsq=rgrid(i)**2
c
cc TFTR DT 66887 mass density profile
c      rhoprf(i)=(1.3+5.3*(1.-rsq)*(1.-0.95*rsq*(1.-rsq)))/6.6
c
cc JT-60U ICRF XP's mass density profile
ckg change
      rhoprf(i)=min((1.0001-alphar*rsq**prho)**arho,1.)
c      rhoprf(i)=min((1.001-alphar*rsq-prho*rsq**2)**arho,1.)
ckg      rhoprf(i)=(1.0001-rsq*prho-rsq**3*arho)
c      rhoprf(i)=1.+rsq*(alphar+ rsq*(prho+ rsq*(arho +
c     &	rsq*(arho1 + rsq*arho2))))
c
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  electron density, temperature & pressure profiles
ccc..  alne is the electron collision Coulomb Logarithm
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cc JT-60U ICRF XP's mass density profile
      denc(i,1)=rhoprf(i)
c      alphan=.99999
ckg      tc(i,1)=(1.-alphan*rsq)
ckg this is start_36484 shot
c      write(*,*) '2nd',i
c(1.003-rsq**1.001)**0.93755
ck change gpprf(i)/pprf(1)/denc(i,1)
c(1.-alphan*rsq)**1.544
      tc(i,1)=(1.001-alphart*rsq**prhot)**arhot
c
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  thermal ion density, temperature & pressure profiles
ccc..  alnd is the thermal ion collision Coulomb Logarithm
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c TFTR DT 66887
c      denc(i,2)=(9.2+5.4*(1.-rsq)+0.32*sin(twopi*(rsq-0.5))
c     &          -10.89*exp(-((rsq-0.14)/2.1)**2))
c     &         /(14.6-10.89*exp(-(0.14/2.1)**2))
c      tc(i,2)=(0.6+17.*(1.-rsq)*(1.-1.85*rsq*(1.-rsq)))/17.6
cc JT-60U ICRF XP's mass density profile
      denc(i,2)=rhoprf(i)
      tc(i,2)=tc(i,1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c TFTR DT 66887
c      denc(i,3)=(9.2+5.4*(1.-rsq)+0.32*sin(twopi*(rsq-0.5))
c     &          -10.89*exp(-((rsq-0.14)/2.1)**2))
c     &         /(14.6-10.89*exp(-(0.14/2.1)**2))
c      tc(i,3)=(0.6+17.*(1.-rsq)*(1.-1.85*rsq*(1.-rsq)))/17.6
cc JT-60U ICRF XP's mass density profile
      denc(i,3)=rhoprf(i)
      tc(i,3)=tc(i,1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc JT-60U ICRF XP's mass density profile
      denc(i,5)=rhoprf(i)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc JT-60U ICRF XP's mass density profile
      denc(i,4)=rhoprf(i)
      tc(i,4)=tc(i,1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do is=1,ict
      denc(i,is)=denc0(is)*denc(i,is)
      tc(i,is)=tc0(is)*tc(i,is)
      pc(i,is)=tc(i,is)*denc(i,is)/(denc0(is)*tc0(is))
      pc(i,is)=pc0(is)*pc(i,is)
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Coulomb Logarithm for computing alpha slowing-down distribution
      alnc(i,1)=24.-alog(sqrt(denc(i,1))/tc(i,1))
      do is=2,ict
      alnc(i,is)=14.2+alog(sqrt(tc(i,1)/denc(i,1))
     &         *rmhp(3)*rmcp(is)*valpha/(rmhp(3)+rmcp(is)))
      end do

ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  beam density, temperature & pressure profiles
ccc..  ebdprf is the beam birth energy radial profile
ccc..  ebdc is the lower beam critical energy profile
ccc..  assume ebdc has same radial profile as teprf
ccc..  ebdcp is the psi derivative of ln(ebdc)
ccc..  alnbd is the beam ion collision Coulomb Logarithm
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      th(i,1)=th0(1)

ckg NN Gorelenkov from orbit2d file
      z1=.5*(1.+denc(i,4)/denc(i,1))-denc(i,3)/denc(i,1)/6.
      xc2=1.294*sqrt(4.*th(i,1)/(3.52e+6*rmhp(1)))
      xc2=z1**0.666666*(.16784/xc2*sqrt(tc(i,1)/1.e3))**2
      ehc(i,1)=xc2*th(i,1)
ckg      ehc(i,1)=ehc0(1)*(tc(i,1)/tc0(1))
ckg      xc2=ehc(i,1)/th(i,1)
      xc=sqrt(xc2)
      xc3=xc2*xc
      if(xc.ne.0.) then
      z1=(1.+xc)/sqrt(1.-xc+xc2)
      z2=atan((2.-xc)/(sqrt(3.)*xc))
      z3=atan(1./sqrt(3.))+z2
      ghprf(i,1)=alog(z1)-sqrt(3.)*z3
      ghprf(i,1)=1.0+(2./3.)*xc2*ghprf(i,1)
      else
      ghprf(i,1)=1.0
      end if

c TFTR DT 66887
c      ph(i,1)=ph0(1)*(1.-rsq)**1.93
c      ph(i,3)=ph0(3)*(1.-rsq)**1.93
c JET 50725c01
c      ph(i,1)=ph0(1)*(1.0001-rsq**0.6)**6.2     !t=67sec
c      ph(i,1)=ph0(1)*(1.0001-rsq**1.1)**10.1    !t=68.6sec
c     ph(i,1)=ph0(1)*(1.0001-rsq**0.9)**3.86         !t=70.5sec
c      ph(i,3)=ph0(3)*(1.-rsq)**1.93
c      ph(i,1)=ph0(1)*(1.00001-alpharh*rsq**prhoh)**arhoh
c      ph(i,3)=ph0(3)*(1.00001-alpharh*rsq**prhoh)**arhoh
c      exp(-abs((sqrt(rsq)-alpharh)/arhoh)**prhoh)
cc
      ph(i,1)=ph0(1)*exp(-abs((rgrid(i)-alpharh)/arhoh)**prhoh)
c12      ph(i,1)=ph0(1)*exp( -(alpharh/prhoh)*tanh( (sqrt(rsq)-arhoh)
c12     & /alpharh ))*0.5
cnng13      ph(i,1)=pc(i,2)
      ph(i,2)=ph(i,1)
c12
      ph(i,3)=ph(i,1)
      if(tip.ne.'q')ph(i,3)=ph0(3)*exp(-abs((rgrid(i)-alpharh)/arhoh)
     &     **prhoh)
c
      cdenbd0=betah0(1)*(bsqav(1)*b0**2*1.0e11)/(4.03*th0(1))
c
ccc    denbd(i)=pbdprf(i)*alog((1.+xc3)/xc3)/gbdprf(i)
c
      if(ehc(i,1).lt.1.e3) then
c   if ebdc is zero the density becomes infinity, therefore set
      ebcutoff=1.e3
      yc3=(ebcutoff/th(i,1))**1.5
      else
      yc3=xc3
      end if
      denh(i,1)=cdenbd0*alog((1.+yc3)/yc3)/ghprf(i,1)
     &         *(ph(i,1)/ph0(1))/(th(i,1)/th0(1))
c
      th(i,2)=th(i,1)
      ehc(i,2)=ehc(i,1)
      ghprf(i,2)=ghprf(i,1)
      denh(i,2)=denh(i,1)

      do is=1,iht-1
      alnh(i,is)=14.2+alog(sqrt(tc(i,1)/denc(i,1))
     &         *rmhp(3)*rmhp(is)*valpha/(rmhp(3)+rmhp(is)))
      end do



      if(iswitch.eq.1) then

ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  hot Alpha particle pressure profile
ccc..  eaprf is the hot Alpha birth energy radial profile
ccc..  ealphc is the lower hot Alpha critical energy profile
ccc..  assume ealphc has same radial profile as teprf
ccc..  ealphcp is the psi derivative of ln(ealphc)
ccc..  alna is the hot ion (alpha) collision Coulomb Logarithm
ccc..  ph = ch * ealpha * ghprf, ealpha is the alpha birth energy (a constant)
ccc..  ch = rkhc * ghprf; 
ccc..  rkhc(r) corresponds to pitch angle integration and poloidal average
ccc..  ghprf(r) corresponds to energy integration 
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
      if(tip.ne.'q')th(i,3)=th0(3)

cccc   assume p-alpha is a polynomial
c2d      rmassdp=2.0
c2d      rmasstp=3.0
c2d      almd=14.2+alog(sqrt(tc(i,1)/denc(i,1))
c2d     &         *rmhp(3)*rmassdp*valpha/(rmhp(3)+rmassdp))
c2d      almt=14.2+alog(sqrt(tc(i,1)/denc(i,1))
c2d     &         *rmhp(3)*rmasstp*valpha/(rmhp(3)+rmasstp))
c2d
c.. assume 50/50 D-T for both thermal and beam ions
c2d
c2d      dz=zc(2)**2*denc(i,2)*(almd/rmassdp+almt/rmasstp)/2.0
c2d     &  +zh(1)**2*denh(i,1)*(almd/rmassdp+almt/rmasstp)/2.0
c2d     &  +zc(5)**2*denc(i,5)*alnc(i,5)/rmcp(5)
c2d	alne=alnc(i,1)
c2d      dz=(0.75*sqrt(pie)*dz/(rmcp(1)*dene(i)*alne(i)))**(2./3.)
c2d      dz=(0.75*sqrt(pie)*dz/(rmcp(1)*denc(i,1)*alne))**(2./3.)
c2d      ehc(i,3)=dz*rmhp(3)*rmcp(1)*tc(i,1)
cc  for comparison purpose with other codes, set   ccccccc
ckg NN Gorelenkov from orbit2d file
      z1=.5*(1.+denc(i,4)/denc(i,1))-denc(i,3)/denc(i,1)/6.
      xc2=1.294*sqrt(4.*th(i,3)/(3.52e+6*rmhp(3)))
      xc2=z1**0.666666*(.16784/xc2*sqrt(tc(i,1)/1.e3))**2
      ehc(i,3)=xc2*th(i,3)
ckg      ehc(i,1)=ehc0(1)*(tc(i,1)/tc0(1))
ckg      xc2=ehc(i,1)/th(i,1)
c
ckg      ehc(i,3)=ehc0(3)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  slowing-down energy integration
ckg      xc2=ehc(i,3)/th(i,3)

      if(xc2.ne.0.) then
c2d	write(*,*) sqrt(xc2)
      xc=sqrt(xc2)
      z1=(1.+xc)/sqrt(1.-xc+xc2)
      z2=atan((2.-xc)/(sqrt(3.)*xc))
      z3=atan(1./sqrt(3.))+z2
      ghprf(i,3)=alog(z1)-sqrt(3.)*z3
      ghprf(i,3)=1.0+(2./3.)*xc2*ghprf(i,3)
      else
      ghprf(i,3)=1.0
      end if

c START comparison for ion Landau damping
c      ph(i,3)=ph0(3)*rhoprf(i)/rhoprf(1)
cph0(3)*(1.-rsq**1.23)**6.5
c TFTR DT 66887
cpprf(1)*(1.-rsq)**3.75
ckg this part is for iswitch != 1, we do not normally use it
      else

cccc for JT-60U ICRF tail protons
cccc ph = n * Th * rkhc for Maxwellian energy distribution
cccc nh = n * rkhc for Maxwellian energy distribution
      if(tip.ne.'q')th(i,3)=th0(3)
      ghprf(i,3)=1.0
c      ph(i,3)=pprf(1)*(1.-rsq)**hscale
      if(tip.ne.'q')ph(i,3)=pprf(1)*exp(-rsq/ascale**2)
c      denh(i,3)=(ph(i,3)*th(1,3))/(th(i,3)*ph(1,3))
      end if


ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccxxxx..  copmpute zeff   ...xxxxxxxxxxxx
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      zeff(i)=0.0
      do 112 is=2,ict
      zeff(i)=zeff(i)+zc(is)**2*denc(i,is)
  112 continue
      do 122 is=1,iht
       zeff(i)=zeff(i)+zh(is)**2*denh(i,is)
  122 continue
      zeff(i)=zeff(i)/denc(i,1)
ccc  for comparison purpose with other codes, set   ccccccc
ckg change      zeff(i)=3.
c1.704
c

   98 continue
      return
      end
ccccc.........
      subroutine inprofa2
c
cccc.....   define all radial profiles
c
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      integer*8 nadres,ierr,lhhh,length
      common/plasmap/rmaj,amin,b0

c
ccc..  read p-profile
      nadres=50
      lhhh=nosurf
c      call zrd(mp1,pprf,nosurf,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', pprf, nadres,lhhh,ierr )
c
      nadres=nadres+nosurf
c      call zrd(mp1,ppprf,nosurf,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', ppprf, nadres,lhhh,ierr )
      if(ishft.eq.0)then
         do is=1,nosurf
            ppprf(is)=ppprf(is)/f
         enddo
      endif
c
      psitot=(psilim-psimin)/twopi
      if(pprf(nosurf).le.0.0) pprf(nosurf)=-pprf(nosurf)
c
ccc   calculate central mass density
      rho0=0.0
      do is=2,ict
         rho0=rho0+denc0(is)*rmcp(is)
      end do
c
      beta0=2.*pprf(1)/bsqav(1)
cc..  Velocities are defined as sqrt(2.*T/m)
      tc0(2)=beta0*rmcp(2)/(2.*4.03*rho0)*(bsqav(1)*b0**2*1.0e11)
c      print *,'test for anlt prfl',bsqav(1),b0**2,rmcp(2),denc0(2),tc0(2
c     &     ),beta0,rho0/(denc0(2)*rmcp(2))
c      tc0(2)=3000.
      do 16 is=1,ict
         tc0(is)=tc0(2)
         vc0(is)=9.79e5*sqrt(2.*tc0(is)/rmcp(is))
         if(is.gt.1) then
            betac0(is)=beta0*denc0(is)*rmcp(is)/(rho0*2.)
         else
            betac0(is)=beta0*0.5
         endif
c4.03*tc0(is)*denc0(is)/(bsqav(1)*b0**2*1.0e11)
         pc0(is)=betac0(is)*pprf(1)/beta0
         wcc(is)=9.58e3*zc(is)*b0/rmcp(is)
   16 continue
c
      do 26 is=1,iht
      vh0(is)=9.79e5*sqrt(2.*th0(is)/rmhp(is))
      ph0(is)=betah0(is)*pprf(1)/beta0
      wch(is)=9.58e3*zh(is)*b0/rmhp(is)
   26 continue

cc..  valpha is the alpha birth velocity or ICRF tail ion thermal speed, etc.
c
      valpha=vh0(3)
c
c
      pie=acos(-1.)
      twopi=2.*pie
c
  678 format("** read itype")
      do 98 i=1, nosurf
      rsq=rgrid(i)**2
c
cc TFTR DT 66887 mass density profile
c      rhoprf(i)=(1.3+5.3*(1.-rsq)*(1.-0.95*rsq*(1.-rsq)))/6.6
c
cc JT-60U ICRF XP's mass density profile
c      rhoprf(i)=1.+rsq*(alphar+ rsq*(prho+ rsq*(arho +
c     &	rsq*(arho1 + rsq*arho2))))
c profile with prescribed parameters from inparam
      rhoprf(i)=(1.0001-alphar*rsq**prho)**arho
c
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  electron density, temperature & pressure profiles
ccc..  alne is the electron collision Coulomb Logarithm
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cc JT-60U ICRF XP's mass density profile
      denc(i,1)=rhoprf(i)
      tc(i,1)=pprf(i)/pprf(1)/denc(i,1)
c      alphan=.9999
c(1.-alphan*rsq)**1.544
c
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  thermal ion density, temperature & pressure profiles
ccc..  alnd is the thermal ion collision Coulomb Logarithm
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c TFTR DT 66887
c      denc(i,2)=(9.2+5.4*(1.-rsq)+0.32*sin(twopi*(rsq-0.5))
c     &          -10.89*exp(-((rsq-0.14)/2.1)**2))
c     &         /(14.6-10.89*exp(-(0.14/2.1)**2))
c      tc(i,2)=(0.6+17.*(1.-rsq)*(1.-1.85*rsq*(1.-rsq)))/17.6
cc JT-60U ICRF XP's mass density profile
      denc(i,2)=rhoprf(i)
      tc(i,2)=tc(i,1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c TFTR DT 66887
c      denc(i,3)=(9.2+5.4*(1.-rsq)+0.32*sin(twopi*(rsq-0.5))
c     &          -10.89*exp(-((rsq-0.14)/2.1)**2))
c     &         /(14.6-10.89*exp(-(0.14/2.1)**2))
c      tc(i,3)=(0.6+17.*(1.-rsq)*(1.-1.85*rsq*(1.-rsq)))/17.6
cc JT-60U ICRF XP's mass density profile
      denc(i,3)=rhoprf(i)
      tc(i,3)=tc(i,1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc JT-60U ICRF XP's mass density profile
      denc(i,4)=rhoprf(i)
      tc(i,4)=tc(i,1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do is=1,ict
      denc(i,is)=denc0(is)*denc(i,is)
      tc(i,is)=tc0(is)*tc(i,is)
      pc(i,is)=tc(i,is)*denc(i,is)/(denc0(is)*tc0(is))
      pc(i,is)=pc0(is)*pc(i,is)
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Coulomb Logarithm for computing alpha slowing-down distribution
      alnc(i,1)=24.-alog(sqrt(denc(i,1))/tc(i,1))
      do is=2,ict
      alnc(i,is)=14.2+alog(sqrt(tc(i,1)/denc(i,1))
     &         *rmhp(3)*rmcp(is)*valpha/(rmhp(3)+rmcp(is)))
      end do


ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  beam density, temperature & pressure profiles
ccc..  ebdprf is the beam birth energy radial profile
ccc..  ebdc is the lower beam critical energy profile
ccc..  assume ebdc has same radial profile as teprf
ccc..  ebdcp is the psi derivative of ln(ebdc)
ccc..  alnbd is the beam ion collision Coulomb Logarithm
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      th(i,1)=th0(1)

ckg      ehc(i,1)=ehc0(1)*(tc(i,1)/tc0(1))
ckg
ckg      xc2=ehc(i,1)/th(i,1)
ckg NN Gorelenkov from orbit2d file
      z1=.5*(1.+denc(i,4)/denc(i,1))-denc(i,3)/denc(i,1)/6.
      xc2=1.294*sqrt(4.*th(i,1)/(3.52e+6*rmhp(1)))
      xc2=z1**0.666666*(.16784/xc2*sqrt(tc(i,1)/1.e3))**2
      ehc(i,1)=xc2*th(i,1)
ckg      ehc(i,1)=ehc0(1)*(tc(i,1)/tc0(1))
ckg      xc2=ehc(i,1)/th(i,1)
      xc=sqrt(xc2)
      xc3=xc2*xc
      if(xc.ne.0.) then
      z1=(1.+xc)/sqrt(1.-xc+xc2)
      z2=atan((2.-xc)/(sqrt(3.)*xc))
      z3=atan(1./sqrt(3.))+z2
      ghprf(i,1)=alog(z1)-sqrt(3.)*z3
      ghprf(i,1)=1.0+(2./3.)*xc2*ghprf(i,1)
      else
      ghprf(i,1)=1.0
      end if

c TFTR DT 66887
c     ph(i,1)=ph0(1)*(1.-rsq)**1.93
ckg     The profile for beam beta for NSTX q_0 = 2.8
c      ph(i,1)=ph0(1)*(1.-rsq**1.33)**3.4
ckg     The profile is for FIRE 2002 snowmass
c      ph(i,1)=ph0(1)*(1.-rsq**0.55)**1.2
ckg example of the profile definition for fast ion species
c D-beams 
c      ph(i,1)=ph0(1)*(pprf(i)/pprf(1))**2
c alphas
c      ph(i,3)=ph0(3)*(pprf(i)/pprf(1))**2
c
      ph(i,1)=ph0(1)*exp(-abs((rgrid(i)-alpharh)/arhoh)**prhoh)
      ph(i,3)=ph0(3)*exp(-abs((rgrid(i)-alpharh)/arhoh)**prhoh)
      ph(i,2)=ph(i,1)
c
      cdenbd0=betah0(1)*(bsqav(1)*b0**2*1.0e11)/(4.03*th0(1))
c
ccc    denbd(i)=pbdprf(i)*alog((1.+xc3)/xc3)/gbdprf(i)
c
      if(ehc(i,1).lt.1.e3) then
c   if ebdc is zero the density becomes infinity, therefore set
      ebcutoff=1.e3
      yc3=(ebcutoff/th(i,1))**1.5
      else
      yc3=xc3
      end if

      denh(i,1)=cdenbd0*alog((1.+yc3)/yc3)/ghprf(i,1)
     &         *(ph(i,1)/ph0(1))/(th(i,1)/th0(1))
c
      th(i,2)=th(i,1)
      ehc(i,2)=ehc(i,1)
      ghprf(i,2)=ghprf(i,1)
      denh(i,2)=denh(i,1)
      ph(i,2)=ph(i,1)
ckg      ph(i,3)=ph(i,1)

      do is=1,iht-1
      alnh(i,is)=14.2+alog(sqrt(tc(i,1)/denc(i,1))
     &         *rmhp(3)*rmhp(is)*valpha/(rmhp(3)+rmhp(is)))
      end do



      if(iswitch.eq.1) then

ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  hot Alpha particle pressure profile
ccc..  eaprf is the hot Alpha birth energy radial profile
ccc..  ealphc is the lower hot Alpha critical energy profile
ccc..  assume ealphc has same radial profile as teprf
ccc..  ealphcp is the psi derivative of ln(ealphc)
ccc..  alna is the hot ion (alpha) collision Coulomb Logarithm
ccc..  ph = ch * ealpha * ghprf, ealpha is the alpha birth energy (a constant)
ccc..  ch = rkhc * ghprf; 
ccc..  rkhc(r) corresponds to pitch angle integration and poloidal average
ccc..  ghprf(r) corresponds to energy integration 
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
cccc   assume p-alpha is a polynomial
      th(i,3)=th0(3)

c2d      rmassdp=2.0
c2d      rmasstp=3.0
c2d      almd=14.2+alog(sqrt(tc(i,1)/denc(i,1))
c2d     &         *rmhp(3)*rmassdp*valpha/(rmhp(3)+rmassdp))
c2d      almt=14.2+alog(sqrt(tc(i,1)/denc(i,1))
c2d     &         *rmhp(3)*rmasstp*valpha/(rmhp(3)+rmasstp))
c2d
c.. assume 50/50 D-T for both thermal and beam ions
c2d
c2d      dz=zc(2)**2*denc(i,2)*(almd/rmassdp+almt/rmasstp)/2.0
c2d     &  +zh(1)**2*denh(i,1)*(almd/rmassdp+almt/rmasstp)/2.0
c2d     &  +zc(5)**2*denc(i,5)*alnc(i,5)/rmcp(5)
c2d	alne=alnc(i,1)
c2d      dz=(0.75*sqrt(pie)*dz/(rmcp(1)*dene(i)*alne(i)))**(2./3.)
c2d      dz=(0.75*sqrt(pie)*dz/(rmcp(1)*denc(i,1)*alne))**(2./3.)
c2d      ehc(i,3)=dz*rmhp(3)*rmcp(1)*tc(i,1)
cc  for comparison purpose with other codes, set   ccccccc
ckg NN Gorelenkov from orbit2d file
      z1=.5*(1.+denc(i,4)/denc(i,1))-denc(i,3)/denc(i,1)/6.
      xc2=1.294*sqrt(4.*th(i,3)/(3.52e+6*rmhp(3)))
      xc2=z1**0.666666*(.16784/xc2*sqrt(tc(i,1)/1.e3))**2
      ehc(i,3)=xc2*th(i,3)
ckg
ckg      ehc(i,3)=ehc0(3)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  slowing-down energy integration
ckg      xc2=ehc(i,3)/th(i,3)

      if(xc2.ne.0.) then
c2d	write(*,*) sqrt(xc2)
      xc=sqrt(xc2)
      z1=(1.+xc)/sqrt(1.-xc+xc2)
      z2=atan((2.-xc)/(sqrt(3.)*xc))
      z3=atan(1./sqrt(3.))+z2
      ghprf(i,3)=alog(z1)-sqrt(3.)*z3
      ghprf(i,3)=1.0+(2./3.)*xc2*ghprf(i,3)
      else
      ghprf(i,3)=1.0
      end if

c TFTR DT 66887
c      ph(i,3)=ph0(3)*(1.-rsq**1.23)**6.5
cpprf(1)*(1.-rsq)**3.75
cben with HINT
c      ph(i,3)=pprf(i)
c      ph(i,1)=pprf(i)
      else

cccc for JT-60U ICRF tail protons
cccc ph = n * Th * rkhc for Maxwellian energy distribution
cccc nh = n * rkhc for Maxwellian energy distribution
      th(i,3)=th0(3)
      ghprf(i,3)=1.0
c      ph(i,3)=pprf(1)*(1.-rsq)**hscale
c      ph(i,3)=pprf(1)*exp(-rsq/ascale**2)
c      denh(i,3)=(ph(i,3)*th(1,3))/(th(i,3)*ph(1,3))
      end if


ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccxxxx..  copmpute zeff   ...xxxxxxxxxxxx
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      zeff(i)=0.0
      do 112 is=2,ict
      zeff(i)=zeff(i)+zc(is)**2*denc(i,is)
  112 continue
ckg      do 122 is=1,iht
ckg      zeff(i)=zeff(i)+zh(is)**2*denh(i,is)
ckg  122 continue
      zeff(i)=zeff(i)/denc(i,1)
ccc  for comparison purpose with other codes, set   ccccccc
c      zeff(i)=2.5
c1.704
c

   98 continue
      return
      end
ccccc.........
ccc
      subroutine inprofe
c
cccc.....   define all radial profiles
c
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      integer*8 nadres,ierr,lhhh,length
      common/plasmap/rmaj,amin,b0

ccc
ccc   central hot particle birth energy eh0 and density, and  beta
ccc   ihdist=1 for singular pitch angle distribution at galh0
ccc   ihdist=0 for nonsingular pitch angle distribution
ccc   ihdir=0 for both co- and counter-going particles
ccc   ihdir=1 for co-going particles
ccc   ihdir=2 for counter-going particles
c
c
ccc..  read p-profile
      nadres=50
      lhhh=nosurf
c      call zrd(mp1,pprf,nosurf,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', pprf, nadres,lhhh,ierr )
      nadres=nadres+nosurf
c      call zrd(mp1,ppprf,nosurf,nadres+ishft,lgivup,ierr)
      call getwa('mpout1', ppprf, nadres,lhhh,ierr )
      if(ishft.eq.0)then
         do is=1,nosurf
            ppprf(is)=ppprf(is)/f
         enddo
      endif
c
      psitot=(psilim-psimin)/twopi
      if(pprf(nosurf).le.0.0) pprf(nosurf)=-pprf(nosurf)
c
      rho0=rhoprf(1)
c
      beta0=2.*pprf(1)/bsqav(1)

cc..  Velocities are defined as sqrt(2.*T/m)
ckg
      bet_trnsp=0.
      do 16 is=1,ict
      denc0(is)=max(denc(1,is),1.)
      tc0(is)=tc(1,is)
      vc0(is)=9.79e5*sqrt(2.*tc0(is)/rmcp(is))
      betac0(is)=4.03*tc0(is)*denc0(is)/(bsqav(1)*b0**2*1.0e11)
      pc0(is)=betac0(is)*pprf(1)/beta0
      wcc(is)=9.58e3*zc(is)*b0/rmcp(is)
c      write(*,*) 'betas= ',is,betac0(is),pc0(is)
      bet_trnsp=bet_trnsp+betac0(is)
   16 continue
c
      do 26 is=1,iht
      denh0(is)=denh(1,is)
      vh0(is)=9.79e5*sqrt(2.*th0(is)/rmhp(is))
      ph0(is)=ph(1,is)/ptot(1)*pprf(1)
      betah0(is)=ph0(is)/pprf(1)*beta0
      wch(is)=9.58e3*zh(is)*b0/rmhp(is)
c      write(*,*) 'betas= ',is,betah0(is)
      bet_trnsp=bet_trnsp+betah0(is)
   26 continue
ckg Note that this factor is for scaling of the total plasma beta which
c     can be different from TRANSP's
c     Also we scale the temperature for the core ions and electrons
c     and density for the hot species.
c     to avoud any normalization just make bet_trnsp=1.
      bet_trnsp=beta0/bet_trnsp
      bet_trnsp=1.
      bet_trnsp1=0.
      do is=1,ict
         tc0(is)=tc0(is)*bet_trnsp
         vc0(is)=9.79e5*sqrt(2.*tc0(is)/rmcp(is))
         betac0(is)=4.03*tc0(is)*denc0(is)/(bsqav(1)*b0**2*1.0e11)
         pc0(is)=betac0(is)*pprf(1)/beta0
         bet_trnsp1=bet_trnsp1+betac0(is)
      enddo
c
      do is=1,iht
ckg         denh0(is)=denh(1,is)*bet_trnsp
ckg         vh0(is)=9.79e5*sqrt(2.*th0(is)/rmhp(is))
ckg         ph0(is)=ph(1,is)/ptot(1)*pprf(1)*bet_trnsp
ckg         betah0(is)=ph0(is)/pprf(1)*beta0
         bet_trnsp1=bet_trnsp1+betah0(is)
      enddo
c      write(*,*) 'bet_new=',bet_trnsp1,'norm=',bet_trnsp
c      write(*,*) ' B_f=',bf2d(1,1),sqrt(bsqav(1))
cc..  valpha is the ICRF tail ion thermal or alpha birth velocity, etc.
c
      valpha=vh0(3)
c
c
      pie=acos(-1.)
      twopi=2.*pie
c
ccc
      hhh=0.
      do i=1,nosurf
         hhh=max(hhh,th(i,3))
      enddo
      do 1 i=1,nosurf
         rsq=rgrid(i)**2
c
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc    density profile
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         rhoprf(i)=rhoprf(i)/rho0

ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  electron density, temperature & pressure profiles
ccc..  alne is the electron collision Coulomb Logarithm
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      do 111 is=1,ict
         tc(i,is)=abs(tc(i,is))*bet_trnsp
      pc(i,is)=tc(i,is)*denc(i,is)/(denc0(is)*tc0(is))
      pc(i,is)=pc0(is)*pc(i,is)
c
  111 continue
      do is=1,iht
         denh(i,is)=denh(i,is)*bet_trnsp
         ph(i,is)=ph(i,is)*bet_trnsp
      enddo
c
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  alnc is the thermal particle collision Coulomb Logarithm
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      alnc(i,1)=24.-alog(sqrt(denc(i,1))/tc(i,1))
      do 113 is=2,ict
      alnc(i,is)=14.2+alog(sqrt(tc(i,1)/denc(i,1))
     &          *rmhp(3)*rmcp(is)*valpha/(rmhp(3)+rmcp(is)))
  113 continue
c
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  alnh is the beam particle collision Coulomb Logarithm
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      do 112 is=1,iht-1
      th(i,is)=th0(is)
      alnh(i,is)=14.2+alog(sqrt(tc(i,1)/denc(i,1))
     &         *rmhp(3)*rmhp(is)*valpha/(rmhp(3)+rmhp(is)))
  112 continue
      if(tip(1:1).ne.'i')then
         th(i,3)=th0(3)
c      else
c         th(i,3)=th(i,3)*th0(3)/hhh
c     this is to test the 'i' df option, that is it should give the same
C     results as 'h' one
c         th(i,3)=th0(3)
      endif
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..  hot particle slowing-down distribution
ccc..  th is the hot particle birth energy radial profile
ccc..  ehc is the lower critical energy
ccc..  assume ehc has same radial profile as teprf for beam particles
ccc..  ehcp is the psi derivative of ln(ehc)
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
ckg NN Gorelenkov from orbit2d file
      z1=.5*(1.+denc(i,4)/denc(i,1))-denc(i,3)/denc(i,1)/6.
      xc2=1.294*sqrt(4.*th(i,1)/(3.52e+6*rmhp(1)))
      xc2=z1**0.666666*(.16784/xc2*sqrt(tc(i,1)/1.e3))**2
      ehc(i,1)=xc2*th(i,1)
      xc2=1.294*sqrt(4.*th(i,2)/(3.52e+6*rmhp(2)))
      xc2=z1**0.666666*(.16784/xc2*sqrt(tc(i,1)/1.e3))**2
      ehc(i,2)=xc2*th(i,2)
ckg      ehc(i,1)=ehc0(1)*(tc(i,1)/tc0(1))
ckg      ehc(i,2)=ehc0(2)*(tc(i,1)/tc0(1))
c
      dz=0.0
      do 110 is=2,ict
      dz=dz+zc(is)**2*denc(i,is)*alnc(i,is)/rmcp(is)
  110 continue
      do 120 is=1,iht-1
      dz=dz+zh(is)**2*denh(i,is)*alnh(i,is)/rmhp(is)
  120 continue
      dz=(0.75*sqrt(pie)*dz/(denc(i,1)/rmcp(1)*alnc(i,1)))**(2./3.)
      ehc(i,3)=dz*rmhp(3)/rmcp(1)*tc(i,1)
ckg NN Gorelenkov from orbit2d file
      z1=.5*(1.+denc(i,4)/denc(i,1))-denc(i,3)/denc(i,1)/6.
      xc2=1.294*sqrt(4.*th(i,3)/(3.52e+6*rmhp(3)))
      xc2=z1**0.666666*(.16784/xc2*sqrt(tc(i,1)/1.e3))**2
      ehc(i,3)=xc2*th(i,3)
ckg
ccc  for comparison purpose with other codes, set   ccccccc
c      ehc(i,3)=ehc0(3)

      do 121 is=1,iht
      ph(i,is)=ph(i,is)/ptot(1)*pprf(1)
c      th(i,is)=th0(is)
c
      xc2=ehc(i,is)/th(i,is)
      xc=sqrt(xc2)
      if(xc.ne.0.) then
      z1=(1.+xc)/sqrt(1.-xc+xc2)
      z2=atan((2.-xc)/(sqrt(3.)*xc))
      z3=atan(1./sqrt(3.))+z2
      ghprf(i,is)=alog(z1)-sqrt(3.)*z3
      ghprf(i,is)=1.0+(2./3.)*xc2*ghprf(i,is)
      else
      ghprf(i,is)=1.0
      end if

  121 continue
c
c  rescale alpha pressure to be same as total pressure at magnetic axis
c      ph(i,3)=pprf(1)*(ph(i,3)/ph0(3))

ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccxxxx..  copmpute zeff   ...xxxxxxxxxxxx
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      zeff(i)=0.0
      do 123 is=2,ict
         zeff(i)=zeff(i)+zc(is)**2*denc(i,is)
  123 continue
c      do 122 is=1,iht
c      zeff(i)=zeff(i)+zh(is)**2*denh(i,is)
c  122 continue
      zeff(i)=zeff(i)/denc(i,1)
ccc  for comparison purpose with other codes, set   ccccccc
c      zeff(i)=2.5
c1.704
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

1     continue
      return
      end
c
      subroutine plotprf
      include 'clich1'
      include 'clich2'
      include 'clich1b'

      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      common/pitchang/fpht(lam,nn,isph),dfpht(lam,nn,isph)
     &,fphc(lamc,nn,isph),dfphc(lamc,nn,isph)
      character*80 buf

      do 11 i=1,nosurf
      do 31 is=1,ict
      pwsc(i,is)=pc(i,is)*wstarc(i,is)
   31 continue

      do 41 is=1,iht
      pwsh(i,is)=ph(i,is)*wstarh(i,is)

      ehc32=ehc(i,is)**1.5
      eh32=th(i,is)**1.5
      ehmix=eh32+ehc32
      driveh(i,is)=-abs(omreal)*(dfphc(lamc,i,is)
     >                      -ehc32/ehmix*fphc(lamc,i,is))
     >         +(th(i,is)*eh32/ehmix-wstarh(i,is))*fphc(lamc,i,is)
   41 continue
   11 continue

      do 32 is=1,ict
         buf='pwsc'
         call posplt1(rgrid,pwsc(1,is),1,nosurf,1,buf)
         buf='tc'
         call posplt1(rgrid,tc(1,is),1,nosurf,2,buf)
         buf='denc'
         call posplt1(rgrid,denc(1,is),1,nosurf,3,buf)
         buf='wstarc'
         call posplt1(rgrid,wstarc(1,is),1,nosurf,4,buf)
         call frame(0)
         buf='tc'
         call posplt1(rgrid,tcp(1,is),1,nosurf,2,buf)
         buf='denc'
         call posplt1(rgrid,dencp(1,is),1,nosurf,3,buf)
         buf='wstarc'
         call posplt1(rgrid,wstartc(1,is),1,nosurf,4,buf)
      call frame(0)
   32 continue

      do 42 is=1,iht
         buf='pwsh'
         call posplt1(rgrid,pwsh(1,is),1,nosurf,1,buf)
         buf='ph'
         call posplt1(rgrid,ph(1,is),1,nosurf,2,buf)
         buf='driveh'
         call posplt1(rgrid,driveh(1,is),1,nosurf,3,buf)
         buf='wstarh'
         call posplt1(rgrid,wstarh(1,is),1,nosurf,4,buf)
         call frame(0)
         buf='ehc'
         call posplt1(rgrid,ehc(1,is),1,nosurf,1,buf)
         buf='ph'
         call posplt1(rgrid,ph(1,is),1,nosurf,2,buf)
         buf='ehc'
         call posplt1(rgrid,ehcp(1,is),1,nosurf,3,buf)
         buf='ph'
         call posplt1(rgrid,php(1,is),1,nosurf,4,buf)
      call frame(0)
   42 continue

       return
      end
c
      subroutine printprf
      include 'clich1'
      include 'clich2'
      include 'clich1b'
c
    2 format(10e12.3)
    3 format(10i5)


      write(17,*) 'RHOPRF profile'
      write(17,2) (rhoprf(isrf),isrf=1,nosurf)
      write(17,*) 'ZEFF profile'
      write(17,2) (zeff(isrf),isrf=1,nosurf)

      do 32 is=1,ict
      write(17,*) 'WSTARC profile'
      write(17,2) (wstarc(isrf,is),isrf=1,nosurf)
      write(17,*) 'WSTARTC profile'
      write(17,2) (wstartc(isrf,is),isrf=1,nosurf)
      write(17,*) 'DenC profile'
      write(17,2) (denc(isrf,is),isrf=1,nosurf)
      write(17,*) 'DenCP profile'
      write(17,2) (dencp(isrf,is),isrf=1,nosurf)
      write(17,*) 'TC profile'
      write(17,2) (tc(isrf,is),isrf=1,nosurf)
      write(17,*) 'TCP profile'
      write(17,2) (tcp(isrf,is),isrf=1,nosurf)
 
   32 continue

      do 42 is=1,iht
      write(17,*) 'WSTARH profile'
      write(17,2) (wstarh(isrf,is),isrf=1,nosurf)
      write(17,*) 'PH profile'
      write(17,2) (ph(isrf,is),isrf=1,nosurf)
      write(17,*) 'PHP profile'
      write(17,2) (php(isrf,is),isrf=1,nosurf)
      write(17,*) 'TH profile'
      write(17,2) (th(isrf,is),isrf=1,nosurf)
      write(17,*) 'THP profile'
      write(17,2) (thp(isrf,is),isrf=1,nosurf)

      write(17,*) 'DenH profile'
      write(17,2) (denh(isrf,is),isrf=1,nosurf)
      write(17,*) 'DenHP profile'
      write(17,2) (denhp(isrf,is),isrf=1,nosurf)
      write(17,*) 'EHC profile'
      write(17,2) (ehc(isrf,is),isrf=1,nosurf)
      write(17,*) 'EHCP profile'
      write(17,2) (ehcp(isrf,is),isrf=1,nosurf)

      write(17,*) 'RKHC profile'
      write(17,2) (rkhc(isrf,is),isrf=1,nosurf)
      write(17,*) 'GHPRF profile'
      write(17,2) (ghprf(isrf,is),isrf=1,nosurf)
      write(17,*) 'GPh profile'
      write(17,2) (gph(isrf,is),isrf=1,nosurf)
      write(17,*) 'CHPRF profile'
      write(17,2) (chprf(isrf,is),isrf=1,nosurf)
      write(17,*) 'CHPPRF profile'
      write(17,2) (chpprf(isrf,is),isrf=1,nosurf)

   42 continue

      return
      end
c
      subroutine taebeta(im1,ihsps,xfow,pk)
      include 'clich1'
      include 'clich1b'
      include 'clich2'
      common/cab/abe0,abe01    

      common/pitch/
     & wdbh(lam,isph),wbounh(lam,isph),wbouh(lam)
     &,wdch(lamc,isph),wcirch(lamc,isph),wcirh(lamc)
     &,wdbc(lam,ispc),wbounc(lam,ispc)
     &,wdcc(lamc,ispc),wcircc(lamc,ispc)

      dimension val(50),beta(50)
     &,chi(nn,isph),cci(nn,ispc)
     &,chit(nn,isph),ccit(nn,ispc)
     &,chic(nn,isph),ccic(nn,ispc)
     &,gabh(nn,isph),gach(nn,isph)
     &,gabh00(nn,isph),gabh01(nn,isph)
     &,gach0p(nn,isph),gach0m(nn,isph)
     &,cecol(nn)
     &,rgrb1(nn),rgrv1(nn),rgrk1(nn),falfven(nn)

      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
c
      common/plasmap/rmaj,amin,b0

      namelist/inp1/b0,rmaj,amin
     &,denc0,tc0,rmcp,zc,rmhp,zh
     &,th0,ehc0,betah0,ihdist,ihdir,galh0,dhlamda
     &,itransp,iswitch
     &,iden,alphar,prho,arho
     &,ascale,hscale

      character*16 nameh(isph),namec(ispc)
      character*80 buf
cc  set defolt input parameters
      call inparam
      call inparam_read(im1,ihsps,xfow,pk)
c
      pie=acos(-1.)
      twopi=2.*pie
c
ccc  choose the sign of om to be same as alpha's wstar to get instability
ccc  for n>0, om=-omreal; for n<0, om=omreal
      om=-abs(omreal)
      omsq=om**2
c      if(ntor.lt.0)then
c         ntor=-ntor
c         om=-om
c         omreal=-omreal
c      endif
c
      write(17,123) omsq
      write(iodat,123) omsq
      write(*,123) omsq
  123 format(1x,'(om/omA)**2 of the TAE mode is ',e12.4)
c
  333 continue
      write(17,55)
      write(iodat,55)
      write(*,55)
   55 format(1x,'**   type in ignited tokamak parameters')
c
c     call talk ("INP1 input desired ignited tokamak parameters.. ")
c2d       write(0,inp1)
c2d       read (5,inp1)

      write(*,556)
      write(iodat,556)
      write(17,556) 
c read the transrun id
      if(itransp.eq.1) then
         call transpin(ihsps)
         call inprofe
      else if(itransp.eq.2) then
         call inprofa2
      else 
         call inprofa
      end if
cnng15
c printing a brute force approach coefficient _xTebf for scaling Te, read mkorb to activate this option
#ifdef _xTebf
      print *,'used scaling factor for Te ',_xTebf
#endif
cnng13
      if(tip.eq.'x')then
         th0(ihsps)=tc(1,2)
      endif
c      print *,'------------- tip=',tip,th0(ihsps),tc(1,2),ihsps
cnng13
      write(*,556)
      write(iodat,556)
      write(17,556)
  556 format(/)

      write(*,inp1)
      write(iodat,inp1)
      write(17,inp1)

cc
cc  define pitch angle distribution for hot particles
      call pitdist
c
cc  note the definition of Alfven speed
      va=2.18e11*b0/sqrt(rho0)
      do 987 i=1,nosurf
      falfven(i)=va*sqrt(bsqav(i))/sqrt(rhoprf(i))
     >          /(2.*twopi*qprf(i)*rmaj*1.e3)
  987 continue

      ve0=vc0(1)
      vd0=vc0(2)
      vbd0=vh0(1)
      vbt0=vh0(2)
c
cc..  Velocities are defined as sqrt(2.*T/m)
cc..  valpha is the alpha birth velocity, or ICRF tail ion thermal velocity  etc.
c
      valphava=valpha/va
c
c Coulomb Logarithm alne=ln(lamda) for e-i collision, lamda=24-ln(sqrt(ne)/te)
c  amu0 is to be multiplied by alne=24.-alog(sqrt(dene)/te) in getting amu
      amu0=1.775e-17*rmaj*q1*sqrt(rmcp(2))/(b0*abs(om))
ckg This constant differs the parallel elctric field term from the
c   drift term, which is the default in Cheng's calculations
      colcoe=abs(om*va/(q1*rmaj*wcc(2)*psitot))
c
c      write(*,'(i4,2e15.6)')(i,qprf(i),rgrid(i),i=1,301)
      do 12 is=1,ict
      cstrc(is)=0.5*vc0(is)*q1*r*r/(wcc(is)*rmaj)*(vc0(is)/valpha)
      cdstc(is)=ntor*cstrc(is)
      cbouc(is)=twopi*q1*r*(vc0(is)/valpha)
      ccirc(is)=cbouc(is)
   12 continue
      write(17,456) (vc0(is),is=1,ict)
      write(iodat,456) (vc0(is),is=1,ict)
      write(*,456) (vc0(is),is=1,ict)
  456 format(' ve0=',e12.4,' vd0=',e12.4,' vt0=',e12.4,
     &       ' vh0=',e12.4,' vz0=',e12.4)
      write(17,57) (betac0(is),is=1,ict)
      write(iodat,57)(betac0(is),is=1,ict)
      write(*,57)(betac0(is),is=1,ict)
   57 format(' betac0(is)=',5e11.3)


      do 22 is=1,iht
      cstrh(is)=0.5*vh0(is)*q1*r*r/(wch(is)*rmaj)*(vh0(is)/valpha)
      cdsth(is)=ntor*cstrh(is)
      cbouh(is)=twopi*q1*r*(vh0(is)/valpha)
      ccirh(is)=cbouh(is)
   22 continue

      ftae=abs(om)*va/(q1*rmaj)/(twopi*1.e3)
      write(17,256) ftae
      write(iodat,256) ftae
      write(*,256) ftae
  256 format('*** The TAE mode frequency is ',f8.1,' kHz ****')

      write(17,56) va,valphava,(vh0(is),is=1,iht)
      write(iodat,56) va,valphava,(vh0(is),is=1,iht)
      write(*,56) va,valphava,(vh0(is),is=1,iht)
   56 format('In cm/sec valf=',e12.4,' valp/va=',e12.4,
     &       ' vbd0=',e12.4,' vbt0=',e12.4,' valpha=',e12.4)
      write(17,58) (betah0(is),is=1,iht)
      write(iodat,58) (betah0(is),is=1,iht)
      write(*,58) (betah0(is),is=1,iht)
   58 format(' betah0(is)=',5e11.3)

      write(17,59) beta0,amu0
      write(iodat,59) beta0,amu0
      write(*,59) beta0,amu0
   59 format(' beta0=',e12.4,' amu0=',e12.4)
c
cc compute psi derivative; rgrid is sqrt(psi)
c
ccc prepare cubic B-spline interpolation
c
      nosf=nosurf
      nosf2=nosf+2
      do 101 i=1,nosurf
      rgsq(i)=rgrid(i)**2
  101 continue
      call splprp(rgsq,nosf)
      call splprp1(nosf,nosf2,rgsq)
c
      call intsplp(nosf,nosf,nosf2,rhoprf,rhopprf,rgsq)
      do 137 i=1,nosurf
      rhopprf(i)=rhopprf(i)/psitot
  137 continue

      do 111 is=1,ict
      call intsplp(nosf,nosf,nosf2,denc(1,is),dencp(1,is),rgsq)
      call intsplp(nosf,nosf,nosf2,tc(1,is),tcp(1,is),rgsq)
      call intsplp(nosf,nosf,nosf2,pc(1,is),pcp(1,is),rgsq)
      do 138 i=1,nosurf
      dencp(i,is)=dencp(i,is)/psitot
      tcp(i,is)=tcp(i,is)/psitot
      pcp(i,is)=pcp(i,is)/psitot
  138 continue
  111 continue

      do 121 is=1,iht
      call intsplp(nosf,nosf,nosf2,denh(1,is),denhp(1,is),rgsq)
      call intsplp(nosf,nosf,nosf2,th(1,is),thp(1,is),rgsq)
      call intsplp(nosf,nosf,nosf2,ph(1,is),php(1,is),rgsq)
      call intsplp(nosf,nosf,nosf2,ehc(1,is),ehcp(1,is),rgsq)
      do 148 i=1,nosurf
      ehcp(i,is)=ehcp(i,is)/psitot
      denhp(i,is)=denhp(i,is)/psitot
      thp(i,is)=thp(i,is)/psitot
      php(i,is)=php(i,is)/psitot
  148 continue
  121 continue


cc obtain psi derivative of ln(ch), ln(cb), ln(ealphc), ln(ebdc) 
c
      rgx1=rgrid(nosurf-2)
      rgx2=rgrid(nosurf-1)
      rgx3=rgrid(nosurf)
      dxx=(rgx3-rgx2)/(rgx2-rgx1)
      dxx1=1.+dxx
c     dxx=1.0

      do 155 is=1,iht
      do 156 i=1,nosurf-1
c
      if(iswitch.ne.1 .and. is.eq.iht) then
         gph(i,is)=(th(i,is)/th0(is))**2.5*rkhc(i,is)
      else
         gph(i,is)=(th(i,is)/th0(is))*ghprf(i,is)*rkhc(i,is)
      endif
      chprf(i,is)=max(ph(i,is),1.e-8)/gph(i,is)
  156 continue
c
      chprf(nosurf,is)=chprf(nosurf-1,is)*dxx1-chprf(nosurf-2,is)*dxx
c
      call intsplp(nosf,nosf,nosf2,chprf(1,is),chpprf(1,is),rgsq)
      do 157 i=1,nosurf-1
         chpprf(i,is)=chpprf(i,is)/(psitot*chprf(i,is))
  157 continue
      chpprf(nosurf,is)=chpprf(nosurf-1,is)*dxx1-chpprf(nosurf-2,is)*dxx
c
  155 continue
c
c
cccccc   compute volume-averaged beta-alpha; <2.*Ph/B**2>
ccc   integrate over sqrt(psi)-grids;  assume uniform grid
      betaav=0.0
      betahav=0.0
      volrad=0.0
      dkeng=0.0

c      betahav1=0.0
c      volrad1=0.0
c      dkeng1=0.0

      do 41 isrf=2,nosurf
      isrf1=isrf-1
      drg=rgrid(isrf)-rgrid(isrf1)

c      rgrb1(isrf)=2.*rgrid(isrf)*ph(isrf,3)*bsqinv(isrf)
c      rgrv1(isrf)=2.*rgrid(isrf)*volr(isrf)
c      rgrk1(isrf)=2.*rgrid(isrf)*delk(isrf)*rhoprf(isrf)
c      betahav1=betahav1+rgrb1(isrf)*drg
c      volrad1=volrad1+rgrv1(isrf)*drg
c      dkeng1=dkeng1+rgrk1(isrf)*drg

      rgrp=rgrid(isrf)*pprf(isrf)*bsqinv(isrf)
     &    +rgrid(isrf1)*pprf(isrf1)*bsqinv(isrf1)
      rgrb=rgrid(isrf)*ph(isrf,3)*bsqinv(isrf)
     &    +rgrid(isrf1)*ph(isrf1,3)*bsqinv(isrf1)
      rgrv=rgrid(isrf)*volr(isrf)
     &    +rgrid(isrf1)*volr(isrf1)
      rgrk=rgrid(isrf)*delk(isrf)*rhoprf(isrf)
     &    +rgrid(isrf1)*delk(isrf1)*rhoprf(isrf1)
      betaav=betaav+rgrp*drg
      betahav=betahav+rgrb*drg
      volrad=volrad+rgrv*drg
      dkeng=dkeng+rgrk*drg

   41 continue
c      dkeng1=dkeng1*omreal**2*oma2
c      betahav1=betahav1/volrad1

      dkeng=dkeng*omreal**2*oma2
      betahav=betahav/volrad
      betaav=betaav/volrad
c      write(iodat,*) (rgrid(isr),isr=1,nosurf)
c      write(iodat,*) (rgrb1(isr),isr=1,nosurf)
c      write(iodat,*) (rgrv1(isr),isr=1,nosurf)
c      write(iodat,*) (rgrk1(isr),isr=1,nosurf)
c      write(*,40) betahav1,dkeng1
c      write(iodat,40) betahav1,dkeng1
c      write(*,444)
      write(*,40) betaav,betahav,dkeng
c      write(iodat,444)
      write(iodat,40) betaav,betahav,dkeng
      write(17,40) betaav,betahav,dkeng
c  444 format('The <beta-h> is computed by assuming that 
c     >the central alpha pressure is same as the central total pressure')
   40 format('<beta> =',e11.4,' <beta-h> =',e11.4,' <del-K> =',e11.4)
c
cccccccc
ccc...  all frequencies are normalized by Alfven frequency
ccc..  particle diamagnetic drift freq. normalized by Alfven freq.
c
      rgx1=rgrid(nosurf-2)
      rgx2=rgrid(nosurf-1)
      rgx3=rgrid(nosurf)
      dxx=(rgx3-rgx2)/(rgx2-rgx1)
      dxx1=1.+dxx
c     dxx=1.0

      do 23 is=1,iht
      do 202 i=1,nosurf-1
ccc  eha is related to the psi derivative of ehc
      eha(i,is)=ntor*cstrh(is)*ehcp(i,is)/(th0(is))
     &           *(ehc(i,is)/th(i,is))**0.5
ccc  wstarh is related to the psi derivative of ch
      wstarh(i,is)=ntor*(2./3.)*cstrh(is)*chpprf(i,is)
     &             *(th(i,is)/th0(is))

      if(iswitch.ne.1 .and. is.eq.iht) then
      wstarh(i,is)=ntor*cstrh(is)*chpprf(i,is)
     &             *(th(i,is)/th0(is))
ccc  wstarh is related to the psi derivative of Th
      wstarth(i,is)=ntor*cstrh(is)*thp(i,is)/th0(is)
      endif

  202 continue
      eha(nosurf,is)=eha(nosurf-1,is)*dxx1-eha(nosurf-2,is)*dxx
      wstarh(nosurf,is)=wstarh(nosurf-1,is)*dxx1
     &                  -wstarh(nosurf-2,is)*dxx
      wstarth(nosurf,is)=wstarth(nosurf-1,is)*dxx1
     &                  -wstarth(nosurf-2,is)*dxx
   23 continue

ckg following defines omegastar for the subsequent calculation of the
ckg Landau damping on different species. If xLandau=0 then not drive
ckg from thermal ions is included
ckg Fu suggests to have it, so we will.
      xLandau=1
      do 14 is=1,ict
      do 105 i=1,nosurf-1
         wstarc(i,is)=xLandau*ntor*cstrc(is)*dencp(i,is)/max(denc(i,is)
     &        ,1.)*(tc(i,is)/tc0(is))
         wstartc(i,is)=xLandau*ntor*cstrc(is)*tcp(i,is)/tc0(is)
  105 continue
      wstarc(nosurf,is)=wstarc(nosurf-1,is)*dxx1
     &                 -wstarc(nosurf-2,is)*dxx
      wstartc(nosurf,is)=wstartc(nosurf-1,is)*dxx1
     &                 -wstartc(nosurf-2,is)*dxx
   14 continue

c  multiply Coulomb Logarithm
ccc amu=(e-collision freq./om)
      do 103 i=1,nosurf-1
      amu(i)=amu0*alnc(i,1)*(denc(i,1)/tc(i,1))**1.5
  103 continue
      amu(nosurf)=amu(nosurf-1)*dxx1-amu(nosurf-2)*dxx
c      write(16,*)'collisional frequency, Coul.Log.'
c      write(16,*)(i,amu(i),alnc(i,1),i=1,nosurf)
c
ccc  store profile information in "otaep" file
c
      write(*,*) 'Read IPLOT, if IPLOT=1 then plot profiles'
      write(iodat,*) 'Read IPLOT, if IPLOT=1 then plot profiles'
c2d
      iplot=0
c     Note that to plot things you should link ncar libraty and to
c     uncomment actual ploting routings in the  posplt1 subroutine
ckg      iplot=1
      call printprf
      write(iodat,*) 'IPLOT= ',iplot
      if(iplot.eq.1) then
      call plotprf
      else
      end if
c
c
cccccccccccccccccccccccccccccccccccccccccccccccc
cc..  For an equilibrium we choose Te0, Tth0(2), dene0, and their profiles.
cc..  We shall compute the threshold of Alpha particle volume averaged beta
cc..  vs. (Valpha/VAlfven).
cc..  Note that only one value of (Valpha/VAlfven) will coincide with the 
cc..  equilibrium condition.
cc..  If the hot particle beta is much smaller than total beta, the calculation
cc..  is consistent with varying Valpha but keeping VAlfven fixed.
cccccccccccccc
cc   val = valpha / valfven
c
      write(17,99)
      write(iodat,99)
      write(*,99)
   99 format(1x,'**** (val/va)min, (val/va)max, idiv')
c2d      read(5,*) valmin, valmax, idiv
      valmin=1.
      valmax=1.
      idiv=1
      if(idiv.eq.1) go to 791
      delval=(valmax-valmin)/(idiv-1.)
      go to 792
  791 valmin=valphava
      delval=0.0
      valmax=valmin
  792 continue
      write(17,83) valmin,valmax,idiv
      write(iodat,83) valmin,valmax,idiv
      write(*,83) valmin,valmax,idiv
   83 format(1x,'(val/va)min=',e12.4,' (val/va)max=',e12.4,' idiv=',i3)
      do 95 id=1,idiv
      val(id)=valmin+(id-1.)*delval
c      if(idiv.eq.1)val(1)=valphava
cc
cc...   radial integration and read in disk files for every surface
ccc
      do 27 is=1,iht
      gr0ch(is)=0.0
      gr0bh(is)=0.0
      grwht(is)=0.0
      grwhc(is)=0.0
   27 continue

      do 17 is=1,ict
      grwct(is)=0.0
      grwcc(is)=0.0
   17 continue

      cei=0.0
c
cc
cc..  Usually there is a large error in computing parallel electric field through surface component of the displacement near magnetic axis.
cc..  This error seriously distorts the trapped electron collisional contribution. 
cc..  Therefore, the radial integration starts at isrf=4
c
c
      sgnomr=om/abs(om)
      istart=4
      if(ishft.eq.0)istart=2
      do 102 irf=istart,nosurf
      rg=rgrid(irf)
      q=qprf(irf)

      do 128 is=1,iht
      phs(is)=chprf(irf,is)
      ech(is)=ehc(irf,is)/th(irf,is)
  128 continue

      do 118 is=1,ict
      pcs(is)=pc(irf,is)
  118 continue
c
c  collision parameters (simplified version)
      amux=amu(irf)
      zeff1=zeff(irf)
c      epsx=epsg(irf)*amin/rmaj
c
      qn=ntor*qprf(irf)
c
cccc...  read disk quantities
      call readsk(irf)
c
      call cheng(irf,q,val(id),colcoe,cecl,om,amux,zeff1
     &,hmin(irf),hmax(irf),qn)
c2d	
      if(ishft.eq.0)return
ccc  multiply 2*rg due to uniform sqrt(psi) grids in the psi integration
      if(irf.eq.1) rg=0.0

      do 29 is=1,iht
      chit(irf,is)=chtr(is)*rg
      chic(irf,is)=chci(is)*rg
      chi(irf,is)=chit(irf,is)+chic(irf,is)

      gabh(irf,is)=wabh0(is)*rg
      gach(irf,is)=wach0(is)*rg
      gabh00(irf,is)=wabh00(is)*rg
      gabh01(irf,is)=wabh01(is)*rg
      gach0p(irf,is)=wach0p(is)*rg
      gach0m(irf,is)=wach0m(is)*rg
   29 continue

      do 19 is=1,ict
      ccit(irf,is)=cctr(is)*rg
      ccic(irf,is)=ccci(is)*rg
      cci(irf,is)=ccit(irf,is)+ccic(irf,is)
   19 continue

      cecol(irf)=cecl*rg

  102 continue

      ibeg=istart+1
c2d	write(*,*) ' sqrt(psi)   trapp   pass   '
      do 151 irf=ibeg,nosurf
      irf1=irf-1
      drg=rgrid(irf)-rgrid(irf1)
c2d       write(*,*) irf,rgrid(irf),chit(irf,3),chic(irf,3),
c     & cecol(irf),cecol(irf1)
      do 28 is=1,iht
      grwht(is)=grwht(is)+(chit(irf,is)+chit(irf1,is))*drg
      grwhc(is)=grwhc(is)+(chic(irf,is)+chic(irf1,is))*drg

      gr0bh0(is)=gr0bh0(is)+(gabh00(irf,is)+gabh00(irf1,is))*drg
      gr0bh1(is)=gr0bh1(is)+(gabh01(irf,is)+gabh01(irf1,is))*drg
      gr0chp(is)=gr0chp(is)+(gach0p(irf,is)+gach0p(irf1,is))*drg
      gr0chm(is)=gr0chm(is)+(gach0m(irf,is)+gach0m(irf1,is))*drg
      gr0bh(is)=gr0bh(is)+(gabh(irf,is)+gabh(irf1,is))*drg
      gr0ch(is)=gr0ch(is)+(gach(irf,is)+gach(irf1,is))*drg
   28 continue
      do 18 is=1,ict
      grwct(is)=grwct(is)+(ccit(irf,is)+ccit(irf1,is))*drg
      grwcc(is)=grwcc(is)+(ccic(irf,is)+ccic(irf1,is))*drg
   18 continue
cfakekg coll      cei=cei-(cecol(irf)+cecol(irf1))*drg
      if(irf.lt.nosurf*.95) cei=cei-(cecol(irf)+cecol(irf1))*drg
cnng14      print *,'irf,cei=', irf, cei
  151 continue

c
      iplot=0

      if(iswitch.eq.1) then
      nameh(1)='D-beam'
      nameh(2)='T-beam'
      nameh(3)='Alpha'
      else
      nameh(1)='D-beam'
      nameh(2)='D-beam'
      nameh(3)='ICRF ion'
      end if
      do is=1,iht
      iplot=iplot+1
      ipl1=iplot-(iplot/4)*4
c2d
      call posplt1(rgrid,chi(1,is),1,nosurf,ipl1,nameh(is))
c2d 
      if(ipl1.eq.0) call frame(0)
      write(17,*) nameh(is)
      write(17,2) (chi(isrf,is),isrf=1,nosurf)
      end do

      iplot=iplot+1
      ipl1=iplot-(iplot/4)*4
c2d
      buf='gap freq'
      call posplt1(rgrid,falfven,1,nosurf,ipl1,buf)
c2d
      if(ipl1.eq.0) call frame(0)

      iplot=iplot+1
      ipl1=iplot-(iplot/4)*4
c2d
      buf='e-collision'
      call posplt1(rgrid,cecol,1,nosurf,ipl1,buf)
c2d
      if(ipl1.eq.0) call frame(0)
      write(17,*) 'trapped electron collision contrib. profile'
      write(17,2) (cecol(isrf),isrf=1,nosurf)

      namec(1)='electron'
      namec(2)='therm-D'
      namec(3)='therm-T'
      namec(4)='hydrogen'
      namec(5)='impurity'
      do is=1,ict
      iplot=iplot+1
      ipl1=iplot-(iplot/4)*4
c2d
      call posplt1(rgrid,cci(1,is),1,nosurf,ipl1,namec(is))
c2d
      if(ipl1.eq.0) call frame(0)
      write(17,*) namec(is)
      write(17,2) (cci(isrf,is),isrf=1,nosurf)
      end do

c2d
      call frame(0)
    2 format(10e12.3)

ccc..  multiply twopi to account for toroidal angle integration in delta-Wk
cc  a factor of 2 is devided due to absorption of 2 into bounce average of rbce
c
cc (growth (or damping) rate/|om|) for each species

      do 49 is=1,iht
      grwht(is)=pi*grwht(is)/(2.0*dkeng)
      grwhc(is)=pi*grwhc(is)/(2.0*dkeng)
      grwh(is)=(grwht(is)+grwhc(is))

      gr0chp(is)=pi*gr0chp(is)/(2.0*dkeng)
      gr0bh0(is)=pi*gr0bh0(is)/(2.0*dkeng)
      gr0chm(is)=pi*gr0chm(is)/(2.0*dkeng)
      gr0bh1(is)=pi*gr0bh1(is)/(2.0*dkeng)
      gr0ch(is)=pi*gr0ch(is)/(2.0*dkeng)
      gr0bh(is)=pi*gr0bh(is)/(2.0*dkeng)
   49 continue

      do 39 is=1,ict
      grwct(is)=pi*grwct(is)/(2.0*dkeng)
      grwcc(is)=pi*grwcc(is)/(2.0*dkeng)
      grwc(is)=(grwct(is)+grwcc(is))
   39 continue

cc    ecoll = collisional damping rate / TAE frequency
      ecoll=pi*cei/(2.0*dkeng)

ccc   critical alpha particle beta(id)
ccc   (grwa*|om|) is the growth rate due alpha at betahav
c grwbd,grwi,grwe,ecoll are grwth rate due beam, ion, trapped-e etc. at the specified plasma condition
ccc..
      grwth1=ecoll
      do is=1,ict
      grwth1=grwth1+grwc(is)
      end do

      if(iht.gt.1) then
      do is=1,iht-1
      grwth1=grwth1+grwh(is)
      end do
      end if

c
      if(grwh(iht).ne.0.) beta(id)=-grwth1*betahav/grwh(iht)
      if(grwth1.ge.0.0) beta(id)=0.0
      if(grwh(iht).eq.0.) write(*,79) val(id)
   79 format('** For val = ',e12.4,'  there is no net alpha
     &  particle resonance')
ccc..  grwth=growth rate / |real frequency| for betah=2.*beta(id)
      if(betahav.eq.0.)then
         grwth=2.0*grwh(iht)+grwth1
      else
         grwth=2.0*grwh(iht)*beta(id)/betahav+grwth1
      endif
      write(17,94)
      write(iodat,94)
      write(*,94)
   94 format('** vh/vA, <beta-h>_crit, grow/|omr|')
      write(17,6) val(id),beta(id),grwth
      write(iodat,6) val(id),beta(id),grwth
      write(*,6) val(id),beta(id),grwth

      write(17,195)
      write(iodat,195)
      write(*,195)
  195 format('** contri. from e-colli., thermal e, D, T, H, C')
      write(17,61) ecoll,(grwc(is),is=1,ict)
      write(iodat,6) ecoll,(grwc(is),is=1,ict)
      write(*,61) ecoll,(grwc(is),is=1,ict)
    6 format(6e13.4)
 61   format('gms_kg',6e13.4)
ckg write simple output
      write(57,*) 'gam_ecoll'
      write(57,*) ecoll
      write(57,*) 'gam_eLandau'
      write(57,*) grwc(1)
      write(57,*) 'gam_DLandau_ZOW'
      write(57,*) grwc(2)
      write(57,*) 'gam_TLandau_ZOW'
      write(57,*) grwc(3)
      write(57,*) 'gam_HLandau_ZOW'
      write(57,*) grwc(4)
      write(57,*) 'gam_CLandau_ZOW'
      write(57,*) grwc(5)
ckg 
      if(iswitch.eq.1) then
      write(17,194)
      write(iodat,194)
      write(*,194)
  194 format('** contribution from D-beam, T-beam, alpha')
      write(17,6) (grwh(is),is=1,iht)
      write(iodat,6) (grwh(is),is=1,iht)
      write(*,6) (grwh(is),is=1,iht)

      write(17,91)
      write(iodat,91)
      write(*,91)
   91 format('* trap. and cir. hot D-beam, T-beam, A contr.=')
      write(17,9)(grwht(is),grwhc(is),is=1,iht)
      write(iodat,9)(grwht(is),grwhc(is),is=1,iht)
      write(*,9) (grwht(is),grwhc(is),is=1,iht)

      else

      write(17,293)
      write(iodat,293)
      write(*,293)
  293 format('** contribution from D-beam, D-beam, ICRF tail ion')
      write(17,6) (grwh(is),is=1,iht)
      write(iodat,6) (grwh(is),is=1,iht)
      write(*,6) (grwh(is),is=1,iht)

      write(17,291)
      write(iodat,291)
      write(*,291)
  291 format('* trap. and cir. hot D-beam, D-beam, ICRF tail ion')
      write(17,9)(grwht(is),grwhc(is),is=1,iht)
      write(iodat,9)(grwht(is),grwhc(is),is=1,iht)
      write(*,9) (grwht(is),grwhc(is),is=1,iht)
      endif

      write(17,191)
      write(iodat,191)
      write(*,191)
  191 format('* trapped thermal e, D, T, H, C contr.')
      write(17,9)(grwct(is),is=1,ict)
      write(iodat,9)(grwct(is),is=1,ict)
      write(*,9) (grwct(is),is=1,ict)

      write(17,192)
      write(iodat,192)
      write(*,192)
  192 format('* circulating thermal e, D, T, H, C contr.')
      write(17,9)(grwcc(is),is=1,ict)
      write(iodat,9)(grwcc(is),is=1,ict)
      write(*,9) (grwcc(is),is=1,ict)

      write(17,92)
      write(iodat,92)
      write(*,92)
   92 format(1x,'** ((gr0bh(is),gr0ch(is)),is=D, T, A)')
      write(17,9) (gr0bh(is),gr0ch(is),is=1,iht)
      write(iodat,9) (gr0bh(is),gr0ch(is),is=1,iht)
      write(*,9) (gr0bh(is),gr0ch(is),is=1,iht)
    9 format(10e11.3)

      write(17,93)
      write(iodat,93)
      write(*,93)
   93 format(1x,'** ((gr0bh0,gr0bh1,gr0chp,gr0chm),is=D, T, A)')
      write(17,16) (gr0bh0(is),gr0bh1(is),gr0chp(is),gr0chm(is)
     &            ,is=1,iht)
      write(iodat,16)(gr0bh0(is),gr0bh1(is),gr0chp(is),gr0chm(is)
     &            ,is=1,iht)
      write(*,16) (gr0bh0(is),gr0bh1(is),gr0chp(is),gr0chm(is)
     &            ,is=1,iht)
   16 format(4e11.3)

   95 continue
      write(17,51)
      write(iodat,51)
      write(*,51)
   51 format(/,1x,'***  type in 1 to restart ')
c2d      read(5,*) irestart
      irestart=0
   50 format(i2)
      if(irestart.eq.1) go to 333
        return
        end
c
      subroutine cheng(isrf,q,val,colcoe,cecl,om,amux,zeff1,hmin,hmax
     &     ,qn)

      include 'clich1'
      include 'clich1b'

      parameter(ndet=12,nnsrf=501)

      common/pitch2d/gal2d(lam,nn),tkbou2d(lam,nn)
     .              ,galc2d(lamc,nn),tkcir2d(lamc,nn)

      common/pitchang/fpht(lam,nn,isph),dfpht(lam,nn,isph)
     &,fphc(lamc,nn,isph),dfphc(lamc,nn,isph)

      common/cab/abe0,abe01

      common/pitch/
     & wdbh(lam,isph),wbounh(lam,isph),wbouh(lam)
     &,wdch(lamc,isph),wcirch(lamc,isph),wcirh(lamc)
     &,wdbc(lam,ispc),wbounc(lam,ispc)
     &,wdcc(lamc,ispc),wcircc(lamc,ispc)

      common/eqqt1/cengh(lam,mt,jb),cengc(lam,mt,jb)
     &,cengcr(lam,mt,jb,4)
     &,cengch(lamc,mt,jcc),cengcc(lamc,mt,jcc)
     &,chb1(mt,jb),chcp1(mt,jcc),chcm1(mt,jcc)
     &,ahb1(mt,jb),ahcp1(mt,jcc),ahcm1(mt,jcc)
     &,chb2(mt,jb),chcp2(mt,jcc),chcm2(mt,jcc)
     &,ahb2(mt,jb),ahcp2(mt,jcc),ahcm2(mt,jcc)
     &,gh1(lam),gh2(lam),ghc1(lamc),ghc2(lamc)
     &,wabh(lam),wach(lamc),wabh1(lam),wach1(lamc)
     &,wabc(lam),wacc(lamc),ab(mt),abc(mt)

      common/ineq/rbce(lam,mt,jb,3),rbcec(lamc,mt,jcc,3)
     .,wd(lam),tkbou(lam),gal(lam),salm(lam)
     .,wdc(lamc),tkcir(lamc),galc(lamc),salmc(lamc)
     .,tdbou(lam),wds(lam)
c
      common/stabin/ntor,lmin,lmax,mmin,mtot,njg,mjx,nosurf
     &  ,igrid,psilim,ltot(4)
	character tip*1
      common/parm1/np,npc,ipbmax,ipcmin,ipcmax,idet
     &,nnsurf,klamb,klamc,ishft,tip

      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      common/coleig/xlam(ndet,nnsrf),cx(nnsrf,ndet,3)
      dimension wabec(ndet),brg(100),ki(100)
ckg
      common/abe/abpe(lam,mt,nn,jb,6)
c
      nsrf=isrf
cc
c calculate the collisonal damping due to trapped electrons
c
cc  for trapped electron collision
      smx=hmax-hmin
      call enginte0(ecollei,pcs(1),smx)
ckg amux is e - coll frequency / omega
      xmu1=amux*2.0/smx
      call ainte1(zeff1,xmu1,nnsrf,nsrf,ndet,aintx,colcoe,cx(1,1,1)
     &     ,xlam(1,nsrf))
cnng14       print *,'isrf, colcoe',isrf,colcoe,aintx,ecollei
c     this is how it should have look like
c      wabec(i)=ecollei*xlam(i,nsrf)*xmu1*cx(nsrf,i,1)
c     >         /(1.0+(xlam(i,nsrf)*xmu1)**2)
c      do 150 i=1,idet
c      wabec(i)=aintx*cx(nsrf,i,1)
c 150  continue
cc..  trapped electron collisional damping
      i1=1
      cecl=ecollei*aintx
cssum1(idet,wabec,i1)
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc   core particle contributions
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
      np=jb
      npc=jcc
ctune      if(nsrf.eq.4) then
c         write(*,*) "Compare the following omega* with om", om
c      endif
c      write(*,'(i4,4e14.5)') nsrf,(wstartc(nsrf,2)
c     &     +wstarc(nsrf,2))*val,wstarc(nsrf,2)*val
      do 811 is=1,ict
      cdstcval=cdstc(is)*val
      cboucval=cbouc(is)*val
      ccircval=ccirc(is)*val
      wstlc1=wstarc(nsrf,is)*val
      wstltc1=wstartc(nsrf,is)*val

      do 75 k=1,lam
      wdbc(k,is)=cdstcval*wd(k)*(tc(nsrf,is)/tc0(is))
      wbounc(k,is)=0.5*cboucval/tkbou(k)*sqrt(tc(nsrf,is)/tc0(is))
   75 continue
c
      do 74 k=1,lamc
      wdcc(k,is)=cdstcval*wdc(k)*(tc(nsrf,is)/tc0(is))
      wcircc(k,is)=ccircval/tkcir(k)*sqrt(tc(nsrf,is)/tc0(is))
   74 continue
c
      do 77 k2=1,mt
cc  l2 is the poloidal mode number
      sl2=1.0
c      sl2=float(k2-1+minm)/qn
      wstlc=wstlc1*sl2
      wstltc=wstltc1*sl2
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..   trapped particle energy integrations
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc (kp-1) is the bounce harmonic number
c
      do 171 kp=1,np
      if(kp.ne.1) go to 18
ccc....   p=0 bounce integral
      do 161 k=1,lam
      call engintc0(rEz,rp2Ez,zengi,om,pcs(is),wstlc,wdbc(k,is),wstltc)
      cengc(k,k2,kp)=zengi*tkbou(k)
      cengcr(k,k2,kp,1)=rp2Ez*tkbou(k)
      cengcr(k,k2,kp,2)=rEz*tkbou(k)
  161 continue
      go to 171
   18 continue
cccxxxxxxxx..  Nonzero p bounce integrals  ..xxxxxxxxxxxxx
      do 172 k=1,lam
       wbouc=(kp-1.0)*wbounc(k,is)
       call engintc(rEz,rp2Ez,zengc,om,pcs(is),wstlc,wbouc,wstltc)
cccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c  note that multiplying 2.0 to account for +p and -p bounce contributions
cccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       cengc(k,k2,kp)=2.0*zengc*tkbou(k)
       cengcr(k,k2,kp,1)=2.0*rp2Ez*tkbou(k)
       cengcr(k,k2,kp,2)=2.0*rEz*tkbou(k)
  172 continue
  171 continue
c
ccccccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc...  circulating particle energy integrations
cccccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc  (kp-1+jcmin) is the circulating harmonic number
      do 67 kp=1,npc
      xkp=(kp-1.0+jcmin-ntor*q)
      do 66 k=1,lamc
       wcirc=xkp*wcircc(k,is)
ccc...  core particle energy integral
       call engintc(rEz,rp2Ez,zengc,om,pcs(is),wstlc,wcirc,wstltc)
       cengcc(k,k2,kp)=zengc*tkcir(k)
       cengcr(k,k2,kp,3)=rp2Ez*tkcir(k)
       cengcr(k,k2,kp,4)=rEz*tkcir(k)
   66 continue
   67 continue
   77 continue
c
      do 32 k=1,lam
      wabc(k)=0.0
      do 33 kp=1,jb
      if(is.eq.1 .and. kp.eq.1)then
c      bhf=0.0
ckg this accounts for trapped electron damping with zeroth bounce resonance
      bhf=1.0
      else
      bhf=1.0
      end if
      if(is.ne.1)then
         do k1=1,mt
            ab(k1)=0.0
            do 34 j1=1,3
               ab(k1)=ab(k1)+rbce(k,k1,kp,j1)*eigfun(isrf,k1,j1)
 34         continue
c            if(is.eq.1.and.k1.eq.2.and.kp.eq.1.and.isrf.eq.100)
c     +         write(*,*)'new rbce',
c     +         (rbce(k,k1,kp,j1),j1=1,3),tpnt1,tpnt2,eigfun(nsrf,k1,1)
c,ab(k1),abpe(k,k1,isrf,kp,2)
c            if(is.eq.1)ab(k1)=abpe(k,k1,isrf,kp,2)
         enddo
      endif
c(abpe(k,k1,isrf,1)*colcoe+abpe(k,k1,isrf,2)
c+abpe(k,k1,isrf,3))
      do 35 k2=1,mt
      do 35 k1=1,mt
ccc...  sum over bounce harmonics
ccc...  the zeroth bounce term of trapped electron is treated elsewhere
      if(is.ne.1)then
         wabc(k)=wabc(k)+bhf*cengc(k,k2,kp)*ab(k1)*ab(k2)
      else
         wabc(k)=wabc(k)+bhf*(cengc(k,k2,kp)*
     +   abpe(k,k1,isrf,kp,2)*abpe(k,k2,isrf,kp,2)
     +   +cengcr(k,k2,kp,1)*
     +    abpe(k,k1,isrf,kp,1)*abpe(k,k2,isrf,kp,1)*colcoe**2
     +   +cengcr(k,k2,kp,2)*
     +    abpe(k,k1,isrf,kp,2)*abpe(k,k2,isrf,kp,1)*colcoe
     +   +cengcr(k,k2,kp,2)*
     +    abpe(k,k2,isrf,kp,2)*abpe(k,k1,isrf,kp,1)*colcoe)
c         if(k2.eq.2.and.k1.eq.2.and.isrf.eq.100)
c     +         write(*,*)'new rbce',kp,
c     +         cengc(k,k2,kp),cengcr(k,k2,kp,2),cengcr(k,k2,kp,1)
c     +         ,abpe(k,k1,isrf,kp,2),abpe(k,k2,isrf,kp,1),colcoe
      endif
   35 continue
   33 continue
c
c..  multiply the scale and normalization factor in the pitch angle integration
c    salm is the numerical weighting factor in pitch angle integration scheme
      wabc(k)=wabc(k)*salm(k)
   32 continue
c
cc..  final integration in pitch angle variable for trapped particles
      i1=1
      cctr(is)=ssum1(lam,wabc,i1)
c
cc
      do 29 k=1,lamc
      wacc(k)=0.0
      do 28 kp=1,jcc
      do 26 k1=1,mt
      abc(k1)=0.0
      do 26 j1=1,3
      abc(k1)=abc(k1)+rbcec(k,k1,kp,j1)*eigfun(isrf,k1,j1)
   26 continue
      do 25 k2=1,mt
      do 25 k1=1,mt
ccc...  sum over transit harmonics
      if(is.ne.1.or.(is.eq.1.and.kp.gt.jb))then
         wacc(k)=wacc(k)+cengcc(k,k2,kp)*abc(k1)*abc(k2)
      else
ckg this part accounts for the parallel electric field in electron Landau damping
         wacc(k)=wacc(k)+cengcc(k,k2,kp)*
     +   abpe(k,k1,isrf,kp,5)*abpe(k,k2,isrf,kp,5)
     +   +cengcr(k,k2,kp,3)*
     +    abpe(k,k1,isrf,kp,4)*abpe(k,k2,isrf,kp,4)*colcoe**2
     +   +cengcr(k,k2,kp,4)*
     +    abpe(k,k1,isrf,kp,5)*abpe(k,k2,isrf,kp,4)*colcoe
     +   +cengcr(k,k2,kp,4)*
     +    abpe(k,k2,isrf,kp,5)*abpe(k,k1,isrf,kp,4)*colcoe
c         if(k2.eq.2.and.k1.eq.2.and.isrf.eq.100)
c     +         write(*,*)'new rbce',kp,
c     +         cengcc(k,k2,kp),cengcr(k,k2,kp,4),cengcr(k,k2,kp,3)
c     +         ,abpe(k,k1,isrf,kp,5),abpe(k,k2,isrf,kp,4),colcoe
      endif
   25 continue
   28 continue
c
c..  multiply the scale and normalization factor in the pitch angle integration
cc  salmc is the numerical weighting factor in pitch angle integration scheme
      wacc(k)=wacc(k)*salmc(k)
   29 continue
c
      i1=1
      ccci(is)=ssum1(lamc,wacc,i1)
c      cc(is)=cctr(is)+ccci(is)
c
  811 continue
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc   slowing-down hot particle contributions
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      iht1=iht-1
      if(iswitch.eq.1) iht1=iht
      do 821 is=1,iht1
      cdsthval=cdsth(is)*val
      cbouhval=cbouh(is)*val
      ccirhval=ccirh(is)*val
      ehalh1=eha(nsrf,is)*val
      wstlh1=wstarh(nsrf,is)*val
c
      wabh00(is)=0.0
      wabh01(is)=0.0
      wach0p(is)=0.0
      wach0m(is)=0.0

      do 5 k=1,lam
      wdbh(k,is)=cdsthval*wd(k)*(th(nsrf,is)/th0(is))
      wbounh(k,is)=0.5*cbouhval/tkbou(k)*sqrt(th(nsrf,is)/th0(is))
      gh1(k)=fpht(k,nsrf,is)
      gh2(k)=fpht(k,nsrf,is)+(2./3.)*dfpht(k,nsrf,is)
    5 continue
      do 4 k=1,lamc
      wdch(k,is)=cdsthval*wdc(k)*(th(nsrf,is)/th0(is))
      wcirch(k,is)=ccirhval/tkcir(k)*sqrt(th(nsrf,is)/th0(is))

cc  this part needs to be reconsidered!!
      ghc1(k)=fphc(k,nsrf,is)
      ghc2(k)=fphc(k,nsrf,is)+(2./3.)*dfphc(k,nsrf,is)
    4 continue
c
      do 7 k2=1,mt
cc  l2 is the poloidal mode number
c      sl2=float(k2-1+minm)/qn
      sl2=1.0
      ehalh=ehalh1*sl2
      wstlh=wstlh1*sl2
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..   trapped particle energy integrations
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc (kp-1) is the bounce harmonic number
c
      do 11 kp=1,np
      if(kp.ne.1) go to 8
ccc....   p=0 bounce integral
      do 44 k=1,lam
ccc...  hot particle energy integral
      call enginth0(zengh,om,phs(is),wstlh,wdbh(k,is)
     &             ,gh1(k),gh2(k),ech(is),ehalh)
      cengh(k,k2,kp)=zengh*tkbou(k)
   44 continue
ccc. integrate over energy and pitch angle for delta function term due to cutoff in energy distribution
      do 41 k=1,lam
      brg(k)=om-wdbh(k,is)
   41 continue
      kcounb0=0
      do 42 k=1,lam-1
      if((brg(k)*brg(k+1)).lt.0.0) then
      kcounb0=kcounb0+1
      ki(kcounb0)=k
      else
      end if
   42 continue
      if(kcounb0.ge.1) then
      do 43 kc=1,kcounb0
      call enghd0(zengh,om,phs(is),nsrf,k2,kp
     &           ,wdbh(1,is),gal,gh1,rbce,tkbou,ech(is),ki(kc))
      wabh00(is)=wabh00(is)+zengh
   43 continue
      else
      end if
      go to 11
    8 continue
cccxxxxxxxx..  Nonzero p bounce integrals  ..xxxxxxxxxxxxx
      do 52 k=1,lam
      wbouh(k)=(kp-1.0)*wbounh(k,is)
cccc...  trapped hot particle energy integral
      call engintht(zengh,om,phs(is),wstlh,wbouh(k),wdbh(k,is)
     &,gh1(k),gh2(k),ech(is),ehalh)
cccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c  note that multiplying 2.0 to account for +p and -p bounce contributions
cccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      cengh(k,k2,kp)=2.0*zengh*tkbou(k)
   52 continue
cc
cc. integrate over energy and pitch angle for alpha's delta function term due to cutoff in Energy distribution
      do 61 k=1,lam
      brg(k)=(om-wdbh(k,is))**2-wbouh(k)**2
   61 continue
      kcounb1=0
      do 62 k=1,lam-1
      if((brg(k)*brg(k+1)).lt.0.0) then
      kcounb1=kcounb1+1
      ki(kcounb1)=k
      else
      end if
   62 continue

      if(kcounb1.ge.1) then
      do 63 kc=1,kcounb1
      call enghdt(zengh,om,phs(is),nsrf,k2,kp
     &  ,wbouh,wdbh(1,is),gal,gh1,rbce,tkbou,ech(is),ki(kc))
      wabh01(is)=wabh01(is)+2.0*zengh
   63 continue
      else
      end if
   11 continue
c
ccccccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc...  circulating particle energy integrations
cccccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc  (kp-1+jcmin) is the circulating harmonic number
      do 27 kp=1,npc
      xkp=(kp-1.0+jcmin-ntor*q)
      do 21 k=1,lamc
       wcirh(k)=xkp*wcirch(k,is)
ccc...  hot particle energy integral
c      call enginthc(zengh,om,phs(is),wstlh,wcirh,wdch(k,is)
c     &             ,ghc1(k),ghc2(k),ech(is),ehalh)
      call engintbm(zengh,zs1,zs2,zs3,zs4,om,phs(is),wstlh
     &,wcirh(k),wdch(k,is)
     &,ghc1(k),ghc2(k),ihdir(is),ech(is),ehalh)
c
      cengch(k,k2,kp)=zengh*tkcir(k)
   21 continue

cc. integrate over energy and pitch angle for delta function term due to cutoff in energy distribution
ccc  delta(E-E_b) term for nonsingular pitch angle distribution
cc. This contribution vanishes for singular pitch distribution
cc  for ihdir=0 both co- and counter- going hot particles
cc  for ihdir=1 co-going hot particles only
cc  for ihdir=2 counter-going hot particles only
cc for co-going particles
      if(ihdir(is).eq.0 .or. ihdir(is).eq.1) then
      do 71 k=1,lamc
      brg(k)=om-wdch(k,is)-wcirh(k)
   71 continue
      kcouncp=0
      do 72 k=1,lamc-1
      if((brg(k)*brg(k+1)).lt.0.0) then
      kcouncp=kcouncp+1
      ki(kcouncp)=k
      else
      end if
   72 continue
      if(kcouncp.ge.1) then
      do 73 kc=1,kcouncp
      call enghmc(zengb,om,phs(is),nsrf,k2,kp
     &,1.,wcirh,wdch(1,is),galc,ghc1,rbcec,tkcir,ech(is),ki(kc))
      wach0p(is)=wach0p(is)+zengb
   73 continue
      else
      end if
      else
      end if
cc for counter going particles
      if(ihdir(is).eq.0 .or. ihdir(is).eq.2) then
      do 81 k=1,lamc
      brg(k)=om-wdch(k,is)+wcirh(k)
   81 continue
      kcouncm=0
      do 82 k=1,lamc-1
      if((brg(k)*brg(k+1)).lt.0.0) then
      kcouncm=kcouncm+1
      ki(kcouncm)=k
      else
      end if
   82 continue
      if(kcouncm.ge.1) then
      do 83 kc=1,kcouncm
      call enghmc(zengb,om,phs(is),nsrf,k2,kp
     &,-1.,wcirh,wdch(1,is),galc,ghc1,rbcec,tkcir,ech(is),ki(kc))
      wach0m(is)=wach0m(is)+zengb
   83 continue
      else
      end if
      else
      end if
   27 continue
    7 continue
c
      do 2 k=1,lam
      wabh(k)=0.0
      do 31 kp=1,jb
      if(kp.eq.1)then
      bhf=0.0
      else
      bhf=1.0
      end if
      do 3 k1=1,mt
      ab(k1)=0.0
      do 3 j1=1,3
      ab(k1)=ab(k1)+rbce(k,k1,kp,j1)*eigfun(isrf,k1,j1)
    3 continue
      do 12 k2=1,mt
      do 12 k1=1,mt
ccc...  sum over bounce harmonics
      wabh(k)=wabh(k)+cengh(k,k2,kp)*ab(k1)*ab(k2)
   12 continue
   31 continue
c..  multiply the scale and normalization factor in the pitch angle integration
cc  salm is the numerical weighting factor in pitch angle integration scheme
      wabh(k)=wabh(k)*salm(k)
    2 continue
cc..  final integration in pitch angle variable
c
c ihdist=0 for slowing-down beam distribution with non-singular pitch angle
c ihdist=1 for slowing-down beam distribution with single pitch angle: galh0
c 
cc  for nonsingular pitch angle distribution of beam
      if(ihdist(is).eq.0)then
      wabh0(is)=wabh00(is)+wabh01(is)
      i1=1
      chtr(is)=ssum1(lam,wabh,i1)+wabh0(is)
	if(isrf.eq.52.and.is.eq.3) then
ckg        write(*,*) 'trapped'
ckg	write(*,'(2e15.6)') (gal(ik),wabh(ik),ik=1,lam)
        endif
      else
cc  for singular pitch angle distribution of beam
ccc There is no contribution from delta(E-E_b) term in beam ion energy in pitch angle integration
      wabh0(is)=0.0
      chtr(is)=0.0
      end if

cc  salm is the numerical weighting factor in pitch angle integration scheme
      if(ihdist(is).eq.1.and.galh0(is).ge.hmin
     &   .and.galh0(is).le.hmax)then
      do 291 k=1,lam
      wabh1(k)=wabh(k)/salm(k)
  291 continue
      call galgrid(galh0(is),iphc(nsrf,is),gal2d(1,nsrf)
     &            ,wabh1,chtr(is),lam)
c     chtr(is)=wabh(iphc(nsrf,is))/salm(iphc(nsrf,is))
      else
      end if
cc
      do 9 k=1,lamc
      wach(k)=0.0
      do 23 kp=1,jcc
      do 6 k1=1,mt
      abc(k1)=0.0
      do 6 j1=1,3
      abc(k1)=abc(k1)+rbcec(k,k1,kp,j1)*eigfun(isrf,k1,j1)
    6 continue
      do 22 k2=1,mt
      do 22 k1=1,mt
ccc...  sum over transit harmonics
      wach(k)=wach(k)+cengch(k,k2,kp)*abc(k1)*abc(k2)
   22 continue
   23 continue
c
c..  multiply the scale and normalization factor in the pitch angle integration
cc  salmc is the numerical weighting factor in pitch angle integration scheme
      wach(k)=wach(k)*salmc(k)
    9 continue
c
cc..  final integration in pitch angle variable
      chci(is)=0.0

cc  For nonsingular pitch angle distribution 
      if(ihdist(is).eq.0)then
ccc  multiply 2.0 to account for the total number of particle weighting in uni-directional beam
      if (ihdir(is).eq.0) wach0(is)=wach0p(is)+wach0m(is)
      if (ihdir(is).eq.1) wach0(is)=2.*wach0p(is)
      if (ihdir(is).eq.2) wach0(is)=2.*wach0m(is)
      i1=1
      chci(is)=ssum1(lamc,wach,i1)+wach0(is)
	if(isrf.eq.52.and.is.eq.3) then
ckg	write(*,*) 'pasiing'
ckg	write(*,'(2e15.6)') (galc(ik),wach(ik),ik=1,lamc)
        endif
      else
      wach0(is)=0.0
      end if
c
cc  For singular pitch angle distribution there is no contribution from delta(E-E_b) term  
      if(ihdist(is).eq.1.and.galh0(is).lt.hmin.and.galh0(is).ge.0.) then
cc  salmc is the numerical weighting factor in pitch angle integration scheme
      do 292 k=1,lamc
      wach1(k)=wach(k)/salmc(k)
  292 continue
      call galcgrid(galh0(is),iphc(nsrf,is),galc2d(1,nsrf)
     &             ,wach1,chci(is),lamc)
c     chci(is)=wach(iphc(nsrf,is))/salmc(iphc(nsrf,is))
      else
      end if

c      if(((nsrf/20)*20).eq.nsrf) then
c      write(*,*) kcounb0,wabh00(is)
c      write(*,*) kcounb1,wabh01(is)
c      write(*,*) kcouncp,wach0p(is)
c      write(*,*) kcouncm,wach0m(is)
c      else
c      end if

c      ch(is)=chtr(is)+chci(is)
  821 continue
c

      if(iswitch.eq.1 .or. iht.lt.3) return
ccc   for ICRF Maxwellian tail ion contribution
      do 521 is=iht,iht
c
      cdsthval=cdsth(is)*val
      cbouhval=cbouh(is)*val
      ccirhval=ccirh(is)*val
      wstlh1=wstarh(nsrf,is)*val
      wstlth1=wstarth(nsrf,is)*val
c
      do  k=1,lam
      wdbh(k,is)=cdsthval*wd(k)*(th(nsrf,is)/th0(is))
      wbounh(k,is)=0.5*cbouhval/tkbou(k)*sqrt(th(nsrf,is)/th0(is))
      gh1(k)=fpht(k,nsrf,is)
      gh2(k)=dfpht(k,nsrf,is)
      end do
c
      do k=1,lamc
      wdch(k,is)=cdsthval*wdc(k)*(th(nsrf,is)/th0(is))
      wcirch(k,is)=ccirhval/tkcir(k)*sqrt(th(nsrf,is)/th0(is))
cc 
      ghc1(k)=fphc(k,nsrf,is)
      ghc2(k)=dfphc(k,nsrf,is)
      end do
c
      do 97 k2=1,mt
cc  l2 is the poloidal mode number
c      sl2=float(k2-1+minm)/qn
      sl2=1.0
      wstlh=wstlh1*sl2
      wstlth=wstlth1*sl2
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc..   trapped particle energy integrations
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc (kp-1) is the bounce harmonic number
c
      do 111 kp=1,np
      if(kp.ne.1) go to 88

ccc....   p=0 bounce integral
      do k=1,lam
ccc...  hot particle energy integral
      call enginthm0(zengh,om,phs(is),wstlh,wdbh(k,is)
     &             ,gh1(k),gh2(k),wstlth)
      cengh(k,k2,kp)=zengh*tkbou(k)
      end do
c
      go to 111
   88 continue

cccxxxxxxxx..  Nonzero p bounce integrals  ..xxxxxxxxxxxxx
      do k=1,lam
      wbouh(k)=(kp-1.0)*wbounh(k,is)
cccc...  trapped hot particle energy integral
      call enginthmt(zengh,om,phs(is),wstlh,wbouh(k),wdbh(k,is)
     &,gh1(k),gh2(k),wstlth)
cccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c  note that multiplying 2.0 to account for +p and -p bounce contributions
cccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      cengh(k,k2,kp)=2.0*zengh*tkbou(k)
      end do

cc
  111 continue
c
ccccccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ccc...  circulating particle energy integrations
cccccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc  (kp-1+jcmin) is the circulating harmonic number
      do 87 kp=1,npc
      xkp=(kp-1.0+jcmin-ntor*q)

      do k=1,lamc
       wcirh(k)=xkp*wcirch(k,is)
ccc...  hot ICRH particle energy integral
ckg
      call enginthc(zengh,om,phs(is),wstlh,wcirh(k),
     &     wdch(k,is),ghc1(k),ghc2(k),ech(is),ehalh)
c      call engintbm(zengh,zs1,zs2,zs3,zs4,om,phs(is),wstlh
c     &,wcirh(k),wdch(k,is)
c     &,ghc1(k),ghc2(k),ihdir(is),ech(is),ehalh)
c
      cengch(k,k2,kp)=zengh*tkcir(k)
      end do

   87 continue
   97 continue
c
c
      do 92 k=1,lam
      wabh(k)=0.0
      do 91 kp=1,jb

      do k1=1,mt
      ab(k1)=0.0
      do j1=1,3
      ab(k1)=ab(k1)+rbce(k,k1,kp,j1)*eigfun(isrf,k1,j1)
      end do
      end do

      do k2=1,mt
      do k1=1,mt
ccc...  sum over bounce harmonics
      wabh(k)=wabh(k)+cengh(k,k2,kp)*ab(k1)*ab(k2)
      end do
      end do

   91 continue
c..  multiply the scale and normalization factor in the pitch angle integration
cc  salm is the numerical weighting factor in pitch angle integration scheme
      wabh(k)=wabh(k)*salm(k)
   92 continue
cc..  final integration in pitch angle variable
c
c ihdist=0 for distribution with non-singular pitch angle
c ihdist=1 for distribution with singlular pitch angle: galh0
c 
      chtr(is)=0.0

      if(ihdist(is).eq.0)then
      i1=1
      chtr(is)=ssum1(lam,wabh,i1)
      end if

cc  salm is the numerical weighting factor in pitch angle integration scheme
      if(ihdist(is).eq.1.and.galh0(is).ge.hmin
     &   .and.galh0(is).le.hmax)then
      do k=1,lam
      wabh1(k)=wabh(k)/salm(k)
      end do
      call galgrid(galh0(is),iphc(nsrf,is),gal2d(1,nsrf)
     &            ,wabh1,chtr(is),lam)
c     chtr(is)=wabh(iphc(nsrf,is))/salm(iphc(nsrf,is))
      end if

cc
      do 390 il=1,lamc
      wach(il)=0.0
      do 393 lc=1,jcc

      do 391 l5=1,mt
      abc(l5)=0.0
      do 391 l4=1,3
      abc(l5)=abc(l5)+rbcec(il,l5,lc,l4)*eigfun(isrf,l5,l4)
  391 continue

      do 392 l2=1,mt
      do 392 l3=1,mt
ccc...  sum over transit harmonics
      wach(il)=wach(il)+cengch(il,l2,lc)*abc(l3)*abc(l2)
  392 continue

  393 continue
c
c..  multiply the scale and normalization factor in the pitch angle integration
cc  salmc is the numerical weighting factor in pitch angle integration scheme
      wach(il)=wach(il)*salmc(il)
  390 continue

c
cc..  final integration in pitch angle variable
cc  For nonsingular pitch angle distribution 
      chci(is)=0.0

      if(ihdist(is).eq.0)then
      i1=1
      chci(is)=ssum1(lamc,wach,i1)
      end if
c
cc  For singular pitch angle distribution 
      if(ihdist(is).eq.1.and.galh0(is).lt.hmin.and.galh0(is).ge.0.) then
cc  salmc is the numerical weighting factor in pitch angle integration scheme
      do k=1,lamc
      wach1(k)=wach(k)/salmc(k)
      end do

      call galcgrid(galh0(is),iphc(nsrf,is),galc2d(1,nsrf)
     &             ,wach1,chci(is),lamc)
c     chci(is)=wach(iphc(nsrf,is))/salmc(iphc(nsrf,is))
      end if

c      ch(is)=chtr(is)+chci(is)


  521 continue

      return
      end
c
      subroutine readsk(nsrf)
      include 'clich1'
      common / chnel / iodat,iomode,iomap1,ioequ1,mp0,mp1,mp2,iotty
      common/ineq/rbce(lam,mt,jb,3),rbcec(lamc,mt,jcc,3)
     .,wd(lam),tkbou(lam),gal(lam),salm(lam)
     .,wdc(lamc),tkcir(lamc),galc(lamc),salmc(lamc)
     .,tdbou(lam),wds(lam)
	character tip*1
      common/parm1/np,npc,ipbmax,ipcmin,ipcmax,idet
     &,nnsurf,klamb,klamc,ishft,tip
      integer*8 nadres,ierr,length
c
      lgivup=1
      length=(3*jb*mt+4)*lam+(3*jcc*mt+4)*lamc
      nadres=(nsrf-1)*length+10+1
c
c      call zrd(mp2,rbce,length,nadres+ishft,lgivup,ierr)
      call getwa('equou1',rbce,nadres,length,ierr)
ccc.....
c     call zch(ioequ,iff,999)
  999 continue
      return
      end
c
      subroutine galgrid(galb0,ipb,gal,tkbou,tg,lam)
      dimension gal(lam),tkbou(lam)
ccc  assign pitch angle grid for a single pitch angle galb0
c
c  gal(1) => hmin, gal(lam) => hmax.
c
      if(galb0.lt.gal(1)) then
      dtg=(tkbou(2)-tkbou(1))/(gal(2)-gal(1))
      dg=gal(1)-galb0
      tg=tkbou(1)-dg*dtg
      return
      else
      end if
      if(galb0.ge.gal(lam)) then
      dtg=(tkbou(lam)-tkbou(lam-1))/(gal(lam)-gal(lam-1))
      dg=galb0-gal(lam)
      tg=tkbou(lam)+dg*dtg
      return
      else
      do 10 i=1,lam-1
      if(galb0.ge.gal(i).and.galb0.lt.gal(i+1)) go to 20
   10 continue
      i=lam-1
   20 ipb=i+1
      dtg=(tkbou(ipb)-tkbou(i))/(gal(ipb)-gal(i))
      dg=galb0-gal(i)
      tg=tkbou(i)+dg*dtg
      end if
      return
      end
c
      subroutine galcgrid(galb0,ipbdc,galc,tkcir,tg,lamc)
      dimension galc(lamc),tkcir(lamc)
ccc  assign pitch angle grid for a single pitch angle galb0
c
c  galc(1) => hmin, galc(lamc) => 0.
c
      if(galb0.lt.galc(lamc)) then
      dtg=(tkcir(lamc)-tkcir(lamc-1))/(galc(lamc)-galc(lamc-1))
      dg=galc(lamc)-galb0
      tg=tkcir(lamc)-dg*dtg
      return
      else
      end if
      if(galb0.ge.galc(1)) then
      dtg=(tkcir(2)-tkcir(1))/(galc(2)-galc(1))
      dg=galb0-galc(1)
      tg=tkcir(1)+dg*dtg
      return
      else
      do 10 i=1,lamc-1
      if(galb0.lt.galc(i).and.galb0.ge.galc(i+1)) go to 20
   10 continue
      i=lamc-1
   20 ipbdc=i+1
      dtg=(tkcir(ipbdc)-tkcir(i))/(galc(ipbdc)-galc(i))
      dg=galb0-galc(ipbdc)
      tg=tkcir(ipbdc)+dg*dtg
      end if
      return
      end
c
      subroutine engintc0(zs1,rp2Ez,zs3,om,p,wstar,wd,wstart)
ckg p in the arguments is not a harmonic number
ccc...  energy integration for p=0 bounce harmonic; ie. magnetic drift resonance term
ccccc... energy integration for Maxwellian distribution function
      if(om.eq.0..or.wd.eq.0.0) go to 1
      z2=om/wd
      if(z2.le.0.0) go to 1
      wstrwd=wstar/wd
      wstrtwd=wstart/wd
      ce0=z2-(wstrwd+wstrtwd*(z2-1.5))
      arg=sqrt(z2)
      if(z2.ge.100.) zs2=0.0
      if(z2.lt.100.) zs2=exp(-z2)
      rp2Ez=-ce0*arg*1.772453851*zs2*p
      zs1=z2*rp2Ez
      zs3=zs1*z2
      return
    1 continue
      rp2Ez=0.
      zs3=0.
      zs1=0.
      return
      end
c
      subroutine enginte0(enge,p,sm)
      enge=2.0*sm*p/sqrt(3.14159)
      return
      end
ckg here xl is e-coll /omega times 2 * eigenvalue 
      subroutine ainte1(zeff,xl,nnsrf,nsrf,ndet,aintx,colcoe,cx,xlam)
      dimension cx(nnsrf,ndet,3),xlam(ndet)
      xend=80.
      xend=20.
      nx=800
      h=xend/float(nx)
      aint5=0.0
c      aint3=0.
      pi2=sqrt(3.1415926)
      xold=0.0
c      mm=10
      do 100 i=1,nx
      x=0.5*h+h*float(i-1)
      xx=(xl*x)**(1./3.)
      xx2=xx**2
      aintx=exp(-xx2)
c      hs=(xx-xold)/float(mm)
c      do 300 j=1,mm
c      y=xold+hs*float(j)-0.5*hs
c      aint3=aint3+hs*exp(-y**2)
c 300  continue
c      xold=xx
      if(xx.le.0.5)then
         hz=4.0*xx/(3.0*pi2)
      else
         hz=aintx/(xx*pi2)+(1.0-0.5/xx2)*erf1(xx)
      end if
      aintx=aintx*x
      hz=zeff+hz
      do idet=1,ndet
         hzz=xlam(idet)*hz
         aint5=aint5+hzz/(hzz**2+x**2)*aintx*(cx(nsrf,idet,1)*colcoe
     &        +cx(nsrf,idet,2)*xx2+cx(nsrf,idet,3)
     &        )**2
      enddo
 100  continue
      aintx=h*2.0*aint5/3.0*xl
c      write(6,*)xend,xl,aint5
c 200  continue
      return
      end
c
      subroutine enginthm0(zs3,om,p,wstar,wd,g1,g2,wstart)
      parameter (pie=3.14159265)
ccc...  energy integration for p=0 bounce harmonic; ie. magnetic drift resonance term
ccccc  for Maxwellian particle distribution with nonsingular pitch angle distribution
      if(wd.eq.0.0) go to 25
      z1=om/wd
      if(z1.le.0.) go to 25
      z1c=z1**1.5
      f01=exp(-z1)
      f01s=z1c*f01
      gz1=g1*z1
      zs1=om*(gz1-g2)-gz1*(wstar+z1*wstart)
      zs1=f01s/abs(wd)*zs1
      zs3=sqrt(pie)*zs1*p
      return
   25 zs3=0.0
      return
      end
c
c
      subroutine enginthmt(zs3,om,p,wstar,wo,wd,g1,g2,wstart)
ccc  trapped particle energy integration for nonzero p bounce terms
cccccccc  for Maxwellian energy distribution
cc  wstart relates to the radial derivative of temperature
ccccccc
      parameter (pie=3.14159265)
      wd2=wd**2
      wo2=wo**2
      ce0=2.0*om*wd
      clo=(wo2+ce0)
      zpl=clo**2-ce0**2
      if(zpl.le.0.0) go to 26
      arg=sqrt(zpl)
      z1=(clo+arg)/(2.0*wd2)
      z2=(clo-arg)/(2.0*wd2)
      zs1=0.0
      zs2=0.0
      if(z1.le.0.0) go to 22
      z1c=z1**1.5
      f01=exp(-z1)
      f01s=z1c*f01
      gz1=g1*z1
      zs1=om*(gz1-g2)-gz1*(wstar+z1*wstart)
      zs1=f01s*(om-wd*z1)*zs1
   22 continue
      if(z2.le.0.0) go to 23
      z2c=z2**1.5
      f02=exp(-z2)
      f02s=z2c*f02
      gz2=g1*z2
      zs2=om*(gz2-g2)-gz2*(wstar+z2*wstart)
      zs2=f02s*(om-wd*z2)*zs2
   23 continue
      zs3=-sqrt(pie)*p*(zs1+zs2)/arg
      return
   26 zs3=0.0
      return
      end
c
      subroutine enginthmc(zs,om,p,wstar,wo,wd,g1,g2,wstart)
ccc  energy integration for nonzero p bounce or all p circulating terms
cccccccc  for Maxwellian energy particle distribution
ccccccc
cc  for both the circulating co- and counter- going particles
      parameter (pie=3.14159265)

      zsp=0.0
      zsm=0.0
      zpl=wo**2+4.0*om*wd
      if(zpl.le.0.0) go to 99
      arg=sqrt(zpl)
      z1=(-wo+arg)/(2.0*wd)
      z2=(-wo-arg)/(2.0*wd)
cc  for co-going particle
      zs1=0.0
      zs2=0.0
      if(z1.le.0.0) go to 12
      z1s=z1*z1
      z1c=z1s*z1
      z1q=z1c*z1
      f01=exp(-z1s)
      f01s=z1q*f01
      gz1=g1*z1s
      zs1=om*(gz1-g2)-gz1*(wstar+z1s*wstart)
      zs1=f01s*(om-wd*z1s)*zs1
   12 continue
      if(z2.le.0.0) go to 13
      z2s=z2*z2
      z2c=z2s*z2
      z2q=z2c*z2
      f02=exp(-z2s)
      f02s=z2q*f02
      gz2=g1*z2s
      zs2=om*(gz2-g2)-gz2*(wstar+z2s*wstart)
      zs2=f02s*(om-wd*z2s)*zs2
   13 continue
      zsp=(zs2+zs1)/arg

cc  for counter-going particle
      z3=-z1
      z4=-z2
      zs3=0.0
      zs4=0.0
      if(z3.le.0.0) go to 22
      z3s=z3*z3
      z3c=z3s*z3
      z3q=z3c*z3
      f03=exp(-z3s)
      f03s=z3q*f03
      gz3=g1*z3s
      zs3=om*(gz3-g2)-gz3*(wstar+z3s*wstart)
      zs3=f03s*(om-wd*z3s)*zs3
   22 continue
      if(z4.le.0.0) go to 23
      z4s=z4*z4
      z4c=z4s*z4
      z4q=z4c*z4
      f04=exp(-z4s)
      f04s=z4q*f04
      gz4=g1*z4s
      zs4=om*(gz4-g2)-gz4*(wstar+z4s*wstart)
      zs4=f04s*(om-wd*z4s)*zs4
   23 continue
      zsm=(zs3+zs4)/arg

   99 continue
c
      zs=sqrt(pie)*(zsp+zsm)*p
c
      return
      end
c
c
c
      subroutine enginth0(zs3,om,p,wstar,wd,g1,g2,ec,ea)
      parameter (pie=3.14159265)
ccc...  energy integration for p=0 bounce harmonic; ie. magnetic drift resonance term
ccccc  for slowing-down particle distribution with nonsingular pitch angle distribution
      if(wd.eq.0.0) go to 25
      xc3=ec**1.5
      eb=om*xc3
      eaa=ea

c      xc3=0.00
c      eb=0.00
c      eaa=0.00

      z1=om/wd
      if(z1.le.0. .or. z1.gt.1.0) go to 25
      wstwd=wstar/wd
      z1c=z1**1.5
      f01=1.0/(z1c+xc3)
      f01s=z1c*f01
      zs1=f01s*(om*(g2-wstwd*g1)+(eaa*z1-eb)*g1*f01)/abs(wd)
      zs3=1.125*pie*zs1*p
      return
   25 zs3=0.0
      return
      end
c
      subroutine enghd0(zs,om,p,isrf,k2,kp,wd,al,g1,rb,tk,ec,kin)
c
      include 'clich1'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      parameter (pie=3.14159265)
      dimension wd(lam),al(lam),g1(lam),rb(lam,mt,jb,3),tk(lam)
c
ccc  magnetic drift resonance contribution
ccc  delta function contribution from slowing-down alpha particles with cutoff energy cutoff
ccc..  integrate over pitch angle space
ccc..  look for pitch angle al0 that satisfies owr=wd resonance
c
      k=kin
      xc3=ec**1.5
c      xc3=0.00
      f0=1.0/(1.0+xc3)

      dal=al(k+1)-al(k)
      dwd=(wd(k+1)-wd(k))/dal
      al0=al(k)+(om-wd(k))/dwd
      gtk1=g1(k+1)*tk(k+1)
      gtk0=g1(k)*tk(k)
      dtk=(gtk1-gtk0)/dal
      gtk=gtk0+(al0-al(k))*dtk
      gbar=(2./3.)*f0*om
cc  this term is always contributing to damping
      zs=1.125*pie*p*gbar*gtk/abs(dwd)
c
      act1=0.0
      do 105 km=1,mt
      ac=0.
      do 106 j1=1,3
      rb0=rb(k,km,kp,j1)
      rb1=rb(k+1,km,kp,j1)
      rbc0=rb0+(al0-al(k))*(rb1-rb0)/dal
      ac=ac+rbc0*eigfun(isrf,km,j1)
  106 continue
      if(km.eq.k2) zc1=ac
      act1=act1+ac
  105 continue
      zs=zs*zc1*act1
c
      return
      end
c
      subroutine enghdt(zs,om,p,isrf,k2,kp,
     &wo,wd,al,g1,rb,tk,ec,kin)
      parameter (pie=3.14159265)
c
ccc  bounce resonance contribution
ccc  delta function contribution from slowing-down alpha particles with cutoff energy cutoff
cc  (w**2-wd**2)-wo**2=0 resonance
ccc  integrate over pitch angle
c 
      include 'clich1'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
cc  assume lam = lamc only
      dimension wo(lam),wd(lam),al(lam),g1(lam),rb(lam,mt,jb,3)
     &         ,tk(lam),wz1(lam),wz2(lam),arg(lam),brg(lam)
ccccccc
      zs1=0.0
      zs2=0.0
      zs=0.0

      xc3=ec**1.5
c      xc3=0.00
      f0=1.0/(1.0+xc3)

      k0=kin
      k1=kin+1
      do 11 k=k0,k1
      wd2=wd(k)**2
      wo2=wo(k)**2
      ce0=2.0*om*wd(k)
      clo=(wo2+ce0)
      zpl=clo**2-ce0**2
      if(zpl.lt.0.0) return
      arg(k)=sqrt(zpl)
      wz1(k)=(clo+arg(k))/(2.0*wd2)
      wz2(k)=(clo-arg(k))/(2.0*wd2)
   11 continue

      if(((1.0-wz1(k0))*(1.0-wz1(k1))).ge.0.0) go to 21
      dal=al(k1)-al(k0)
      dwz1=(wz1(k1)-wz1(k0))/dal
      al1=al(k0)+(1.-wz1(k0))/dwz1
      gtk1=g1(k1)*tk(k1)*(om-wd(k1))/arg(k1)
      gtk0=g1(k0)*tk(k0)*(om-wd(k0))/arg(k0)
      dgtk=(gtk1-gtk0)/dal
      gtk=gtk0+(al1-al(k0))*dgtk

      gbar=(2./3.)*f0*om
      zs1=1.125*pie*p*gbar*gtk/abs(dwz1)

c
      act1=0.0
      do 105 km=1,mt
      ac1=0.0
      do 106 j1=1,3
      rb0=rb(k0,km,kp,j1)
      rb1=rb(k1,km,kp,j1)
      rbc0=rb0+(al1-al(k0))*(rb1-rb0)/dal
      ac1=ac1+rbc0*eigfun(isrf,km,j1)
  106 continue
      if(km.eq.k2) zc1=ac1
      act1=act1+ac1
  105 continue
      zs1=zs1*zc1*act1
   21 continue

      if(((1.0-wz2(k0))*(1.0-wz2(k1))).ge.0.0) go to 22
      dal=al(k1)-al(k0)
      dwz2=(wz2(k1)-wz2(k0))/dal
      al2=al(k0)+(1.-wz2(k0))/dwz2
      gtk1=g1(k1)*tk(k1)*(om-wd(k1))/arg(k1)
      gtk0=g1(k0)*tk(k0)*(om-wd(k0))/arg(k0)
      dgtk=(gtk1-gtk0)/dal
      gtk=gtk0+(al2-al(k0))*dgtk

      gbar=(2./3.)*f0*om
      zs2=1.125*pie*p*gbar*gtk/abs(dwz2)
c
      act2=0.0
      do 108 km=1,mt
      ac2=0.0
      do 107 j1=1,3
      rb0=rb(k0,km,kp,j1)
      rb1=rb(k1,km,kp,j1)
      rbc0=rb0+(al2-al(k0))*(rb1-rb0)/dal
      ac2=ac2+rbc0*eigfun(isrf,km,j1)
  107 continue
      if(km.eq.k2) zc2=ac2
      act2=act2+ac2
  108 continue
      zs2=zs2*zc2*act2
   22 continue
c
      zs=-(zs1+zs2)
      return
      end
c
      subroutine enghmc(zs,om,p,isrf,k2,kp,
     &sigma,wo,wd,al,g1,rb,tk,ec,kin)
c
ccc  transit resonance contribution
ccc  delta(E-E_b) function contribution from slowing-down alpha particles
cc   for uni-directional circulating beam (w-wd)-sigma*wo=0 resonance
cc   integrate over pitch angle space
c 
cc  ec is the ratio of the lower critical energy to beam birth energy
c
      include 'clich1'
      parameter (pie=3.14159265)
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      dimension wo(lamc),wd(lamc),al(lamc),g1(lamc)
     &,rb(lamc,mt,jcc,3),tk(lamc)
      dimension wz1(lamc),wz2(lamc),arg(lamc),brg(lamc)
ccccccc
      zs1=0.0
      zs2=0.0
      zs=0.0

      xc3=ec**1.5
c      xc3=0.00
      f0=1.0/(1.0+xc3)

      k0=kin
      k1=k0+1

      do 11 k=k0,k1
      wo2=wo(k)**2
      ce0=4.0*om*wd(k)
      zpl=wo2+ce0
      if(zpl.le.0.0) return
      arg(k)=sqrt(zpl)
      wz1(k)=(-sigma*wo(k)+arg(k))/(2.0*wd(k))
      wz2(k)=(-sigma*wo(k)-arg(k))/(2.0*wd(k))
   11 continue

      if(((1.0-wz1(k0))*(1.0-wz1(k1))).ge.0.0) go to 21
      dal=al(k1)-al(k0)
      dwz1=(wz1(k1)-wz1(k0))/dal
      al1=al(k0)+(1.-wz1(k0))/dwz1
      gtk1=g1(k1)*tk(k1)/arg(k1)
      gtk0=g1(k0)*tk(k0)/arg(k0)
      dgtk=(gtk1-gtk0)/dal
      gtk=gtk0+(al1-al(k0))*dgtk

      gbar=(2./3.)*f0*om
      zs1=0.5625*pie*p*gbar*gtk/abs(dwz1)
c
      act1=0.0
      do 105 km=1,mt
      ac1=0.0
      do 106 j1=1,3
      rb0=rb(k0,km,kp,j1)
      rb1=rb(k1,km,kp,j1)
      rbc0=rb0+(al1-al(k0))*(rb1-rb0)/dal
      ac1=ac1+rbc0*eigfun(isrf,km,j1)
  106 continue
      if(km.eq.k2) zc1=ac1
      act1=act1+ac1
  105 continue
      zs1=zs1*zc1*act1
c
   21 continue
      if(((1.0-wz2(k0))*(1.0-wz2(k1))).ge.0.0) go to 22
      dal=al(k1)-al(k0)
      dwz2=(wz2(k1)-wz2(k0))/dal
      al2=al(k0)+(1.-wz2(k0))/dwz2
      gtk1=g1(k1)*tk(k1)/arg(k1)
      gtk0=g1(k0)*tk(k0)/arg(k0)
      dgtk=(gtk1-gtk0)/dal
      gtk=gtk0+(al2-al(k0))*dgtk

      gbar=(2./3.)*f0*om
      zs2=0.5625*pie*p*gbar*gtk/abs(dwz2)
c
      act2=0.0
      do 108 km=1,mt
      ac2=0.0
      do 107 j1=1,3
      rb0=rb(k0,km,kp,j1)
      rb1=rb(k1,km,kp,j1)
      rbc0=rb0+(al2-al(k0))*(rb1-rb0)/dal
      ac2=ac2+rbc0*eigfun(isrf,km,j1)
  107 continue
      if(km.eq.k2) zc2=ac2
      act2=act2+ac2
  108 continue
      zs2=zs2*zc2*act2
   22 continue
c
      zs=zs1+zs2
      return
      end
c
      subroutine engintc(zs1,rp2Ez,zs3,om,p,wstar,wo,wstart)
ckg p in the arguments is not a harmonic number
ccc  energy integration for nonzero p bounce or all p circulating terms
ccc..   magnetic drift frequencies are neglected
ccccc... energy integration for maxwellian distribution function
      if(om.eq.0.0.or.wo.eq.0.0) go to 1
      z2=(om/wo)**2
      zstar2=om*wstar/wo**2
      zstart2=om*wstart/wo**2
      ce0=z2-(zstar2+zstart2*(z2-1.5))
      arg=sqrt(z2)
      if(z2.ge.100.) zs2=0.0
      if(z2.lt.100.) zs2=exp(-z2)
      rp2Ez=-ce0*arg*1.772453851*zs2*p
      zs1=z2*rp2Ez
      zs3=z2*zs1
      return
    1 continue
      rp2Ez=0.
      zs3=0.
      zs1=0.
      return
      end
c
      subroutine engintht(zs3,om,p,wstar,wo,wd,g1,g2,ec,ea)
ccc  trapped Alpha energy integration for nonzero p bounce terms
cccccccc  for slowing-down particle distribution
cc  ea relates to the radial derivative of lower critical Alpha energy
cc  ec is the ratio of the lower critical energy to Alpha birth energy
ccccccc
      parameter (pie=3.14159265)
      xc3=ec**1.5
      eb=om*xc3
      eaa=ea

c      xc3=0.00
c      eb=0.00
c      eaa=0.00

      wd2=wd**2
      wo2=wo**2
      ce0=2.0*om*wd
      clo=(wo2+ce0)
      zpl=clo**2-ce0**2
      if(zpl.le.0.0) go to 26
      arg=sqrt(zpl)
      z1=(clo+arg)/(2.0*wd2)
      z2=(clo-arg)/(2.0*wd2)
      zs1=0.0
      zs2=0.0
      if(z1.le.0.0.or.z1.gt.1.0) go to 22
      z1c=z1**1.5
      f01=1.0/(z1c+xc3)
      f01s=z1c*f01
      zs1=f01s*(om-wd*z1)*(om*g2-wstar*z1*g1+(eaa*z1-eb)*g1*f01)
   22 continue
      if(z2.le.0..or.z2.gt.1.0) go to 23
      z2c=z2**1.5
      f02=1.0/(z2c+xc3)
      f02s=z2c*f02
      zs2=f02s*(om-wd*z2)*(om*g2-wstar*z2*g1+(eaa*z2-eb)*g1*f02)
   23 continue
      zs3=-1.125*pie*p*(zs1+zs2)/arg
      return
   26 zs3=0.0
      return
      end
c
      subroutine enginthc(zs,om,p,wstar,wo,wd,g1,g2,ech,ea)
ccc  energy integration for nonzero p bounce or all p circulating terms
cccccccc  for slowing-down Alpha particle distribution
ccccccc
cc  for both the circulating co- and counter- going Alpha particles
cc  ea relates to the radial derivative of lower critical Alpha energy
cc  ech is the ratio of the lower critical energy to Alpha birth energy
cc  velocity integration
      parameter (pie=3.14159265)
      xc3=ech**1.5
      eb=om*xc3
      eaa=ea

c      xc3=0.00
c      eaa=0.00
c      eb=0.00

      zsp=0.0
      zsm=0.0
      zpl=wo**2+4.0*om*wd
      if(zpl.le.0.0) go to 99
      arg=sqrt(zpl)
      z1=(-wo+arg)/(2.0*wd)
      z2=(-wo-arg)/(2.0*wd)
cc  for co-going Alpha
      zs1=0.0
      zs2=0.0
      if(z1.le.0.0 .or. z1.gt.1.0) go to 12
      z1s=z1*z1
      z1c=z1s*z1
      z1q=z1c*z1
      f01=1.0/(z1c+xc3)
      zs1=z1q*(om*g2-wstar*z1s*g1+(eaa*z1s-eb)*g1*f01)*f01
   12 continue
      if(z2.le.0.0 .or. z2.gt.1.0) go to 13
      z2s=z2*z2
      z2c=z2s*z2
      z2q=z2c*z2
      f02=1.0/(z2c+xc3)
      zs2=z2q*(om*g2-wstar*z2s*g1+(eaa*z2s-eb)*g1*f02)*f02
   13 continue
      zsp=(zs2+zs1)/arg

cc  for counter-going Alpha
      z3=-z1
      z4=-z2
      zs3=0.0
      zs4=0.0
      if(z3.le.0.0 .or. z3.gt.1.0) go to 22
      z3s=z3*z3
      z3c=z3s*z3
      z3q=z3c*z3
      f03=1.0/(z3c+xc3)
      zs3=z3q*(om*g2-wstar*z3s*g1+(eaa*z3s-eb)*g1*f03)*f03
   22 continue
      if(z4.le.0.0 .or. z4.gt.1.0) go to 23
      z4s=z4*z4
      z4c=z4s*z4
      z4q=z4c*z4
      f04=1.0/(z4c+xc3)
      zs4=z4q*(om*g2-wstar*z4s*g1+(eaa*z4s-eb)*g1*f04)*f04
   23 continue
      zsm=(zs3+zs4)/arg

   99 continue
c
      zs=1.125*pie*(zsp+zsm)*p
c
      return
      end
c
c
      subroutine engintbm(zs,zs1,zs2,zs3,zs4,om,p,wstar
     >,wo,wd,g1,g2,ihdir,ech,eba)
ccc  energy integration for nonzero p bounce or all p circulating terms
cccccccc  for slowing-down hot particle distribution
ccccccc
cc  for co-going particle (circulating particle only)
cc  for ihdir=0 both co- and counter- going particles
cc  for ihdir=1 co-going particle only
cc  for ihdir=2 counter-going particle only
cc  eba relates to the radial derivative of lower critical particle energy
cc  ech is the ratio of the lower critical energy to particle birth energy
cc  velocity integration
      parameter (pie=3.14159265)
      xc3=ech**1.5
      ebb=om*xc3
c      xc3=0.00
c      eba=0.00
c      ebb=0.00
      zsp=0.0
      zsm=0.0
      zpl=wo**2+4.0*om*wd
      if(zpl.le.0.0) go to 99
      arg=sqrt(zpl)
      z1=(-wo+arg)/(2.0*wd)
      z2=(-wo-arg)/(2.0*wd)
cc  for co-going particles
      if(ihdir.eq.0 .or. ihdir.eq.1)then
      zs1=0.0
      zs2=0.0
      if(z1.le.0.0 .or. z1.gt.1.0) go to 12
      z1s=z1*z1
      z1c=z1s*z1
      z1q=z1c*z1
      f01=1.0/(z1c+xc3)
      zs1=z1q*(om*g2-wstar*z1s*g1+(eba*z1s-ebb)*g1*f01)*f01
   12 continue
      if(z2.le.0.0 .or. z2.gt.1.0) go to 13
      z2s=z2*z2
      z2c=z2s*z2
      z2q=z2c*z2
      f02=1.0/(z2c+xc3)
      zs2=z2q*(om*g2-wstar*z2s*g1+(eba*z2s-ebb)*g1*f02)*f02
   13 continue
      zsp=(zs2+zs1)/arg
      else
      end if

cc  for counter-going particles
      if(ihdir.eq.0 .or. ihdir.eq.2)then
      z3=-z1
      z4=-z2
      zs3=0.0
      zs4=0.0
      if(z3.le.0.0 .or. z3.gt.1.0) go to 22
      z3s=z3*z3
      z3c=z3s*z3
      z3q=z3c*z3
      f03=1.0/(z3c+xc3)
      zs3=z3q*(om*g2-wstar*z3s*g1+(eba*z3s-ebb)*g1*f03)*f03
   22 continue
      if(z4.le.0.0 .or. z4.gt.1.0) go to 23
      z4s=z4*z4
      z4c=z4s*z4
      z4q=z4c*z4
      f04=1.0/(z4c+xc3)
      zs4=z4q*(om*g2-wstar*z4s*g1+(eba*z4s-ebb)*g1*f04)*f04
   23 continue
      zsm=(zs3+zs4)/arg
      else
      end if

   99 continue
c
      if(ihdir.eq.0) zs=1.125*pie*(zsp+zsm)*p
ccc  multiply 2.0 to account for the total number of particle weighting in one-directional going particles
      if(ihdir.eq.1) zs=2.25*pie*zsp*p
      if(ihdir.eq.2) zs=2.25*pie*zsm*p
c
      return
      end
c
      subroutine zeta(z,zetaoz)
      complex z,zetaoz,dzetaz,term,fmult,terme,an1,bn1,zsquar,
     1  hold,temp1,temp2,ddzeta,dddzet,cexp
      real imagte,imagmu,imagse,imagsu
      erro=1.e-15
      zsquar=z*z
      x=real(z)
      y=aimag(z)
      fn=real(zsquar)
      if (y.gt.0.) go to 99
      if (abs(fn).lt.174..and.abs(aimag(zsquar)).lt.5.e4) go to 98
      if (fn.gt.0.) go to 97
      write(16,1) z
1     format (" argument wp of subroutine zeta has too large a negative
     1imaginary part, wp = ",1pe14.7," + ",e14.7," i")
97    hold=(0.,0.)
      go to 99
98    hold=(0.,1.77245385090551603)*cexp(-zsquar)
99    if (x*x+y*y.gt.16.) go to 200
      if (abs(y).ge.1.) go to 300
      realte=-2.*x
      imagte=-2.*y
      realmu=.5*(imagte*imagte-realte*realte)
      imagmu=-imagte*realte
      realsu=realte
      imagsu=imagte
      if (x.eq.0..and.y.eq.0.) go to 103
      fn=3.
100   realse=realte
      imagse=imagte
      realte=(realse*realmu-imagse*imagmu)/fn
      imagte=(realse*imagmu+imagse*realmu)/fn
      realse=realsu
      imagse=imagsu
      realsu=realsu+realte
      imagsu=imagsu+imagte
      fn=fn+2.
      if (abs(realse-realsu).gt.erro.or.abs(imagse-imagsu).gt.
     1  erro) go to 100
103   x=realsu
      fn=imagsu
      if (y.gt.0.)hold=(0.,1.77245385090551603)*cexp(-zsquar)
      zetaoz=cmplx(x,fn)+hold
      go to 401
200   fn=5.
      dddzet=6.
      term=dddzet
      fmult=.5/zsquar
201   terme=term
      term=term*fmult*fn*(fn-1.)/(fn-3.)
      zetaoz=term/terme
      if (abs(real(zetaoz))+abs(aimag(zetaoz)).gt.1.) go to 250
      zetaoz=dddzet
      dddzet=dddzet+term
      fn=fn+2.
      if (cabs(zetaoz-dddzet).gt.erro) go to 201
250   dddzet=dddzet/(zsquar*zsquar)
      if (y.gt.0.) go to 260
      fn=1.
      if (y.lt.0.) fn=2.
      dddzet=dddzet-4.*fn*hold*z*(2.*zsquar-3.)
260   ddzeta=-(4.+(zsquar-.5)*dddzet)/(z*(2.*zsquar-3.))
      dzetaz=(2.-z*ddzeta)/(2.*zsquar-1.)
      zetaoz=-(1.+.5*dzetaz)/z
      go to 401
300   if (y.lt.0.) z=conjg(z)
      terme=(1.,0.)
      term=(0.,0.)
      dzetaz=term
      fmult=terme
      n=0
      an1=z
      bn1=-z*z+.5
301   temp1=bn1*term+an1*terme
      temp2=bn1*fmult+an1*dzetaz
      zetaoz=temp1/temp2
      dzetaz=(zetaoz-term/fmult)/zetaoz
      if (abs(real(dzetaz)).lt.erro.and.abs(aimag(dzetaz)).lt.erro)
     1  go to 302
      bn1=bn1+2.
      n=n+1
      an1=-.5*float(n*(n+n-1))
      terme=term
      dzetaz=fmult
      term=temp1
      fmult=temp2
      if (n.lt.30) go to 301
302   if (y.ge.0.) go to 401
      zetaoz=conjg(zetaoz)+2.*hold
401   return
      end
c
      subroutine plotr
      include 'clich1'
      include 'clich2'
c
c      call frame(0)
c      call posplt1(rgrid,rshear,1,nosurf,1,'shear')
c
c      call posplt1(rgrid,rcurvn,1,nosurf,2,'curvn')
c
c      call posplt1(rgrid,rsigma,1,nosurf,3,'sigma')
c
c      call posplt1(rgrid,rgrgrpp,1,nosurf,4,'grgrpp')
c     read(59,101) idummy
c
c      call frame(0)
c      call posplt1(rgrid,rxjprxj,1,nosurf,1,'xjprxj')
c
c      call posplt1(rgrid,rxjsq,1,nosurf,2,'xjsq')
c
c      call posplt1(rgrid,rxsdpxs,1,nosurf,3,'xsdpxs')
c
c      call posplt1(rgrid,rgptdtgp,1,nosurf,4,'gptdtgp')
c     read(59,101) idummy
  101 format(i2)
c
c      call frame(0)
c      call posplt1(rgrid,rxjdtxj,1,nosurf,1,'xjdtxj')
c
c      call posplt1(rgrid,rcurvs,1,nosurf,2,'curvs')
c
c      call posplt1(rgrid,rgptgpp,1,nosurf,3,'gptgpp')
c
c      call posplt1(rgrid,rgpdtgp,1,nosurf,4,'gpdtgp')
c
c      call frame(0)
c      call posplt1(rgrid,rqp,1,nosurf,1,'qp')
c      call posplt1(rgrid,agptdtgp,1,nosurf,2,'agptdtgp')
c      call posplt1(rgrid,agrgrpp,1,nosurf,3,'agrgrpp')
c      call posplt1(rgrid,agptdth,1,nosurf,4,'agptdth')
ccc  plot trapped particle pitch angle domain
c      call frame(0)
c      call posplt3(rgrid,hmin,hmax,1,nosurf,1,'hmin&hmax')
c      call posplt1(rgrid,wstara,1,nosurf,2,'wstara')
      return
      end
      subroutine posplt1( x,yr,lmax,jd, ipos, wlabel )
c
      dimension x(jd), yr(jd,lmax)
      character*80 wlabel(1)
c      dimension wlabel(2)
      character*80 string(2)
c
      ymax = -1.0e33
      ymin =  1.0e33
c
      do 40 l=1,lmax
      do 40 j = 2, jd
c
      if ( yr(j,l) .gt. ymax ) ymax = yr(j,l)
      if ( yr(j,l) .lt. ymin ) ymin = yr(j,l)
c 
   40 continue
c
      xmin = x(1)
      xmax = x(jd)
      if ( ymax .eq. ymin ) ymax = ymin + 1.0
      dx = xmax - xmin
      dy = ymax - ymin
      xw = xmin + dx/3.0
      yw = ymax + dy/20.0
c
      go to (110,120,130,140), ipos
c
  110 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.67,0.97)
      go to 200
  120 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.67,0.97)
      go to 200
  130 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.3,0.6)
      go to 200
  140 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.3,0.6)
c
  200 continue
c
      jdm1=jd-1
      do 50 l=1,lmax
ckg      call trace(x(2),yr(2,l), jdm1,-1,-1,0.0,0.0)
   50 continue
ckg      call setlch ( xw,yw, 0,1,0,-1 )
c
      write(string,300 ) wlabel
ckg      call gtext(string,-1,-1)
  300 format ( 2a8 )
c
      return
      end
c
      subroutine posplt2( x,yr,lmax,jd, ipos, isrf, wlabel )
c
      dimension x(jd), yr(jd,lmax)
      dimension wlabel(1)
      character*80 string(10)
c
      ymax = -1.0e33
      ymin =  1.0e33
c
      do 40 l=1,lmax
      do 40 j = 2, jd
c
      if ( yr(j,l) .gt. ymax ) ymax = yr(j,l)
      if ( yr(j,l) .lt. ymin ) ymin = yr(j,l)
c
   40 continue
c
      xmin = x(1)
      xmax = x(jd)
      if ( ymax .eq. ymin ) ymax = ymin + 1.0
      dx = xmax - xmin
      dy = ymax - ymin
      xw = xmin + dx/3.0
      yw = ymax + dy/20.0
c
      go to (110,120,130,140), ipos
c
  110 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.67,0.97)
      go to 200
  120 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.67,0.97)
      go to 200
  130 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.3,0.6)
      go to 200
  140 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.3,0.6)
c
  200 continue
c
      jdm1=jd-1
      do 50 l=1,lmax
ckg      call trace(x(2),yr(2,l), jdm1,-1,-1,0.0,0.0)
   50 continue
ckg      call setlch ( xw,yw, 0,1,0,-1 )
c
      write(string,300 ) wlabel, isrf
ckg      call wrtstr(string,2)
  300 format ( a8, i5 )
c
      return
      end
      subroutine posplt3( x,yr,zr,lmax,jd, ipos, isrf, wlabel )
c
      dimension x(jd), yr(jd,lmax), zr(jd,lmax)
      dimension wlabel(1)
      character*80 string(10)
c
      ymax = -1.0e33
      ymin =  1.0e33
c
      do 40 l=1,lmax
      do 40 j = 2, jd
c
      if ( yr(j,l) .gt. ymax ) ymax = yr(j,l)
      if ( yr(j,l) .lt. ymin ) ymin = yr(j,l)
      if ( zr(j,l) .gt. ymax ) ymax = zr(j,l)
      if ( zr(j,l) .lt. ymin ) ymin = zr(j,l)
c
   40 continue
c
      xmin = x(1)
      xmax = x(jd)
      if ( ymax .eq. ymin ) ymax = ymin + 1.0
      dx = xmax - xmin
      dy = ymax - ymin
      xw = xmin + dx/3.0
      yw = ymax + dy/20.0
c
      go to (110,120,130,140), ipos
c
  110 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.67,0.97)
      go to 200
  120 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.67,0.97)
      go to 200
  130 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.3,0.6)
      go to 200
  140 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.3,0.6)
c
  200 continue
c
      jdm1=jd-1
      do 50 l=1,lmax
ckg      call trace(x(2),yr(2,l), jdm1,-1,-1,0.0,0.0)
ckg      call trace(x(2),zr(2,l), jdm1,-1,-1,0.0,0.0)
   50 continue
ckg      call setlch ( xw,yw, 0,1,0,-1 )
c
      write(string,300 ) wlabel, isrf
ckg      call wrtstr(string,2)
  300 format ( a8, i5 )
c
      return
      end
c
      subroutine talk ( words )
c
      include 'clich1'
      include 'clich2'
c
      character*80 words( 8 )
c
      write(*,10 ) words
   10 format ( 1x, 10a8 )
c
      return
      end
c
c
      subroutine extrap1(y,x,n)
      dimension y(n),x(n),z(100)
ccc  extrapolate y values at x=0 and 1; x is the square 
ccc root of the  normalized toroidal flux
      nm1=n-1
      dy1=(y(2)-y(1))/(x(2)-x(1))
      z(1)=y(1)-dy1*x(1)
      dyn=(y(n)-y(nm1))/(x(n)-x(nm1))
      z(n)=y(n)+dyn*(1.-x(n))
      y(1)=z(1)
      y(n)=z(n)
      return
      end
c

      subroutine transpin(ihsps)
      include 'clich1'
      include 'clich1b'
      include 'clich2'
       parameter(ndat2=ndat+2)
c      parameter(ndat=40,ndat2=ndat+2)
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
     >,betath(ndat),betaadat(ndat),pthdat(ndat)

ccc read radial profiles from TRANSP data
c moved to datin_pest for consistency with reading taefl.txt file open(10,file='transp.dat',status='old')
c
c   read data files from TRANSP output
c   the data arrays are in normalized poloidal flux (psi) coordinate
ckg      if(ishft.eq.-1)then
ckg         call datin
ckg      else
         call datin_pest(r,zh(ihsps),tip,rmcp,ispc)
cng         print *,r,rmhp,tip
cng         stop
ckg      endif
c
ccc prepare cubic spline interpolation
ccc rgrid is uniform in the square root of psi grid
ccc rgsq is the psi grid
cccc
      do 9 i=1,nn
      rgsq(i)=rgrid(i)**2
    9 continue
      call splprp(rsdat,ndat)
      call splprp1(nosurf,ndat2,rgsq(1))
c
cc interpolate the quantities to the rgsq grid points
c
      call intspl(nosurf,ndat,ndat2,ptotdat,ptot(1),rgsq(1))

      call intspl(nosurf,ndat,ndat2,rhodat,rhoprf(1),rgsq(1))
      call intspl(nosurf,ndat,ndat2,dedat,denc(1,1),rgsq(1))
      call intspl(nosurf,ndat,ndat2,dddat,denc(1,2),rgsq(1))
      call intspl(nosurf,ndat,ndat2,dtdat,denc(1,3),rgsq(1))
      call intspl(nosurf,ndat,ndat2,dzdat,denc(1,5),rgsq(1))
      call intspl(nosurf,ndat,ndat2,dhdat,denc(1,4),rgsq(1))

      call intspl(nosurf,ndat,ndat2,dbddat,denh(1,1),rgsq(1))
      call intspl(nosurf,ndat,ndat2,dbtdat,denh(1,2),rgsq(1))
      call intspl(nosurf,ndat,ndat2,dadat,denh(1,3),rgsq(1))

c      call intspl(nosurf,ndat,ndat2,drdat,denr(1),rgsq(1))

      call intspl(nosurf,ndat,ndat2,tedat,tc(1,1),rgsq(1))
      call intspl(nosurf,ndat,ndat2,tddat,tc(1,2),rgsq(1))
      call intspl(nosurf,ndat,ndat2,ttdat,tc(1,3),rgsq(1))
      call intspl(nosurf,ndat,ndat2,thdat,tc(1,4),rgsq(1))
      call intspl(nosurf,ndat,ndat2,tzdat,tc(1,5),rgsq(1))
c H minority temperature, valid only if tip='i' and set in inprofe
      call intspl(nosurf,ndat,ndat2,thmindat,th(1,3),rgsq(1))

      call intspl(nosurf,ndat,ndat2,pbddat,ph(1,1),rgsq(1))
      call intspl(nosurf,ndat,ndat2,pbtdat,ph(1,2),rgsq(1))
      call intspl(nosurf,ndat,ndat2,padat,ph(1,3),rgsq(1))
ckg this array is not actually used
      call intspl(nosurf,ndat,ndat2,betaadat,betaa(1),rgsq(1))
c           rewind 10
c
   99 format(5e14.5)

      return
      end
c
      subroutine smooth(y,n)
      dimension y(n),z(400)
      if(n.gt.400) stop 'Check n comp to its limit'
      nm1=n-1
      do 1 i=2,nm1
      z(i)=(y(i-1)+y(i)+y(i+1))/3.0
    1 continue
      z(1)=(y(1)+y(2))/2.0
      z(n)=(y(nm1)+y(n))/2.0
      do 2 i=1,n
      y(i)=z(i)
    2 continue
      return
      end
c
c
cccccc****************************************
c  spline package
ccccccc****************************************
c
      subroutine intspl(n,nsp,nsp2,f,ff,xg)
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2)
      common/splinc/gspl(nnp2,nn),gsplp(nnp2,nn)
      dimension ff(n),xg(n)
      dimension f(nsp),cf(nsp2)
c
      call depose(f,cf)
c
      do 1 i=1,n
      ff(i)=0.0
      do 2 isp=1,nsp2
      ff(i)=ff(i)+cf(isp)*gspl(isp,i)
    2 continue
    1 continue
c
      return
      end
      subroutine intsplp(n,nsp,nsp2,f,ff,xg)
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2)
      common/splinc/gspl(nnp2,nn),gsplp(nnp2,nn)
      dimension ff(n),xg(n)
      dimension f(nsp),cf(nsp2)
c
      call depose(f,cf)
c
      do 1 i=1,n
      ff(i)=0.0
      do 2 isp=1,nsp2
      ff(i)=ff(i)+cf(isp)*gsplp(isp,i)
    2 continue
    1 continue
c
      return
      end
      subroutine splprp(rg,nsp)
      include 'clich1a'
      parameter(n=nogrid,np2=n+2)
      common/rcl/r(n)
      common/sacl/sa(4,4,np2)
      common/phi01cl/p01(n)/phi02cl/p02(n)/phi03cl/p03(n)
      common/phi11cl/p11(n)/phi12cl/p12(n)/phi13cl/p13(n)
      common/phi21cl/p21(n)/phi22cl/p22(n)/phi23cl/p23(n)
      common/phii1cl/pi1(n)/phii2cl/pi2(n)
      common/phii3cl/pi3(n)/phii4cl/pi4(n)
      common/ncl/nm2,nm1,nnn,np1,nnnp2
      common/phi4cl/dum4(n)
      common/phi5cl/dum5(n)
      dimension rg(nsp)
ccc
ccc  set up spline grids; nsp must be less than or equal to nogrid
c
      nm1=nsp-1
      nm2=nsp-2
      nnn=nsp
      np1=nsp+1
      nnnp2=nsp+2
      do 1 i=1,nsp
      r(i)=rg(i)
    1 continue
c
      call bsplcof
      call deset(0)
c
      return
      end
      subroutine splprp1(n,nsp2,xg)
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2)
      common/splinc/gspl(nnp2,nn),gsplp(nnp2,nn)
      dimension xg(n)
c
      do 1 i=1,n
      do 2 isp=1,nsp2
      call spval(isp,xg(i),spl,1)
      call spvalp(isp,xg(i),splp,1)
      gspl(isp,i)=spl
      gsplp(isp,i)=splp
    2 continue
    1 continue
      return
      end
      subroutine spval(i,x,y,m)
c---computes m values of the i-th b-spline at the m x-values and
c   returns them in y. note that the x-values are assumed to be
c---ordered.
c---
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2)
      dimension x(m),y(m)
      common/rcl/r(nn)
      common/sacl/sa(4,4,nnp2)
      common/ncl/nm2,nm1,n,np1,np2
      kk=1
      do 10 j=1,m
      y(j)=0.0
      do 11 k=kk,n
      if(x(j).le.r(k)) go to 12
   11 continue
   12 l=max0(2,k)-i+3
      kk=k
      if(l.lt.1.or.l.gt.4) go to 10
      y(j)=(sa(1,l,i)+x(j)*(sa(2,l,i)+x(j)*(sa(3,l,i)+x(j)*sa(4,l,i))))
   10 continue
      return
      entry spvalp(i,x,y,m)
      kk=1
      do 20 j=1,m
      y(j)=0.0
      do 21 k=kk,n
      if(x(j).le.r(k)) go to 22
   21 continue
   22 l=max0(2,k)-i+3
      kk=k
      if(l.lt.1.or.l.gt.4) go to 20
      y(j)=sa(2,l,i)+x(j)*(2.0*sa(3,l,i)+3.0*x(j)*sa(4,l,i))
   20 continue
      return
      entry spvalpp(i,x,y,m)
      kk=1
      do 30 j=1,m
      y(j)=0.0
      do 31 k=kk,n
      if(x(j).le.r(k)) go to 32
   31 continue
   32 l=max0(2,k)-i+3
      kk=k
      if(l.lt.1.or.l.gt.4) go to 30
      y(j)=2.0*sa(3,l,i)+6.0*sa(4,l,i)*x(j)
   30 continue
      return
      end
      subroutine depose(fn,c)
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2)
      dimension fn(1000),c(1)
      common/rcl/r(nn)
      common/ncl/nm2,nm1,n,np1,np2
      common/phi01cl/p01(nn)/phi02cl/p02(nn)/phi03cl/p03(nn)
      common/phi11cl/p11(nn)/phi12cl/p12(nn)/phi13cl/p13(nn)
      common/phi4cl/e(nn)/phi5cl/f(nn)
      nm3=n-3
      if(abs(fn(1)).gt.abs(10.*fn(2))) go to 20
c---find derivative of fcn at end pts.
       d11=(fn(2)-fn(1))/(r(2)-r(1))
       d12=(fn(3)-fn(2))/(r(3)-r(2))
       d13=(fn(4)-fn(3))/(r(4)-r(3))
       d21=(d12-d11)/(r(3)-r(1))
       d22=(d13-d12)/(r(4)-r(2))
       d31=(d22-d21)/(r(4)-r(1))
      fp0=d11-d21*r(2)+d31*r(2)*r(3)
      c(1)=fn(1)/p01(1)
   22 continue
c---
       d11=(fn(nm2)-fn(nm3))/(r(nm2)-r(nm3))
       d12=(fn(nm1)-fn(nm2))/(r(nm1)-r(nm2))
       d13=(fn(n)-fn(nm1))/(r(n)-r(nm1))
       d21=(d12-d11)/(r(nm1)-r(nm3))
       d22=(d13-d12)/(r(n)-r(nm2))
       d31=(d22-d21)/(r(n)-r(nm3))
       fp1=+r(n)*(d13+(r(n)-r(nm1))*(d22+d31*(r(n)-r(nm2))))
c---compute e and f.
      f(2)=(fp0-p11(1)*c(1))/p12(1)
      e(2)=0.0
      e(3)=-p03(2)/p02(2)
      f(3)=(fn(2)-p01(2)*f(2))/p02(2)
      do 10 i=4,n
       e(i)=1.0/(p01(i-1)*e(i-1)+p02(i-1))
       f(i)=(fn(i-1)-p01(i-1)*f(i-1))*e(i)
       e(i)=-p03(i-1)*e(i)
   10 continue
c---compute c.
      c(np2)=fn(n)/p03(n)
      c(np1)=(fp1-p13(n)*c(np2))/p12(n)
       do 12 i=2,n
       j=np2-i
       c(j)=e(j)*c(j+1)+f(j)
   12 continue
      return
   20 continue
       d11=(fn(3)-fn(2))/(r(3)-r(2))
       d12=(fn(4)-fn(3))/(r(4)-r(3))
       d13=(fn(5)-fn(4))/(r(5)-r(4))
       d21=(d12-d11)/(r(4)-r(2))
       d22=(d13-d12)/(r(5)-r(3))
       d31=(d22-d21)/(r(5)-r(2))
      fn0=fn(2)+(r(1)-r(2))*(d11+(r(1)-r(3))*(d21
     2                          +(r(1)-r(4))* d31))
      fp0=d11+d21*((r(1)-r(3))+(r(1)-r(2)))+d31*((r(1)-r(3))*
     2    (r(1)-r(4))+(r(1)-r(2))*(r(1)-r(4))+
     3 (r(1)-r(2))*(r(1)-r(3)))
      c(1)=fn0/p01(1)
      go to 22
      end
      subroutine repose(c,v)
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2)
      dimension c(1),v(1)
      common/rcl/r(nn)
      common/ncl/nm2,nm1,n,np1,np2
      common/phi01cl/p01(nn)/phi02cl/p02(nn)/phi03cl/p03(nn)
      common/phi11cl/p11(nn)/phi12cl/p12(nn)/phi13cl/p13(nn)
      common/phi21cl/p21(nn)/phi22cl/p22(nn)/phi23cl/p23(nn)
      common/phii1cl/pi1(nn)/phii2cl/pi2(nn)/phii3cl/pi3(nn)
      common/phii4cl/pi4(nn)
c---computes the function values at the test points from the spline coef
      do 10 i=1,n
      v(i)=c(i)*p01(i)+c(i+1)*p02(i)+c(i+2)*p03(i)
   10 continue
      return
      entry reposp(c,v)
      do 11 i=1,n
      v(i)=c(i)*p11(i)+c(i+1)*p12(i)+c(i+2)*p13(i)
   11 continue
      return
      entry repospp(c,v)
      do 12 i=1,n
      v(i)=c(i)*p21(i)+c(i+1)*p22(i)+c(i+2)*p23(i)
   12 continue
      return
      entry reposi(c,v)
      v(1)=0.0
      do 13 i=2,n
      v(i)=v(i-1)+c(i-1)*pi1(i)+c(i)*pi2(i)+c(i+1)*pi3(i)+
     2            c(i+2)*pi4(i)
   13 continue
      return
      end
      subroutine splerr(c,er)
      dimension er(1),c(1)
c---relative error estimate of spline fit.
c---note# the routine only returns a value in er(i) if
c---      the error is greater than the initial value of er(i).
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2)
      common/ncl/nm2,nm1,n,np1,np2
      common/sacl/sa(4,4,nnp2)
      common/rcl/r(nn)
      fcn(i,l,x)=sa(1,l,i)+x*(sa(2,l,i)+x*(sa(3,l,i)+x*sa(4,l,i)))
      fval=0.0
      do 11 i=1,n
      fval=fval+abs(c(i)*fcn(i,4,r(i))+c(i+1)*fcn(i+1,3,r(i))+
     2 c(i+2)*fcn(i+2,2,r(i)))
   11 continue
      fval=fval/n
      if(fval.eq.0.0) fval=1.
      do 10 i=2,nm1
      err=0.02625*(c(i-1)*(           -sa(4,4,i-1))
     2              +c(i)  *(sa(4,4,i)  -sa(4,3,i)  )
     3              +c(i+1)*(sa(4,3,i+1)-sa(4,2,i+1))
     4              +c(i+2)*(sa(4,2,i+2)-sa(4,1,i+2))
     5              +c(i+3)*(sa(4,1,i+3)))
     6         *(amax1(r(i+1)-r(i),r(i)-r(i-1))**3)
      err=abs(err/fval)
      if(err.gt.er(i))er(i)=err
   10 continue
      return
      end
      subroutine bsplcof
c---routine computes the coefficients of the b-splines which form
c   a basis over the set of knots, r, with repeated knots
c   at the end points.  s(j,l,i) j-power, l-segment, i-spline.
c---
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2,nsa=16*nnp2)
      common/ncl/nm2,nm1,n,np1,np2
      common/rcl/r(nn)
      common/sacl/sa(nsa)
c---i=1
      c4=4.00/((r(2)-r(1))**4)
      idx=16
      sa(idx)   =              -c4
      sa(idx-1) =           3.0*c4*r(2)
      sa(idx-2) =          -3.0*c4*r(2)*r(2)
      sa(idx-3) =               c4*r(2)*r(2)*r(2)
c---i=2
      c4=4.0/((r(3)-r(2))*(r(3)-r(1))**3)
      c3=4.0/((r(2)-r(3))*(r(2)-r(1))**3)
      idx=32
      sa(idx)   =              -c4
      sa(idx-1) =           3.0*c4*r(3)
      sa(idx-2) =          -3.0*c4*r(3)*r(3)
      sa(idx-3) =               c4*r(3)*r(3)*r(3)
      sa(idx-4) =sa(idx)       -c3
      sa(idx-5) =sa(idx-1) +3.0*c3*r(2)
      sa(idx-6) =sa(idx-2) -3.0*c3*r(2)*r(2)
      sa(idx-7) =sa(idx-3)     +c3*r(2)*r(2)*r(2)
c---i=3
      c4=4.0/((r(4)-r(2))*(r(4)-r(3))*(r(4)-r(1))**2)
      c3=4.0/((r(3)-r(2))*(r(3)-r(4))*(r(3)-r(1))**2)
      c2=4.0/((r(2)-r(3))*(r(2)-r(4))*(r(2)-r(1))**2)
      idx=48
      sa(idx)   =              -c4
      sa(idx-1) =          +3.0*c4*r(4)
      sa(idx-2) =          -3.0*c4*r(4)*r(4)
      sa(idx-3) =               c4*r(4)*r(4)*r(4)
      sa(idx-4) =sa(idx)       -c3
      sa(idx-5) =sa(idx-1) +3.0*c3*r(3)
      sa(idx-6) =sa(idx-2) -3.0*c3*r(3)*r(3)
      sa(idx-7) =sa(idx-3)     +c3*r(3)*r(3)*r(3)
      sa(idx-8) =sa(idx-4)     -c2
      sa(idx-9) =sa(idx-5) +3.0*c2*r(2)
      sa(idx-10)=sa(idx-6) -3.0*c2*r(2)*r(2)
      sa(idx-11)=sa(idx-7)     +c2*r(2)*r(2)*r(2)
c---i=4,n-1
      do 10 i=4,nm1
      m1=i-3
      m2=i-2
      m3=i-1
      m4=i
      m5=i+1
      c4=4.0/((r(m5)-r(m1))*(r(m5)-r(m2))*(r(m5)-r(m3))*(r(m5)-r(m4)))
      c3=4.0/((r(m4)-r(m1))*(r(m4)-r(m2))*(r(m4)-r(m3))*(r(m4)-r(m5)))
      c2=4.0/((r(m3)-r(m1))*(r(m3)-r(m2))*(r(m3)-r(m4))*(r(m3)-r(m5)))
      c1=4.0/((r(m2)-r(m1))*(r(m2)-r(m3))*(r(m2)-r(m4))*(r(m2)-r(m5)))
      idx=16*i
      sa(idx)   =              -c4
      sa(idx-1) =          +3.0*c4*r(m5)
      sa(idx-2) =          -3.0*c4*r(m5)*r(m5)
      sa(idx-3) =               c4*r(m5)*r(m5)*r(m5)
      sa(idx-4) =sa(idx)       -c3
      sa(idx-5) =sa(idx-1) +3.0*c3*r(m4)
      sa(idx-6) =sa(idx-2) -3.0*c3*r(m4)*r(m4)
      sa(idx-7) =sa(idx-3)     +c3*r(m4)*r(m4)*r(m4)
      sa(idx-8) =sa(idx-4)     -c2
      sa(idx-9) =sa(idx-5) +3.0*c2*r(m3)
      sa(idx-10)=sa(idx-6) -3.0*c2*r(m3)*r(m3)
      sa(idx-11)=sa(idx-7)     +c2*r(m3)*r(m3)*r(m3)
      sa(idx-12)=sa(idx-8)     -c1
      sa(idx-13)=sa(idx-9) +3.0*c1*r(m2)
      sa(idx-14)=sa(idx-10)-3.0*c1*r(m2)*r(m2)
      sa(idx-15)=sa(idx-11)    +c1*r(m2)*r(m2)*r(m2)
   10 continue
c---i=n
      m1=n-3
      m2=n-2
      m3=n-1
      m4=n
      c4=12.0/((r(m4)-r(m1))*(r(m4)-r(m2))*(r(m4)-r(m3)))
      c3=-4.0/((r(m4)-r(m2))*(r(m4)-r(m3))*(r(m4)-r(m1))**2)
     2   -4.0/((r(m4)-r(m1))*(r(m4)-r(m3))*(r(m4)-r(m2))**2)
     3   -4.0/((r(m4)-r(m1))*(r(m4)-r(m2))*(r(m4)-r(m3))**2)
      c2= 4.0/((r(m3)-r(m1))*(r(m3)-r(m2))*(r(m3)-r(m4))**2)
      c1= 4.0/((r(m2)-r(m1))*(r(m2)-r(m3))*(r(m2)-r(m4))**2)
      idx=16*n
      sa(idx-4) =              -c3
      sa(idx-5) =c4        +3.0*c3*r(m4)
      sa(idx-6) = -2.0*c4  -3.0*c3*r(m4)*r(m4)
      sa(idx-7) =c4*r(m4)*r(m4)+c3*r(m4)*r(m4)*r(m4)
      sa(idx-8) =sa(idx-4)     -c2
      sa(idx-9) =sa(idx-5) +3.0*c2*r(m3)
      sa(idx-10)=sa(idx-6) -3.0*c2*r(m3)*r(m3)
      sa(idx-11)=sa(idx-7)     +c2*r(m3)*r(m3)*r(m3)
      sa(idx-12)=sa(idx-8)     -c1
      sa(idx-13)=sa(idx-9) +3.0*c1*r(m2)
      sa(idx-14)=sa(idx-10)-3.0*c1*r(m2)*r(m2)
      sa(idx-15)=sa(idx-11)    +c1*r(m2)*r(m2)*r(m2)
c---i=n+1
      c4= 12.0/((r(m4)-r(m2))    *(r(m4)-r(m3)))
      c3=-12.0/((r(m4)-r(m3))    *(r(m4)-r(m2))**2)
     2   -12.0/((r(m4)-r(m2))    *(r(m4)-r(m3))**2)
      c2=  4.0/((r(m4)-r(m3))    *(r(m4)-r(m2))**3)
     2    +4.0/((r(m4)-r(m2))**2 *(r(m4)-r(m3))**2)
     3    +4.0/((r(m4)-r(m2))    *(r(m4)-r(m3))**3)
      c1=  4.0/((r(m3)-r(m2))    *(r(m3)-r(m4))**3)
      idx=16*(n+1)
      sa(idx-8) =                        -c2
      sa(idx-9) =          c3        +3.0*c2*r(m4)
      sa(idx-10)=-c4  -2.0*c3*r(m4)  -3.0*c2*r(m4)*r(m4)
      sa(idx-11)= c4*r(m4)+c3*r(m4)*r(m4)+c2*r(m4)*r(m4)*r(m4)
      sa(idx-12)=sa(idx-8)     -c1
      sa(idx-13)=sa(idx-9) +3.0*c1*r(m3)
      sa(idx-14)=sa(idx-10)-3.0*c1*r(m3)*r(m3)
      sa(idx-15)=sa(idx-11)    +c1*r(m3)*r(m3)*r(m3)
c---i=n+2
      c4= +4.0/(r(m4)-r(m3))
      c3=-12.0/(r(m4)-r(m3))**2
      c2=+12.0/(r(m4)-r(m3))**3
      c1= -4.0/(r(m4)-r(m3))**4
      idx=16*(n+2)
      sa(idx-12)=                          -c1
      sa(idx-13)=            c2        +3.0*c1*r(m4)
      sa(idx-14)=  -c3  -2.0*c2*r(m4)  -3.0*c1*r(m4)*r(m4)
      sa(idx-15)=c4+c3*r(m4)+c2*r(m4)*r(m4)+c1*r(m4)*r(m4)*r(m4)
      return
      end
      subroutine deset(jh)
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2)
      common/rcl/r(nn)
      common/ncl/nm2,nm1,n,np1,np2
      common/phi01cl/p01(nn)/phi02cl/p02(nn)/phi03cl/p03(nn)
      common/phi11cl/p11(nn)/phi12cl/p12(nn)/phi13cl/p13(nn)
      common/phi21cl/p21(nn)/phi22cl/p22(nn)/phi23cl/p23(nn)
      common/phii1cl/pi1(nn)/phii2cl/pi2(nn)/phii3cl/pi3(nn)
      common/phii4cl/pi4(nn)
      common/bvpcl/d1,d2,d3,d4,d5,d6
ckg      data pi1/nn*0.0/,pi2/nn*0.0/,pi3/nn*0.0/,pi4/nn*0.0/
c---subroutine sets up values of splines at the knots.
      do 10 i=1,n
      pi1(i)=0.
      pi2(i)=0.
      pi3(i)=0.
      pi4(i)=0.
      call spval(i,r(i),p01(i),1)
      call spval(i+1,r(i),p02(i),1)
      call spval(i+2,r(i),p03(i),1)
      call spvalp(i,r(i),p11(i),1)
      call spvalp(i+1,r(i),p12(i),1)
      call spvalp(i+2,r(i),p13(i),1)
      call spvalpp(i,r(i),p21(i),1)
      call spvalpp(i+1,r(i),p22(i),1)
      call spvalpp(i+2,r(i),p23(i),1)
   10 continue
      p02(1)=0.
      p03(1)=0.
      p01(n)=0.
      p02(n)=0.
      p13(1)=0.
      p11(n)=0.
      do 12 i=2,n
      pi1(i)=gauss6(i-1,i,jh)
      pi2(i)=gauss6(i,i,jh)
      pi3(i)=gauss6(i+1,i,jh)
      pi4(i)=gauss6(i+2,i,jh)
   12 continue
      d1=p11(1)
      d2=p12(1)
      d3=p12(n)
      d4=p13(n)
      d5=p03(n)
      d6=p01(1)
      return
      end
      function gauss6(k,i,jh)
      include 'clich1a'
      parameter(nn=nogrid,nnp2=nn+2)
      dimension f(6),x(6)
      common/rcl/r(nn)
      data w1,w2,w3/.467913934572691,.360761573048139,.171324492379170/
      data x1,x2,x3/.238619186083197,.661209386466265,.932469514203152/
      rd=0.5*(r(i)-r(i-1))
      rp=0.5*(r(i)+r(i-1))
      x(1)=rp-rd*x3
      x(2)=rp-rd*x2
      x(3)=rp-rd*x1
      x(4)=rp+rd*x1
      x(5)=rp+rd*x2
      x(6)=rp+rd*x3
      call spval(k,x(1),f(1),6)
      gauss6=rd*(w3*(x(1)**jh*f(1)+x(6)**jh*f(6))
     2          +w2*(x(2)**jh*f(2)+x(5)**jh*f(5))
     4          +w1*(x(3)**jh*f(3)+x(4)**jh*f(4)))
      return
      end
      subroutine mapg(xmin,xmax, ymin,ymax, x,y, a,b)
      return
      end
c
c**************************************************************************
c Here we transform the Pest eigenvector to Nova format
      subroutine pest2nova_egv
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      dimension wrk(nn,mt,3),wrk_s(nn,mt,3)
      data hhh/-1000./
      save emax
      if(emax.ne.0..and.hhh.eq.-1000.) emax=0.
      call equiv(wrk,eigfun,nn*mt*3)
      call equiv(wrk_s,eigfun_s,nn*mt*3)
      htheta=twopi/real(mth)
      isrf=nsrf
      idelta=1
      do imm=1,mt
         eigfun(isrf,imm,1)=0.
         eigfun(isrf,imm,3)=0.
         eigfun_s(isrf,imm,1)=0.
         eigfun_s(isrf,imm,3)=0.
         immi=max(imm-idelta,1)
         imma=min(imm+idelta,mt)
         theta=-htheta
         do itheta=1,mth
            theta=theta+htheta
            hhh_s=q*(-grpsth(itheta)**2+grthsq(itheta)*grpssq(itheta))
     &           *f*xjacob(itheta)/r/g/r
            hhh=hhh_s*cos((minm+imm-1)*theta)
            hhh_s=hhh_s*sin((minm+imm-1)*theta)
            hh1=xjacob(itheta)*bsq(itheta)/r/r/g*sin((minm+imm-1)
     &           *theta)
            hh1_s=xjacob(itheta)*bsq(itheta)/r/r/g*cos((minm+imm-1)
     &           *theta)
            hcs=f*grpsth(itheta)/grpssq(itheta)
            do imh=immi,imma
               ss=wrk_s(isrf,imh,1)*sin((minm+imh-1)*theta)
               cs=wrk(isrf,imh,1)*cos((minm+imh-1)*theta)
c
               eigfun(isrf,imm,1)=eigfun(isrf,imm,1)+cs*hhh
               eigfun_s(isrf,imm,1)=eigfun_s(isrf,imm,1)+ss*hhh_s
c
               eigfun(isrf,imm,3)=eigfun(isrf,imm,3)-(wrk(isrf,imh,3)
     &              *sin((minm+imh-1)*theta)+cs*hcs)*hh1
               eigfun_s(isrf,imm,3)=eigfun_s(isrf,imm,3)+(wrk_s(isrf
     &              ,imh,3)*cos((minm+imh-1)*theta)+ss*hcs)*hh1_s
            enddo
         enddo
         eigfun(isrf,imm,1)=-eigfun(isrf,imm,1)*htheta/pi
         eigfun(isrf,imm,3)=-eigfun(isrf,imm,3)*htheta/pi
         eigfun_s(isrf,imm,1)=-eigfun_s(isrf,imm,1)*htheta/pi
         eigfun_s(isrf,imm,3)=-eigfun_s(isrf,imm,3)*htheta/pi
         eigfun(isrf,imm,2)=pp*eigfun(isrf,imm,1)
         eigfun_s(isrf,imm,2)=pp*eigfun_s(isrf,imm,1)
         emax=max(emax,abs(eigfun(isrf,imm,1)))
      enddo
c      write(*,*) eigfun(isrf,2-minm,1),eigfun(isrf,2-minm,3)
c     &     ,eigfun_s(isrf,2-minm,1),eigfun_s(isrf,2-minm,3)
c      if(isrf.ge.nosurf)then
c         do i=1,nosurf
c            do j=1,3
c               do imm=1,mt
c                  eigfun(i,imm,j)=eigfun(i,imm,j)/emax
c               enddo
c            enddo
c         enddo
c      endif
      return
      end
c**************************************************************************
c is actually used for nova not pest
c running with ASTRA needs care how R0, amin, b0 are read, as well as nastra is set = ndat
      subroutine datin_pest(r,szh,tip,rmcp,ispc)
      include 'gridparam'
c      parameter(ndat=20)
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
     >,betath(ndat),betaadat(ndat),pthdat(ndat)
c     >,ufastpa(ndat),ufastpp(ndat)
      dimension betaH(ndat),bdry(ndat),rmcp(ispc)
c2d     >,dvol(ndat),darea(ndat),surfa(ndat)
      common / chnel / iodat,iomode,iomap1,ioequ1,mp0,mp1,mp2,iotty
      common/plasmap/rmaj,amin,b0
      character*80 title,sc
      character*(*) tip
      parameter(ntaefl=101,nprfls=15)
      dimension ptaefl(ntaefl,nprfls)
c      parameter(nastra=96,nprfastr=66)
      parameter(nastra=105,nprfastr=66)
      dimension prastra(nastra,nprfastr)
      logical ftransp/.false./, ftaefl/.false./,fastra19/.false./
      logical fvzee/.false./
      integer io/10/
c
      open(io,file='transp.dat',err=334,status='old')
      ftransp=.true.
      goto 335
 334  open(io,file='astra.dat',err=337,status='old')
      if(nastra/=ndat) stop 'change ndat to = ntaefl '
      fastra19=.true.
      goto 335
 337  open(io,file='taefl.txt',err=338,status='old')
      if(ntaefl-1/=ndat) stop 'change ndat to = ntaefl-1 '
      ftaefl=.true.
c here we are having taefl file reading option
      write(*,*)'No transp.dat file was found, found ',
     &     'taefl.txt to read profiles'
      read(io,*)
      read(io,*)
      read(io,*) b0
      read(io,*)
      read(io,*) rmaj
      read(io,*)
      read(io,*) amin
c      amin=amin
      b0=b0*1.e4
      do i=1,9
         read(io,*)
      enddo
      do i=1,ntaefl
         read(io,*) (ptaefl(i,j),j=1,nprfls)
      enddo
      print *,'Read the profiles',b0,rmaj,amin
      goto 335
 338  open(io,file='cat_307002_00006_walphas.dat',err=333,status='old')
      if(ntaefl-1/=ndat) stop 'change ndat to = ntaefl-1 '
      fvzee=.true.
      rmcp(2)=2.5
c here we are having taefl file reading option
      write(*,*)'No transp.dat file was found, found ',
     & 'cat_307002_00006_walphas.dat to read profiles'
      read(io,*)
      do i=1,ntaefl
         read(io,*) (ptaefl(i,j),j=1,14)
c         print *,(ptaefl(i,j),j=1,14)
      enddo
      print *,'Read the Vanzee profiles and from NOVA_param',
     &     b0,rmaj,amin
 335  continue
      if(ftransp) then
         read(10,8) sc
c      read(10,'(a26,f11.4)') sc(1:26),b0
         read(10,'(a16,e10.3)') sc(1:16),b0
         write(*,*) sc(1:16),b0
 32      continue 
         read(10,8) sc
         if(sc(2:9).ne.'Profiles')goto 32
c taefl branch already known true b0
c--------------------------------------------
         read(10,8) sc
c x=zone centr grid = sq_root(normlz toroid flux) (~r/a) @X BepHo
         read(10,99,ERR=33) (rtsdat(k),k=1,ndat)
c     
         read(10,8) sc
c     x=zone boundary grid = sq_root(normlz toroid flux) (~r/a) @XB BepHo
         read(10,99) (bdry(k),k=1,ndat)
      elseif(ftaefl)then
         do i=1,ndat
            bdry(i)=ptaefl(i+1,1) ! sqrt(tor flux)
            rtsdat(i)=bdry(i)-.5/ndat
            if(i==1) then
               fluxp(1)=0.
            else
               fluxp(i)=fluxp(i-1)+(ptaefl(i,1)**2-ptaefl(i-1,1)
     *              **2)*2./(ptaefl(i,2)+ptaefl(i-1,2))
            endif
         enddo
         do i=1,ndat
            fluxp(i)=fluxp(i)/fluxp(ndat)
c            print *,fluxp(i),ptaefl(i+1,1)
         enddo
      elseif(fvzee)then
         do i=1,ndat
            bdry(i)=ptaefl(i+1,1) ! sqrt(tor flux)
            rtsdat(i)=bdry(i)-.5/ndat
            if(i==1) then
               fluxp(1)=0.
            else
               fluxp(i)=ptaefl(i+1,14)
            endif
         enddo
      else
         write(*,*)'Found astra.dat file only.'
c         ter He Tst R=6.2 a=2 B=2.65 I=7.5 q=3.69 <n>=3.64 Time=95.126 dt=100            
c         read(io,'(a14,f4.2,a2,f2.0,a2,f5.2)') sc,rmaj,sc,amin,sc,b0
c        iter HeH+ Be R=6.2 a=2 B=2.65 I=7.5 q=4.34 <n>=3.02 Time=****** dt=200 
         read(io,'(a16,f4.2,a2,f2.0,a2,f5.2)') sc,rmaj,sc,amin,sc,b0
         write(*,*) 'ASTRA is read with R0,a,B0',rmaj,amin,b0
c         write(*,*) 'Skiping 19 lines'
         write(*,*) 'Skiping 21 lines'
         b0=b0*1.e4
         do i=1,20
            read(io,*)
         enddo
         do i=1,nastra
            read(io,*) (prastra(i,j),j=1,nprfastr)
            fluxp(i)=prastra(i,54)
            qdat(i)=prastra(i,7)
            if(i==1) then
               fluxt(i)=fluxp(i)*qdat(i)
            else
               fluxt(i)=fluxt(i-1)+(fluxp(i)-fluxp(i-1))*qdat(i) ! tor flux 
            endif
         enddo
         print *,'Read Astra profiles',b0,rmaj,amin
         bdry(1)=sqrt(fluxt(1)/fluxt(ndat))
         rtsdat(1)=bdry(1)/2.
         do i=2,ndat
            bdry(i)=sqrt(fluxt(i)/fluxt(ndat)) ! norm sqrt(tor flux)
            rtsdat(i)= bdry(i)-.5/ndat !sqrt(fluxt(i)-fluxt(i-1)) ! bdry(i)-.5/ndat
c            print *,fluxp(i),qdat(i),bdry(i)
         enddo
c         stop
      endif
c
      if(ftransp) then
         do iread=1,4
            read(10,8) sc
c     Midplane major radii (cm) @XB      RPLOT abrev = RMJMP   
            read(10,99) (shift(k),k=1,ndat)
         enddo
         rmaj=shift(ndat)

         read(10,8) sc
c     Midplane minor radii (cm) @XB      RPLOT abrev = RMNMP 
         read(10,99) (shift(k),k=1,ndat)
         amin=shift(ndat)
         b0=b0/rmaj*10000.

         read(10,8) sc
c     Shafranov shift (cm) of zone bdry
         read(10,99) (shift(k),k=1,ndat)

         do iread=1,2
            read(10,8) sc
c Toroidal flux (Webers) at zone bdry
            read(10,99) (fluxt(k),k=1,ndat)
         enddo

         read(10,8) sc
c Poloidal Flux (Webers) at zone bdry
         read(10,99) (fluxp(k),k=1,ndat)
      endif

      if(ftransp) then
c
         read(10,8) sc
c Electron density at x (/cm3)
         read(10,99) (dedat(k),k=1,ndat)

         read(10,8) sc
c Hydrogen density at x (/cm3)
         read(10,99) (dhdat(k),k=1,ndat)

         read(10,8) sc
c Thermal deuterium density at x (/cm3)
         read(10,99) (dddat(k),k=1,ndat)

         read(10,8) sc
c Thermal tritium density at x (/cm3)
         read(10,99) (dtdat(k),k=1,ndat)

         read(10,8) sc
c Impurity density at x (/cm3)
         read(10,99) (dzdat(k),k=1,ndat)

         read(10,8) sc
c Total ion density (NH+ND+NT+NIMP) at x (/cm3)
         read(10,99) (diondat(k),k=1,ndat)

         read(10,8) sc
c Fast alpha density NFI at x (/cm3)
         read(10,99) (dadat(k),k=1,ndat)

         do iread=1,2
            read(10,8) sc
c     Total beam ion density BDENS at x (/cm3) - not used
            read(10,99) (dbeamdat(k),k=1,ndat)
         enddo

         read(10,8) sc
c  D-Beam ion density BDENS_D at x (/cm3)
         read(10,99) (dbddat(k),k=1,ndat)

         read(10,8) sc
c T-Beam ion density BDENS_T at x (/cm3)
         read(10,99) (dbtdat(k),k=1,ndat)

         read(10,8) sc
c Impurity Temperature at x (keV)
         read(10,99) (tzdat(k),k=1,ndat)

         read(10,8) sc
c Hydrogenic (ions) Temperature at x (keV)
         read(10,99) (tddat(k),k=1,ndat)

         read(10,8) sc
c Electron Temperature at x (keV)
         read(10,99) (tedat(k),k=1,ndat)

         do iread=1,6
            read(10,8) sc
c     q_psi at zone bdry
            read(10,99) (qdat(k),k=1,ndat)
         enddo
      elseif(fastra19)then
         do i=1,ndat
c            bdry(i) sqrt(tor flux)
            dedat(i)=1.e13*prastra(i,4)
            dhdat(i)=1.e13*prastra(i,8)
            dddat(i)=1.e13*prastra(i,12)
            dtdat(i)=0.
            dzdat(i)=1.e13*prastra(i,13)
            diondat(i)=dddat(i)+dzdat(i)+dhdat(i)+dtdat(i)
            dadat(i)=0.
            dbddat(i)=1.e13*prastra(i,16)
            dbtdat(i)=0.
            tzdat(i)=prastra(i,11)*1.e3
            tddat(i)=prastra(i,11)*1.e3
            tedat(i)=prastra(i,3)*1.e3
         enddo
      elseif(fvzee)then
         do i=1,ndat
            dedat(i)=1.e13*ptaefl(i+1,4)
            dhdat(i)=0.
            dddat(i)=1.e13*ptaefl(i+1,4)/ptaefl(i+1,5)
            dtdat(i)=0.
            dzdat(i)=0.
            diondat(i)=dddat(i)+dzdat(i)+dhdat(i)+dtdat(i)
            dadat(i)=0.
            dbddat(i)=1.e13*ptaefl(i+1,8)
            dbtdat(i)=0.
            tzdat(i)=ptaefl(i+1,3)*1.e3
            tddat(i)=ptaefl(i+1,3)*1.e3
            tedat(i)=ptaefl(i+1,2)*1.e3
         enddo
c         continue
      else
         do i=1,ndat
c            bdry(i)=ptaefl(i+1,1) ! sqrt(tor flux)
            dedat(i)=ptaefl(i+1,5)*1.e13
            dhdat(i)=0.
            dddat(i)=ptaefl(i+1,4)*1.e13
            dtdat(i)=0.
            dzdat(i)=ptaefl(i+1,6)*1.e13
            diondat(i)=dddat(i)+dzdat(i)+dhdat(i)+dtdat(i)
            dadat(i)=0.
            dbddat(i)=ptaefl(i+1,3)*1.e13
            dbtdat(i)=0.
            tzdat(i)=ptaefl(i+1,8)*1.e3
            tddat(i)=ptaefl(i+1,8)*1.e3
            tedat(i)=ptaefl(i+1,9)*1.e3
            qdat(i)=ptaefl(i+1,2)
         enddo
      endif


      if(ftransp) then
         read(10,8) sc
c Electron energy density 3*NE*TE/2 (Joules/cm3)
         read(10,99) (ue(k),k=1,ndat)

         read(10,8) sc
c Ion energy density 3*NI*TI/2 (Joules/cm3)
         read(10,99) (ui(k),k=1,ndat)

         read(10,8) sc
c D-beam parallel energy density (Joules/cm3) UBPAR_D
         read(10,99) (ubdpar(k),k=1,ndat)

         read(10,8) sc
c D-beam perpendicular energy density (Joules/cm3) UBPRP_D
         read(10,99) (ubdper(k),k=1,ndat)

         read(10,8) sc
c T-beam parallel energy density (Joules/cm3)
         read(10,99) (ubtpar(k),k=1,ndat)

         read(10,8) sc
c T-beam perpendicular energy density (Joules/cm3)
         read(10,99) (ubtper(k),k=1,ndat)

         read(10,8) sc
c Fusion ion parallel energy density (Joules/cm3)
         read(10,99) (ufipa(k),k=1,ndat)

         read(10,8) sc
c Fusion ion perpendicular energy density (Joules/cm3)
         read(10,99) (ufipp(k),k=1,ndat)

         read(10,8) sc
c  Rotation energy density (Joules/cm3)
         read(10,99) (uphi(k),k=1,ndat)
c      if(int(szh).eq.1)then
c         print *, 'reading fast ion energies'
         read(10,8) sc
c Fast ion parallel energy density (Joules/cm3)
         read(10,99) (ufastpa(k),k=1,ndat)

         read(10,8) sc
c Fast ion perpendicular energy density (Joules/cm3)
         read(10,99) (ufastpp(k),k=1,ndat)
c      else
c         print *, 'read fusion ion energies'
c         do iread=1,2
c            read(10,8) sc
c            read(10,99) (utot(k),k=1,ndat)
c         enddo
c      endif
         read(10,8) sc
c     Total energy den (UE+UI+UPHI+UFASTPA+UFASTPP)
         read(10,99) (utot(k),k=1,ndat)
         
         read(10,8) sc
c DT neutron emission rate (#/cm3/sec)
         read(10,99) (ftotdt(k),k=1,ndat)
         
         read(10,8) sc
c Total neutron emission rate (#/cm3/sec)
         read(10,99) (ttntx(k),k=1,ndat)
 
         do iread=1,3
            read(10,8) sc
c     Thermal Toroidal Beta (=BTE+BTI) at x
            read(10,99) (betath(k),k=1,ndat)
         enddo
         
         do iread=1,2
            read(10,8) sc
c     Minority beta  (TOROIDAL)  BTMIN
            read(10,99) (betaH(k),k=1,ndat)
         enddo

         if(int(szh).eq.1)then
            read(10,8) sc
c     dummy reading
            read(10,99) (betadat(k),k=1,ndat)
         else
            read(10,8) sc
c     Alpha Toroidal Beta at x BTFI
            read(10,99) (betaadat(k),k=1,ndat)
         endif

         read(10,8) sc
c Total Toroidal Beta (=BTPL+BTBM+BTMIN+BTFI) at x BTTOT
         read(10,99) (betadat(k),k=1,ndat)

         read(10,8) sc
c Thermal Pressure 2(UI+UE)/3 (Pascals)   PPLAS
         read(10,99) (pthdat(k),k=1,ndat)
         
         read(10,8) sc
c Total Pressure (thermal+3*UPHI/2+UFASTPA+UFASTPP/2) PTOWB
         read(10,99) (ptotdat(k),k=1,ndat)

c
         read(10,8) sc
c H minority temperature TMINI_H
         read(10,99) (thmindat(k),k=1,ndat)
      
         if(tip(1:1).eq.'i'.or.tip(1:1).eq.'h')then
            read(10,8) sc
c H minority density  NMINI_H
            if(sc(7:13).eq.'NMINI_H')then
               read(10,99) (dhmindat(k),k=1,ndat)
            else
               write(*,*) 'WARNING'
               write(*,*) 'i did not find H minority T prof.'
               write(*,*) 'switching to `h` distr.func.'
               tip='h'
            endif
         endif

c      read(10,8) sc
c Z_eff profile
c      read(10,99) (zeffdat(k),k=1,ndat)
         rewind 10              ! finished reading transp.dat
      endif

      if(ftransp) then
         do 21 k=1,ndat
cc  thermal Tritium and Hydrogen have same temperature as thermal Deuterium 
            ttdat(k)=tddat(k)
            thdat(k)=tddat(k)
cc  define beam-D and beam-T pressure
            pbddat(k)=(ubdpar(k)+ubdper(k)*0.5)*1.e6
            pbtdat(k)=(ubtpar(k)+ubtper(k)*0.5)*1.e6
cc  H-min.press.=tot.f.ion - D and T beam energy - fusion ion density
            if(int(szh).eq.1)then
               padat(k)=((ufastpa(k)-ufipa(k)-ubdpar(k)-ubtpar(k))+(
     &            ufastpp(k)-ufipp(k)-ubdper(k)-ubtper(k))*0.5)*1.e6
            elseif(tip(1:1).eq.'i'.or.tip(1:1).eq.'h')then
cc  minority pressure used for benchmark with Irvine ?? gives the same equiv beta of minorities as the transp beta of beam ions
               padat(k)=dhmindat(k)*thmindat(k)*0.25e-12
            else
cc  alpha pressure
               padat(k)=(ufipa(k)+ufipp(k)*0.5)*1.e6
            endif
ckg this part is to include cherez zhopu ICRF ions (don't need it)
c      padat(k)=padat(1)*betaadat(k)/betaadat(1)
c      print *,'----------====== padat is',padat(k)

   21    continue
      elseif(fastra19)then
         do k=1,ndat
c            bdry(k)=ptaefl(k+1,1) ! sqrt(tor flux)
cc  thermal Tritium and Hydrogen have same temperature as thermal Deuterium 
            ttdat(k)=tddat(k)
            thdat(k)=tddat(k)
cc  define beam-D and beam-T pressure from transp %
            pbddat(k)=prastra(k,51)*1.e-2
            pbtdat(k)=0.
            ptotdat(k)=(prastra(k,51)+prastra(k,55)+prastra(k,59))*1.e-2
cc  H-min.press.=tot.f.ion - D and T beam energy - fusion ion density
            if(int(szh).eq.1)then
               padat(k)=prastra(k,51)*1.e-2*0.5    !cnng19 /2 to make beta as from astra
            elseif(tip(1:1).eq.'i')then
cc  minority pressure used for benchmark with Irvine ?? gives the same equiv beta of minorities as the transp beta of beam ions
               padat(k)=0.
            else
cc  alpha pressure
               padat(k)=prastra(k,55)*1.e-2
            endif
c            print *,fluxp(k),qdat(k),bdry(k),padat(k)
         enddo
c         stop
c         continue
      elseif(fvzee)then
         do k=1,ndat
cc  thermal Tritium and Hydrogen have same temperature as thermal Deuterium
            ttdat(k)=tddat(k)
            thdat(k)=tddat(k)
cc  define beam-D and beam-T pressure from Pa
            pbddat(k)=ptaefl(k+1,9)*1.e-6*1.e6*1.e3
            pbtdat(k)=0.
c ptaefl(k+1,12)*1.e3
cc  H-min.press.=tot.f.ion - D and T beam energy - fusion ion density
            if(int(szh).eq.1)then
               padat(k)=ptaefl(k+1,13)*1.e3
            elseif(tip(1:1).eq.'s')then
cc  minority pressure used for benchmark with Irvine ?? gives the same equiv beta of minorities as the transp beta of beam ions
               padat(k)=ptaefl(k+1,13)*1.e3
            else
cc  alpha pressure
               padat(k)=0.
            endif
            ptotdat(k)=padat(k)+pbddat(k)+(ptaefl(k+1,2)*ptaefl(k+1,4)
     &       +ptaefl(k+1,3)*ptaefl(k+1,4)/ptaefl(k+1,5))*1.602e-3*1.e3          ! 1.602x10^-19x10^3eVx10^13cm^-3
         enddo
      else                      ! do equiv "reading" of taefl.txt file
         do k=1,ndat
c            bdry(k)=ptaefl(k+1,1) ! sqrt(tor flux)
cc  thermal Tritium and Hydrogen have same temperature as thermal Deuterium 
            ttdat(k)=tddat(k)
            thdat(k)=tddat(k)
cc  define beam-D and beam-T pressure from transp JLES/CM3 to kPa
            pbddat(k)=ptaefl(k+1,10)*1.e3
            pbtdat(k)=0.
            ptotdat(k)=ptaefl(k+1,12)*1.e3
cc  H-min.press.=tot.f.ion - D and T beam energy - fusion ion density
            if(int(szh).eq.1)then
               padat(k)=0.
            elseif(tip(1:1).eq.'i')then
cc  minority pressure used for benchmark with Irvine ?? gives the same equiv beta of minorities as the transp beta of beam ions
               padat(k)=0.25e-12
            else
cc  alpha pressure
               padat(k)=0.
            endif
         enddo
      endif
 336  continue
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
      call smooth(thmindat,ndat)
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
     &         +1.*dhmindat(i)
         if(tip(1:1).eq.'i')dadat(i)=dhmindat(i)
   13 continue
c
c
cc extrapolate the quantities to normalized sqrt(tor_flux)=0 and 1 
      if(ftransp) call boundtoc(bdry,fluxp,rtsdat,ndat) ! go from boundary grid to central one for the grid
      if(ftransp) call extrap1(fluxp,rtsdat,ndat)
c
      call extrap1(dedat,rtsdat,ndat)
      call extrap1(rhodat,rtsdat,ndat)
      call extrap1(ptotdat,rtsdat,ndat)
      call extrap1(betadat,rtsdat,ndat)
ckg
      call boundtoc(bdry,qdat,rtsdat,ndat)
      call extrap1(qdat,rtsdat,ndat)
c
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
      call extrap1(thmindat,rtsdat,ndat)

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
 33   write(*,*) sc
      write(*,*) '--> above line should match ndat, which is:',
     &     ndat
      write(*,*)'--> but it does not'
      write(*,*)'--> modify ndat in gridparam and recompile gotae'
      stop ' --> Stopping for now.'
 333  stop 'no transp.dat or taefl.txt or astra.dat present, stopping'
      end
c
c**************************************************************************
c     Here we calculate TAE growth rate using analytical formula for
C     simple benchmark in simple tokamak approximation
c     [see Gorelenkov et.al. NF 2005 on ITER.]
c     Added are the estimates for the effective scatt. and drag frequencies
c     as per Duarte writeup, /u/ngorelen/Ql2dtae/Ql2dtaeetex/averagingDuarte.pdf
c
      subroutine gg_ae_theory(ihsps,immax,isrfmaxin,wtae,omstar
     &        ,b0_vd,valfvkvi,ai_vd,zi_vd,vi_vd)
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      include 'orb2d.par'
      common/functn/omreal1,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
ckg limit how close you can get to the boundaries
      isrfmax=max(isrfmaxin,11)
      isrfmax=min(isrfmax,nosurf-11)
      p2gc=max(chidelt(1)**2,0.1)
      print *, '--------------------------------------'
      print *, 'Analytical estimate of TAE growth rate'
      write(*,*)'At peak Xir surface index =',isrfmaxin
      print *, 'For beams use Chi0 = ',sqrt(p2gc)
ckg evaluate derivative in the vicinity using ten points
      sum1=0
      sum2=0
      do i=isrfmax-10,isrfmax+10
         sum1=sum1+(ph(i,ihsps)-ph(isrfmax,ihsps))*(ppsi(i)-ppsi(isrfmax
     &        ))
         sum2=sum2+(ppsi(i)-ppsi(isrfmax))**2
      enddo
      sum2=sum2*ph(isrfmax,ihsps)
      gwskgw=-omstar*0.5
c     this was attempt to use analytical formula for the omega*, but
c     mine gwskgw is good (npoBepeHo 3 pa3a)
c      gwst_anlt=sum1*(minm+immax-1)*qprf(nosurf)*rr(1,1)*vi/(rr(nosurf,1
c     &     )-rr(1,1))**2/(sum2*valfvkvi*omreal1*gwci)
c      stop
c beam zow gg, no fow effects
      ggthkgw=-3.*betah0(ihsps)*ph(isrfmax,ihsps)/ph(1,ihsps)*
     .     3.1415926*0.125*qprf(isrfmax)**2*
     .     valfvkvi**2*(1.+p2gc)**2/p2gc*(3.
     .     -(3.-2.*p2gc+3.*p2gc**2)/(p2gc*(1.+p2gc))
     .     +2.*valfvkvi**2*gwskgw*sum1/sum2/p2gc)
c slowing down gg with fow, i.e. at the plato regime
      ggthkgw_se=5.*3.1415926*qprf(isrfmax)**2*ppsi(isrfmax)
     &     *betah0(ihsps)*ph(isrfmax,ihsps)/ph(1,ihsps)*(-sum1/sum2-1.
     &     /gwskgw)*valfvkvi*(1.-valfvkvi**2)
c 
      ggthkgw_m=0.
      ggthkgw_sezow=0.
      ggthkgw_sezowd=0.
      v_crit=0.58
c sum over fundamental resonances
      do i=1,2
         gcr=1.-(i*2-3)*qprf(nosurf)/abs(qprf(isrfmax)*omreal1)
         gcr=abs(gcr)
         gcr=valfvkvi/max(gcr,0.001)
         gcrv=gcr/v_crit
         gcrv=min(gcrv,1.)
         gcrv=max(gcrv,gcr)
c                       check if gcr is not big enough for passing particles
         if(gcr.le.1..and.(gcr.gt.sqrt(rgrid(isrfmax)*(rr(nosurf,1)
     &        -rr(1,1))/rr(1,1))))then
            write(*,*) 'v||res/vA    =',gcr,' < Chi_separatrix/2 ='
     &           ,sqrt(rgrid(isrfmax)*(rr(nosurf,1)-rr(1,1))/rr(1,1))
            write(*,*) 'v||res/vA/vcr=',gcrv
c maxwellian with isotropic pitch distribution without fow effects
c  for two fundamental v||=vA and v||=vA/3 resonance
            gcr2=gcr**2
            gcr4=gcr2**2
            gcrv2=gcrv**2
            gcrv4=gcrv2**2
            ggthkgw_m=ggthkgw_m+sqrt(3.1415926)*0.25*(qprf(nosurf)
     &           /omreal1)**2*betah0(ihsps)*ph(isrfmax,ihsps)/ph(1,ihsps
     &           )*(-gwskgw*sum1/sum2-1.)*gcr*(0.5+gcr2+gcr4)*exp(
     &           -gcr2)
c slowing down with isotropic pitch distribution without fow effects
c  for two fundamental v||=vA and v||=vA/3 resonance
            ggthkgw_sezow=ggthkgw_sezow+(
c omega* part
     &           -gwskgw*sum1/sum2*(gcrv-gcr)*
     &           gcr4*(1.+0.3333*(6.+1./gcrv2+1./gcr2+1./(gcr*gcrv))
     &           /(gcr*gcrv))
c omega* part with v_crit
     &           -gwskgw*sum1/sum2*gcr**7*(0.5/v_crit**3)*
     &           (0.3333/(gcrv2*gcrv4)+1./gcrv4+1./gcrv2-2.33333)
c Landau damping part
     &           -1.5*gcr2*(gcrv-gcr)*(2.+0.3333*(gcr2+gcrv2+gcr*gcrv)+
     &           1./(gcr*gcrv) )
     &           )
c damping at the cutoff
            ggthkgw_sezowd=ggthkgw_sezowd-gcr*(0.5+gcr2+0.5*gcr4)
         endif
      enddo
      ggthkgw_sezow =ggthkgw_sezow *3.*3.1415926*(qprf(nosurf)/omreal1)
     &     **2*0.03125*betah0(ihsps)*ph(isrfmax,ihsps)/ph(1,ihsps)
     &     /(1.-0.6*v_crit**2)
      ggthkgw_sezowd=ggthkgw_sezowd*3.*3.1415926*(qprf(nosurf)/omreal1)
     &     **2*0.03125*betah0(ihsps)*ph(isrfmax,ihsps)/ph(1,ihsps)
     &     /(1.-0.6*v_crit**2)
      write(*,*) 'r/a,dlnfkdps=',rgrid(isrfmax),sum1/sum2
      write(*,*) 'betah0      =',betah0(ihsps)
      write(*,*) 'gg_beam_zow =',ggthkgw
      write(*,*) 'gg_sldwn_zow=',ggthkgw_sezow,'cutoff damping'
     &     ,ggthkgw_sezowd,'<- gg_sldwn_zow is without cuttoff'
      write(*,*) 'gg_sldwn_fow=',ggthkgw_se
      write(*,*) 'gg_maxw_zow =',ggthkgw_m
c this is place reserved for Duarte contribution
c
      a_vd=0.5*abs(rr(nosurf,1)-rr(nosurf,mth/2))
      R0_vd=0.5*abs(rr(nosurf,1)+rr(nosurf,mth/2))
      dne_vd=denc(isrfmax,1)/1.e14
      dni_vd=(denc(isrfmax,2)+denc(isrfmax,3))/1.e14
      Te_vd=tc(isrfmax,1)/1.e3
      z1_vd=0.5 - denc(isrfmax,3)/(1.e14*6.*dne_vd)
      r_vd=rgrid(isrfmax)*sqrt(qprf(isrfmax))*a_vd
      q_vd=qprf(isrfmax)
      dqdr_vd=rqp(isrfmax)*psitot/(sqrt(q_vd)*a_vd)
      epsilon_vd=r_vd/r0_vd
      tauSmod_vd=0.01*(Te_vd**(3./2.))*Ai_vd/(((Zi_vd)**2)*dne_vd) !tauS is the characteristic time for slowing-down
      xmplusl_vd=ntor*q_vd
      n_vd=ntor
      vstar_vd=0.16*((z1_vd)**(1./3.))*
     &     (Te_vd**(1./2.))*10000000.  !vstar is the critical velocity at which ion drag starts to be more important than electron drag
      vA_vd=valfvkvi*vi_vd*1.e6        !Alfvén velocity
      xkparal_vd=n_vd /(2. * R0_vd)
      vparal_vd=q_vd * R0_vd * vA_vd *xkparal_vd/((1.+2.*
     &     (q_vd**2.))**(1./2))        !this expression is for BAAE gap, modified Alfven branch
      tauB_vd = (2.*3.14* q_vd * R0_vd)/vA_vd     !tau bounce
      averageScatt_vd=6.28*((r_vd**2.)*R0_vd)/(vparal_vd*q_vd*
     &     tauB_vd)
      averageDrag_vd=6.28*(q_vd*(R0_vd**2.))/((2.**(1./2))* tauB_vd)
c Fokker-Planck coefficients
      xnuPerp_vd=(1./2)*((vstar_vd)**3.)/((vparal_vd**3.)*tauSmod_vd)
      xnuDrag_vd=1./tauSmod_vd
      dOmegadPphi_vd=xmplusl_vd*(vparal_vd/(R0_vd*q_vd**2.))*dqdr_vd* 
     &     (Ai_vd*1.04*(1.e-8)*q_vd/(B0_vd * r_vd)) !derivative of Omega with respect to P_phi
      print *, "dOmegadPphi=", dOmegadPphi_vd
c Effective collisional coefficients

      xnuEffScatt_vd = (averageScatt_vd*xnuPerp_vd*(vparal_vd**2.)*
     &     (dOmegadPphi_vd**2.))**(1./3)

      xnuEffDrag_vd=(averageDrag_vd*(2.**(1./2))*(tauSmod_vd**(-1.))*
     &   (1.+(vstar_vd**3.)/(vparal_vd**3.))* dOmegadPphi_vd)**(1./2)

      print *, "nuEffScatt_vd=",xnuEffScatt_vd
      print *, "nuEffDrag_vd=",xnuEffDrag_vd
      
c
      return
      end
c**************************************************************************
ckg this is subroutine used to calculate radiative damping accroding to paper Liming Yu, 2009
ckg Written By Yahui Wang
ckg date: Oct. 9, 2019
ckg POP, Kinetic damping of Alfven eigenmodes in general tokamak geometry.
ckg clindrial coordinate is needed.

      subroutine RD_RSAE(immax)
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      include 'orb2d.par'
      integer (kind=4)jx,jy
      integer (kind=4) gqmin		! the grid number in radius direction where q_min exist.
      dimension f_alfven(nn)
      dimension dq2dr(nn)
      dimension dqdr(nn)
      real    (kind=8)v_alfven
      real (kind=8) rad_damp
      real (kind=8) r_min	! minor radius at q_min
      real (kind=8) rho_i	! the gyro radius of ion paritcles at q_min
      real (kind=8) rho_s,rho_sp
      real (kind=8) lambda_eign
      real (kind=8) S_eign
      dimension coeff_rad(2)
      parameter(nrad=400)
      dimension w1(nrad),w2(nrad)
      dimension rmin(nn,nts)
      real rr00,zz00
      dimension deltar(nn,nts)
      dimension deltarsq(nn)
      real (kind=8) rb0
      real (kind=8) temp1,temp2
!for test
      dimension xx(nn)
      real (kind=8) qmin
      dimension qx(nn)
      logical test
!-----------------------------------------------------------------------
!! variables to estimate the continue damping accroding to paper (Zonca, POP 2002)
      real (kind=8) S_c, gamma_com

      common/plasmap/ rmaj,amin,b0
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
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
     >,betath(ndat),betaadat(ndat),pthdat(ndat)
c      common/dat/fluxp(ndat),rdat(ndat),rhodat(ndat)
      write(*,*) "call RD_RSAE"
      open(137,file='datcon',status='old')
c here om is TAE eigenvalue from NOVA code (omega**2)
ckg      write(6,*)"input om"
ckg      read(5,*)om 
c     here w1 and w2 are first and second continuum curves as
c     calculated by continuum code (in omega**2) 
c     (w1 and w2 should be non-zero and finite)
c
      read(137,*)n1,n2
      do i=n1,n2
      read(137,*)w1(i),w2(i)
!      write(*,*)i,w1(i),w2(i)
      enddo
      close(137)
      coeff_rad(1)=0.1125
      coeff_rad(2)=1.2206
      write(*,*) 'ntor=', ntor
      write(*,*) 'the polodial mode number is', minm+immax-1
!      write(*,*) "nn in radiative_damping_w is", nn
!      write(*,*)  'major radius is', rmaj	! in unit of (cm)
      va=2.18e11*b0/sqrt(rho0)			! in unit of (cm/s)
      do jx=1,nn,1
         f_alfven(jx)=va*sqrt(bsqav(jx))/sqrt(rhoprf(jx))
     >          /(2.*twopi*qprf(jx)*rmaj)		! in unit of (Hz)
      enddo
      write(*,*) 'the alfven speed and center density are:'
      write(*,'(2e16.7)') va,rho0
      write(*,*)'omreal in unit of (Hz) is',omreal*va/rmaj/qprf(nn)
     &   /twopi

      do jx=1,nn,1
!        write(*,*)jx,rgrid(jx),grpssq2d(jx,1),grpssq2d(jx,mth/2)

      enddo 

      rr00=(rr(1,1)+rr(1,mth/2+1))*0.5
      zz00=(zz(1,1)+zz(1,mth/2+1))*0.5
      do jx=1,nn
         do i=1,nts
            rmin(jx,i)=sqrt((rr(jx,i)-rr00)**2.0+(zz(jx,i)-zz00)**2.0)
         enddo
      enddo
      do jx=2,nn
        do i=1,nts
          deltar(jx,i)=rmin(jx,i)-rmin(jx-1,i)
        enddo
        deltarsq(jx)=rgrid(jx)**2.0-rgrid(jx-1)**2.0
      enddo
c      open(unit=1001,file='fluxp.dat')
c      do jx=1,ndat,1
c        read(1001,*)fluxp(jx)
c      enddo
c      close(1001)


      rb0=b0/10000.0*rmaj/100.0/((rr(1,1)+rr(1,mth/2+1))/2.0)      !in unit of (1)
!      do jx=1,mth/2+2
!        write(*,*)jx,rr(1,jx),rb0
!      enddo 


      do jx=2,nn-1,1
        dqdr(jx)=(qprf(jx+1)-qprf(jx-1))/(rgrid(jx+1)**2.0-
     &    rgrid(jx-1)**2.0)/(fluxp(ndat)-fluxp(1))
     &    *sqrt(grpssq2d(jx,1))*rb0
        dqdr(jx)=(deltarsq(jx)**2.0*qprf(jx+1)-deltarsq(jx+1)**2.0*
     &   qprf(jx-1)-(deltarsq(jx)**2.0-deltarsq(jx+1)**2.0)*qprf(jx))
     &  /(deltarsq(jx+1)**2.0*deltarsq(jx)+deltarsq(jx)**2.0*
     &  deltarsq(jx+1))/(fluxp(ndat)-fluxp(1))*sqrt(grpssq2d(jx,1))
     &  *rb0
        temp1=(rgrid(jx+1)**2.0-rgrid(jx-1)**2.0)/
     &   (rr(jx+1,1)-rr(jx-1,1))*(fluxp(ndat)-fluxp(1))
        temp1=(deltar(jx,1)**2.0*rgrid(jx+1)**2.0-deltar(jx+1,1)**2.0
     &    *rgrid(jx-1)**2.0-(deltar(jx,1)**2.0-deltar(jx+1,1)**2.0)
     &    *rgrid(jx)**2.0)/(deltar(jx+1,1)**2.0*deltar(jx,1)
     &   +deltar(jx,1)**2.0*deltar(jx+1,1))*(fluxp(ndat)-fluxp(1))


        temp2=(deltar(jx,1)**2.0*qprf(jx+1)-deltar(jx+1,1)**2.0
     &    *qprf(jx-1)-(deltar(jx,1)**2.0-deltar(jx+1,1)**2.0)*qprf(jx))
     &    /(deltar(jx+1,1)**2.0*deltar(jx,1)+deltar(jx,1)**2.0
     &    *deltar(jx+1,1))

!        write(*,*)jx,temp1
!     &   ,sqrt(grpssq2d(jx,1))*rb0,dqdr(jx),temp2


!        dq2dr(jx)=(qprf(jx+1)-2.0*qprf(jx)+qprf
!     &  (jx-1))*4.0/(rgrid(jx+1)**2.0-rgrid(jx-1)**2.0)**2.0*
!     &  grpssq2d(jx,1) +0.5*(qprf(jx+1)-qprf(jx-1))/(rgrid(jx+1)**2.0
!     &  -rgrid(jx-1)**2.0)**2.0*(grpssq2d(jx+1,1)-grpssq2d(jx-1,1))

        dq2dr(jx)=2.0*(deltarsq(jx)*qprf(jx+1)+deltarsq(jx+1)*qprf
     &   (jx-1)-(deltarsq(jx)+deltarsq(jx+1))*qprf(jx))/((deltarsq
     &   (jx)+deltarsq(jx+1))*deltarsq(jx)*deltarsq(jx+1))*
     &  grpssq2d(jx,1) +0.5*(deltarsq(jx)**2.0*qprf(jx+1)-deltarsq(jx+1)
     &   **2.0*qprf(jx-1)-(deltarsq(jx)**2.0-deltarsq(jx+1)**2.0)*
     &  qprf(jx))*(deltarsq(jx)**2.0*grpssq2d(jx+1,1)-deltarsq(jx+1)
     &   **2.0*grpssq2d(jx-1,1)-(deltarsq(jx)**2.0-deltarsq(jx+1)**2.0)*
     &  grpssq2d(jx,1))/(deltarsq(jx+1)**2.0*deltarsq(jx)+deltarsq(jx)
     &  **2.0*deltarsq(jx+1))
        dq2dr(jx)=dq2dr(jx)/(fluxp(ndat)-fluxp(1))**2.0
     &    *(rb0)**2.0
       temp1=dq2dr(jx)
        dq2dr(jx)=2.0*(deltar(jx,1)*qprf(jx+1)+deltar
     &   (jx+1,1)*qprf(jx-1)-(deltar(jx,1)+deltar(jx+1,1))*qprf(jx))
     &   /((deltar(jx,1)+deltar(jx+1,1))*deltar(jx,1)*deltar(jx+1,1))
       temp2=dq2dr(jx) 
!       write(*,*)jx,temp1,temp2

        dq2dr(jx)=2.0*(deltarsq(jx)*qprf(jx+1)+deltarsq(jx+1)*qprf
     &   (jx-1)-(deltarsq(jx)+deltarsq(jx+1))*qprf(jx))/((deltarsq
     &   (jx)+deltarsq(jx+1))*deltarsq(jx)*deltarsq(jx+1))*grpssq2d
     &  (jx,mth/2+1)+0.5*(deltarsq(jx)**2.0*qprf(jx+1)-deltarsq(jx+1)
     &   **2.0*qprf(jx-1)-(deltarsq(jx)**2.0-deltarsq(jx+1)**2.0)*
     &  qprf(jx))*(deltarsq(jx)**2.0*grpssq2d(jx+1,mth/2+1)-deltarsq
     &   (jx+1)**2.0*grpssq2d(jx-1,mth/2+1)-(deltarsq(jx)**2.0-
     &   deltarsq(jx+1)**2.0)*grpssq2d(jx,mth/2+1))/(deltarsq(jx+1)**2.0
     &   *deltarsq(jx)+deltarsq(jx)**2.0*deltarsq(jx+1))
        dq2dr(jx)=dq2dr(jx)/(fluxp(ndat)-fluxp(1))**2.0
     &    *(rb0)**2.0
       temp1=dq2dr(jx)
        dq2dr(jx)=2.0*(deltar(jx,mth/2+1)*qprf(jx+1)+deltar
     &   (jx+1,mth/2+1)*qprf(jx-1)-(deltar(jx,mth/2+1)+deltar
     &   (jx+1,mth/2+1))*qprf(jx))/((deltar(jx,mth/2+1)+deltar
     &   (jx+1,mth/2+1))*deltar(jx,mth/2+1)*deltar(jx+1,mth/2+1))
       temp2=dq2dr(jx)       

      enddo

      do jx=2,nn-1,1
!        dq2dr(jx)=(qprf(jx+1)-2.0*qprf(jx)+qprf
!     &  (jx-1))*4.0/(rmin(jx+1,mth/2+1)-rmin(jx-1,mth/2+1))**2.0
!        dq2dr(jx)=(qprf(jx+1)-2.0*qprf(jx)+qprf
!     &  (jx-1))*4.0/(rmin(jx+1,1)-rmin(jx-1,1))**2.0
        dq2dr(jx)=2.0*(deltar(jx,mth/2+1)*qprf(jx+1)+deltar
     &   (jx+1,mth/2+1)*qprf(jx-1)-(deltar(jx,mth/2+1)+deltar
     &   (jx+1,mth/2+1))*qprf(jx))/((deltar(jx,mth/2+1)+deltar
     &   (jx+1,mth/2+1))*deltar(jx,mth/2+1)*deltar(jx+1,mth/2+1))
       temp1=dq2dr(jx)
        dq2dr(jx)=2.0*(deltar(jx,1)*qprf(jx+1)+deltar
     &   (jx+1,1)*qprf(jx-1)-(deltar(jx,1)+deltar(jx+1,1))*qprf(jx))
     &   /((deltar(jx,1)+deltar(jx+1,1))*deltar(jx,1)*deltar(jx+1,1))
       temp2=dq2dr(jx)
       dq2dr(jx)=(temp1+temp2)/2.0 
! dq2dr in unit of (1/m^2)
! rr in unit of (m)
        dq2dr(jx)=dq2dr(jx)/10000.0		! in unit of (1/cm^2)
      enddo

      do jx=2,nn-1,1
        dq2dr(jx)=2.0*(deltarsq(jx)*qprf(jx+1)+deltarsq(jx+1)*qprf
     &   (jx-1)-(deltarsq(jx)+deltarsq(jx+1))*qprf(jx))/((deltarsq
     &   (jx)+deltarsq(jx+1))*deltarsq(jx)*deltarsq(jx+1))*
     &  grpssq2d(jx,1) +0.5*(deltarsq(jx)**2.0*qprf(jx+1)-deltarsq(jx+1)
     &   **2.0*qprf(jx-1)-(deltarsq(jx)**2.0-deltarsq(jx+1)**2.0)*
     &  qprf(jx))*(deltarsq(jx)**2.0*grpssq2d(jx+1,1)-deltarsq(jx+1)
     &   **2.0*grpssq2d(jx-1,1)-(deltarsq(jx)**2.0-deltarsq(jx+1)**2.0)*
     &  grpssq2d(jx,1))/(deltarsq(jx+1)**2.0*deltarsq(jx)+deltarsq(jx)
     &  **2.0*deltarsq(jx+1))
        dq2dr(jx)=dq2dr(jx)/(fluxp(ndat)-fluxp(1))**2.0
     &    *(rb0)**2.0
       temp1=dq2dr(jx)
        dq2dr(jx)=2.0*(deltarsq(jx)*qprf(jx+1)+deltarsq(jx+1)*qprf
     &   (jx-1)-(deltarsq(jx)+deltarsq(jx+1))*qprf(jx))/((deltarsq
     &   (jx)+deltarsq(jx+1))*deltarsq(jx)*deltarsq(jx+1))*grpssq2d
     &  (jx,mth/2+1)+0.5*(deltarsq(jx)**2.0*qprf(jx+1)-deltarsq(jx+1)
     &   **2.0*qprf(jx-1)-(deltarsq(jx)**2.0-deltarsq(jx+1)**2.0)*
     &  qprf(jx))*(deltarsq(jx)**2.0*grpssq2d(jx+1,mth/2+1)-deltarsq
     &   (jx+1)**2.0*grpssq2d(jx-1,mth/2+1)-(deltarsq(jx)**2.0-
     &   deltarsq(jx+1)**2.0)*grpssq2d(jx,mth/2+1))/(deltarsq(jx+1)**2.0
     &   *deltarsq(jx)+deltarsq(jx)**2.0*deltarsq(jx+1))
        dq2dr(jx)=dq2dr(jx)/(fluxp(ndat)-fluxp(1))**2.0
     &    *(rb0)**2.0
       temp2=dq2dr(jx) 
       dq2dr(jx)=(temp1+temp2)/2.0
       dq2dr(jx)=dq2dr(jx)/10000.0		! in unit of (1/cm^2)    
      enddo

!      do jx=2,nn-1,1
!        dq2dr(jx)=(qprf(jx+1)-2.0*qprf(jx)+qprf
!     &  (jx-1))*4.0/(rmin(jx+1,mth/2+1)-rmin(jx-1,mth/2+1))**2.0
!       temp1=dq2dr(jx)
!        dq2dr(jx)=(qprf(jx+1)-2.0*qprf(jx)+qprf
!     &  (jx-1))*4.0/(rmin(jx+1,1)-rmin(jx-1,1))**2.0
!       temp2=dq2dr(jx)
!       dq2dr(jx)=(temp1+temp2)/2.0 
!! dq2dr in unit of (1/m^2)
!! rr in unit of (m)
!        dq2dr(jx)=dq2dr(jx)/10000.0		! in unit of (1/cm^2)
!      enddo

      dq2dr(1)=dq2dr(2)
      dq2dr(nn)=dq2dr(nn-1)
      gqmin=minloc(qprf,1)
      gqmin=min(gqmin,nn-10)
      gqmin=max(11,gqmin)
      write(*,*)"gqmin=",gqmin
      r_min=0.5*(rr(gqmin,1)-rr(gqmin,mth/2+1))*100.0	! in unit of (cm)
      rho_i=(2.0*tc(gqmin,2)*1.602*10.0**(-19.0)*2.0*1.6726*10.0
     &   **(-27.0))**0.5/(rb0*sqrt(bsqav(gqmin))*
     &   1.602*10.0**(-19.0))     ! in unit of (m), tc in unit of (eV) 
      rho_i=rho_i*100.0			! in unit of (cm)
      rho_sp=(2.0*th(gqmin,3)*1.602*10.0**(-19.0)*2.0*1.6726*10.0
     &   **(-27.0))**0.5/(rb0*sqrt(bsqav(gqmin))*
     &   1.602*10.0**(-19.0))     ! in unit of (m), tc in unit of (eV)
      rho_sp=rho_sp*100.0 		! in unit of (cm)

      rho_s=sqrt((tc(gqmin,1)/tc(gqmin,2)+0.75)*rho_i**2.0+denh(gqmin,1)
     &  /denc(gqmin,1)*rho_sp**2.0/(1.0+rho_sp**2.0*
     &  ((minm+immax-1)/r_min)**2.0))

!      rho_s=sqrt((tc(gqmin,1)/tc(gqmin,2)+0.75)*rho_i**2.0)
      write(*,*)tc(gqmin,1),tc(gqmin,2),denc(gqmin,1),denc(gqmin,2)
      write(*,*)th(gqmin,1),th(gqmin,2),th(gqmin,3),denh(gqmin,1)
     &  ,denh(gqmin,2),denh(gqmin,3)
      write(*,*)rho_i,rho_s,rho_sp
      write(*,*)tc(gqmin,1)/tc(gqmin,2)*rho_i**2.0,+0.75*rho_i**2.0,
     & denh(gqmin,1)/denc(gqmin,1)*rho_sp**2.0/(1.0+rho_sp**2.0*
     &  ((minm+immax-1)/r_min)**2.0)         

      lambda_eign=qprf(gqmin)*rho_s**2.0*((minm+immax-1.0)/r_min)**4.0
     &  *(ntor*qprf(gqmin)/(minm+immax-1.0)-1.0)/dq2dr(gqmin)
      s_eign=-(minm+immax-1)*rmaj*qprf(gqmin)**2.0*(ntor-(minm+immax-1)
     &  /qprf(gqmin))/rmaj *((omreal)**2.0
     &  -w1(gqmin))/(w1(gqmin)*(r_min)**2.0
     &  *dq2dr(gqmin))
      rad_damp=coeff_rad(1)*s_eign*exp(-coeff_rad(2)*s_eign/
     &  sqrt(abs(lambda_eign)))

      write(*,*)"ntor=",ntor
      write(*,*)"m=",minm+immax-1
      write(*,*)"rmaj=",rmaj
      write(*,*)"r_min=",r_min
      write(*,*)"qprf(min)=",qprf(gqmin)
      write(*,*)"(n-m/q)/rmaj=",(ntor-(minm+immax-1)/qprf(gqmin))/rmaj
      write(*,*)"omreal",omreal,omreal*va/rmaj/qprf(nn)/twopi
      write(*,*)"omega_A",sqrt(w1(gqmin)),
     &           sqrt(w1(gqmin))*va/rmaj/qprf(nn)/twopi 
      write(*,*)"dq2dr,rho_s",dq2dr(gqmin),rho_s
      write(*,*)"lambda_eign=",lambda_eign
      write(*,*)"S_eign=",S_eign
      write(*,*)"The ratio of Radative damping to wave frequency is:"
     & ,rad_damp


      S_c=(r_min**2.0*dq2dr(gqmin)/qprf(gqmin)**2.0)**0.5
      gamma_con=-2.0*omreal*exp(-4.0*(-ntor*(ntor*qprf(gqmin)
     &          -(minm+immax-1))*qprf(gqmin)/qprf(nn))**0.5/S_c)
      write(*,*)'gamma_con=',gamma_con,S_c,-ntor*(ntor*qprf(gqmin)
     &         -(minm+immax-1))*qprf(gqmin)/qprf(nn),qprf(gqmin)


      open(unit=1001,file='rsae_rd.txt')
      write(1001,*)"ntor=",ntor
      write(1001,*)"m=",minm+immax-1
      write(1001,*)"rmaj=",rmaj
      write(1001,*)"r_min=",r_min
      write(1001,*)"qprf(min)=",qprf(gqmin)
      write(1001,*)"(n-m/q)/rmaj=",
     &  (ntor-(minm+immax-1)/qprf(gqmin))/rmaj
      write(1001,*)"omreal",omreal,omreal*va/rmaj/qprf(nn)/twopi
      write(1001,*)"omega_A",sqrt(w1(gqmin)),
     &           sqrt(w1(gqmin))*va/rmaj/qprf(nn)/twopi 
      write(1001,*)"dq2dr,rho_s",dq2dr(gqmin),rho_s
      write(1001,*)"contribution to rho_s:electon,thermal ion,fast ion",
     & tc(gqmin,1)/tc(gqmin,2)*rho_i**2.0,+0.75*rho_i**2.0,
     & denh(gqmin,1)/denc(gqmin,1)*rho_sp**2.0/(1.0+rho_sp**2.0*
     &  ((minm+immax-1)/r_min)**2.0) 
      write(1001,*)"lambda_eign=",lambda_eign
      write(1001,*)"S_eign=",S_eign
      write(1001,*)"The ratio of Radative damping to wave frequency is:"
     & ,rad_damp
      write(1001,*)'gamma_con=',gamma_con
      close(1001)
      open(unit=1002,file='eq_pro.dat')
      write(1002,'(8a18)')"r","q","te","ti","th","ne","ni","nb"
      do jx=1,nn
        write(1002,'(8e18.9)')rgrid(jx),qprf(jx),tc(jx,1),tc(jx,2)
     &    ,th(jx,3),denc(jx,1),denc(jx,2),denh(jx,1)      
      enddo
      close(1002)
      open(unit=1003,file='eigfun.dat')
      do jx=1,nn
        do jy=1,mt
          write(1003,*)eigfun(jx,jy,1)
        enddo
      enddo
      close(1003)
      return
      end
c**************************************************************************


