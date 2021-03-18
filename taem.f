      include 'clich1'
      include 'clich2'
      include 'clich1b'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      common/plasmap/rmaj,amin,b0
ckg_rad
      dimension anq(nn),alamx(nn),w(nn),tfx(nn)
ckg_rad
      character pest,legname*6
      integer*8 nmap1,lhhh,igivup,ierr,nadres,len
c
ckg      istat=ishell('rm -f otae')
      open(16,file='otae',status='unknown')
ckg      istat=ishell('rm -f otaep')
      open(17,file='otaep',status='unknown')
ckg      istat=ishell('rm -f otaep')
      open(57,file='smpl_out',status='unknown')
c      open(7,file='../Src_NOVAK/leg32',status='old')
      if(lam.gt.99)then
         write(legname(1:3),'(i3)')lam
         legname=legname(1:3)//'leg'
      else
         write(legname(1:2),'(i2)')lam
         legname=legname(1:2)//'leg'
      endif
      iread=8
      open(iread,file=legname,status='old',err=101)
      goto 102
 101  call makeleg(iread,legname)
ckg      stop 'FF'
c 
ckg      istat=ishell('rm -f plottae')
c2d      call ncarcgm(1,'plottae')
c Default values are defined here
 102  call defolt(im1,ihsps,xfow,pk)
c You can redefine toroidal mode number here
ckg change      ntor = 3
      pest='n'
      if(pest.eq.'y')then
         ishft=0
c
         mp1=42
c
         mp0=41
      else
         ishft=-1
      endif
c
      call input
c
c
      nmap1 = 50+9*nogrid+6*(nsf+ishft+1)*nts+5
      nmpdsk  = 12*(nsf+ishft+1)*nts+nogrid+2*nts+50+5
      is = 0
      ndsk = 1
c replacing      call zop ( outmap1, mp1, nmap1, ndsk, is, 100 )
c      call zop ( mp1, 'mpout1', nmap1, ndsk, is, 100 )
      nmap1 = nmap1/512
      lhhh = 0
      call wopen('mpout1',nmap1,lhhh,ierr)
      ndsk = 3
c replacing      call zop ( iomode, mp0, nmpdsk, ndsk, is, 100 )
c      call zop ( mp0, 'mapdsk', nmpdsk, ndsk, is, 100 )
      nmap1 = nmpdsk/512
      call wopen('mapdsk',nmap1,lhhh,ierr)
ccc..    read equilibrium quantities
c replacing      call getwa('mapdsk',ntitle(1),41,1,1,100)
c
      nadres=1
      len=41
      call getwa('mapdsk',ntitle(1),nadres,len,ierr)
      nosurf=nosurf1
c      call zrd(mp0,ntitle(1),41,1+ishft,igivup,ierr)
      write(*,'(20a)') ntitle
c      write(*,*) dat, nx, nz, nosurf1, mth, ngrid,
c     .                alx, alz, xzero, xma, r, p0, gp00,
c     .                psimin, psilim, psipls,
c     .                betag, betap, upsiln, nnj, mmj
      njg = nnj
      mjx = mmj
      igrid=ngrid
ccccccccccccccccccccccccccccccccccccccccccccccccc
ccc..  find out gconst to rescale q
      lgivup=1
      nadres=50+nosurf*2
c replacing      call zrd(outmap1,q,1,nadres,ierr)
c      call zrd(mp1,q,1,nadres+ishft,igivup,ierr)
c
      len=1
      call getwa('mpout1',q,nadres,len,ierr)
      nadres=nadres+nosurf*2
c replacing      call zrd(outmap1,g,1,nadres,ierr)
c      call zrd(mp1,g,1,nadres+ishft,igivup,ierr)
c
      call getwa('mpout1',g,nadres,len,ierr)
      nadres=nadres+nosurf*2
c replacing      call zrd(outmap1,f,1,nadres,ierr)
c      call zrd(mp1,f,1,nadres+ishft,igivup,ierr)
c
      call getwa('mpout1',f,nadres,len,ierr)
ccc   redefine q
c      q=r*g/f
c      gconst=(qscale*f/r)**2-g**2
c      if(qscale.eq.1.0) gconst=0.0
ccccccccccc
      nadres = 50 + 3*nosurf - 1
c replacing      call zrd ( outmap1, q1, 1, nadres, lgivup, 100 )
c      call zrd(mp1,q1,1,nadres+ishft,igivup,ierr)
c
      call getwa('mpout1',q1,nadres,len,ierr)
      nadres = nadres + nosurf * 2
c replacing      call zrd ( outmap1, g1, 1, nadres, lgivup, 100 )
c      call zrd(mp1,g1,1,nadres+ishft,igivup,ierr)
c
      call getwa('mpout1',g1,nadres,len,ierr)
      nadres = nadres + nosurf * 2
c replacing      call zrd ( outmap1, f1, 1, nadres, lgivup, 100 )
c      call zrd(mp1,f1,1,nadres+ishft,igivup,ierr)
c
      call getwa('mpout1',f1,nadres,len,ierr)
c      write(*,*) 'q,g,f,q1,g1,f1',q,g,f,q1,g1,f1
c      q1=r*g1/f1
c      g2=g1*g1
c      gnew2=g2+gconst
c      gnew=sqrt(gnew2)
c      q1=gnew/g1*q1
cccc...  oma2=(normalized alfven freq.)**2
      r2=r*r
      oma2=1./(q1*q1*r2)
cccccccccccccccccccccccccccccccccccccccccccccccccc
cccc
      nnx=nosurf
      xr=(psilim-psimin)/twopi
ccc.. define size of output file ' ioequ1 '
      length=(3*jb*mt+4)*lam+(3*jcc*mt+4)*lamc
      len0=length*nosurf+10
      len=len0+12*nosurf+5
      ndsk=1
c replacing      call zowc(outequ1,mp2,len,ndsk,iiff,100)
ckg      istat=ishell('rm -f equou1')
ckg cdeb
c      call zop(mp2,'equou1',500*512, ndsk, 0, 100 )
c
      lhhh = 0
ckg      nmap1 = 500
      nmap1 = 300000
      nmap1 = length
      call wopen('equou1',nmap1,lhhh,ierr)
c replacing      call zwr(outequ1,ntor,10,1,1,100)
cdeb        call zwr(mp2,ntor,10,1+ishft,1,100)
      lhhh = 1
c+ishft
      nadres=1
      lhhh=10
c      print *, ntor,lmin,lmax,kmin,ktot,njg,mjx,nnx,igrid,xr
c     &,ltot,ltotsq,ltots2,ltots8
      call putwa('equou1',ntor,nadres,lhhh,ierr)
c
ccc
ccc...  read TAE mode eigenfunctions and eigenfrequency ( zeroth order)
c
      lens=nosurf*(mt*3+4.*nts)+4
      is = 0
      ndsk = 1
      symm_as=1.
      len=1
      if(pest.ne.'y')then
         ioc=98
c         call zop(ioc, 'eigenf', lens+512, ndsk, 0, 100 )
c
         nmap1 = lens/512
         lhhh=0
         call wopen('eigenf', nmap1, lhhh, ierr)
ccc   ..    read equilibrium quantities
c         call zrd(ioc,omreal,1,1+ishft,igivup,ierr)
c
         lhhh=1
         call getwa ( 'eigenf', omreal, lhhh, len, ierr )
c         call zrd(ioc,eigfun,lens,2+ishft,igivup,ierr)
c         call zcl(ioc,ierr)
cdeb
c         lens=nosurf*(3*mt+4*nts)+1
c         length=lens/512+2
c         call wopen('eigenf',length,0,ierr)
         lhhh=lhhh+1
         len=nosurf*mt*3
         call getwa ( 'eigenf', eigfun, lhhh, len, ierr )
         lhhh=lhhh+len+1
         len=1
         call getwa ( 'eigenf', gamcont, lhhh, len, ierr )
         write(*,*) 'Continuum gamma/omega=',gamcont
         write(57,*) 'Continuum damping'
         write(57,*) gamcont
c         call zcl(ioc,ierr)

         call zeroo(eigfun_s(1,1,1),nn*mt*3)
cnng14 only third component is made equal for derivative calculations for E||
         call equiv(eigfun_s(1,1,3),eigfun(1,1,3),nn*mt)
c
         call wclose('eigenf',ierr)
      else
c read the eigenmode structure from PEST code
         open(unit=45,file='eigvector',status='old')
         read(45,*)omreal
         read(45,*)ntor
c         omreal=sign(sqrt(abs(omreal))*q1,omreal)
         omreal=sign(sqrt(abs(omreal)),omreal)
         read(45,*)nosurf_egv
         read(45,*)npolhrm
         do i_egv=1,npolhrm
            read(45,*)m_egv
            do i_surf=1,nosurf_egv
               read(45,*)(eigfun((i_surf-1)+1,m_egv-minm+1,j_egv),j_egv
     &              =1,3,2),eigfun((i_surf-1)+1,m_egv-minm+1,2)
            enddo
         enddo
         do i_egv=1,npolhrm
            read(45,*,end=111,err=111)m_egv
            do i_surf=1,nosurf_egv
               read(45,*,end=111,err=111)(eigfun_s((i_surf-1)+1,m_egv
     &              -minm+1,j_egv),j_egv=1,3,2),eigfun_s((i_surf-1)+1
     &              ,m_egv-minm+1,2)
            enddo
         enddo
         symm_as=0.
         goto 112
 111     call zeroo(eigfun_s(1,1,1),nn*mt*3)
         symm_as=1.
 112     continue
c            do j_egv=1,3
c This subroutine is filling the gaps in the eigenvectors as read from
c Pest output file eigenvector, which has the values for the eigenvector
c at i_r = 1,3,5...
c               call fill(eigfun(1,m_egv-minm+1,j_egv),eigfun(nosurf_egv
c     &              *2-1,m_egv-minm+1,j_egv),nosurf_egv-1)
c            enddo
      endif
ckg find the peak of the radial mode structure among all the harmonics
ckg to compute the omega_star
      isrfmax=0
      immax=0
      hhh=0.
      do i=1,nn
         do im=1,mt
            if(hhh.lt.abs(eigfun(i,im,1)))then
               hhh=abs(eigfun(i,im,1))
               isrfmax=i
               immax=im
            endif
         enddo
      enddo
      write (*,*) 'Eigf peaks at isrf=',isrfmax,' for i_m=',immax
c      omrealstar=2.4*abs(ntor)/10.
      omrealstar=0.
c  0.23*abs(ntor)/25.
ckg that was for iter
cthis 3.1 is for a17run of the shot 107001 and 3.525 is for the etai=1 imposed
c3.525*abs(ntor)/10.
c used omst correction for d3d shot 107001, typically should be zero
c precession rotation
c2.4*abs(ntor)/10.
      write(*,*)'read omreal',omreal,'consider correction omega_*pi=',
     &     omrealstar
      omreal=0.5*(omrealstar+sqrt(omrealstar**2+omreal**2*4.))
      write(*,*)"will use omreal",omreal
c      stop
ckg_pess this part is to make calculations pessimistic, 
c        which is to take only two main harmonics.
      ipess=0
      if(ipess.eq.1) call keep2mainhrmnks(eigfun(1,1,1),mt,nn)
ckg_pess
c
c        add nonuniform grid
c
      eps=1.0
      do 11 i=1,nthetaj
      rx1=(1.0/nthetaj)*i
      rx2=(2.0-eps)*rx1+(eps-1.0)*rx1**2
      rx3=(1.0/nthetaj)*(i-1.0)
      rx4=(2.0-eps)*rx3+(eps-1.0)*rx3**2
      rx(i)=(rx4+rx2)*0.25*pi
      drx(i)=(rx2-rx4)*0.5*pi
  11  continue
c
      icount=0
cnng12      nosurf11=nosurf-4
      nosurf11=nosurf
ckg      do 80 isrf=1,nosurf
c     Note that it has to start from 1 as it may interfare with reading
c     RGRID in dskin
ckg_rad
      tfxa=0.0
ckg_rad
      do 80 isrf=1,nosurf
      nsrf = isrf
c
c         write(*,*) 'before dskin'
      call dskin
c         write(*,*) 'after dskin'
      if(nsrf.eq.1)q0=q
c Here we transform the Pest eigenvector to the Nova format
      if(pest.eq.'y') call pest2nova_egv
c2d      write(*,*) ' isrf,f,fp,rrg=',isrf,f,fp,rrg,curvn(1),curvs(16)
c2d      write(*,*)  rgrid(isrf),eigfun(isrf,4-minm+1,1)
c don't use      call matrix
c2d

      call deltak(engk)
      delk(nsrf)=engk
ckg_rad
      anq(nsrf) = abs(ntor) * qprf(nsrf)
      if(nsrf.eq.1)then
         tfx(nsrf)=0.0
         w(nsrf)=0.0
      else
         w(nsrf)=delk(nsrf)*rg*(rgrid(nsrf)-rgrid(nsrf-1))
         tfxa=tfxa+qprf(nsrf)*(rgrid(nsrf)**2-rgrid(nsrf-1)**2)
         tfx(nsrf)=sqrt(tfxa)
      endif
ckg_rad
      
c
      if(isrf.gt.3) call kinetic
c      write(*, *)'Step==',5,pest
      call interfa
      call dskout
c
cccc  plot trapped-particle terms
c      if((nsrf/10)*10.ne.nsrf) go to 81
c      call frame(0)
c      call posplt2(tgrid,bf,1,mth,1,nsrf,'bf')
c      call posplt1(dtb,tboun,1,nthetaj,3,'tboun')
c      call posplt1(dtc,tcirc,1,nthetaj,4,'tcirc')
c     call frame(0)
c     call posplt2(gal,rbce(1,mt+jmin+1,1,1),1,lam,2,1,'rbce1')
c     call posplt1(gal,rbce(1,mt+jmin+1,2,1),1,lam,4,'rbce2')
c     call frame(0)
c     call posplt2(gal,rbce(1,mt+jmin+2,1,1),1,lam,2,2,'rbce1')
c     call posplt1(gal,rbce(1,mt+jmin+2,2,1),1,lam,4,'rbce2')
c     call frame(0)
c     call posplt2(galc,rbcec(1,mt+jmin+1,3,1),1,lamc,2,1,'rbcec0')
c     call posplt1(galc,rbcec(1,mt+jmin+1,4,1),1,lamc,4,'rbcec1')
c      call frame(0)
c      call posplt2(galc,tkcir,1,lamc,1,nsrf,'tkcir')
c      call posplt1(gal,tkbou,1,lam,2,'tkbou')
c      call posplt1(galc,wdc,1,lamc,3,'wdc')
c      call posplt1(gal,wd,1,lam,4,'wdt')
c     call posplt1(galc,rbcec(1,mt+jmin+1,5,1),1,lamc,4,'rbcec2')
c     call posplt1(dtb,deltb,1,nthetaj,3,'deltb')
c     call posplt1(dtc,deltc,1,nthetaj,4,'deltc')
c   81 continue
c
 80   continue
      write(*,*) 'q0,q1,omreal=',q0,q1,omreal
cccc
c
cc..  store equilibrium profiles
c
      nadres=len0+1
      lgivup=1
      nlen=10*nosurf
ctem       call putwa('equou1',qprf,nadres,nosurf,ierr)
c      call zwr(mp2,qprf,nosurf,nadres+ishft,1,100)
      nadres=nadres+nosurf
ctem      call putwa('equou1',gprf,nadres,nosurf,ierr)
c      call zwr(mp2,gprf,nosurf,nadres+ishft,1,100)
      nadres=nadres+nosurf
ctem      call putwa('equou1',pprf,nadres,nosurf,ierr)
c      call zwr(mp2,pprf,nosurf,nadres+ishft,1,100)
      nadres=nadres+nosurf
ctem      call putwa('equou1',rhoprf,nadres,nosurf,ierr)
c      call zwr(mp2,rhoprf,nosurf,nadres+ishft,1,100)
      nadres=nadres+nosurf
ctem      call putwa('equou1',rgrid,nadres,nosurf,ierr)
c      call zwr(mp2,rgrid,nosurf,nadres+ishft,1,100)
      nadres=nadres+nosurf
ctem      call putwa('equou1',hmin,nadres,nosurf,ierr)
c      call zwr(mp2,hmin,nosurf,nadres+ishft,1,100)
      nadres=nadres+nosurf
ctem      call putwa('equou1',hmax,nadres,nosurf,ierr)
c      call zwr(mp2,hmax,nosurf,nadres+ishft,1,100)
c      call zcl(mp2,ierr)
c
c      call wclose('equou1',ierr)
c
cc....
ccc...  perturbative method for computing stability criteria for TAE modes
ccc
c use im1 according to notes on DF below
c2d also uncomment KINETICS above and never comment it to avoid /0 as
ckg it defines all the defaults
      call taebeta(im1,ihsps,xfow,pk)
ckg_rad
c need to put better effective mass and elongation, but for burning plasma is OK.
      write(*,*)'   im1==================',im1
      if(im1.ne.1)then
ckg      aeff=2.5
      aeff=rho0/(denc0(2)+denc0(3)+denc0(4)+denc0(4))
      write(*,*) ' aeff=',aeff,rho0,denc0(2)+denc0(3)+denc0(4)+denc0(4)
      elong=1.8
      call elongs(elong,nosurf/2)
      vs0=9.79e5*sqrt(2./aeff)
      wcix=9.58e3*b0/aeff
      rhos0=vs0/wcix
ckg      write(*,*)tc0(1),b0,rhos0,amin,omreal
      do nsrf=1,nosurf
         w(nsrf)=rhoprf(nsrf)*w(nsrf)
         if(nsrf.eq.1)then
            alamx(nsrf)=0.0
         else
            shearx=tfx(nsrf)*(qprf(nsrf)-qprf(nsrf-1))/qprf(nsrf)
     >           /(tfx(nsrf)-tfx(nsrf-1))
            rhos=rhos0*sqrt(tc(nsrf,1))*sqrt(0.75*tc0(2)/tc0(1)+1.)
            alamx(nsrf)=2.82*(rhos*sqrt(tfxa)/(amin*tfx(nsrf)))
     &           *qprf(nsrf)*abs(shearx*ntor)
         endif
      enddo
ckg now we compute the radiative damping: --> gam_rad
      n0=2
      call radiative_damping(alamx(1),elong,gam_rad,anq(1),nn,n0,nosurf
     &     -n0,omreal**2,w(1))
ckg_rad
      endif
      close(iodat)
      close(17)
c2d       call plote
ckg **************************************************************************
ckg change See F's comments on the convention for particle indexes 
c      ihsps=1
c
      iter=0
c
      ntgr=20
c these are for NSTX or other ST
c      iter=2
c      ntgr=30
c
      b0oo=b0/10000.
      ai=rmhp(ihsps)
      scalvavf=1.
      symm=1.
c      symm=0.
      vi=1.298
      vi=sqrt(th0(ihsps)/(1.e+6)*4./3.52/rmhp(ihsps))*vi
      zi=zh(ihsps)
c The variable tip is related to the distribution function type.
c It is used in the subroutine distrifun (file orbit2d.f)
c!!!!!!!it is read two times from NOVAK_param file since some variables 
c       are not in the common blocks
c
c     's' for Slowing down distribution, also
c     'd' for Slowing down distribution with v_cr = v_a0 / 3.
c     'e' Maxwellian for test with jet equilibrium 
c         (It'll switch to im1= 0, e.g.TAE's automatically)
c     't' slowing down for test with tftr equilibrium
c     'j' slowing down for test with jet equilibrium
c     'm' Maxwellian with standard input
c     'x' Maxwellian with thermal ions for proper Landau damping calculations
c   Below is for Gaussian pitch distribution centered at p=P0GA
c   and having width delta_p = DPGA, p=mu*Baxis/E (sometimes in the
c   literature refered as =lambda). Add im1=3 for better treatment 
c   of the plasma rotation. 
c
c     'c' co passing only with slowing down distribution
c     'h' ICRF created minority with Maxwell velocity
c         f ~ exp(-(p-chidelt(1))^2/chidelt(2)^2), so that
c         trapped ion with bounce tip at the magnetic axis has p=1. For
c         example for ICRH at res.layer at R, chidelt(1)=R/Raxis, 
c         chidelt(2)=dR/Raxis, where dR is width of the res.layer
c         chidelt(i) array is from the NOVA_param file under the chi0 name
c	'h' He3 minority related input works if in NOVA_param 
c	    (ihsps=2 is not important, any #, 1,2,3 will do) ai=3, zi=2, and tip='h'.
c          checked for on sparc_12345I64 with the information read from transp.dat file
c          by variables
c           tmini_3     HE3 ICRF MINORITY 2/3<E>         EV              :X 
c           nmini_3     HE3 ICRF MINORITY DENSITY        EV              :X 
c     'i' is similar to 'h' option for ICRF but with H temperature read 
c         from transp.dat file: if tip(1:1).eq.'i'.and.int(szh).eq.1
c           tmini_h     H ICRF MINORITY 2/3<E>  
c           nmini_h     H ICRF MINORITY DENSITY 
c         Those two option 'h','i' work so far with transp input, i.e. itransp=1 in NOVA_param.
c     'r' this symbol is for special distribution function for START shot
c         #35305 at 26.3 ms for the form see subroutine DISTRFUN
c     'g' (i recommend to use 'l' option, which is better. Lorentz model
c         'g' option is not good for really wide distributions. )
c         is for general form in Gaussian pitch angle with velocity
c         dependent width, so that parameters are read from NOVAK_param
c         file with slowing down in velocity. 
c         f ~ exp(-(chi-chidelt(1))^2/dchi^2), where chi = v_||/v on the LFS
c         dchi^2=chidelt(2)+sqrt(psi)*chidelt(3)-
c                        |((v/v0)^3/chidelt(5))*(1.+chidelt(6)**3)|
c          -chidelt(4)*ln|----------------------------------------|
c                        |  (v/v0)^3/chidelt(5)+chidelt(6)**3     |
c         Simple form with constant width is if chidelt(3)=chidelt(4)=0. 
c     'l' slowing down in velocity, Gaussian in pitch angle, chi=v_||/v
c         initially, i.e. at the injection energy. Later, i.e. at lower
c         energies, the spread in pitch angle is calculated according 
c         to Gorelenkov-Berk distribution function model from NF'05. 
c         by expanding the pitch angle distribution in Legendre
c         polynomials and evolving it as velocity decreases. 
c         For input use chidelt(1) as central value, sqrt(chidelt(2)) as the width,
c         which are taken at the injection energy. 
c         In this case f ~ exp(-(chi-chidelt(1))^2/chidelt(2))
c     'b' is like 'l' but with the fractional contributions from full, half 
c         and third energy sources
c     'q' is the distribution function to be read from the CQL produced
c         file. It is called cql.dat and is selfcontained. In case new file
c         is to be read it is advised to make a soft link to it.
c         In this case the temperature of this specie is read from cql3d
c         data file and supersedes the transp or analytical values. 
c         It is also important that its maximum value is normalized to 
c         the NOVAK_param file value for particle energy to be more flexible 
c         in physics studies. 
c         Particle mass is read from cql3d file too - set up in cql3dinit.f90.
c         In case of itransp=0(analytical profiles) betah0 is read from NOVAK_param
c         Particle charge is set to 1 in the same subroutine.
c
c     'o' is to skip fast ions with FOW effects
c
c     fast ion beta profile is generated from perpendicular and parallel
C     energy densities if the data is read from TRANSP (itransp=1). 
c       It is determined by (UBPRP_D/2+UBPAR_D) for deuterium beams
c                           (UBPRP_T/2+UBPAR_T) for tritium beams
c                           (UFIPP/2+UFIPA) for fusion ions
c                           ((UFASTPP-UBPRP_D-UBPRP_T-UFIPP)/2+(UFASTPA-UBPAR_D-UBPAR_T-UFIPA))
c                           for H-minority fast ions if ihsps=3 and its mass set to 1
c     Type of distribution function ('tip') can be changed in NOVAK_param
c     file, which is generated automatically at first start
c     Variable im1 is 0 for normal stability runs of *AE
c                     1     ideal kink modes with w=0
c                     3 is like 0 but with the rotation included
c                     10is a special run for EP redistribution due to TAE activity
      p0ga=chidelt(1)
      dpga=chidelt(2)
c
      if(symm_as.ne.symm) then
         write(*,*) ' --------- symmetry warning'
         write(*,*) ' Your eigenmode structure file has symm=',symm_as
         write(*,*) ' While you set symm=',symm,'! Make it consistent!'
         write(*,*) ' --------- symmetry warning'
      endif
      nosurf=nosurf11
      if(tip.eq.'o')stop 'skipping fast ion FOW part'
c      tip='b'
      istart=1
      call orbit2d(ai,b0oo,deltk,dpga,gamkom,gamkomflr,im1,iter,immax
     &     ,isrfmax,ntgr,omkom,omkomflr,p0ga,pk,scalvavf,symm,tip,vi
     &     ,xfow,zi,istart)
      write(*,*) '1 gam/om,gamflr/om,omr/om,omrflr/om',gamkom,gamkomflr
     &     ,omkom,omkomflr
c      call system("rm gam100.100d")
c      tip='b'
c      chidelt(1)=0.2205
c      chidelt(2)=0.1484
c      istart=1
c      call orbit2d(ai,b0oo,deltk,dpga,gamkom,gamkomflr,im1,iter,immax
c     &     ,isrfmax,ntgr,omkom,omkomflr,p0ga,pk,scalvavf,symm,tip,vi
c     &     ,xfow,zi,istart)
c      write(*,*) '2 gam/om,gamflr/om,omr/om,omrflr/om',gamkom,gamkomflr
c     &     ,omkom,omkomflr
ckg19 radiative damping of RSAEs by Y.Wang
ckg
      if(itransp.eq.1) call RD_RSAE(immax)
      if(im1.eq.1) then
ckg analytical estimates of fast ion contribution to m=1/n=1 stabilization
ckg see w_tw_s program in tae.f file for comments on how it is done
c
         call w_tw_s(b0,deltk,ihsps,omreal,omkom,omkomflr,w2_t,w2_s
     &        ,rdw_k,rdw_kflr)
         write(*,*) ' For M=1 stability analysis used w2_t,w2_s=',w2_t
     &        ,w2_s
         write(*,*) '------------------------'
      endif
close(57)
      stop
  100 call errmes ( iodat, 4hmain )
  300 format ( i1 )
c
      end
c************************************************************************
c calculate element for collisional trapped electron damping
c Here AL is pitch angle
      subroutine bouncee(csss,al,eta,etap,k,tpnt1,tpnt2,tkbouc)
c
      include 'clich1'
      include 'clich2'
      common/functn/omreal,eigfun(nn,mt,3),eigfun_s(nn,mt,3)
      common/rbccol/abe(nn,lam),abp(nn,lam),abe_av(nn,lam),abe_avav(lam
     &     ,nn)
      common/coefs/c1(nts),c2(nts),c3(nts),c4(nts),igm1(nts)
     &,ig(nts),igp1(nts),igp2(nts)
      common/cossin/cost1(nts),sint1(nts),cost2(nts),costb(nts)
     &,sintb(nts),cost(nts),sint(nts),cost3(nts),sint3(nts)
     &,cost1p(2*nts),sint1p(2*nts),cost2p(2*nts)
     &,cost1m(2*nts),sint1m(2*nts),cost2m(2*nts)
     &,costp(nts),costm(nts),sintp(nts),sintm(nts)
      dimension tkbouc(lam),rbcep(lam,mt,3)
c
cc..  xx2 is the radial scale factor
      real nq
      nq=abs(ntor*q)
      xx2=psitot*2.*rrg
c
      bf0=0.
      bf0d=0.0
c      wd(k)=0.
      tpp=tpnt2+tpnt1
      tpm=tpnt2-tpnt1
      wj=pi/nthetaj
      wj0=2.0/twopi
ccc.. bounce average of magnetic drift freq. wd(k)
ccc   integrate from 0 to tpnt2 due to up-down symmetry over nthetaj points
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
c This is xjacob(i)*bf(i)
      at(i)=c1(i)*ss3(igm1(i))+c2(i)*ss3(ig(i))
     &     +c3(i)*ss3(igp1(i))+c4(i)*ss3(igp2(i))
      ss2(i)=c1(i)*ss1(igm1(i))+c2(i)*ss1(ig(i))
     &     +c3(i)*ss1(igp1(i))+c4(i)*ss1(igp2(i))
      if(ss2(i).lt.0.) ss2(i)=-ss2(i)
c      ss4(i)=c1(i)*bf(igm1(i))+c2(i)*bf(ig(i))
c     &     +c3(i)*bf(igp1(i))+c4(i)*bf(igp2(i))
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
      grpp(i)=c1(i)*grpssq(igm1(i))+c2(i)*grpssq(ig(i))
     &     +c3(i)*grpssq(igp1(i))+c4(i)*grpssq(igp2(i))
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
ckg Note, that ss2 here is v_||^2/v^2
      do 25 i=1,nthetaj
      gt1x=sqrt((dt(i)-tpnt1)*(tpnt2-dt(i))*ss2(i))
c      gt1(i)=sqrt((dt(i)-tpnt1)*(tpnt2-dt(i))/ss2(i))
      gt1(i)=gt1x/ss2(i)
      gt1i(i)=gt1(i)*at(i)
      gt1id=gt1x*at(i)
c      wd(k)=wd(k)+gt1i(i)*(fwk(i)+fwb(i)*al)*drx(i)
      bf0d=bf0d+gt1id*drx(i)
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
      tkbouc(k)=bf0*2.0
      tdbou(k)=bf0d*2.0
      wjb=1.0/bf0
c      wd(k)=wd(k)*wjb
      wjb1=wjb
c
ccc..  bounce average of quantities involving exp(i*m*th)
c
ccc..  j=1 is for p=0 bounce term
c	j=1
c      jm1=j-1
c      if(nsrf.eq.50)then
c      write(16,*)"grpp profile"
c      write(16,*)(i,grpp(i),i=1,nthetaj)
c      end if
      do 32 l=1,mt
         l2=l+minm-1
         alq=l2-ntor*q
c     
         do i=1,3
            rbce(k,l,1,i)=0.0
            rbcee(k,l,i)=0.0
         enddo
         rbce(k,l,2,1)=0.0
         rbcep(k,l,1)=0.0
c     
         do 33 i=1,nthetaj
            cost(i)=cos(alq*dt(i))
            sint(i)=sin(alq*dt(i))
            cost1(i)=cost(i)*drx(i)
            sint1(i)=sint(i)*drx(i)
            cost2(i)=gt1i(i)/ds1(i)*cost1(i)
ckg Utilization
            sint3(i)=cost(i)*gt1(i)*rg*drx(i)*real(l2)/grpp(i)
            sint(i)=cost2(i)/ds1(i)
            costb(i)=gt1i(i)*cost1(i)
            cost1(i)=al*ap(i)*costb(i)
ckg Utilization
            costb(i)=gd(i)*costb(i)
            sintb(i)=gt1i(i)*sint1(i)
            sint1(i)=al*as(i)*sintb(i)
ckg Utilization
            sintb(i)=gs(i)*sintb(i)
c            sint1(i)=gt1i(i)*(gs(i)+al*as(i))*sint1(i)
c            cost1(i)=gt1i(i)*(gd(i)+al*ap(i))*cost1(i)
c            sint1(i)=gt1i(i)*(gs(i)+al*as(i))*sint1(i)
c     kg here ds1 = |B|
c     drx  delta theta
c     grpp is |grad psi|^2
c     gt1i is dt
c     wjb1 is 1/tau_b
            cost3(i)=grpp(i)*gt1i(i)*cost(i)/ds1(i)*drx(i)
 33      continue
         do 34 i=1,nthetaj
            rbce(k,l,1,1)=rbce(k,l,1,1)+cost1(i)+costb(i)
            rbce(k,l,1,3)=rbce(k,l,1,3)+sint1(i)+sintb(i)
            rbce(k,l,1,2)=rbce(k,l,1,2)+cost2(i)
c
            rbcee(k,l,1)=rbcee(k,l,1)+costb(i)
            rbcee(k,l,3)=rbcee(k,l,3)+sintb(i)
            rbcee(k,l,2)=rbcee(k,l,2)+sint(i)
c
            rbce(k,l,2,1)=rbce(k,l,2,1)+sint3(i)
            rbcep(k,l,1)=rbcep(k,l,1)+cost3(i)
 34      continue
         do i=1,3
            rbce(k,l,1,i)=rbce(k,l,1,i)*wjb1
            rbcee(k,l,i)=rbcee(k,l,i)*wjb1
         enddo
         rbce(k,l,1,2)=rbce(k,l,1,2)*al
c
         rbce(k,l,2,1)=rbce(k,l,2,1)*wjb1
         rbcep(k,l,1)=rbcep(k,l,1)*wjb1
c      
c            if(l.eq.2.and.abs(q-1.3).lt.0.01)
c     +         write(*,*)'new rbce',
c     +         (rbce(k,l,1,jj),jj=1,1),tpnt1,tpnt2,al
c
c     
 32   continue
c     
      abe(nsrf,k)=0.0
      abp(nsrf,k)=0.0
      abe_av(nsrf,k)=0.0
      if(nsrf.eq.1) go to 98
      isf0=max(2,nsrf-1)
      isf0=min(nosurf-3,isf0)
      if(nsrf.gt.1) then
         sss=0.5/rgrid(nsrf)
      else
         sss=0.0
      end if
      do 40 l=1,mt
cnng14 smooth an extra saved array for further derivations
        if(nsrf.eq.4)then
           call smooth(eigfun_s(1,l,3),nn)
           call smooth(eigfun_s(1,l,3),nn)
           call smooth(eigfun_s(1,l,3),nn)
        endif
        call der4(rgrid(isf0),eigfun_s(isf0,l,3),rgrid(nsrf),dyn)
        do ijb=1,jb
        abpe(k,l,nsrf,ijb,1)=csss*rbce(k,l,ijb,1)*eigfun(nsrf,l,1)
ckg this is parallel electric field term, may give large contribution to
ckg the damping if there is singularity
     &       +rbcep(k,l,1)*dyn*sss
c
c define v.E terms for electrons
        enddo
        abp(nsrf,k)=abp(nsrf,k)+abpe(k,l,nsrf,2,1)
c        if(nsrf.eq.46.or.nsrf.eq.45.or.nsrf.eq.47) then
c           print *,'nsrf,l,dyn',nsrf,l,dyn,abpe(k,l,nsrf,2,1)
c        endif
        do ijb=1,jb
        habp=0
        habe=0
	do j=1,3
           habp=habp+rbce(k,l,ijb,j)*eigfun(nsrf,l,j)
           habe=habe+rbcee(k,l,j)*eigfun(nsrf,l,j)
        enddo
        abpe(k,l,nsrf,ijb,2)=habp
        abpe(k,l,nsrf,ijb,3)=habe
        enddo
        abe(nsrf,k)=abe(nsrf,k)+abpe(k,l,nsrf,1,2)
        abe_av(nsrf,k)=abe_av(nsrf,k)+abpe(k,l,nsrf,1,3)
 40   continue
ckg note that abpe is not used right after this subroutine. It is defined in the bounce 
      if(k.ne.lam)then
         abe_avav(k,nsrf)=abe_avav(k+1,nsrf)+abe(nsrf,k)*tkbouc(k)*eta
     &        +abe(nsrf,k+1)*tkbouc(k+1)*etap
      else
         abe_avav(k,nsrf)=0.
      endif
c
c abp is the parallel electric field term
c abe is the curvature term
c
c      if(nsrf.gt.1) then
c         abp(nsrf,k)=abp(nsrf,k)*0.5/(rgrid(nsrf))
c      else
c         abp(nsrf,k)=0.0
c      end if
   98 continue
c
      return
      end
c************************************************************************
c Calculate total parallel electric field includes FLR and distribution averaging
      subroutine E_parallel(abe_avav,abe_av,abe,abp,abpe,csss,galx,gale
     &     ,jb,lam,mt,nn,nsrf,salmc,sm,smc,tkbouc,tkcir,tpnt1c,tpnt2c)
      common/colli2/rbcee(100),sst1(100)
      dimension abe_avav(lam,nn),tkbouc(lam),tkcir(lam)
     &     ,salmc(lam),galx(lam),gale(lam),tpnt1c(lam),tpnt2c(lam)
     &     ,abe(nn,lam),abe_av(nn,lam),abp(nn,lam),abpe(lam,mt,nn,jb,3)
      do i=1,lam
         sst1(i)=0.
      enddo
      sst=0.
      sss=0.
c      write(*,*)'So far1',nsrf
      do i=lam,1,-1
         sss=sss+salmc(i)*tkcir(i)
         if(i.lt.lam)then
            sst=sst+(galx(i)+galx(i+1))
            sst1(i)=sst1(i+1)+(galx(i)*tkbouc(i)+galx(i+1)*tkbouc(i
     &           +1))
         else
            sst1(i)=0.
         endif
      enddo
c     Number of trapped particles at \theta = 0
      do i=1,lam
         sst1(i)=sst1(i)*sm/sst
      enddo
c      write(*,*) sst1(1)/sss,sqrt(sm/smc),sst1(lam/2)/sss
c      sst=sst*sm/float(lam)
      hhh=sm*0.75/(float(lam))
c     Number of passing particles is SSS
      do i=1,lam
         abe_avav(i,nsrf)=abe_avav(i,nsrf)*hhh/sss
      enddo
      do i=1,lam
         call avav(gale(i),galx(i),i,tpnt1c,tpnt2c,tkbouc)
      enddo
      do i=1,lam
         if(csss.lt.1.e-5)then
            abe_avav(i,nsrf)=0.
            rabp=0.
         else
            abe_avav(i,nsrf)=rbcee(i)/abe_av(nsrf,i)
            rabp=-(sss+sst1(i))/sss
         endif
         abe_av(nsrf,i)=(rabp+abe_avav(i,nsrf))*abe_av(nsrf,i)
         abp(nsrf,i)=abp(nsrf,i)*rabp
c         do imt=1,mt
c            do ijb=1,jb
c            abpe(i,imt,nsrf,ijb,1)=rabp*abpe(i,imt,nsrf,ijb,1)
c            abpe(i,imt,nsrf,ijb,3)=(rabp+abe_avav(i,nsrf))
c     +         *abpe(i,imt,nsrf,ijb,3)
c            enddo
c         enddo
      enddo
      return
      end
c************************************************************************
c Here we compute the third averaged term for parallel electric field
c which enters through the perturbed trap.elec. density
      subroutine avav(al,eta,k,tpnt1c,tpnt2c,tkbouc)
c
      include 'clich1'
      include 'clich2'
      common/rbccol/abe(nn,lam),abp(nn,lam),abe_av(nn,lam),abe_avav(lam
     &     ,nn)
      common/coefs/c1(nts),c2(nts),c3(nts),c4(nts),igm1(nts)
     &,ig(nts),igp1(nts),igp2(nts)
      dimension tkbouc(lam),tpnt1c(lam),tpnt2c(lam)
c
      bf0=0.
c
      tpp=tpnt2c(k)+tpnt1c(k)
      tpm=tpnt2c(k)-tpnt1c(k)
ccc.. bounce average of magnetic drift freq. wd(k)
ccc   integrate from 0 to tpnt2 due to up-down symmetry over nthetaj points
cccc  add tboun(theta), the bounce orbit time for trapped particle
      do 22 i=1,nthetaj
         dt(i)=0.5*(tpp+tpm*sin(rx(i)))
   22 continue
c
ccc  set up coefficients for 4-point Lagrange interpolation
      do 11 i=1,nthetaj
         call fx4p(dt(i),i,mth2,dth)
   11 continue
c
      do 1 is=1,mth2
         ss1(is)=1.-al*bf(is)
    1 continue

      do 23 i=1,nthetaj
         at(i)=c1(i)*ss3(igm1(i))+c2(i)*ss3(ig(i))
     &        +c3(i)*ss3(igp1(i))+c4(i)*ss3(igp2(i))
         ss2(i)=c1(i)*ss1(igm1(i))+c2(i)*ss1(ig(i))
     &        +c3(i)*ss1(igp1(i))+c4(i)*ss1(igp2(i))
ckg Note, that ss2 here is v_||^2/v^2
         if(ss2(i).lt.0.) ss2(i)=-ss2(i)
         gt1i(i)=at(i)*sqrt((dt(i)-tpnt1c(k))*(tpnt2c(k)-dt(i))/ss2(i))
c     
 23   continue
      do 26 i=1,nthetaj
         bf0=bf0+gt1i(i)*drx(i)*gimme_abe(abe_avav(1,nsrf),dt(i),lam
     &        ,tpnt2c)
 26   continue
      rbcee2(k)=bf0*2.0/tkbouc(k)
c      write(*,*) (abe_avav(i,nsrf),i=1,lam)
c
c      tkbouc(k)=bf0*2.0
      return
      end
c************************************************************************
c Note plus minus symmetriya predpolagaetsya
      function gimme_abe(abe_avav,dt,lam,tpnt2c)
      dimension abe_avav(lam),tpnt2c(lam)
      i=indx(tpnt2c,lam,dt)-1
      i=max(i,1)
      i=min(lam-3,i)
      call fun4(tpnt2c(i),abe_avav(i),dt,gimme_abe)
      return
      end
c************************************************************************
      subroutine inparam_read(im1,ihsps,xfow,pk)
c
cccc.....   Create data file with parameters to be used in the run
ckg         If file does not exist create it with deafault parameters
ckg         from previously executed subroutine inparam
ckg         To add new variable you can put it in any order. The default
ckg         value will be used, while old variables would get values
ckg         from old NOVAK_param file
c
      include 'clich1'
      include 'clich2'
      include 'clich1b'
      common/plasmap/rmaj,amin,b0
      character defin*72,sample*40
      logical klyu(2)
      dimension wrk(3)
      data io/23/,defin/' '/,sample/' '/,klyu/2*.true./
      open(io,file='NOVAK_param',status='unknown')
 1    i1=2
c-----------------------------
      defin='Toroidal mode number'
      sample='ntor'
      i2=i1+3
      hhh=real(ntor)
      nwrk=1
      call read_vec(defin,i1,i2,io,klyu,sample,hhh,nwrk)
      ntor=int(hhh)
c-----------------------------
      defin='Key=1=>TRANSP.dat,=0 =>anal.prfs,=2=>anal+ T(r)<=p(r)/n(r)'
      sample='itransp'
      i2=i1+6
      hhh=real(itransp)
      call read_vec(defin,i1,i2,io,klyu,sample,hhh,1)
      itransp=int(hhh)
c-----------------------------
      defin='Major rad. of geom. center [cm] if(itransp|=1)'
      sample='rmaj'
      i2=i1+3
      call read_vec(defin,i1,i2,io,klyu,sample,rmaj,1)
c-----------------------------
      defin='Minor rad. of last surface [cm]'
      sample='amin'
      i2=i1+3
      call read_vec(defin,i1,i2,io,klyu,sample,amin,1)
c-----------------------------
      defin='Vacuum mag. field at geom. center [Gauss]'
      sample='B0'
      i2=i1+1
      call read_vec(defin,i1,i2,io,klyu,sample,b0,1)
c-----------------------------
      defin='Electron density [cm^-3] at magnetic axis'
      sample='dn_e'
      i2=i1+3
      call read_vec(defin,i1,i2,io,klyu,sample,denc0(1),1)
c-----------------------------
      defin='Electron temperature [eV] at magnetic axis'
      sample='T_e'
      i2=i1+2
      call read_vec(defin,i1,i2,io,klyu,sample,tc0(1),1)
c-----------------------------
      defin='ICTs central density [cm^-3] for D,T,H,C'
      sample='dn_i'
      i2=i1+3
      call read_vec(defin,i1,i2,io,klyu,sample,denc0(2),ict-1)
c-----------------------------
      defin='ICTs central temperature [eV] for D,T,H,C'
      sample='T_i'
      i2=i1+2
      call read_vec(defin,i1,i2,io,klyu,sample,tc0(2),ict-1)
c-----------------------------
      defin='Central fast particle beta: D,T,alphas'
      sample='betah0'
      i2=i1+5
      call read_vec(defin,i1,i2,io,klyu,sample,betah0(1),3)
c-----------------------------
      defin='Fast particle energy [eV]: D,T,alphas'
      sample='energyh'
      i2=i1+6
      call read_vec(defin,i1,i2,io,klyu,sample,th0(1),3)
c      do i=1,3
c         ehc0(i)=th0(i)/9.
c      enddo
c-----------------------------
      defin='1st parameter for plasma density'
      sample='alphar'
      i2=i1+5
      call read_vec(defin,i1,i2,io,klyu,sample,alphar,1)
c-----------------------------
      defin='2nd parameter for plasma density'
      sample='prho'
      i2=i1+3
      call read_vec(defin,i1,i2,io,klyu,sample,prho,1)
c-----------------------------
      defin='3rd parameter for plasma density'
      sample='arho'
      i2=i1+3
      call read_vec(defin,i1,i2,io,klyu,sample,arho,1)
c-----------------------------
      defin='Fast particle mass'
      sample='rmhp'
      i2=i1+3
      call read_vec(defin,i1,i2,io,klyu,sample,rmhp(1),3)
c-----------------------------
      defin='Fast particle charge'
      sample='zh'
      i2=i1+1
      call read_vec(defin,i1,i2,io,klyu,sample,zh(1),3)
c-----------------------------
      defin='Fast particle index 1 to 3'
      sample='ihsps'
      i2=i1+4
      hhh=real(ihsps)
      call read_vec(defin,i1,i2,io,klyu,sample,hhh,1)
      ihsps=int(hhh)
c-----------------------------
      defin=
     &'alphart,prhot,arhot:Tthermal~(1-alphart*rsq**prhot)**arhot'
      sample='alphart'
      i2=i1+6
      wrk(1)=alphart
      wrk(2)=prhot
      wrk(3)=arhot
      call read_vec(defin,i1,i2,io,klyu,sample,wrk,3)
      alphart=wrk(1)
      prhot=wrk(2)
      arhot=wrk(3)
c-----------------------------
c      defin='alpharh,prhoh,arhoh:bet_h~(1.-alpharh*Psi^prhoh)^arhoh'
      defin=
     &'alpharh,prhoh,arhoh:bet_h~exp(-|(rsq^1/2-alpharh)/arhoh|^prhoh)'
      sample='alpharh'
      i2=i1+6
      wrk(1)=alpharh
      wrk(2)=prhoh
      wrk(3)=arhoh
      call read_vec(defin,i1,i2,io,klyu,sample,wrk,3)
      alpharh=wrk(1)
      prhoh=wrk(2)
      arhoh=wrk(3)
c-----------------------------
      defin='3D distr. funct. params'
      sample='chi0'
      i2=i1+3
      call read_vec(defin,i1,i2,io,klyu,sample,chidelt(1),6)
c-----------------------------
      defin='Finite Orbit Width: 0< xfow <1'
      sample='xfow'
      i2=i1+3
      call read_vec(defin,i1,i2,io,klyu,sample,xfow,1)
c-----------------------------
      defin='Key to m1 stability calculations Im1'
      sample='im1'
      i2=i1+2
      hhh=real(im1)
      call read_vec(defin,i1,i2,io,klyu,sample,hhh,1)
      im1=int(hhh)
c-----------------------------
      defin=
     & 'Velocity grid exponent: pk=1 equidistant in v or 0.5 in energy'
      sample='pk'
      pk=1.
      i2=i1+1
      call read_vec(defin,i1,i2,io,klyu,sample,pk,1)
c-----------------------------
      defin='Distr.function type s-slow.down,m -mxwll.(see taem.f)'
      sample='_dtype'
      i2=i1+5
cnng15 tip below can be made a constant like 100.
      call read_vec(defin,i1,i2,io,klyu,sample,tip,1)
c
      if(.not.klyu(2).and.klyu(1))then
         rewind(io)
         klyu(1)=.false.
         goto1
      endif
      close(io)
      if(.not.klyu(1)) stop 'check file NOVAK_param and run again'
      return
      end
c**************************************************
c this subr. read data from IO unit number if klyu(1)=.true.
c or writes in it if klyu(1)=.false. 
c it returns klyu(2)=.false. if somethings is wrong 
      subroutine read_vec(dumb0,i1,i2,io,klyu,sample,work,nwork)
      dimension work(nwork)
      character dumb*72,dumb0*72,sample*40
      logical klyu(2)
      character tip*1
      common/parm1/np,npc,ipbmax,ipcmin,ipcmax,idet
     &,nnsurf,klamb,klamc,ishft,tip
c
      if(.not.klyu(1))goto10
      rewind(io)
      i=i2-i1+1
 1    read(io,'(72a)',end=9) dumb
      if(dumb(i1:i2).eq.sample(1:i)) then
         if(sample(1:1).eq.'n'.or.sample(1:1).eq.
     &        'k'.or.sample(1:1).eq.'i') then
c            read(io,'(7i7)',end=9)(work(i),i=1,nwork)
            read(io,*,end=9)(work(i),i=1,nwork)
            write(*,'(40a)')(' ',i=1,i1-1),sample(1:i2-i1+1)
c            write(*,'(7i7)')(work(i),i=1,nwork)
            write(*,*)(work(i),i=1,nwork)
         else if(sample(1:1).eq.'_')then
            read(io,'(7a)',end=9)(work(i),i=1,nwork)
            write(*,'(40a)')(' ',i=1,i1-1),sample(1:i2-i1+1)
            write(*,'(7a)')(work(i),i=1,nwork)
         else
c            read(io,'(7e10.3)',end=9) (work(i),i=1,nwork)
            read(io,*,end=9) (work(i),i=1,nwork)
            write(*,'(40a)')(' ',i=1,i1-1),sample(1:i2-i1+1)
            write(*,'(7e10.3)') (work(i),i=1,nwork)
         endif
      else
         goto 1
      endif
      return
 9    if(klyu(1))then
         klyu(2)=.false.
         return
      endif
 10   write(*,*)'!!default|old value of ',sample(1:i2-i1+1)
     &     ,' is written'
      write(io,*)'  ',dumb0
      write(io,'(40a)')(' ',i=1,i1-1),sample(1:i2-i1+1)
      if(sample(1:1).eq.'n'.or.sample(1:1).eq.'k'.or.sample(1:1).eq.'i')
     &     then
         write(io,'(7i4)')(int(work(i)),i=1,nwork)
c         write(io,'(7e10.3)')(work(i),i=1,nwork)
      else if(sample(1:1).eq.'_')then
         write(io,'(7a)')(work(i),i=1,1)
      else
         write(io,'(7e10.3)') (work(i),i=1,nwork)
      endif
      return
      end
c************************************************************************
ckg_rad this is subroutine written by G. Fu.
      subroutine radiative_damping(al,elong,gam_rad,anq,nn,n1a,n2a,om,w
     &     )
      parameter(nrad=400)
      dimension eps(100),g(100),wint(100),alam(100),eps1(100)
      dimension anq(nn),al(nn),px(nrad),w(nn),psi(100),qx(nrad)
      dimension w1(nrad),w2(nrad),imax(100),igap(100)
ckg      open(99,file='radfile',status='old')
      if(nn.gt.nrad)stop'increase dimension of nrad in radiat._damp.'
      rave1=0.0
      open(37,file='datcon',status='old',err=400)
c
c here om is TAE eigenvalue from NOVA code (omega**2)
c
ckg      write(6,*)"input om"
ckg      read(5,*)om
c
c     here w1 and w2 are first and second continuum curves as
c     calculated by continuum code (in omega**2) 
c     (w1 and w2 should be non-zero and finite)
c
      read(37,*)n1,n2
      do i=n1,n2
      read(37,*)w1(i),w2(i)
      enddo
      close(37)
c
cc
c  take output from NOVA-K code
c where anq(i) = ntor * qprf(i)
c
c w(i) is a local weigh (xi**2) calculated like this in NOVA-K:
c
c      call deltak(engk)
c      delk(nsrf)=engk
c      w(isrf)=rhoprf(isrf)*delk(isrf)*rg*drg
c
c al=alamx(is) is calculated like this in NOVA-K code:
c
c       vs0=9.79e5*sqrt(2.*te0/aeff)
c	rhos0=vs0/wcix
c       rhos=rhos0*sqrt(teprf(is))*sqrt(0.75*ti0/te0+1.)
c	amin=0.32*rmaj
c       if(is.eq.1)then
c       shearx=0.0
c       alamx(1)=0.0
c	else
c       shearx=tfx(is)*(qprf(is)-qprf(is-1))/qprf(is)
c     >            /(tfx(is)-tfx(is-1))
c  alamx(is)=2.82*(rhos/(amin*tfx(is)))*qprf(is)*abs(shearx)*ntor
c      end if
c
c where amin is the minor radius assuming a/R=0.32 !!!
c where shearx is magnetic shear, tfx is the normalized minor radius
c calculated from the toroidal flux tfxa:
c	tfxa=0.0
c      do 84 isrf=2,nosurf
c      tfxa=tfxa+qprf(isrf)*(rgrid(isrf)**2-rgrid(isrf-1)**2)
c      tfx(isrf)=sqrt(tfxa)
c 84   continue
c      do 85 isrf=2,nosurf
c      tfx(isrf)=tfx(isrf)/sqrt(tfxa)
c 85   continue
c     tfx(1)=0.0
ckg
ckg       read(99,*)n1a,n2a
ckg       do i=n1a,n2a
ckg       read(99,*)anq(i),w(i),al(i)
ckg       enddo
c
c
c
c determine local maximum of the first continuum curve
c >> i.e., the lower bound of each gap
      ii=0
      do i=n1,n2-2
      if(w1(i).lt.w1(i+1).and.w1(i+2).lt.w1(i+1))then
      ii=ii+1
      imax(ii)=i+1
ckg      write(*,*)"ii,i",ii,i,"w1=",w1(i+1),"w2=",w2(i+1)
      endif
      enddo

c     determine if it is a gap and if om is inside the gap
c
      iii=0
      do i=1,ii
      im=imax(i)
      iq=anq(imax(i))
      dq=anq(imax(i))-float(iq)-0.5
      fom=(om-w1(im))*(w2(im)-om)
ckg      write(*,*)"im,iq,dq,fom",im,iq,dq,fom
      if(abs(dq).lt.0.2.and.fom.gt.0.)then
      iii=iii+1
      igap(iii)=im
c
      om1=w1(im)
      om2=w2(im)
      eps(iii)=(om2-om1)/(om1+om2)
      eps1(iii)=(om2-om1)/(2.0*om)
      g(iii)=-1.0+2.0*(om-om1)/(om2-om1)
      alam(iii)=al(im)/eps(iii)**1.5
ckg      write(*,*)eps(iii),g(iii),alam(iii)
      endif
      enddo
c
c     calculate weight for each gap
c
      do i=1,iii
      im=igap(i)
      iq=anq(im)
      wint(i)=0.0
      do j=n1a,n2a
      fs=(anq(j)-float(iq))*(float(iq+1)-anq(j))
      if(fs.gt.0.0)then
      wint(i)=wint(i)+w(j)
      endif
      enddo
      enddo
c
c     
ckg      write(*,*)"eps(i),alamd,xg1,wint(i)"
      rave=0.0
      rave1=0.0
ckg      write(6,*)"input elongation parameter"
ckg      read(5,*)elong
c
c     xg1 is the local damping at each gap
c

      do i=1,iii
         delt=1.0-g(i)**2
         alamd=alam(i)/sqrt(elong)
         axi=ax(g(i))
         xg1=delt*exp(-2.0*axi/alamd)*0.5*eps1(i)
         rave=rave+wint(i)*xg1
         rave1=rave1+wint(i)
ckg      write(*,*)eps(i),alamd,xg1,wint(i)
 900     format(1x,4(2x,1e12.6))
      enddo
 400  continue
      if(rave1.gt.0.00000001)then
         rave=rave/rave1
      else
         rave=0.
      endif
      write(*,*)"global radiative damping gam/om =",-rave
      write(*,*)"at elongation",elong
      write(57,*)"radiative damping for modes in the TAE gap"
      write(57,*) -rave
      return
      end

	function ax(g0)
	pi=3.14159
	n=100
	ax=0.0
	dh=pi/float(n)/2.0
	do 100 i=1,n
	th=dh*float(i-1)
	y=(sin(th))**2
	ax=ax+cos(th)*sqrt(1.-g0+2.0*g0*y
     >    -(1.0+g0)*y**2)
 100	continue
	ax=ax*dh*(1.0+g0)
	return
	end
ckg_rad
ckg_pess
      subroutine keep2mainhrmnks(eigfun,mt,nn)
      dimension eigfun(nn,mt,3)
      data jm2/1/,jm3/2/
      emax3=-1000.
      emax2=-1000.
      do j=1,mt
         emax1=-1000.
         do i=2,nn
            if(emax1.lt.abs(eigfun(i,j,1)))then
               emax1=abs(eigfun(i,j,1))
            endif
         enddo
         if(emax1.gt.emax2)then
            jm3=jm2
            emax3=emax2
            jm2=j
            emax2=emax1
         elseif(emax1.gt.emax3)then
            jm3=j
            emax3=emax1
         endif
      enddo
      write(*,*) 'maximum harmonics are', jm2,jm3
      do j=1,mt
         if(j.ne.jm2.and.j.ne.jm3)then
            do i=1,nn
               eigfun(i,j,1)=eigfun(i,j,1)/100.
               eigfun(i,j,2)=eigfun(i,j,2)/100.
               eigfun(i,j,3)=eigfun(i,j,3)/100.
            enddo
         endif
      enddo
      return
      end
ckg_pess
