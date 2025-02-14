c*********************************************************
c It defines the grid for Berk formula for the estimate of
c TAE amplitude
      subroutine berkgrid(dbob,iberk,dbobmax,ximax)
      dimension dbob(iberk)
      dbob(1)=1.e-7
      dbob(iberk)=1.e-2
      dbob(iberk-1)=(dbob(iberk)/dbob(1))**(1./real(iberk-1))
      do i=2,iberk-1
         dbob(i)=dbob(i-1)*dbob(iberk-1)
      enddo
      return
      end
c*********************************************************
      function q_psi(psi_av)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      iind=indx(ppsi(1),nosurf,psi_av)
      q_psi=ynn(ppsi(iind),qoo(iind),psi_av)
      return
      end
c*********************************************************
c psi using index of the surface
      function psi_indx(index)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      psi_indx=ppsi(index)
      return
      end
c*********************************************************
      subroutine psi_res(w_rotres)
      include 'clich1'
      include 'clich2'
      include 'orb2d.par'
      dimension xrt(nn),ind(nn)
c assume no more then two roots
      call roots(qoo(1),ppsi(1),ind(1),nosurf,1.,1.e-6,0.,xrt(1)
     &     ,nrt)
      w_rotres=xrt(nrt)
      return
      end
c***************************************************************************
      subroutine read_vectr(i1,i2,ier,io,sample,work,nwork)
      dimension work(nwork)
      character dumb*72,sample*40,lblprf*40
      rewind(io)
      i=i2-i1+1
 1    read(io,'(72a)',err=9) dumb
      if(dumb(i1:i2).eq.sample(1:i)) then
         read(io,'(5e13.4)') (work(i),i=1,nwork)
      else
         goto 1
      endif
      ier=0
      return
 9    ier=1
      return
      end
c***************************************************************************
c  this subroutine approximates the rotation by the exponential dependence:
c       w_rot0-w_rotres**p_rot*(w_rot0-w_rot1)
      subroutine readrot(w_rot0,w_rot1,p_rot,ierr)
      include 'gridparam'
      parameter (n=2,iw=n+1)
c      parameter (ndat=20,n=2,iw=n+1)
      dimension psipol(ndat),w_rot(ndat),x(3),w(iw,5),w6(iw,n)
      character sample*40
      common /rotat/w_rot,psipol,p_rot0
      external funct1,monit
      data io/10/,tol/1.e-7/,maxcal/100/,ifail/0/
      open(unit=io,file='transp.dat',err=1,status='old')
      ierr=0
c---------
ckg      sample='thermal plasma toroidal rotation'
ckg      i1=2
ckg      i2=18
      ishft=0
      if(ndat.gt.99) ishft=1
      sample='0 OMEGA'
      i1=5+ishft
      i2=11+ishft
      call read_vectr(i1,i2,ierr,io,sample,w_rot(1),ndat)
      if(abs(w_rot(1)).gt.tol.and.abs(w_rot(ndat)).gt.tol) then
c take not-exactly the first point in the rotation profile
         w_rot0=(w_rot(1)+w_rot(2))*0.5
         x(1)=w_rot0/w_rot(1)
c take not-exactly the last point in the rotation profile
         w_rot1=(w_rot(ndat)+w_rot(ndat-1))*0.5
         x(2)=w_rot1/w_rot(ndat)
c---------
ckg      sample='Poloidal Flux'
ckg      i1=2
ckg      i2=13
         sample='0 PLFLX'
         i1=5+ishft
         i2=11+ishft
         call read_vectr(i1,i2,ierr,io,sample,psipol(1),ndat)
         do i=1,ndat
            psipol(i)=psipol(i)/psipol(ndat)
            w_rot(i)=w_rot(i)/w_rot0
         enddo
         p_rot=alog((w_rot(ndat/2)-w_rot(1))/(w_rot(ndat)-w_rot(1)))
     &        /alog(psipol(ndat/2))
         p_rot0=p_rot
         x(3)=0.
c---------
         write(*,*) 'Rotn ',n,x,f,tol,iw,MAXCAL, IFAIL
         call e04ccf(n, X, F, TOL, IW, W(1,1),W(1,2),W(1,3),W(1,4),W(1
     &        ,5),W6(1,1), FUNCT1, MONIT, MAXCAL, IFAIL)
         p_rot=x(3)+p_rot0
         w_rot1=x(2)*w_rot(ndat)*w_rot0*1.e-6
         w_rot0=x(1)*w_rot(1)*w_rot0*1.e-6
cnng13
c      write(*,'(3f13.5)') ((psipol(i),w_rot(i),w_rot0+psipol(i)**p_rot
c     &     *(w_rot1-w_rot0)),i=1,ndat)
cnng13
c      write(*,*) x,f
      endif
      close(io)
      return
 1    ierr=1
      close(io)
      return
      end
c***************************************************************************
c looks like it integrates the rotation profile
c but only two parameters are found x(1) and x(2)
      subroutine funct1(n,x,f)
      include 'gridparam'
      integer n
c      parameter(ndat=20,nwrk=33,ntry=10)
      parameter(nwrk=33,ntry=15)
      common /rotat/w_rot(ndat),psipol(ndat),p_rot0
      dimension work(nwrk),x(n)
      p_rot=p_rot0
      if(n.eq.3)p_rot=x(3)+p_rot
      do i=1,ntry
         work(i)=(w_rot(i)-x(1)*w_rot(1)-psipol(i)**(p_rot)
     &        *(x(2)*w_rot(ndat)-x(1)*w_rot(1)))**2
      enddo
      call asimp(psipol(1),work(1),ntry,f)
c      write(*,'(4f13.5)') f,x
      return
      end
c***************************************************************************
      subroutine monit(fmin,fmax,sim,n,n1,ncall)
      real sim(n1,n)
c      write(*,*)'ncall',ncall,'fmin',fmin
c      write(*,*)((sim(i,j),j=1,n),i=1,n1)
      return
      end
c***************************************************************************
      subroutine monit1(fmin,fmax,sim,n,n1,ncall)
      real sim(n1,n)
c      write(*,*)'ncall',ncall,'fmin',fmin
c      write(*,*)((sim(i,j),j=1,n),i=1,n1)
      return
      end
c***************************************************************************
c  this subroutine approximates the resonance frequency along the P_phi direction
c     or along the slanted direction 
c     here dwp is array of initial values of the coefficient for the Taylor-like expansion
c          vres is the value of the resonance velocity in K units
c          the idea is to approximate the slanted path from the resonance point
      subroutine wresminimize(wresPfiV,dwp,dwpo,d2wpo,nmuc
     ^     ,imu,npfc,nvc,qpf,v,qpfres,vres,omstar,wb,ifail)
      include 'orbit.grid'
      include 'gridparam'
      parameter(nxw=2,iw=nxw+1)
c      parameter (ndat=20,nxw=2,iw=nxw+1)
      dimension wresPfiV(npfc,nvc),qpf(npfc),v(nvc)
     &     ,wb(nvc,nmuc,npfc),wbc(nv,nmu,npf)
      dimension psipol(ndat),w_rot(ndat),x(3),w(iw,5),w6(iw,nxw)
      external fungw0,monit1
      common /lbqc/wresPfiVc(npf,nv),qpfc(npf),vc(nv),qpfresc,
     &     vresc,omstarc,wbc,imuc
c      data tol/1.e-7/,maxcal/1000/
      data tol/1.e-6/,maxcal/1000/,x/3*0./
c---------
c this "if" is to check the dimensions consistency and stop if they are not
c  they should be OK always
      if(npfc.ne.npf.or.nvc.ne.nv.or.nmu.ne.nmuc) 
     &     stop 'stopping in wresminimize'
      x(1)=-dwp*0.3
      x(1)=0.
      x(2)=0.
cnxw      x(3)=0.
      wresPfiVc=wresPfiV
      qpfc=qpf
      vc=v
      qpfresc=qpfres
      vresc=vres
      wbc=wb
      imuc=imu
      omstarc=omstar
c---------
cnng21 keep checking if a fitting of DF to 2nd order polyn works better than Taylor
c      write(*,*) 'Wresminimize',nxw,x,vres,omstar,qpfres,dwpo
      call e04ccf(nxw, X, F, TOL, IW, W(1,1), W(1,2), W(1,3), W(1,4)
     &     , W(1,5),W6(1,1), fungw0, MONIT1, MAXCAL, IFAIL)
      dwpo=x(1)
      d2wpo=x(2)
cnng13
c      write(*,'(3f13.5)') ((psipol(i),w_rot(i),w_rot0+psipol(i)**p_rot
c     &     *(w_rot1-w_rot0)),i=1,ndat)
cnng13
c      write(*,*) x,f
      return
      end
c***************************************************************************
c provides the integrand for resonance frequency Omega
      subroutine fungw0(n,x,f)
c dimensions of sparse grids
      include 'orbit.grid'
      include 'gridparam'
      integer n
      common /lbqc/wresPfiVc(npf,nv),qpfc(npf),vc(nv),qpfresc,vresc
     ^     ,omstarc,wbc(nv,nmu,npf),imuc
      dimension work(npf),x(n)
c these limit the range of rectangle where the rsnce freq is approximated
c good choice for d3d_153071Duarte was gdp2v=0.004 gdp2pf=0.04
c or for tftr_103101 at gdp2v=0.001(or 0.004) gdp2pf=0.1 (& 0.08)
c for d3d_159243 mostly gdp2v=0.004, gdp2pf=0.04; 0.01 0.08 for n=4
c  following two const are used just below for gaus widths and they are **2
      gdp2v=0.125                !0.05; 0.004 !0.001   could be up to ~0.25
      gdp2pf=0.25               !0.04; 0.08 !0.1                     ~0.5 
      i=1
      do while ( i <= npf )
         iv=1
         do while ( iv <= nv )
c add to the sum only if particle/orbit exist, i.e. if its bounce time !=0
            if(abs(wbc(iv,imuc,i)).gt.1.e-10) then
c if derivatives are the ones which they should be in the considered case
c which is tftr_103101 mode and dOm/dv>0, and dOm/dPf>0
c we get the following expansion
c the underlying comparison m.b.done for Omega_ij and we are looking at the 
c approximation with up to second polynomial in Pfi ignoring v means that this
c is only a 1D problem or more sophisticated 2D 
               work(i)=(wresPfiVc(i,iv)-x(1)*( (qpfc(i)-qpfresc)
     ^              -(vc(iv)-vresc)*(omstarc*vc(iv))
     ^              )-(qpfc(i)-qpfresc)**2*x(2)*0.5
     ^              )**2*exp(-(vc(iv)-vresc)**2 /gdp2v -
     ^              (qpfresc-qpfc(i))**2/gdp2pf)
            endif
            iv=iv+1
         enddo
         i=i+1
      enddo
      call asimp(qpfc(1),work(1),npf,f)
c      write(*,'(4f13.5)') f,x
c      p_rot=p_rot0
c      if(n.eq.3)p_rot=x(3)+p_rot
      return
      end
c***************************************************************************
