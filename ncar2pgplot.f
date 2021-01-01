      subroutine grafinit(i)
      if(i.eq.1)then
         write(*,*) 'Initiated PGPLOT'
      else if(i.eq.2)then
         write(*,*) 'closed PGPLOT'         
      endif
      return
      end
      subroutine frame
      return
      end
      subroutine contgraf
      return
      end
      subroutine plchhq
      return
      end
      subroutine twodgraf(x,y,n,npl,col,pat,colbkg
     $     ,logd,lbl,lblx,lbly,lblrun,itty,itty2,klog,yerr,ierb)
      parameter(nin=501,nplin=10)
c      integer*8 n,npl,logd,itty,itty2,klog,ierb
      integer n,npl,logd,itty,itty2,klog,ierb
      character*(*) col(npl),pat(npl),colbkg,lbl,lblx,lbly,lblrun
      character*40 lblv
      double precision x(n),y(n,npl),yerr(n)
      dimension xin(nin),yin(nin,nplin),in(nplin),lblv(nplin)
      data lblv/nplin*'                                  '/
      lblv(1)='m='
      if(lbl(1:3).eq.'gbf')lblv(1)='i='
      if(n.gt.nin.or.npl.gt.nplin)then
         stop 'Check dimensions in NCAR2PGPLOT'
      endif
      do i=1,npl
         in(i)=n
         if(i.ne.1) lblv(i)='               '
         if(i.eq.1.and.abs(klog).lt.10)then
            write(lblv(1)(3:3),'(i1)')abs(klog)
         elseif(i.eq.1.and.abs(klog).lt.100)then
            write(lblv(1)(3:4),'(i2)')abs(klog)
         elseif(i.eq.1.and.abs(klog).lt.1000)then
            write(lblv(1)(3:5),'(i3)')abs(klog)
         elseif(i.eq.1)then
            stop 'Are you out of your mind: NOVA with m > 1000?'
         elseif(abs(klog+i-1).lt.10)then
            write(lblv(i)(1:1),'(i1)')abs(klog+i-1)
         elseif(abs(klog+i-1).lt.100)then
            write(lblv(i)(1:2),'(i2)')abs(klog+i-1)
         elseif(abs(klog+i-1).lt.1000)then
            write(lblv(i)(1:3),'(i3)')abs(klog+i-1)
         endif
      enddo
      xmax=-1.e30
      xmin=1.e30
      ymax=-1.e30
      ymin=1.e30
      do i=1,n
         xin(i)=real(x(i))
         xmax=max(xmax,xin(i))
         xmin=min(xmin,xin(i))
         do j=1,npl
            yin(i,j)=real(y(i,j))
            ymax=max(ymax,yin(i,j))
            ymin=min(ymin,yin(i,j))
         enddo
c         write(*,*)yin(i,2),yin(i,5)
      enddo
c      write(*,*)int4(n),int4(n),int4(npl),xin(1),xmin,xmax,yin(1,1),ymin
c     &     ,ymax,lbl,lblx,lbly
ckg   if junk=0 show interactive window, if =-10, plot the first
c     graph and exit
      junk=-10
      call show_31(lbl,lblx,lbly,lblv(1),1,in(1),nin,npl ! int4(npl)
     &     ,nin,xin(1),xmin,xmax,yin(1,1),ymin,ymax,1,7,junk)
      return
      end
      function stf1(x)
      double precision xd
      xd=x
      stf1=stf(xd)
      return
      end
      function ssum1(lam,wabh,i1)
      integer*8 lam,i1
      double precision wabh(lam),ssum1
c      dimension wabhs(100)
c      lams=lam
      ssum1=0.
      do i=i1,lam
         ssum1=ssum1+wabh(i)
      enddo
c      i1s=i1
c      ssum1=ssum(lams,wabhs(1),i1s)
c      do i=1,lams
c         wabh(i)=wabhs(i)
c      enddo
c      i1=i1s
      return
      end
