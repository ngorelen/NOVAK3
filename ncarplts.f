c ncaplts v1.1 \\ by derek fox
c
c                      ********************
c                      *                  *
c                      *     ncarplts     *
c                      *                  *
c                      ********************
c......................................................................
c
c    The purpose of this plot package is to make the process of 
c  graphing with the ncar graphics library easier and the output of
c  the process more elegant.  You should find that the routines in 
c  this package offer far more power than ncargks readily and wil-
c  lingly yields forth on its own.
c
c    Bugs/complaints/suggestions to derekfox@phoenix.princeton.edu.
c  I can't promise to fix it, but I do promise to reply.
c......................................................................
c
c  USE:
c
c       (1)  CALL GRAFINIT(1) to initialize graphics;
c  
c  Then, for each plot:
c  
c       (2) call one of the comprehensive graphing routines
c                TWODGRAF, CONTGRAF, or THRDGRAF, or
c
c       [2'] call regular ncar routines to draw a basic plot   
c                (this package provides no comprehensive routine
c                 for plotting in THREED space, for example)
c
c       [3]  call subsidiary routines from ncarplts or the
c                regular ncar library to add some finishing
c                touches, if you wish; and then
c
c       (4)  CALL FRAME to terminate the plot.
c
c  Finally, 
c
c       (5)  CALL GRAFINIT(2) to terminate graphics.
c
c
c  COMPILATION:  requires ncar and ncargks libraries, v3.0.  Compile
c     ncarplts.f and link to your program and to the libraries.
c
c        An easy way to do this is to use ncar's special version of
c     the fortran compiler, ncargf77:
c
c             ncargf77 ncarplts.f program.f
c
c     where the name of your program replaces "program.f".  This 
c     will produce the executable "a.out".  
c
c        In addition, GNU Make Makefiles are included for the suns
c     and crays which compile the example program and are readily
c     modified.  
c
c        Duplication Cautions should be expected because of the 
c     replacement of NCAR's default TRN32S and CPMPXY routines.  
c
c
c  ROUTINES:
c
c  Comprehensive graphing routines (for details, see the 
c    full argument listings later in this documentation):  
c
c       contgraf(zplot,m1dim,m2dim,xmn,xmx,ymn,ymx,kscale,
c         kchz,chzcon,krect,xij,yij,nodec,kcolor,labelplt,
c         labelx,labely,labrun,ntty,niorec)
c
c                                 framed color contour plot
c
c       thrdgraf(xra,yra,zplot,m1dim,m2dim,labelplt,labelx,
c         labely,labelz,labrun,ntty,niorec,x3dview,kscale)
c
c                                        boxed surface plot
c
c       twodgraf(xplot,yplot,kpts,kcurves,kcrvclr,kcrvtype,
c         kbackclr,delymax,labelplt,labelx,labely,labrun,
c         ntty,niorec,kscale,yploterr,kerrplot)
c
c                         color 2d plot with multiple lines
c
c  Subsidiary routines (see the routines themselves for 
c    further information):
c
c       colram:   slave routine for contgraf (coloring
c                 of areas between contour lines).
c       cpmpxy:   slave routine for contgraf (remapping
c                 of data for irregular grids).
c       curv3d:   draws dashed curves in THREED space.
c       findlen:  returns length of a string ending in 
c                 '$'.
c       grafinit: initializes & finalizes graphics (see 
c                 above, under USE).
c       hlin3s:   allows (somewhat erratic) hidden-line 
c                 drawing in SRFACE space.
c       intchr:   reads integer number into a string.
c       lblpos:   figures out where to draw the labels
c                 on a 3d box, so that they face the 
c                 viewer.
c       line3d:   draws dashed lines in THREED space.
c       line3s:   draws dashed lines in SRFACE space.
c       nicify:   "rounds off" numbers for plot edges
c                 and chooses an interval size.
c       realchr:  reads a real number into a string.
c       setclr:   sets color for all ncar purposes, by
c                 name.  
c       setpat:   sets pattern for all ncar line-
c                 drawing purposes, by name.  
c       trn32s:   slave routine for thrdgraf, hlin3s,
c                 line3s (carries out perspective 
c                 transformation).  
c
c
c  ARGUMENT LISTINGS:
c
c  Contgraf:
c *********************************************************************
c *  zplot(i,j)     the m1dim by m2dim array to be plotted            *
c *  m1dim          first dimension of zplot                          *
c *  m2dim          second dimension of zplot                         *
c *  xmn            value of x at left boundary of frame              *
c *  xmx            value of x at right boundary of frame             *
c *  ymn            value of y at lower boundary of frame             *
c *  ymx            value of y at upper boundary of frame             *
c *  kscale         = 1: plot zplot                                   *
c *                 = 2: plot log10(zplot)                            *
c *  kchz           = 0: conpack chooses contour levels; chzcon ir-   *
c *                      relevant                                     *
c *                 =-1: chzcon contains number of contour levels     *
c *                      desired (as a real number)                   *
c *                 =-2: chzcon contains desired spacing of contour   *
c *                      levels (the interval as a real number)       *
c *                 > 0: chzcon is a real array with kchz elements    *
c *                      giving the values at which contour lines     *
c *                      should be drawn                              *
c *  chzcon         = real number or real array; meaning depends on   *
c *                   value of "kchz" (see above)                     *
c *  krect          = 0: data lies on evenly-spaced rectangular grid; *
c *                      xij and yij ignored                          *
c *                 = 1: NCAR ezmap mapping (x=longitude, y=latitude) *
c *                      xij(1) and xij(2) are min and max longitude; *
c *                      yij(1) and yij(2) are min and max latitude;  *
c *                      both coordinates measured in degrees.        *
c *                 = 2: Polar coordinate mapping (x=radius, y=theta) *
c *                      xij(1) and xij(2) are min and max radius;    *
c *                      yij(1) and yij(2) are min and max theta;     *
c *                      angle is measured in degrees.                *
c *                 = 3: x=x(i) and y=y(j); xij and yij contain one-  *
c *                      dimensional arrays for the mapping           *
c *                 = 4: x=x(i,j) and y=y(i,j); xij and yij contain   *
c *                      two-dimensional arrays of the mapping        *
c *  xij            =(real array):  meaning depends on value of krect *
c *  yij            =(real array):  meaning depends on value of krect *
c *  nodec          number of decimal places to use for legend labels *
c *                  (contgraf decides whether exponential or not)    *
c     *  kcolor         = 0: for output to mono printer (color lines)     *
c *                 <>0: for output to screen only or color printer   *
c *                      (color fill)                                 *
c *  labelplt(i)    title of plot                                     *
c *  labelx(i)      x-axis label                                      *
c *  labely(i)      y-axis label                                      *
c *  labrun(i)      run label                                         *
c *  ntty           logical unit number for writing error messages    *
c *  niorec         2nd logical unit number for writing error messages*
c *                                                                   *
c *  note:                                                            *
c *                                                                   *
c *  (1) labelplt, labelx, labely, labrun must end with the character *
c *      "$".                                                         *
c *********************************************************************
c
c  Thrdgraf:
c *********************************************************************
c *  xra(i)         array of x-coordinates (m1dim long)               *
c *  yra(j)         array of y-coordinates (m2dim long)               *
c *  zplot(i,j)     the m1dim by m2dim array to be plotted            *
c *  m1dim          first dimension of zplot                          *
c *  m2dim          second dimension of zplot                         *
c *  labelplt(i)    title of plot                                     *
c *  labelx(i)      x-axis label                                      *
c *  labely(i)      y-axis label                                      *
c *  labelz(i)      z-axis label                                      *
c *  labrun(i)      run label                                         *
c *  ntty           logical unit number for writing error messages    *
c *  niorec         2nd logical unit number for writing error messages*
c *  x3dview(i)     viewing position coordinates (x,y,z), in units of *
c *                 the plot frame box.  The box is 1 unit on a side, *
c *                 with its (min x, min y, min z) corner at the      *
c *                 origin.                                           *
c *  kscale         = 1: plot zplot                                   *
c *                 = 2: plot log10(zplot)                            *
c *                                                                   *
c *  note:                                                            *
c *                                                                   *
c *  (1) labelplt, labelx, labely, labelz and labrun must end with    *
c *       the character "$" for proper placement.                     *
c *  (2) labelx, labely, and labelz must be all-caps, standard ANSI   *
c *      fortran 77 characters.                                       *
c *  (3) xra and yra must be monotonically increasing functions of    *
c *      their indices.                                               *
c *********************************************************************
c
c Twodgraf:
c *********************************************************************
c *  xplot(i)       array of (kpts) x-coordinates for the points      *
c *  yplot(i,j)     kpts by kcurves array of y-coordinates            *
c *  kpts           number of points per curve                        *
c *  kcurves        number of curves                                  *
c *  kcrvclr(i)     (kcurves) array of color names for the curves     *
c *  kcrvtype(i)    (kcurves) array of pattern names for the curves   *
c *  kbackclr       color name for grid, labels, background           *
c *  delymax        maximum number of decades allowed for log scale   *
c *  labelplt       title for the plot (*)                            *
c *  labelx         x-axis label (*)                                  *
c *  labely         y-axis label (*)                                  *
c *  labrun         "run label" -- printed in small type at bottom (*)*
c *  ntty           logical unit number for writing error messages    *
c *  niorec         2nd logical unit number for writing error messages*
c *  kscale         = 1: linear scale                                 *
c *                 = 2: semi-log scale (log y-axis)                  *
c *  yploterr(i,j)  kpts by kcurves array of y error values           *
c *  kerrplt        = 0: no error bars                                *
c *                 = 1: error bars                                   *
c *                                                                   *
c *  note:                                                            *
c *                                                                   *
c *  (*)'d strings must end in '$' to be centered properly.           *
c *********************************************************************
c
c  The code for all routines can be found below, in alphabetical order.   
c......................................................................
cend documentation
c*dk colram
c**********************************************************************
cbeg                        colram                                    *
c**********************************************************************
      subroutine colram(xcra,ycra,ncra,iaia,igia,naia)
c
      dimension xcra(*),ycra(*),iaia(*),igia(*)
c......................................................................
c        colram fills in the area between contour lines in color.
c......................................................................
c
      ifll=1
cncar
      call cpgeti('NCL - No. of Cnt Lines',nocl)
      step = 45./nocl
cncar
c
c  If any of the area identifiers is negative, then don't fill
c
      do 10  i=1,naia
        if(iaia(i).lt.0) ifll=0
   10 continue
c
c  Otherwise, fill the area with the color defined by its area
c    identifier relative to group 3, the contour line group
c
      if(ifll.ne.0)then
        ifll=0
        do 20  i=1,naia
          if(igia(i).eq.3)ifll=iaia(i)
   20   continue
        if(ifll.gt.0.and.ifll.le.nocl+1)then
cncar
          call gsfaci(nint( (ifll-1)*step ) + 1)
          call gfa(ncra-1,xcra,ycra)
cncar
        endif
      endif
c
      return
      end
c
cend colram
c*dk contgraf
c**********************************************************************
cbeg                        contgraf                                  *
c**********************************************************************
      subroutine contgraf(zplot,m1dim,m2dim,xmn,xmx,ymn,ymx,kscale,
     &  kchz,chzcon,krect,xij,yij,nodec,kcolor,labelplt,labelx,labely,
     &  labrun,ntty,niorec)
c
        character*(*) labelx,labely,labelplt,labrun
	character*50 chclv, llbs(50)
        character*8 kbackclr
	real lgnd(4),grph(4),grid(4),xlplot(2),
     &   ylplot(2)
        dimension zplot(*),chzcon(*),xij(*), yij(*)
        external cpdrpl
        external colram
c                                                      Workspace arrays
        dimension zzplot(500000)
        dimension rwrk(5000), iwrk(1000), iama(1000000)
        dimension xcra(5000), ycra(5000), iaia(10), igia(10), 
     &       iasf(13), lclr(50), tksa(10), tksb(10)
c       
        common /cpmpcm1/ comxi(1000),comyj(1000)
        common /cpmpcm2/ comxij(1000,1000),comyij(1000,1000)
        common /cpmpinf/ mxdim,nydim
c
        data tksa / 1.e36, 1.e36, 6., 1.e36, 0.010, 0., 
     &    1.e36, 1.e36, 0., 0. /
        data tksb / 1.e36, 1.e36, 6., 1.e36, 0., 0.015, 
     &    1.e36, 1.e36, 0., 0.010 /
c
        data grph / 0.0, 1.0, 0.0, 1.0 /
        data grid / 0.1, 0.9, 0.1, 0.9 /
	data lgnd / 0.905, 1., 0., 0.9 /
        data iasf / 13*1 /
c......................................................................
c       contgraf plots a contour plot of zplot using the conpack
c    utility of the ncar graphics library v3.0, frames the plot,
c    and writes a legend.  it was cannibalized from the contgraf
c    routine of degraf 60.0 written by derek fox under supervision 
c    of daren stotler, 8/92.  further modified and improved by 
c    derek fox 1/93.
c
c    calling arguments
c    -----------------
c
c    zplot(i,j)     the m1dim by m2dim array to be plotted
c    m1dim          first dimension of zplot
c    m2dim          second dimension of zplot
c    xmn            value of x at left boundary of frame
c    xmx            value of x at right boundary of frame
c    ymn            value of y at lower boundary of frame
c    ymx            value of y at upper boundary of frame
c    kscale         = 1: plot zplot
c                   = 2: plot log10(zplot)
c    kchz           = 0: conpack chooses contour levels; chzcon ir-
c                        relevant
c                   =-1: chzcon contains number of contour levels 
c                        desired (as a real number)
c                   =-2: chzcon contains desired spacing of contour
c                        levels (the interval as a real number)
c                   > 0: chzcon is a real array with kchz elements
c                        giving the values at which contour lines 
c                        should be drawn
c    chzcon         = real number or real array; meaning depends on 
c                     value of "kchz" (see above)
c    krect          = 0: data lies on evenly-spaced rectangular grid;
c                        xij and yij ignored
c                   = 1: NCAR ezmap mapping (x=longitude, y=latitude)
c                        xij(1) and xij(2) are min and max longitude;
c                        yij(1) and yij(2) are min and max latitude;
c                        both coordinates measured in degrees.
c                   = 2: Polar coordinate mapping (x=radius, y=theta)
c                        xij(1) and xij(2) are min and max radius;
c                        yij(1) and yij(2) are min and max theta; 
c                        angle is measured in degrees.  
c                   = 3: x=x(i) and y=y(j); xij and yij contain one-
c                        dimensional arrays for the mapping
c                   = 4: x=x(i,j) and y=y(i,j); xij and yij contain
c                        two-dimensional arrays of the mapping
c    xij            =(real array):  meaning depends on value of krect
c    yij            =(real array):  meaning depends on value of krect
c    nodec          number of decimal places to use for legend labels
c                     (contgraf decides whether exponential or not)
c    kcolor         = 0: for output to mono printer (color lines)
c                   <>0: for output to screen only or color printer
c                        (color fill)
c    labelplt(i)    title of plot
c    labelx(i)      x-axis label
c    labely(i)      y-axis label
c    labrun(i)      run label
c    ntty           logical unit number for writing error messages
c    niorec         2nd logical unit number for writing error messages
c
c    note:
c
c    (1) labelplt, labelx, labely, labrun must end with the character 
c        "$".  
c......................................................................
c    
c                       Set all GKS aspect source flags to 'individual'
cncar                                              and force solid fill
        call gsasf(iasf)
        call gsfais(1)
cncar                                    Find the length of all strings
c                                         <nlblplt,nlblx,nlbly,nlabrun>
        call findlen(labelplt,nlblplt)
        call findlen(labelx,nlblx)
        call findlen(labely,nlbly)
        call findlen(labrun,nlabrun)
c
c                                Set color for grid, background, labels
c                      (also see the legend plotting routines, however)
        kbackclr='cyan'
c                                 Set the dash patterns for solid, dash
c                                                  <isolidpat,idashpat>
        isolidpat = 65535
        idashpat  = 52428
c                                        Read arrays into common blocks
c                                          for non-rectangular mappings
        if(krect.eq.1 .or. krect.eq.2)then
          xc1 = xij(1)
          xcm = xij(2)
          yc1 = yij(1)
          ycn = yij(2)
        endif
c
        if(krect.eq.3 .or. krect.eq.4)then
          mxdim = m1dim
          nydim = m2dim
          xc1 = 1.
          xcm = real(m1dim)
          yc1 = 1.
          ycn = real(m2dim)
        endif
c
        if(krect.eq.3)then
          do 2  i=1,m1dim
            comxi(i) = xij(i)
    2     continue
c
          do 4  j=1,m2dim
            comyj(j) = yij(j)
    4     continue
        endif
c
        if(krect.eq.4)then
          do 8 j=1,m2dim
            do 6 i=1,m1dim
              index = (j-1)*m1dim + i
              comxij(i,j) = xij(index)
              comyij(i,j) = yij(index)
    6       continue
    8     continue
        endif
c
c
 10     j12dim=m1dim*m2dim
        if(j12dim.eq.0 .or. j12dim.gt.500000)go to 1010
        if(xmn.ge.xmx)go to 1012
        if(ymn.ge.ymx)go to 1016
c                                                          clear zzplot
        do 12 i=1,j12dim
 12       zzplot(i)=0.0
c
        zmin=+1.e30
        zmax=-1.e30
        go to(15,20)kscale
c                                                 kscale=1: linear plot
 15     do 16 i=1,j12dim
          zzplot(i)=zplot(i)
          zmin=min(zmin,zzplot(i))
          zmax=max(zmax,zzplot(i))
 16     continue
        go to 50
c                                                    kscale=2: log plot
 20     do 25 i=1,j12dim
          zploti=abs(zplot(i))
          if(zploti.eq.0.0)zzplot(i)=0.0
          if(zploti.ne.0.0)then
            zzplot(i)=alog10(zploti)
            zmin=min(zmin,zzplot(i))
            zmax=max(zmax,zzplot(i))
          endif
 25     continue
c
        if(zmin.ge.zmax)go to 1035
c
        do 40 i=1,j12dim
 40       if(zzplot(i).eq.0.0)zzplot(i)=zmin
c                                                           plot zzplot
cncar
 50     call set(0.,1.,0.,1., 0.,1.,0.,1., 1)
        call setclr(kbackclr)
c                                                       Draw the labels
c
        call plchhq(0.5,0.95,labelplt(1:nlblplt), 0.015, 0., 0.)
        call plchhq(0.5,0.05,labelx(1:nlblx), 0.01, 0., 0.)
        call plchhq(0.02,0.5,labely(1:nlbly), 0.01, 90., 0.)
c
c                                                   Initialize the plot
        call cprset
        call cpsetc('HLT - Hi/Lo lable Text string','''')
        call cpsetc('ILT - Info Lable Text string','''')
        call cpseti('LLP - Line Label Positioning',3)
        call cpsetr('T2D - 2D Smoothing',2.5)
c                                                 Turn on special value
        call cpsetr('SPV - SPecial Value',1.e36)
        call cpseti('PAI - Param Array Index',-2)
        call cpseti('CLU - Contour Level Use flag',1)
        call cpsetr('CLL - Contour Level Line width',2.)
c
c                         For irregular mappings, outline the grid edge
        if(krect.eq.1 .or. krect.eq.2 .or. krect.eq.4)then
          call cpseti('PAI - Param Array Index',-1)
          call cpseti('CLU - Contour Level Use flag',1)
          call cpsetr('CLL - Contour Level Line width',2.)
        endif
c
        call cpseti('MAP - Type of Mapping',krect)
c
        if(krect.ne.0)then
          call cpsetr('XC1 - X Coord at i=1',xc1)
          call cpsetr('XCM - X Coord at i=M',xcm)
          call cpsetr('YC1 - Y Coord at j=1',yc1)
          call cpsetr('YCN - Y Coord at j=N',ycn)
          call cpseti('SET - do SET call?',0)
          call set(grid(1),grid(2),grid(3),grid(4),
     &             xmn,xmx,ymn,ymx,1)
        else
          call cpsetr('VPL - View Portal Left',grid(1))
          call cpsetr('VPR - View Portal Right',grid(2))
          call cpsetr('VPB - View Portal Bottom',grid(3))
          call cpsetr('VPT - View Portal Top',grid(4))
          call cpsetr('VPS - View Portal Scale',0.0)
        endif
c
        call cprect(zzplot,m1dim,m1dim,m2dim, rwrk,5000, iwrk,1000)
c
        if(kchz.gt.0)then
          call cpseti('CLS - Cont Lvl Selection',0)
          nocl = kchz
          call cpseti('NCL - Number of Cont Lvls',nocl)
          do 52  i=1,nocl
            call cpseti('PAI - Param Array Index',i)
            call cpsetr('CLV - Cont Level Value',chzcon(i))
   52     continue
        else
          if(kchz.eq.-1)
     &      call cpseti('CLS - Cont Lvl Selection',nint(-chzcon(1)))
          if(kchz.eq.-2)then
            call cpsetr('CIS - Cont Interval Specifier',chzcon(1))
            call cpsetr('CMN - Cont Minimum',zmin)
            call cpsetr('CMX - Cont Maximum',zmax)
            call cpseti('CLS - Cont Lvl Selection',1)
          endif
          call cppkcl(zzplot,rwrk,iwrk)
          call cpgeti('NCL - Number of Cont Lvls',nocl)
        endif
c                                               Adjust colors, patterns
        if(kcolor.eq.0)then
          clrstp = 45./(nocl-1)
          do 55  i=1,nocl
            call cpseti('PAI - Param Array Index',i)
            call cpseti('CLU - Cont Lvl Use',3)
	    call cpseti('CLC - Cont Lvl Color',nint((i-1)*clrstp)+1)
            if(mod(i,2).eq.1) then
              call cpseti('CLD - Cont Lvl Dash Ptn',idashpat)
            else
              call cpseti('CLD - Cont Lvl Dash Ptn',isolidpat)
            endif
   55     continue
        endif
c
c
c                                       Prepare autograph to frame plot
        call set(0.,1.,0.,1., 0.,1.,0.,1., 1)
        call agseti('FRAME.',2)
        call anotat('$','$',0,0,0,0)
        call agsetp('GRAPH.',grph,4)
        call agsetp('GRID.',grid,4)
        call agseti('X/NICE.',0)
	call agsetf('X/MIN.',xmn)
	call agsetf('X/MAX.',xmx)
        call agseti('Y/NICE.',0)
	call agsetf('Y/MIN.',ymn)
	call agsetf('Y/MAX.',ymx)
        if(kcolor.ne.0)then
          call agsetp('LEFT/TICKS.',tksa,10)
          call agsetp('BOTTOM/TICKS.',tksa,10)
          call agseti('RIGHT/CONTROL.',-1)
          call agseti('TOP/CONTROL.',-1)
        endif
c                                                        Frame the plot
	call agstup(0.,0,0,0,0,0.,0,0,0,0)
	call agback
c                                            Restore autograph defaults
        if(kcolor.ne.0)then
          call agsetp('LEFT/TICKS.',tksb,10)
          call agsetp('BOTTOM/TICKS.',tksb,10)
          call agseti('RIGHT/CONTROL.',4)
          call agseti('TOP/CONTROL.',4)
        endif
	call agseti('X/NICE.',-1)
	call agsetf('X/MIN.',1.e36)
	call agsetf('X/MAX.',1.e36)
	call agseti('Y/NICE.',-1)
	call agsetf('Y/MIN.',1.e36)
	call agsetf('Y/MAX.',1.e36)
	call agseti('FRAME.',1)
c
c                                          Set up for line labels/color
        call arinam(iama,400000)
c       
        if(krect.ne.0)then
          call set(grid(1),grid(2),grid(3),grid(4),
     &             xmn,xmx,ymn,ymx,1)
        endif
c
        if(kcolor.eq.0)then
c                                                Add labels to area map
          call cplbam(zzplot,rwrk,iwrk,iama)
c                                                     Draw the contours
          call cpcldm(zzplot,rwrk,iwrk,iama,cpdrpl)
c                                                      Draw line labels
          call cplbdr(zzplot,rwrk,iwrk)
c                                                        (kcolor.ne.0):
        else
c                                         Add contour lines to area map
          call cpclam(zzplot,rwrk,iwrk,iama)
c                                                         Color the map
          call arscam(iama,xcra,ycra,5000,iaia,igia,10,colram)
c
        endif
c
c                                                         Set up legend
   75   call set(0.,1.,0.,1., 0.,1.,0.,1., 1)
        lgnd(3) = lgnd(4) - 0.02*nocl - 0.015
c
        if(kcolor.ne.0)clrstp=45./nocl
c
c                                   Find proper decimal form for levels
c
        if(kchz.le.0)then
          call cpgetr('CIU - Cont Interval Used',ciu)
        else
          ciu = (zmax-zmin)/real(nocl)
        endif
c
        if(ciu.gt.100 .or. ciu.lt.0.01)nodec = -nodec
c
c                                       Set legend label & color arrays
        do 80  i=1,nocl
          call cpseti('PAI - Param Array Index',i)
          call cpgetr('CLV - Cont Lvl Value',clv)
          call realchr(clv,nodec,chclv,nchclv)
          llbs(i) = chclv(1:nchclv)
          lclr(i) = nint( (i-1)*clrstp ) + 1
   80   continue
          lclr(nocl+1) = 46
c                                                     Legend with lines
        if(kcolor.eq.0)then
          xlplot(1) = lgnd(1) + 0.005
          xlplot(2) = xlplot(1) + 0.02
c                                                       text for legend
          do 85  i=1,nocl
            yloff = lgnd(4) - 0.02*(nocl-i+1)
	    call setclr(kbackclr)
            call plchhq(xlplot(2)+0.004,yloff, llbs(i), 0.008, 0., -1.)
c
c                                                 draw lines for legend
            ylplot(1) = yloff+0.007
            ylplot(2) = yloff-0.007
            if(mod(i,2).eq.0) call setpat('solid')
            if(mod(i,2).eq.1) call setpat('dash')
            call gsplci(lclr(i))
            call curved(xlplot,ylplot,2)
   85     continue
        endif
c                                                Legend with color fill
        if(kcolor.ne.0)then
c                          (color 16 is cyan -- see setclr for details)
          call lbseti('CLB - Color of LaBels',16)
          call lbseti('CBL - Color of Box Lines',16)
          call lblbar(1, lgnd(1),lgnd(2),lgnd(3),lgnd(4), nocl+1, 
     &      0.33,1., lclr,1, llbs,nocl,1) 
        endif
c                                                 Restore some defaults
        call cpseti('CLS - Cont Lvl Selection',16)
        call dashdb(ior(isolidpat,0))
c                                                        Plot run label
 90     call set(0.,1.,0.,1., 0.,1.,0.,1., 1)
        call setclr(kbackclr)
        call plchhq(0.1,0.01, labrun(1:nlabrun), 0.008, 0., -1.)
cncar
c
 99     return
c----------------------------------------------------------------------
c                                                                errors
 1010   write(ntty,10010)labelplt(1:40),j12dim
        write(niorec,10010)labelplt(1:40),j12dim
        go to 1099
 1012   write(ntty,10012)
        write(niorec,10012)
        go to 1099
 1014   write(ntty,10014)
        write(niorec,10014)
        go to 1099
 1016   write(ntty,10016)
        write(niorec,10016)
        go to 1099
 1018   write(ntty,10018)
        write(niorec,10018)
        go to 1099
 1035   write(ntty,10035)
        write(niorec,10035)
 1099   return
c----------------------------------------------------------------------
c                                                     format statements
10010   format(x,a,x,'error in contgraf: m1dim*m2dim=',x,i7,x,
     &   '.eq.0.or.gt.500000')
10012   format(' error in contgraf: xmn.ge.xmx')
10014   format(' error in contgraf: delx.le.0')
10016   format(' error in contgraf: ymn.ge.ymx')
10018   format(' error in contgraf: dely.le.0')
10035   format(' no contour plot made: zmin.ge.zmax')
        end
cend contgraf
c*dk cpmpxy
c**********************************************************************
c                           cpmpxy                                    *
c**********************************************************************
        subroutine cpmpxy (imap,xin,yin,xout,yout)
c......................................................................
c  cpmpxy is the utility used by CONPACK to remap data points prior to 
c      contouring.  This version allows for four types of special 
c      plots:
c
c      imap=1:  NCAR EZMAP mapping.  x is longitude in degrees, and y
c               is latitude.  
c      imap=2:  Polar coordinates.  x is radius and y is polar angle,
c               in degrees.
c      imap=3:  Orthogonal, unevenly-spaced.  x=comxi(i) and 
c               y=comyj(j).  Calling routine must fill these arrays.
c      imap=4:  Generalized distortion.  x=comxij(i,j) and 
c               y=comyij(i,j).  Calling routine must fill the arrays.
c......................................................................
c
        common /cpmpcm1/ comxi(1000),comyj(1000)
        common /cpmpcm2/ comxij(1000,1000),comyij(1000,1000)
        common /cpmpinf/ mxdim,nydim
c
        xout = xin
        yout = yin
c
        if(imap.eq.1)then
          call maptrn(yin,xin,yout,xout)
        endif
c
        if(imap.eq.2)then
          xout = xin * cos(.017453292519943*yin)
          yout = xin * sin(.017453292519943*yin)
        endif
c
        if(imap.eq.3)then
          i = max(1,min(mxdim-1,int(xin)))
          j = max(1,min(nydim-1,int(yin)))
          xout = (real(i+1)-xin)*comxi(i) + (xin-real(i))*comxi(i+1)
          yout = (real(j+1)-yin)*comyj(j) + (yin-real(j))*comyj(j+1)
        endif
c
        if(imap.eq.4)then
          i = max(1,min(mxdim-1,int(xin)))
          j = max(1,min(nydim-1,int(yin)))
          xout = (real(j+1)-yin)*( (real(i+1)-xin)*comxij(i,j) + 
     &      (xin-real(i))*comxij(i+1,j) )  +  (yin-real(j))*
     &    ( (real(i+1)-xin)*comxij(i,j+1) + (xin-real(i))*
     &      comxij(i+1,j+1) )
          yout = (real(j+1)-yin)*( (real(i+1)-xin)*comyij(i,j) + 
     &      (xin-real(i))*comyij(i+1,j) )  +  (yin-real(j))*
     &    ( (real(i+1)-xin)*comyij(i,j+1) + (xin-real(i))*
     &      comyij(i+1,j+1) )
        endif
c
        return
        end
c
cend cpmpxy
c*dk curv3d
c**********************************************************************
cbeg                        curv3d                                    *
c**********************************************************************
	subroutine curv3d(x3ra,y3ra,z3ra,nopts)
c......................................................................
c  curv3d plots a (colored, dashed) line through the projections of 
c      the nopts number of points defined by x3ra,y3ra,z3ra.  It is a
c      THREED utility, and must be preceded by a SET3 call.
c......................................................................
	dimension x3ra(*),y3ra(*),z3ra(*)
	dimension xra(100),yra(100)
c
	if(nopts.lt.2 .or. nopts.gt.100)go to 100
cncar
	do 10  i=1,nopts
	  call trn32t(x3ra(i),y3ra(i),z3ra(i),xra(i),yra(i),dummy,2)
	  xra(i) = cpux(ifix(xra(i)))
	  yra(i) = cpuy(ifix(yra(i)))
 10     continue
c
	call curved(xra,yra,nopts)
cncar
100	return
	end
cend curv3d
c*dk findlen
c**********************************************************************
cbeg                         findlen                                  *
c**********************************************************************
      subroutine findlen(string,lstring)
c......................................................................
c Finds the `$' that signals termination of strings and returns its
c    position, so that strings can be centered appropriately.
c......................................................................
c
      character*(*) string
c
      length = len(string)
c
      do  10 i=1,length
        if (string(i:i).eq.'$') go to 20
   10 continue
c
   20 lstring = i-1
C
      return
      end
cend findlen
c*dk grafinit
c**********************************************************************
cbeg                       grafinit                                   *
c**********************************************************************
        subroutine grafinit(kinit)
c
c......................................................................
c     grafinit initializes (kinit=1) and finalizes (kinit=2) the
c   graphical output.  
c......................................................................
c
        integer kinit
    5   go to (10,100)kinit
c
c                                              kinit=1: initialize ncar
cncar
   10   call opngks
   28   call gstxfp(12,2)
c                                                    Set up color table
c
c  (0=blck,1=violet,6=blue,16=cyan,26=green,36=yellow,41=orng,46=red;
c   51=magenta,61=white,and other numbers 1-60 are in-between colors)
c
        do 20 i=1,60
          if(i.le.46)then
            hue = 270. - (i-1)*6.
          else
            hue = 360. - (i-46)*6.
          endif
          call hsvrgb( hue, 1., 1., red, green, blue )
          call gscr( 1, i, red, green, blue )
   20   continue
c
        call gscr( 1,  0, 0., 0., 0.)
        call gscr( 1, 61, 1., 1., 1.)
c                                           Set function code delimiter
c                                           character for plchhq to ''
c
	call pcsetc('FC - function code delimiter','''')
c
c                                            A little self-promotion...
        call set(0.,1.,0.,1., 0.,1.,0.,1., 1)
        call setclr('white')
cccc        call plchhq(0.5,0.54, 'ncarplts 1.1', 0.07, 0., 0.)
c
        call setclr('violet')
cccc        call plchhq(0.5,0.41, 'by derek fox', 0.025, 0., 0.)
cccc        call plchhq(0.5,0.37, 'princeton plasma physics laboratory',
cccc     &       0.01, 0., 0.)
c
cccc        call frame
cncar
c
   99   return
c
c                                              kinit=2: finish graphics
cncar
 100    call clsgks
cncar
 199    return
        end
cend grafinit
c*dk hlin3s
c**********************************************************************
cbeg                         hlin3s                                   *
c**********************************************************************
        subroutine hlin3s(x1,y1,z1,x2,y2,z2)
c......................................................................
c     hlin3s draws a line between the 2d projections of the points
c   (x1,y1,z1) and (x2,y2,z2).  it is a SRFACE utility, and must be
c   preceded by the appropriate SRFACE or TRN32S call.  
c
c     the line will (most of the time) not draw over a surface that
c   has been drawn by SRFACE.  the line may be colored but not 
c   dashed.
c......................................................................
c
cncar
        call trn32s(x1,y1,z1,xt1,yt1,dummy,1)
        call trn32s(x2,y2,z2,xt2,yt2,dummy,1)
c
	ixt1 = xt1
	iyt1 = yt1
	ixt2 = xt2
	iyt2 = yt2
c
        call draws(ixt1,iyt1, ixt2,iyt2, 1,0)
cncar
 100    return
        end
cend hlin3s
c*dk intchr
c**********************************************************************
cbeg                         intchr                                   *
c**********************************************************************
	subroutine intchr(intno,chint,lchint)
c......................................................................
c intchr takes an integer and returns a character string chint with 
c     intno taking up the first lchint places.
c......................................................................
	character*(*) chint
c
	nintno=intno
c
	if(nintno.eq.0)then
	  chint = '0'
	  lchint = 1
	  go to 20
	endif
c
	nascii0 = ichar('0')
c
	chint(1:1) = ' '
	if(nintno.lt.0)chint(1:1) = '-'
	  nintno = iabs(nintno)
c
	lchint = int(log10(real(nintno)))+1
c
	do 10  i=lchint,1,-1
	  j = lchint - i + 2
	  chint(j:j) = char(nascii0+mod(nintno,10**i)/(10**(i-1)))
  10    continue
c
	if(chint(1:1).eq.'-')then
	  lchint=lchint+1
	else
	  chint = chint(2:lchint+1)
	endif
c
  20    return
	end
cend intchr
c*dk lblpos
c**********************************************************************
cbeg                          lblpos                                  *
c**********************************************************************
	subroutine lblpos(vuncar,ixydis,iyxdis,izxdis,izydis)
c
	dimension vuncar(6)
c......................................................................
c       lblpos determines which side of the plot volume each label
c     should appear on, and sets its arguments appropriately.
c
c	vuncar:  contains coordinates of eye position, followed by
c                coordinates of point looked at.  for the routine 
c                to work correctly, this second point should be the 
c                center of the plot volume.
c     returns
c
c	iabdis:	 .eq.1 for a-axis label on max-b side of box
c                .eq.-1 for a-axis label on min-b side of box
c......................................................................
c
	dispx = vuncar(1) - vuncar(4)
	dispy = vuncar(2) - vuncar(5)
c
	ixydis = nint(sign(1.,dispy))
	iyxdis = nint(sign(1.,dispx))
	izxdis = ixydis
	izydis = -iyxdis
c
	return
	end
cend lblpos
c*dk line3d
c**********************************************************************
cbeg                         line3d                                   *
c**********************************************************************
        subroutine line3d(x1,y1,z1,x2,y2,z2)
c......................................................................
c  line3d draws a (colored, dashed) line between the projections of 
c      the points (x1,y1,z1) and (x2,y2,z2).  It is a THREED utility,
c      and must be preceded by a SET3 call.
c......................................................................
c
        real xra(2),yra(2)
cncar
        call trn32t(x1,y1,z1,xra(1),yra(1),dummy,2)
        call trn32t(x2,y2,z2,xra(2),yra(2),dummy,2)
c
	do 10  i=1,2
	  xra(i) = cpux(ifix(xra(i)))
	  yra(i) = cpuy(ifix(yra(i)))
  10    continue
c
        call curved(xra,yra,2)
cncar
100	return
        end
cend line3d
c*dk line3s
c**********************************************************************
cbeg                         line3s                                   *
c**********************************************************************
        subroutine line3s(x1,y1,z1,x2,y2,z2)
c......................................................................
c  line3s draws a (colored, dashed) line between the projections of 
c      the points (x1,y1,z1) and (x2,y2,z2).  It is a SRFACE utility,
c      and must be preceded by a TRN32S or SRFACE call.  The line 
c      will not be hidden.  
c......................................................................
c
        real xra(2),yra(2)
cncar
        call trn32s(x1,y1,z1,xra(1),yra(1),dummy,1)
        call trn32s(x2,y2,z2,xra(2),yra(2),dummy,1)
c
        call curved(xra,yra,2)
cncar
100	return
        end
cend line3s
c*dk nicify
c**********************************************************************
cbeg                         nicify                                   *
c**********************************************************************
        subroutine nicify(xmin,xmax,delx)
c
	dimension nicenum(6)
	data nicenum / 20, 10, 5, 4, 2, 1 /
	data zfuzz / 1.e-10 /
c......................................................................
c Resets the boundaries of a graph so that it looks nice.
c......................................................................
c
	if(abs(xmax-xmin).lt.zfuzz)then
	  if(xmax.ne.0.)then
	    if(abs(xmax-xmin)/xmax.gt.zfuzz)go to 10
       	    xmin = 0.
	    xmax = 2.*xmax
	  else
	    xmin = 0.
	    xmax = 1.
	  endif
	endif
c
   10	if(xmax-xmin.lt.0.)then
	  temp = xmin
	  xmin = xmax
	  xmax = temp
	endif
c
        delx = 10**( aint( alog10(xmax-xmin) - 1. ) )
	if(xmax-xmin.lt.10.)delx = delx/10.
c
        if ( abs(amod( xmax-xmin, delx )) .lt. zfuzz ) go to 20
c
        if (abs(amod(xmin,delx)).gt.zfuzz) then
          xmin = aint( xmin/delx )*delx
          if (xmin.lt.0.) xmin = xmin - delx
        endif
c
        if (abs(amod(xmax,delx)).gt.zfuzz) then
          xmax = aint( xmax/delx + 1. )*delx
          if (xmax.lt.0.) xmax = xmax - delx
        endif
c
   20   mjrx = int( (xmax-xmin)/delx )
c
	do 30  i=1,5
	  if(mjrx.ge.4*nicenum(i))go to 100
   30   continue
c
  100   delx = delx * nicenum(i)
	return
        end
cend nicify
c*dk realchr
c**********************************************************************
cbeg                         realchr                                  *
c**********************************************************************
        subroutine realchr(realno,nopl,chreal,lchreal)
c......................................................................
c  realchr takes a real number realno and returns a string chreal with
c     the number in character form in the first lchreal positions.
c     nopl is the number of decimal places desired; if nopl.lt.0 the 
c     number will be in exponential form.
c......................................................................
	character*(*) chreal
	character*3 chexp
c
	inopl = nopl
	zrealno = realno
	iexpflag = 0
	nascii0 = ichar('0')
	zfuzz = 1.e-10
c
	if(inopl.lt.0)then
          inopl=iabs(inopl)
	  if(zrealno.eq.0.)go to 5
	  if(abs(zrealno).lt.1.)zfuzz=-zfuzz
	  iexp = int(log10(abs(zrealno))+zfuzz)
	  if(abs(zrealno).lt.1.)iexp=iexp-1
	  zrealno = zrealno / (10.**iexp)
	  iexpflag = 1
	endif
c
	if(inopl.eq.0)zrealno=nint(zrealno)
c                                              digits to left of dec pt
   5    call intchr(int(zrealno),chreal,lint)
	lchreal = lint
c
	if(inopl.eq.0)go to 100
c                                  (because intchr never returns '-0':)
c
	if(zrealno.lt.0. .and. zrealno.gt.-1.) then
	  chreal(2:lint+1) = chreal(1:lint)
	  chreal(1:1) = '-'
	  lint = lint + 1
	endif
c                                                        decimal places
	zrealno=abs(zrealno)
	chreal(lint+1:lint+1) = '.'
	if(inopl.ne.1)then
	  do 10  i=1,inopl-1
	    j = lint + i + 1
	    chreal(j:j) = char(nascii0+mod(int(zrealno*10.**i),10))
  10      continue
	endif
        j = lint + inopl + 1
	chreal(j:j) = char(nascii0+mod(nint(zrealno*10.**inopl),10))
        lchreal = j
c
	if(iexpflag.eq.0)go to 100
c                                               exponent for exp'l form
	call intchr(iexp,chexp,lexp)
c                                          (this assumes use of PLCHHQ)
cncar
	chreal(lchreal+1:lchreal+6) = '*10''S'''
cncar
	do 20  i=1,lexp
	  j = lchreal + 6 + i
	  chreal(j:j) = chexp(i:i)
  20    continue
	lchreal = j
c
 100	return
	end
cend realchr
c*dk setclr
c**********************************************************************
cbeg                      setclr                                      *
c**********************************************************************
	subroutine setclr(color)
c......................................................................
c     Emulates DISSPLA routine for setting color for lines & text.
c......................................................................
	common /thrint/ ithrmj, ithrmn, ithrtx
        character*(*) color
        character*3 clr
c
        ishade = 1
        nascii0 = ichar('0')
        do 10  i=4,12
          index = ichar(color(i:i))
          if( index.gt.nascii0 .and. index.le.nascii0 + 5) then
            ishade = index - nascii0
            go to 15
          endif
   10   continue
   15   continue
c
        clr = color(1:3)
        ncolor = 61
c
        if (clr.eq.'bla' .or. clr.eq.'blk')   ncolor = 0
        if (clr.eq.'whi' .or. clr.eq.'wht')   ncolor = 61
c
        if (clr.eq.'pur' .or. clr.eq.'prp')   ncolor =  0 + ishade
        if (clr.eq.'blu')                     ncolor =  5 + ishade
        if (clr.eq.'tur' .or. clr.eq.'trq')   ncolor = 10 + ishade
        if (clr.eq.'cya' .or. clr.eq.'cyn')   ncolor = 15 + ishade
        if (clr.eq.'bgr')                     ncolor = 20 + ishade
        if (clr.eq.'gre' .or. clr.eq.'grn')   ncolor = 25 + ishade
        if (clr.eq.'ygr')                     ncolor = 30 + ishade
        if (clr.eq.'yel' .or. clr.eq.'ylw')   ncolor = 35 + ishade
        if (clr.eq.'ora' .or. clr.eq.'orn')   ncolor = 40 + ishade
        if (clr.eq.'red')                     ncolor = 45 + ishade
        if (clr.eq.'mag' .or. clr.eq.'mgn')   ncolor = 50 + ishade
        if (clr.eq.'vio' .or. clr.eq.'vlt')   ncolor = 55 + ishade
cncar
	call sflush
	call gsplci(ncolor)
	call gsfaci(ncolor)
	call gstxci(ncolor)
	ithrmj = ncolor
	ithrmn = ncolor
	ithrtx = ncolor
cncar
        return
        end
cend setclr
c*dk setpat
c**********************************************************************
cbeg                      setpat                                      *
c**********************************************************************
	subroutine setpat(pattern)
c
        character*(*) pattern
c......................................................................
c     Emulates DISSPLA routine for setting dashed-line patterns.
c......................................................................
c
        ipat = 65535
c
	if (pattern.eq.'blank')     ipat =     0
        if (pattern.eq.'solid')     ipat = 65535
	if (pattern.eq.'longdash')  ipat = 61680
	if (pattern.eq.'dash')      ipat = 52428
	if (pattern.eq.'dot')       ipat = 43690
	if (pattern.eq.'chndot')    ipat = 58596
cncar
	call dashdb(ior(ipat,0))
	call agseti('DASH/PAT/1.',ipat)
cncar
        return
        end
c
cend setpat
c*dk thrdgraf
c**********************************************************************
cbeg                         thrdgraf                                 *
c**********************************************************************
        subroutine thrdgraf(xra,yra,zplot,m1dim,m2dim,
     &   labelplt,labelx,labely,labelz,labrun,
     &   ntty,niorec,x3dview,kscale)
c
        character*(*) labelplt,labelx,labely,labelz,labrun
	character*50 chlab,dumlblz,chzscale
        dimension zplot(*)
        dimension x3dview(*),xra(*),yra(*)
        dimension vuncar(6)
c                                                      Workspace arrays
        dimension work(100000),zzplot(50000),sxra(10000),syra(10000)
c
        common / srfip1 / ifr,   istp,  irots, idrx,  idry,
     &                    idrz,  iupper,iskirt,ncla,  theta,
     &                    hskirt,chi,   clo,   cinc,  ispval
c......................................................................
c       thrdgraf plots a three-dimensional plot of zplot, using
c    the ncar graphics library.  it was written by daniel heifetz,
c    pppl, x2620, in disspla for the program degraf and translated 
c    & modified by derek fox under supervision of daren stotler, 
c    pppl, 7/92.  this version is for more general use.
c
c    calling arguments
c    -----------------
c    xra(i)         array of x-coordinates (m1dim long)
c    yra(j)         array of y-coordinates (m2dim long)
c    zplot(i,j)     the m1dim by m2dim array to be plotted
c    m1dim          first dimension of zplot
c    m2dim          second dimension of zplot
c    labelplt(i)    title of plot
c    labelx(i)      x-axis label
c    labely(i)      y-axis label
c    labelz(i)      z-axis label
c    labrun(i)      run label
c    ntty           logical unit number for writing error messages
c    niorec         2nd logical unit number for writing error messages
c    x3dview(i)     viewing position coordinates (x,y,z), in units of 
c                   the plot frame box.  The box is 1 unit on a side, 
c                   with its (min x, min y, min z) corner at the 
c                   origin.
c    kscale         = 1: plot zplot
c                   = 2: plot log10(zplot)
c
c    note:
c
c    (1) labelplt, labelx, labely, labelz and labrun must end with 
c         the character "$" for proper placement.  
c    (2) labelx, labely, and labelz must be all-caps, standard ANSI
c        fortran 77 characters to be plotted.
c    (3) xra and yra must be monotonically increasing functions of
c        their indices.  
c......................................................................
c
c                                               Find lengths of strings
c                                   <nlblplt,nlblx,nlbly,nlblz,nlabrun>
        call findlen(labelplt,nlblplt)
        call findlen(labelx,nlblx)
        call findlen(labely,nlbly)
        call findlen(labelz,nlblz)
        call findlen(labrun,nlabrun)
c
c                                               Set xmin,xmax,ymin,ymax
        xmin = xra(1)
        xmax = xra(m1dim)
        ymin = yra(1)
        ymax = yra(m2dim)
c
 10     j12dim=m1dim*m2dim
        if((j12dim.eq.0).or.(j12dim.gt.50000))go to 1010
        do 12 i=1,50000
 12       zzplot(i)=0.0
c                              find minimum and maximum values of zplot
c                                                           <zmin,zmax>
        go to(20,30)kscale
c                                                 kscale=1: linear plot
 20     zmin=10000.
        zmax=-10000.
        do 25 i=1,j12dim
          zzplot(i)=zplot(i)
	  zmin=min(zmin,zplot(i))
 25       zmax=max(zmax,zplot(i))
        if(zmax.eq.0.)go to 1058
	call nicify(zmin,zmax,zdelz)
        go to 50
c                                                       log/linear plot
c                                                             take logs
 30     zmin=1000.
        zmax=-1000.
        do 34 i=1,j12dim
          zploti=abs(zplot(i))
          if(zploti.eq.0.)go to 32
          zzplot(i)=alog10(zploti)
          zmin=min(zmin,zzplot(i))
          zmax=max(zmax,zzplot(i))
 32       continue
 34     continue
        if(zmax.lt.zmin)go to 1058
c                                                     set zmax and zmin
        zmax=aint(zmax)
        zmin=aint(zmin)
        if (zmax.gt.0.) zmax=zmax+1.
        if (zmin.lt.0.) zmin=zmin-1.
	call nicify(zmin,zmax,zdelz)
        do 36 i=1,j12dim
 36       if((zzplot(i).eq.0.0).or.(zzplot(i).lt.zmin))zzplot(i)=zmin
c
c                                                     Make a cube shape
c                                                      <xmag,ymag,zmag>     
 50     zxtnt = zmax - zmin
        xxtnt = xmax - xmin
        yxtnt = ymax - ymin
        xmag = 1./xxtnt
        ymag = 1./yxtnt
        zmag = 1./zxtnt
c                                                    Remap data to grid
	do 52  i=1,m1dim
          sxra(i) = xra(i) * xmag
 52     continue
c
        do 54  j=1,m2dim
          syra(j) = yra(j) * ymag
 54     continue
c
        do 56  i=1,j12dim
          zzplot(i) = zzplot(i) * zmag
 56     continue
c
	xmn = xmin * xmag
	xmd = xmn + 0.5
	xmx = xmn + 1.
	ymn = ymin * ymag
	ymd = ymn + 0.5
	ymx = ymn + 1.
        zmn = zmin * zmag
	zmd = zmn + 0.5
        zmx = zmn + 1.
c
	vuncar(1) = xmn + x3dview(1)
	vuncar(2) = ymn + x3dview(2)
	vuncar(3) = zmn + x3dview(3)
	vuncar(4) = xmd
	vuncar(5) = ymd
	vuncar(6) = zmd
c                                                            Plot title
cncar
        call set(0.,1.,0.,1., 0.,1.,0.,1., 1)
        call setclr('cyan')
        call plchhq(0.5,0.95, labelplt(1:nlblplt), 0.015, 0., 0.)
c
c                                             Check for standard errors
c
        zview = x3dview(1)**2 + x3dview(2)**2 + x3dview(3)**2
        if(zview.le.0.0)go to 1050
        if(ymx.le.0.0)go to 1056
        if(kscale.eq.2)then
          if(zmax.le.zmin)go to 1058
          if(zdelz.le.0.0)go to 1059
        endif
c                                                           Plot zzplot
        call setclr('blue')
        ifr = 0
	call set(0.,1.,0.,1., 1.,1024.,1.,1024., 1)
        call srfgk(sxra,syra,zzplot,work,m1dim,m1dim,m2dim,vuncar,0.)
c
c                                                      Plot axis labels
        call setclr('cyan')
	call lblpos(vuncar,ixydis,iyxdis,izxdis,izydis)
	xlaby = ymd + 0.50*ixydis
	ylabx = xmd + 0.50*iyxdis
c
	if(izxdis*izydis.gt.0)then
	  izxflag = 0
	  izyflag = izydis
	  iztop = izxdis + izydis
	  zlabx = xmd + 0.50*izxdis
	  zlaby = ymd + 0.64*izydis
	else
	  izxflag = izxdis
	  izyflag = 0
	  iztop = 2*izxdis + izydis
	  zlabx = xmd + 0.64*izxdis
	  zlaby = ymd + 0.50*izydis
	endif
c
c                                                        Box the graph:
c                                         ...first always-visible lines
c
	izdis = nint(sign(1.,vuncar(3)-vuncar(6)))
c
        call line3s(xmx,ymn,zmd+izdis/2.,xmx,ymx,zmd+izdis/2.)
        call line3s(xmn,ymn,zmd+izdis/2.,xmx,ymn,zmd+izdis/2.)
        call line3s(xmn,ymn,zmd+izdis/2.,xmn,ymx,zmd+izdis/2.)
        call line3s(xmn,ymx,zmd+izdis/2.,xmx,ymx,zmd+izdis/2.)
c
	call line3s(xmd+iyxdis/2.,ymd+ixydis/2.,zmd+izdis/2.,
     &     xmd+iyxdis/2.,ymd+ixydis/2.,zmd-izdis/2.)
	call line3s(xmd+iyxdis/2.,ymd-ixydis/2.,zmd+izdis/2.,
     &     xmd+iyxdis/2.,ymd-ixydis/2.,zmd-izdis/2.)
	call line3s(xmd-iyxdis/2.,ymd+ixydis/2.,zmd+izdis/2.,
     &     xmd-iyxdis/2.,ymd+ixydis/2.,zmd-izdis/2.)
c
	call line3s(xmd+iyxdis/2.,ymd+ixydis/2.,zmd-izdis/2.,
     &     xmd+iyxdis/2.,ymd-ixydis/2.,zmd-izdis/2.)
	call line3s(xmd+iyxdis/2.,ymd+ixydis/2.,zmd-izdis/2.,
     &     xmd-iyxdis/2.,ymd+ixydis/2.,zmd-izdis/2.)
c
c                                                 ...then hidden lines
c
	call hlin3s(xmd-iyxdis/2.,ymd-ixydis/2.,zmd-izdis/2.,
     &     xmd-iyxdis/2.,ymd-ixydis/2.,zmd+izdis/2.)
	call hlin3s(xmd-iyxdis/2.,ymd-ixydis/2.,zmd-izdis/2.,
     &     xmd-iyxdis/2.,ymd+ixydis/2.,zmd-izdis/2.)
	call hlin3s(xmd-iyxdis/2.,ymd-ixydis/2.,zmd-izdis/2.,
     &     xmd+iyxdis/2.,ymd-ixydis/2.,zmd-izdis/2.)
c
c                                            Plot ticks, numeric labels
	tklen = 0.02
	axslab = 2. * tklen
	znlabx = xmd + 0.50*izxdis
	znlaby = ymd + 0.50*izydis
c                                                                z-axis
c
c                              (a. find proper decimal form for labels)
c
	nodec = 1
	dumlblz(1:nlblz) = labelz(1:nlblz)
c
        if(zdelz.gt.100. .or. zdelz.lt.0.01)then
          izexp = int(alog10(zdelz))
	  call intchr(iabs(izexp),chzscale,lchzscale)
	  if(zdelz.gt.1.)dumlblz(nlblz+1:nlblz+7) = ' (/10**'
	  if(zdelz.lt.1.)dumlblz(nlblz+1:nlblz+7) = ' (*10**'
	  dumlblz(nlblz+8:nlblz+lchzscale+7) = chzscale(1:lchzscale)
	  nlblz = nlblz+lchzscale+8
	  dumlblz(nlblz:nlblz) = ')'
	else
	  izexp = 0
	  if(zdelz.ge.5.0)nodec=0
          if(zdelz.lt.0.5)nodec=2
        endif
c                                                 (b. write the labels)
	do 65  iztk=0,int((zmax-zmin)/zdelz)
	  tkz = zmin + float(iztk)*zdelz 
	  zlbl = tkz / 10.**float(izexp)
	  call line3s(znlabx,znlaby,zmag*tkz,
     &      znlabx+tklen*izxflag,znlaby+tklen*izyflag,zmag*tkz)
	  call realchr(zlbl,nodec,chlab,lchlab)
	  call pwrzs(znlabx+axslab*izxflag,znlaby+axslab*izyflag,
     &      zmag*tkz, chlab(1:lchlab), lchlab, 15, 3, iztop, 0)
 65     continue
c
        call pwrzs(zlabx,zlaby,zmd,
     &       dumlblz(1:nlblz), nlblz, 20, 3, iztop, 0)
c                                                                x-axis
        call pwrzs(xmd,xlaby,zmn-0.14, 
     &       labelx(1:nlblx), nlblx, 20, -ixydis, 3, 0)
	call line3s(xmn,xlaby,zmn, xmn,xlaby,zmn-tklen)
	call intchr(nint(xmin),chlab,lchlab)
	call pwrzs(xmn,xlaby,zmn-axslab, ' '//chlab(1:lchlab)//' ',
     &    lchlab+2, 15, -1*ixydis, 3, ixydis)
c
	call line3s(xmx,xlaby,zmn, xmx,xlaby,zmn-tklen)
	call intchr(nint(xmax),chlab,lchlab)
	call pwrzs(xmx,xlaby,zmn-axslab, ' '//chlab(1:lchlab)//' ',
     &    lchlab+2, 15, -1*ixydis, 3, -ixydis)
c                                                                y-axis
        call pwrzs(ylabx,ymd,zmn-0.14, 
     &       labely(1:nlbly), nlbly, 20, 2*iyxdis, 3, 0)
	call line3s(ylabx,ymn,zmn, ylabx,ymn,zmn-tklen)
	call intchr(nint(ymin),chlab,lchlab)
	call pwrzs(ylabx,ymn,zmn-axslab, ' '//chlab(1:lchlab)//' ',
     &    lchlab+2, 15, 2*iyxdis, 3, -iyxdis)
c
	call line3s(ylabx,ymx,zmn, ylabx,ymx,zmn-tklen)
	call intchr(nint(ymax),chlab,lchlab)
	call pwrzs(ylabx,ymx,zmn-axslab, ' '//chlab(1:lchlab)//' ',
     &    lchlab+2, 15, 2*iyxdis, 3, iyxdis)
c
c
        call setpat('solid')
c                                                        Plot run label
cncar
 90     call set(0.,1.,0.,1., 0.,1.,0.,1., 1)
        call setclr('cyan')
        call plchhq(0.5,0.92, labrun(1:nlabrun), 0.008, 0., 0.)
c                                                  reset srface param's
        ifr = 1
cncar
 99     return
c----------------------------------------------------------------------
c                                                                 error
 1010   write(ntty,10010)labelplt(1:40),j12dim
        write(niorec,10010)labelplt(1:40),j12dim
        go to 1099
 1050   write(ntty,10050)labelplt(1:40)
        write(niorec,10050)labelplt(1:40)
        go to 1099
 1055   write(ntty,10055)labelplt(1:40),xmx,xmn
        write(niorec,10055)labelplt(1:40),xmx,xmn
        go to 1099
 1056   write(ntty,10056)labelplt(1:40),ymx
        write(niorec,10056)labelplt(1:40),ymx
        go to 1099
 1058   write(ntty,10058)labelplt(1:40),zmax,zmin
        write(niorec,10058)labelplt(1:40),zmax,zmin
        go to 1099
 1059   write(ntty,10059)labelplt(1:40),zdelz
        write(niorec,10059)labelplt(1:40),zdelz
 1099   return
c----------------------------------------------------------------------
c                                                     format statements
10010   format(x,a40,x,'error in thrdgraf: m1dim*m2dim =',x,i7,
     &   x,'.eq.0.or.gt.50000')
10050   format(x,a40,x,'error in thrdgraf: abs(x3dview).le.0')
10055   format(x,a40,x,'error in thrdgraf: xmx=',1pg9.2,'.le.xmn=',
     &   g9.2)
10056   format(x,a40,x,'error in thrdgraf: ymx=',1pg9.2,'.le.0.0')
10058   format(x,a40,x,'error in thrdgraf: zmax=',1pg9.2,'.le.zmin=',
     &   g9.2)
10059   format(x,a40,x,'error in thrdgraf: zdelz=',1pg9.2,'.le.0.0')
        end
cend thrdgraf
cncar
c*dk trn32s
c**********************************************************************
cbeg                           trn32s                                 *
c**********************************************************************
      SUBROUTINE TRN32S (X,Y,Z,XT,YT,ZT,IFLAG)
cmod
c......................................................................
c This routine is a slight modification of the ncar routine, designed 
c so that SRFACE uses a reduced area for its plot, and labels and lines
c drawn with line3s don't run up against a boundary.
c
c Modifications are bracketed by 'cmod' lines, and are in lower case.
c......................................................................
cmod
C
C THIS ROUTINE IMPLEMENTS THE 3-SPACE TO 2-SPACE TRANSFOR-
C MATION BY KUBER, SZABO AND GIULIERI, THE PERSPECTIVE
C REPRESENTATION OF FUNCTIONS OF TWO VARIABLES. J. ACM 15,
C 2, 193-204,1968.
C IFLAG=0 ARGUMENTS
C X,Y,Z    ARE THE 3-SPACE COORDINATES OF THE INTERSECTION
C          OF THE LINE OF SIGHT AND THE IMAGE PLANE.  THIS
C          POINT CAN BE THOUGHT OF AS THE POINT LOOKED AT.
C XT,YT,ZT ARE THE 3-SPACE COORDINATES OF THE EYE POSITION.
C
C IFLAG=1 ARGUMENTS
C X,Y,Z    ARE THE 3-SPACE COORDINATES OF A POINT TO BE
C          TRANSFORMED.
C XT,YT    THE RESULTS OF THE 3-SPACE TO 2-SPACE TRANSFOR-
C          MATION.
C          USE IFIX(XT) AND IFIX(YT) IN GPL CALLS.
C ZT       NOT USED.
C IF LL (IN COMMON) =0 XT AND YT ARE IN THE SAME SCALE AS X, Y, AND Z.
C
      COMMON /PWRZ1S/ XXMIN      ,XXMAX      ,YYMIN      ,YYMAX      ,
     1                ZZMIN      ,ZZMAX      ,DELCRT     ,EYEX       ,
     2                EYEY       ,EYEZ
      COMMON /SRFBLK/ LIMU(1024) ,LIML(1024) ,CL(41)     ,NCL        ,
     1                LL         ,FACT       ,IROT       ,NDRZ       ,
     2                NUPPER     ,NRSWT      ,BIGD       ,UMIN       ,
     3                UMAX       ,VMIN       ,VMAX       ,RZERO      ,
     4                IOFFP      ,NSPVAL     ,SPVAL      ,BIGEST
      DIMENSION       NLU(7)     ,NRU(7)     ,NBV(7)     ,NTV(7)
C
C SAVE INSERTED BY BEN DOMENICO 9/8/85 BECAUSE OF ASSUMPTION THAT
C   JUMP, JUMP2, AND JUMP3 ARE PRESERVED BETWEEN CALLS.
C   THERE MAY BE OTHER SUCH ASSUMPTIONS AS WELL.
C
      SAVE
C
C PICTURE CORNER COORDINATES FOR LL=1
C
cmod
      data nlu(1),nru(1),nbv(1),ntv(1)/  200,824,  140,764/
c      DATA NLU(1),NRU(1),NBV(1),NTV(1)/  10,1014,  10,1014/
cmod
C
C PICTURE CORNER COORDINATES FOR LL=2
C
      DATA NLU(2),NRU(2),NBV(2),NTV(2)/  10, 924,  50, 964/
C
C PICTURE CORNER COORDINATES FOR LL=3
C
      DATA NLU(3),NRU(3),NBV(3),NTV(3)/ 100,1014,  50, 964/
C
C PICTURE CORNER COORDINATES FOR LL=4
C
      DATA NLU(4),NRU(4),NBV(4),NTV(4)/  10,1014,  10,1014/
C
C PICTURE CORNER COORDINATES FOR LL=5
C
      DATA NLU(5),NRU(5),NBV(5),NTV(5)/  10,1014,  10,1014/
C
C PICTURE CORNER COORDINATES FOR LL=6
C
      DATA NLU(6),NRU(6),NBV(6),NTV(6)/  10, 512, 256, 758/
C
C PICTURE CORNER COORDINATES FOR LL=7
C
      DATA NLU(7),NRU(7),NBV(7),NTV(7)/ 512,1014, 256, 758/
C
C STORE THE PARAMETERS OF THE SET32 CALL FOR USE WHEN
C TRN32 IS CALLED.
C
      IF (IFLAG)  40, 10, 40
   10 CONTINUE
      ASSIGN  60 TO JUMP3
      IF (IOFFP .EQ. 1) ASSIGN  50 TO JUMP3
      AX = X
      AY = Y
      AZ = Z
      EX = XT
      EY = YT
      EZ = ZT
C
C AS MUCH COMPUTATION AS POSSIBLE IS DONE DURING EXECUTION
C THIS ROUTINE WHEN IFLAG=0 BECAUSE CALLS IN THAT MODE ARE INFREQUENT.
C
      DX = AX-EX
      DY = AY-EY
      DZ = AZ-EZ
      D = SQRT(DX*DX+DY*DY+DZ*DZ)
      COSAL = DX/D
      COSBE = DY/D
      COSGA = DZ/D
      SINGA = SQRT(1.-COSGA*COSGA)
      ASSIGN 120 TO JUMP2
      IF (LL .EQ. 0) GO TO  20
      ASSIGN 100 TO JUMP2
      DELCRT = NRU(LL)-NLU(LL)
      U0 = UMIN
      V0 = VMIN
      U1 = NLU(LL)
      V1 = NBV(LL)
      U2 = NRU(LL)-NLU(LL)
      V2 = NTV(LL)-NBV(LL)
      U3 = U2/(UMAX-UMIN)
      V3 = V2/(VMAX-VMIN)
      U4 = NRU(LL)
      V4 = NTV(LL)
      IF (NRSWT .EQ. 0) GO TO  20
      U0 = -BIGD
      V0 = -BIGD
      U3 = U2/(2.*BIGD)
      V3 = V2/(2.*BIGD)
C
C THE 3-SPACE POINT LOOKED AT IS TRANSFORMED INTO (0,0) OF
C THE 2-SPACE.  THE 3-SPACE Z AXIS IS TRANSFORMED INTO THE
C 2-SPACE Y AXIS.  IF THE LINE OF SIGHT IS CLOSE TO PARALLEL
C TO THE 3-SPACE Z AXIS, THE 3-SPACE Y AXIS IS CHOSEN (IN-
C STEAD OF THE 3-SPACE Z AXIS) TO BE TRANSFORMED INTO THE
C 2-SPACE Y AXIS.
C
   20 IF (SINGA .LT. 0.0001) GO TO  30
      R = 1./SINGA
      ASSIGN  70 TO JUMP
      RETURN
   30 SINBE = SQRT(1.-COSBE*COSBE)
      R = 1./SINBE
      ASSIGN  80 TO JUMP
      RETURN
   40 CONTINUE
      XX = X
      YY = Y
      ZZ = Z
      GO TO JUMP3,( 50, 60)
   50 IF (ZZ .EQ. SPVAL) GO TO 110
   60 Q = D/((XX-EX)*COSAL+(YY-EY)*COSBE+(ZZ-EZ)*COSGA)
      GO TO JUMP,( 70, 80)
   70 XX = ((EX+Q*(XX-EX)-AX)*COSBE-(EY+Q*(YY-EY)-AY)*COSAL)*R
      YY = (EZ+Q*(ZZ-EZ)-AZ)*R
      GO TO  90
   80 XX = ((EZ+Q*(ZZ-EZ)-AZ)*COSAL-(EX+Q*(XX-EX)-AX)*COSGA)*R
      YY = (EY+Q*(YY-EY)-AY)*R
   90 GO TO JUMP2,(100,120)
cmod
  100 xx = amin1(1009.,amax1(15.,u1+u3*(fact*xx-u0)))
      yy = amin1(1009.,amax1(15.,v1+v3*(fact*yy-v0)))
c      XX = AMIN1(U4,AMAX1(U1,U1+U3*(FACT*XX-U0)))
c      YY = AMIN1(V4,AMAX1(V1,V1+V3*(FACT*YY-V0)))
      GO TO 120
cmod
  110 XX = NSPVAL
      YY = NSPVAL
C
  120 XT = XX
      YT = YY
      RETURN
      END
cend trn32s
cncar
c*dk twodgraf
c**********************************************************************
cbeg                       twodgraf                                   *
c**********************************************************************
        subroutine twodgraf(xplot,yplot,kpts,kcurves,
     &   kcrvclr,kcrvtype,kbackclr,delymax,labelplt,
     &   labelx,labely,labrun,ntty,niorec,kscale,yploterr,kerrplot)
c
        character*(*) kcrvclr(*),kcrvtype(*),kbackclr,labelplt,
     &   labelx,labely,labrun
        character*9 chjcrv
        real lgnd(4),grph(4),grid(4),ugrid(4)
        dimension yplot(kpts*kcurves),xplot(kpts),yploterr(kpts),
     &    xeplot(2),yeplot(2),xlplot(2),ylplot(2)
c                                                      Workspace arrays
        dimension yyplot(5000)
c
c        data grph / 0.0, 1.0, 0.11, 0.88 /
        data grph / 0.03,.97, 0.35, 0.98 /
        data grid / 0.12, 0.9, 0.12, 0.9 /
        data lgnd / 0.91, 1., 0., 0.80 /
c......................................................................
c     twodgraf plots kcurves curves from yplot vs. xplot on a 
c   linear/linear (kscale=1) or log/linear (kscale=2) graph,
c   using the NCAR library's AUTOGRAPH package.  
c
c   calling arguments
c   -----------------
c
c   xplot(i)       array of (kpts) x-coordinates for the points
c   yplot(i,j)     kpts by kcurves array of y-coordinates
c   kpts           number of points per curve
c   kcurves        number of curves
c   kcrvclr(i)     (kcurves) array of color names for the curves
c   kcrvtype(i)    (kcurves) array of pattern names for the curves
c   kbackclr       color name for grid, labels, background
c   delymax        maximum number of decades allowed for log scale
c   labelplt       title for the plot (*)
c   labelx         x-axis label (*)
c   labely         y-axis label (*)
c   labrun         "run label" -- printed in small type at bottom (*)
c   ntty           logical unit number for writing error messages
c   niorec         2nd logical unit number for writing error messages
c   kscale         = 1: linear scale
c                  = 2: semi-log scale (log y-axis)
c   yploterr(i,j)  kpts by kcurves array of y error values
c   kerrplt        = 0: no error bars
c                  = 1: error bars
c
c   note:
c
c     (*)'d strings must end in '$' to be centered properly.  
c......................................................................
c
c                                        Find the length of all strings
c                                         <nlblplt,nlblx,nlbly,nlabrun>
        call findlen(labelplt,nlblplt)
        call findlen(labelx,nlblx)
        call findlen(labely,nlbly)
        call findlen(labrun,nlabrun)
c                                  Find value of autograph "null" param
cncar
	call aggetf('NUL/1.',spval)
cncar
c
   10   zfuzz = 1.e-4
        jtotpts=kpts*kcurves
        if(jtotpts.gt.5000)go to 1010
c                                                       Find ymin, ymax
	ymax=-1.e30
	ymin=+1.e30
c
        do 20  i=1,jtotpts
	  if(kscale.eq.2 .and. yplot(i).eq.0)then
	    yyplot(i) = spval
	  else
            yyplot(i) = yplot(i)
            if(yyplot(i).eq.spval)go to 20
	    if(kscale.eq.2)yyplot(i)=abs(yyplot(i))
   15       if(kerrplot.eq.0)then
	      ymin=min(ymin,yyplot(i))
	      ymax=max(ymax,yyplot(i))
	    else
	      ymin=min(ymin,yyplot(i)-yploterr(i))
	      ymax=max(ymax,yyplot(i)+yploterr(i))
	    endif
	  endif
   20   continue
c
        if(ymin.ge.ymax)go to 1020
        if(abs(ymax-ymin)/ymax.lt.zfuzz)go to 1020
c
	if(kscale.eq.2 .and. (ymin*10.**delymax).gt.ymax)then
	  ymin = 10.**(int(log10(ymax)-delymax))
	  if(ymax.gt.1.)ymin=ymin*10.
	  do 40  i=1,jtotpts
	    if(yyplot(i).eq.spval)go to 40
	    yyplot(i)=max(ymin,yyplot(i))
   40     continue
	endif
c                                                        Plot yplot...
cncar
   50   call set(0.,1.,0.,1., 0.,1.,0.,1., 1)
        call setclr(kbackclr)
c                                                   Write title, labels
c
        titley = grph(4) - 0.015
        call plchhq(0.5, titley, labelplt(1:nlblplt), 0.015, 0., 0.)
        xlblx = grph(1) + (grph(2)-grph(1))*(grid(1)+grid(2))/2
        xlbly = grph(3) + 0.55*(grph(4)-grph(3))*grid(3)
        call plchhq(xlblx, xlbly, labelx(1:nlblx), 0.01, 0., 0.)
        ylblx = grph(1) + 0.40*(grph(2)-grph(1))*grid(1)
        ylbly = grph(3) + (grph(4)-grph(3))*(grid(3)+grid(4))/2
        call plchhq(ylblx, ylbly, labely(1:nlbly), 0.01, 90., 0.)
c
c                                            Graph grid, numeric labels
c
        call agseti('FRAME.',2)
        call agsetp('GRAPH.',grph,4)
c        call agsetp('GRID.',grid,4)
	call agsetf('Y/MIN.',ymin)
	call agsetf('Y/MAX.',ymax)
        if(iabs(kscale).eq.1)call agseti('Y/LOG.',0)
        if(kscale.eq.2)call agseti('Y/LOG.',1)
c
        call anotat('$','$',0,0,0,0)
        call agseti('BACKGND.',2)
	call agstup(xplot,1,0,kpts,1,yplot,0,0,0,0)
	call agback
c                                               Increase polyline width
 60     call gslwsc(2.)
c                                                        Graph curve(s)
	nocrvs=kcurves
	do 65 jcurve=1,kcurves
          indexx=kpts*(jcurve-1)+1
	  if(kcrvtype(jcurve).eq.'blank')then
	    nocrvs = nocrvs - 1
	    go to 65
	  endif
          call setpat(kcrvtype(jcurve))
          call setclr(kcrvclr(jcurve))
	  call agcurv(xplot,1,yyplot(indexx),1,kpts,1)
 65     continue
c
	call setpat('solid')
c                                                      graph error bars
        if(kerrplot.ne.0)then
          plxmn = grph(1) + (grph(2)-grph(1))*grid(1)
          plxmx = grph(1) + (grph(2)-grph(1))*grid(2)
          plymn = grph(3) + (grph(4)-grph(3))*grid(3)
          plymx = grph(3) + (grph(4)-grph(3))*grid(4)
          call aggetp('SEC/USER.',ugrid,4)
          call set(plxmn,plxmx,plymn,plymx, ugrid(1),ugrid(2),
     &             ugrid(3),ugrid(4), kscale)
	  zdxbar=0.007*(ugrid(2)-ugrid(1))
c
          do 75 jcurve=1,kcurves
            call setclr(kcrvclr(jcurve))
            do 70 jpt=1,kpts
              indexx=kpts*(jcurve-1)+jpt
	      if(yyplot(indexx).eq.spval)go to 70
c                                                          vertical bar
              xeplot(1)=xplot(jpt)
              xeplot(2)=xeplot(1)
              yeplot(1)=yyplot(indexx)-yploterr(indexx)
              yeplot(2)=yyplot(indexx)+yploterr(indexx)
	      if(kscale.eq.2 .and. yeplot(1).le.0.)yeplot(1)=ymin
              call curve(xeplot,yeplot,2)
c                                                    top horizontal bar
              xeplot(1)=xplot(jpt)-zdxbar
              xeplot(2)=xeplot(1)+2.0*zdxbar
              ztemp=yeplot(1)
              yeplot(1)=yeplot(2)
              call curve(xeplot,yeplot,2)
c                                                 bottom horizontal bar
              yeplot(2)=ztemp
              yeplot(1)=yeplot(2)
              call curve(xeplot,yeplot,2)
 70         continue
 75       continue
        endif
c                                              return to thin polylines
        call gslwsc(1.)
c                                                         set up legend
        call set(0.,1.,0.,1., 0.,1.,0.,1., 1)
        if(kcurves.eq.1)go to 80
c
          lgnd(3) = lgnd(4) - 0.02*nocrvs - 0.005
c 
          xlplot(1) = lgnd(1) + 0.005
          xlplot(2) = xlplot(1) + 0.025
c                                                       text for legend
	  yloff = lgnd(4) + 0.01
          do 78 jcurve=1,kcurves
	    if(kcrvtype(jcurve).eq.'blank')go to 78
	    call intchr(jcurve,chjcrv,njcrv)
            yloff = yloff - 0.02
	    call setclr(kbackclr)
            call plchhq(xlplot(2)+0.005,yloff, chjcrv(1:njcrv),
cccc            call plchhq(xlplot(2)+0.005,yloff, kcrvtype(jcurve),
cccc     &       0.008, 0., -1.)
     &       0.008, 0., -1.)
c                            draw (wide) lines for legend
            call gslwsc(2.)
            ylplot(1) = yloff+0.007
            ylplot(2) = yloff-0.007
	    call setpat(kcrvtype(jcurve))
cccc	    call setpat('solid')
            call setclr(kcrvclr(jcurve))
            call curved(xlplot,ylplot,2)
            call gslwsc(1.)
 78     continue
c                                                       print run label
 80     call setclr(kbackclr)
	call setpat('solid')
        rnlblx = grph(1) + (grph(2)-grph(1))*grid(1)
        rnlbly = grph(3)+0.01
        call plchhq(rnlblx,rnlbly,labrun(1:nlabrun),0.008,0.,-1.)
c
c                                                              clean up
        call agseti('BACKGND.',1)
	call agseti('Y/LOG.',0)
	call agsetf('Y/MIN.',1.e36)
	call agsetf('Y/MAX.',1.e36)
	call agseti('FRAME.',1)
cncar
c
 99     return
c----------------------------------------------------------------------
c                                                               error
 1010   write(ntty,10010)labelplt(1:40),kpts
        write(niorec,10010)labelplt(1:40),kpts
 1011 	return
c
 1020   write(ntty,10020)labelplt(1:40),ymin,ymax
	write(niorec,10020)labelplt(1:40),ymin,ymax
 1021   return
c
 1030   write(ntty,10030)labelplt(1:40),xmin,xmax
	write(niorec,10020)labelplt(1:40),xmin,xmax
 1031   return
c----------------------------------------------------------------------
c                                                     format statements
10010   format(x,a40,x,'error in twodgraf: kpts=',x,i5,
     &   x,'.gt.5000')
10020	format(x,a40,x,'error in twodgraf: ymin=',1pg9.2,
     &   ' and ymax=',g9.2)
10030	format(x,a40,x,'error in twodgraf: xmin=',1pg9.2,
     &   ' and xmax=',g9.2)
c
        end
cend twodgraf
