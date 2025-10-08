Simple instructions on how to run NOVA-C (former NOVA-K) postprocessor which is version 3.
Some typical problems with NOVA are attempted to be described in the
file uFAQ in the source directory.
From v.3.0 FLR are accounted for through k_theta (dominant) and k_radial 
contributions to the Bessel function. 

(A PREFERRABLE VERSION OF THIS README EXISTS IN GOOGLE DRIVE, 
MY Drive/Epresearch/NOVAReadme TO BE SYNCED WITH THIS FILE.
Please use that Readme for all edits.)

1. compiling the NOVA-C code (skip to 2. if you don't need it)
   from your working directory run this:
   gmake -f /u/ngorelen/SRC/NOVAK3/mkorb
   if you need to use some capabilities of NOVA-C which are under 
   development add a extra flag 
   gmake -f /u/ngorelen/SRC/NOVAK3/mkorb DB="-D _dev" 
   but some other capabilities can be turn on, read "mkorb" comments for details

   In building it minimum external libraries are required. I use pgplot
   graphics, so that pgplot library is recommended for some
   plot. However the calls to it can be turned off safely without
   effecting the results. 

2. I assume that the mode structure file exists. In my case this is
   egn03w.1438E+01. This is the default case from NOVA run. Look at
   Readme file of NOVA for details on how to get this mode structure
   file. Link the mode structure file to the standard file name NOVA-C
   uses. !!! dimensions are not recorded into eigenf! Must be consistent 
   with the production lengths! such as in clichlmax: minm_1=-3,maxm_1=20 !!!
   
   ln -sf egn03w.1438E+01 eigenf

3. Run the executable file gotae: 

  ./gotae

4. To see the extended chirping criterium contour map do this:
	after running gotae run the python script like this
	/u/ngorelen/work/exe/contmap4u.py Out/contmap.u4m

after this it generates the default parameter file for NOVA-C, which is
NOVAK_param

other parameter file is clich1b; it contains mass, charge of EPs; in
support of RBQ we keep there CONDEF, CONDE, CONDI, one fluid eff,
electr, ion thermal ion conductivities.

It contains all the information required to run the the NOVA-C. Check
this file change the appropriate parameters to what you need (it has
comments) and rerun ./gotae to get the stability results.

In particular make sure that n number is what you need. Default case
should be ntor=3. 

If transp.dat (taefl.txt) is read then it is used to read rmaj,amin,B0

Comments on the types of the distribution function (dtype in NOVAK_param
file) are in the taem.f of the source file. They are also at the end of
this document. 

The results should be similar/same as in the file

/u/ngorelen/work/nova/zma/smpl_out, which follows

 gam_ecoll
 -1.215388686420203E-03
 gam_eLandau
 -2.396208215541913E-04
 gam_DLandau
 -8.399222380960909E-09
 gam_TLandau
 -3.392252017357413E-16
 gam_HLandau
 -9.404220052485388E-13
 gam_CLandau
 0.000000000000000E+00
 radiative damping for modes in the TAE gap
 -1.683368854538708E-02
 Fast ion beta, beta_h
 8.000000000000000E-02
 >Finite Orbit Width res. (FOW) for fast ion specie:
  Mass/m_p
 2.000000000000000
  Charge/e_p
 1.000000000000000
  Velocity in 10^9 cm/sec
 0.9784043131548431
  Va/Valfv=
 1.008564299313750
 Fast ion growth rate w/o FLR, gam_h
 -5.378821842684051E-02
 Fast ion growth rate withFLR, gam_hFLR
 -3.885640119857010E-02
 Critical beta is to be calculated
 beta_hcrit=-beta_h*(all dampings above)/gam_h or 
 beta_hcrit=-beta_h*(all dampings above)/gam_hFLR or 

5. ultimate result can be checked with the transp run for TFTR case,
   which is published in ~/work/PAPERS/Gorelenkov_PoP99_nlin.pdf

   On the cluster all the results are in ~/work/nova/tftr_103101
 --- or checked it with the d3d_159243 case published in ~/work/PAPERS/Taimourzadeh_NF18aev2.pdf

********************************************************************************
Appendix I. comments of the use of different distribution functions of
fast ions. 

c The argument tip is related to the distribution function type.
c It is used in the subroutine distrifun (file orbit2d.f)
c The following explanation is given in taem.f file.
c
c     's' for Slowing down distribution, also
c     'd' for Slowing down distribution with v_cr = v_a0 / 3.
c     'e' Maxwellian for test with jet equilibrium 
c         (It'll switch to im1= 0, e.g.TAE's automatically)
c     't' slowing down for test with tftr equilibrium
c     'j' slowing down for test with jet equilibrium
c     'm' Maxwellian with standard input taken from transp.dat, which is 
c     	  determined by (UBPRP_D/2+UBPAR_D) for deuterium beams +++
c                           (UBPRP_T/2+UBPAR_T) for tritium beams ???
c                           (UFIPP/2+UFIPA) for fusion ions 	  ???
c                           ((UFASTPP-UBPRP_D-UBPRP_T-UFIPP)/2+(UFASTPA-UBPAR_D-UBPAR_T-UFIPA)) ???
c                           for H-minority fast ions if ihsps=3 and its mass set to 1		
c     'x' Maxwellian with thermal ions for proper D/T Landau damping calculations 
c         (add im1=3 for rot, which maybe important for resonance)
c         if ihsps=1, i.e. deuterium specie is prescribed. 
c   Below is for Gaussian pitch distribution centered at p=P0GA
c   and having width delta_p = DPGA, p=mu*Baxis/E (sometimes in the
c   literature refered as =lambda). 
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
c        extra input is used here: chidelt(4)=cbot (see BOTfeature.dvi for explanations and derivation)
c	       	     	     	   		    = fraction of EP pressure in the BOT feature
c				   chidelt(5)=delta = width of BOT feature in norm. velocity
c				   chidelt(6)=x0    = BOT centr. velocity normalized to EP birth velocity
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
c                     3 is like 0 but with the rotation included in the parametric w_rot dependence
c                     10is a special run for EP redistribution due to TAE activity


5.a To see DF plot do this:
    gmake -f ~/SRC/NOVAK3/mkorb DB="-D dumpdf"
    run gotae3
    gmake -f ~/SRC/DATA/mkplot FILE=plt_df_nova && ./xplt_df_nova
    it will plot DF vs Pphi, mu coutour map.

-------------=================  Production runs  ================---------------

10. cp /u/ngorelen/work/exe/runeqAnis.pyt .
    ln -sf transp_t153.dat transp.dat
    gmake -f ~/SRC/Equil/mkeqsolAnis
    make mnths0=158 inside gridparam
    ./runeqAnis.pyt
11. cp /u/ngorelen/work/exe/runmap.scr .
    ./runmap.scr
12. cp /u/ngorelen/work/exe/runcont.scr .
    ./runcont.scr
13. cp /u/ngorelen/work/exe/runInitStabN.py ./runInitStab.py
    Look for: 
    Nmin = 1?
    Nmax = 10?
    qa = 11?
    nmbrgrd=151? inside gridparam
    line=line.replace('= 3', '= 10')? make it 40 for om2 upper range from plotcont
    line=line.replace('= 1', '= 3')?  --//--
14. cp /u/ngorelen/work/exe/runstab_Li_Ep.scr .
    &inp iden=4 ntor=${i} lmax=32 plmax=50. 15. 5. ns = 3 151 1 / ? 151 is max gridparam dims
    @ l0=54 ? 79 for st40 to make chi0=0.79
    @ l1=55 ? ....
15. cp /u/ngorelen/work/exe/runStabNovak.py .
    Nmin = 1?
    Nmax = 2?
    qa = 11? m.b. as in 13. for the same eigenf length
    if iline == 45 : line=' 0.200E+01 0.300E+01 0.100E+01' ? add for H beams if needed
    if iline == 48 : line=' 0.100E+01 0.100E+01 0.100E+01' ? for H beams --//--
    if iline == 51 : line=line.replace('   1', '   1')? or make it 3 for alphas or H beams
#        if iline == 24 : line=line.replace('SPS = d', 'SPS = a')? only for alphas
#        if iline == 31 : line=line.replace('TYPE = l', 'TYPE = s')? to make default, which already is
#        if iline == 43 : line=line.replace(' z0=10', ' z0=20')? only for alphas and ST40 when H beams
#        if iline == 44 : line=line.replace(' z1=11', ' z1=21')? only for alphas -//-
    if iline == 33 : line=' 0.800E+05 0.800E+05 0.352E+07\n' ? look for 0.800E+05 if ihsps=1 or 2
    line=line.replace('s', 'l')  ? m.b. in the default set up
    
    ./runStabNovak.py

16. cp /u/ngorelen/work/exe/runGrowthNovak.py .
    Nmin = 3  ?
    Nmax = 30 ?
    fls=glob.glob("N"+str(i)+"/outa*") ? for alphas, change it to d for beam ions
    list1[3] = "i" ? use i for thermal ion Landau damping
#    list1[32] = "1"

17. Creating a file like Xin04w.2140E+01 for ORBIT to run is possible
    only for Boozer equilibrium coordinates. After you find the mode that
    file will be created. And after that and NOVA-C run this file will be
    updated with the additional information.
