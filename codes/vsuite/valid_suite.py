__author__ = "joanne cohn"
__email__  = "jcohn@berkeley.edu"
__version__= "1.1"  #updated BWC M*(Mh) from newer version of paper


import numpy as N
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab

def comments():
    """
Please do email me if you have questions!
Appendix of http://arxiv.org/abs/1609.03956 has information on how to run this program if more documentation is needed.
    
generate 7 plots using runsuite(): below.
4 are stellar mass functions:
all, quiescent, star forming, and all on one page, compared to several observations described below
(quiescent/star forming division at log sfr = -0.49 + (0.65+slopeval) (logM* - 10) +1.07 *(z-0.1)+shiftval, for slopeval=shiftval=0,
Moustakas et al eq 2), although many of papers listed use UVJ.

1 is stellar mass-sfr diagram [can be compared with e.g., Moustakas et al 2013, but not overplotted with it]

1 is ssfr in 4 stellar mass bins* (no cut on ra, dec for this)

1 is stellar mass to halo mass diagram for central galaxies, compared to Behroozi, Wechsler, Conroy 2013 and Moster,Naab, White 2013 fits
Behroozi,Wechsler,Conroy 2013 use Mvir
Moster, Naab,White 2013 use M200

If you use this program, please reference the papers and people who measured
all of these data!!
They are listed below "%%%"


USAGE:
runsuite(zcen, "inputfile.dat",hval,omm,slopeval,shiftval, boxside,runname,delz,ramin,ramax,decmin,decmax):
zcen is central redshift
fname = "inputfile.dat" described below, can call it something else
 if you want.  ascii text.
hval  = hubble constant
omm = omega_matter (e.g. 0.31)
slopeval = in sfr-M* bimodal diagram, **change in** slope of line to
separate star-forming and quiescent from PRIMUS
shiftval = change in shift of line between star forming and quiescent
 from PRIMUS

PRIMUS starforming and quiescent split by: 
log SFR = log sfrmin -0.49 + (0.65+slopeval) (logM* - 10) +1.07 *(z-0.1) + shiftval 

boxside = in Mpc/h for fixed time, any negative number if light cone
runname = string, such as "run0"
if lightcone, delz,ramin,ramax,decmin,decmax listed next.
if fixed time these arguments (delz, ramin,ramax,decmin,decmax)
are ignored and are not needed.

 files needed:
 from your simulation: requires "inputfile.dat" in the form of
 log10 m* (0) sfr (1), ra (2), dec (3), zred(4), ifsat (5) log10 m_halo (6)
 units:
log10 M* [M_o]
sfr units are per yr (not gyr)
 ra, dec the usual
 zred = redshift
ifsat = 0 for central, 1 for sat
 m_halo = halo mass (Mvir, [M_o])

Comparisons are made with data files, listed below and in this directory: 
--note that aside from (1),(5), (6), these were copied from tables and plots,
please let me know if you find errors!  thank you.

1. moustakas_z%s.smf,
provided at www.peterbehroozi.com/data.html,observational-data.tar.gz

2. fig3_bwc_12.dat
points from fig 3, Behroozi, Peter S., Wechsler, Risa H., Conroy, Charlie
    The Average Star Formation Histories of Galaxies in Dark Matter Halos from z=0-8


3. smf_all_supergrid01.txt, smf_starforming_supergrid01.txt,
smf_quiescent_supergrid01.txt
 Moustakas et al, 2013, table 4.

4. smf_zfourge_%s_supergrid.txt, %s= all, starforming, quiescent,
 Tomczak et al 2014, table 1
 

5. Vmax_NN*.dat, MaxLik**.dat, 
provided at
http://cosmos.phy.tufts.edu/~danilo/MuzzinEtal2013/Muzzin_et_al._(2013).html
paper is

6.  henriques_all.dat, henriques_quiescent.dat,henriques_starforming.dat:
 points from figs 2 and 7 of Henriques et al
http://arxiv.org/abs/1410.0365, data provided at
 http://galformod.mpa-garching.mpg.de/public/LGalaxies/figures_and_data.php

7. viper_sch_*.dat from Moutard et al,    Moutard et al, 2016
    The VIPERS Multi-Lambda Survey. II
Diving with massive galaxies in 22 square degrees since z = 1.5
    arxiv: 1602.05917 v3
    table 2.
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Full references for papers:
    Behroozi,Wechsler,Conroy
    from papers http://arxiv.org/abs/1207.6105 , http://arxiv.org/abs/1209.3013
        The Average Star Formation Histories of Galaxies in Dark Matter Halos from z = 0-8
     2013,ApJ,770,57 and
     On the Lack of Evolution in Galaxy Star Formation Efficiency,
     2013 ApJ, 762, L31

       Henriques, Bruno M. B.; White, Simon D. M.; Thomas, Peter A.; Angulo, Raul; Guo, Qi; Lemson, Gerard; Springel, Volker; Overzier, Roderik
Galaxy formation in the Planck cosmology - I. Matching the observed evolution of star formation rates, colours and stellar masses
 2015
http://arxiv.org/abs/1410.0365
 MNRAS,451,2663 
data tables: http://galformod.mpa-garching.mpg.de/public/LGalaxies/figures_and_data.php

    Moster, Naab & White, 2013
Galactic star formation and accretion histories from matching galaxies to dark matter haloes
http://arxiv.org/abs/1205.5807
 MNRAS, 428, 3121
      
    Moustakas, John, et al,
    PRIMUS: Constraints on Star Formation Quenching and Galaxy Merging, and the Evolution of the Stellar Mass Function from z = 0-1
    http://arxiv.org/abs/1301.1688
    ApJ, 2013, 767, 50
     
   Moutard et al, 2016
    The VIPERS Multi-Lambda Survey. II
Diving with massive galaxies in 22 square degrees since z = 1.5
    arxiv: 1602.05917 v3
    
Adam Muzzin, Danilo Marchesini, Mauro Stefanon, Marijn Franx, Henry J. McCracken, Bo Milvang-Jensen, James S. Dunlop, J. P. U. Fynbo, Gabriel Brammer, Ivo Labbe, Pieter van Dokkum, 2013,
          The Evolution of the Stellar Mass Functions of Star-Forming and Quiescent Galaxies to z = 4 from the COSMOS/UltraVISTA Survey
http://arxiv.org/abs/1303.4409
2013, ApJ, 777, 18

    Tomczak et al, 2014,
    Galaxy Stellar Mass Functions from ZFOURGE/CANDELS: An Excess of Low-mass Galaxies since z = 2 and the Rapid Buildup of Quiescent Galaxies
    arXiv:1309.5972
	2014,ApJ,783,85

    A companion reference to Tomczak et al is the description of how the data/catalogs were put
    together for the survey, to appear in Straatman et al, submitted.

%%%%%%%%%%%%%%%%%
Many thanks to P. Behroozi, M. van Daalen, A. Gonzalez, L. Guzzo, B. Henriques, J. Moustakas, M. White for help in putting the data and
plots together.  
    
    
%%%%%%%%%%%%    
model choices: BC03, Maraston (->BC03 0.14dexM*),Pegase, FSPS for sps model
type of sf: ssp,burst (exponential), constant
type of dust: Blanton-Roweis, Charlot-Fall, Calzetti
 A comprehensive analysis of uncertainties affecting the stellar mass-halo mass
 relation for 0<z<4,
http://arxiv.org/abs/1001.0015
imf: salpeter,kroupa,chabrier
To rescale stellar masses from Chabrier or Kroupa to Salpeter IMF,
we divide by constant factors 0.61 and 0.66, respectively. Madau and Dickinson, 2014,Cosmic Star Formation History, http://arxiv.org/abs/1403.0007, page 14
    """



def chiofz(zval=0.45,omm=0.31):
  """
  comoving distance to redshift zval
  in Mpc/h
  omega matter = omm
  use for volume if specify region with ra/dec
  """
  Nint = 300000
  zp1int = N.linspace(1,zval+1,Nint)
  ez = N.sqrt(omm*zp1int*zp1int*zp1int + (1-omm))
  tmp = 2997.925* N.trapz(1/ez, dx = zp1int[1]-zp1int[0])
  return(tmp)


def getsimstellar(zcen=0.45,addcolor=0,fname="galshort.dat",hval=0.67,omm=0.31,slopeval=0,shiftval=0,boxside=100,delz=0.02,ramin=-2,ramax=-2,decmin=2,decmax=-2,scatterval=0):
   """
   usage:
   getsimstellar(zval,addcolor,inputfil,hval,omm,slopeval,shiftval,boxside,delz,ra_min,ra_max,dec_min,dec_max)

   zcen: central redshift for simulation

   addcolor=0  all galaxies
   addcolor=1  red
   addcolor=2  blue
   fname:  galaxy data file, more below
   hval   = hubble constant (e.g. 0.67)
   shiftval,slopeval =changes from PRIMUS SFR-M* active quiescent classification, set to 0 if want to be as simple as possible
   omm    = omega_matter (e.g. 0.31)   
   boxside: positive for periodic box, negative for light cone
   boxside = side of box when periodic box
   **or**
   boxside<0 (e.g. -1) if light cone
   delz = use z range zcen+delz > z < zcen-delz for light cone [ignored for per]   
   ra_min,ra_max,dec_min,dec_max :min/max ra and dec for light cone [ignored for per]

   this is n(M) not N(>M)
   gal units M_o
   smf_*supergrid*01,03,10 used.
   color log sfrmin -0.49 + (0.65+slopeval) (logM* - 10) +1.07 *(z-0.1) + shiftval

   galaxy data file fname entries on each line, one per galaxy
   example:
#a = 0.9947 
# M*(0) [M_o], sfr(1) [M_o/yr],ra(2),dec(3),zred(4),ifsat(5),logmh(6)[M_o] 
   1.146e+01   4.831e-01    1.   1.  0.01  0  14.4675 
   9.696e+00   7.124e-03    1.   1.  0.01  1  11.2347 
   1.142e+01   1.355e-01    1.   1.  0.01  0  14.4215 
   8.386e+00   2.415e-03    1.   1.  0.01  1  9.5894 
   etc...

   [boxside > 0, i.e. fixed time, give boxside in units of Mpc/h]
   log10 mstellar (no h), sfr (M_o/Gyr)

   [boxside < 0,use any boxside value < 0, lightcone]
   log10 mstellar (no h), sfr (M_o/Gyr), ra, dec, redshift 
   """
   ff = open(fname)
   gals = N.loadtxt(ff)
   ff.close()
   
   logstell = gals[:,0]
   sfr      = gals[:,1]
   if (boxside < 0):
       print "using light cone"
       ra = gals[:,2]
       dec = gals[:,3]
       redz = gals[:,4]
       #need ra, dec, redshift,delz
       chimax = chiofz(zcen+delz,omm)  #[Mpc/h]
       chimin = chiofz(zcen-delz,omm)  #[Mpc/h]
       print "ramin,ramax, decmin,decmax %5.4f    %5.4f   %5.4f  %5.4f \n"%(ramin,ramax,decmin,decmax)
       angvol = -(N.cos((90-decmin)*N.pi/180) - N.cos((90-decmax)*N.pi/180))*(N.pi*(ramax-ramin)/180.)
       chivol =(chimax*chimax*chimax - chimin*chimin*chimin)/3.
       vol = chivol*angvol  # in [Mpc/h]^3
       # truncate galaxy sample to light cone
       jj = N.nonzero((ra>ramin)&(ra<ramax)&(dec>decmin)&(dec<decmax)&
        (redz<zcen+delz)&(redz>zcen-delz))[0]
       sfr = sfr[jj]
       logstell = logstell[jj]
       redz     = redz[jj]
   if (boxside>0):
       print "using periodic box, side %8.2f Mpc/h"%(boxside)
       vol = boxside*boxside*boxside
       redz = zcen
   #units:
   # want mpc not mpc/h
   vol = vol/(hval*hval*hval)

   ## add random scatter as function of z
   if (scatterval==1):
     sigval = 0.07 +0.04*redz
     logstell += N.random.normal(0,sigval,logstell.size)

   
   jj = N.arange(logstell.size)
   #note color cut is assuming h70's in units
   if (addcolor>0):
       #moustakas 2013 units to compare
       #Moustakas box is in units of [Mpc/h70]^3, we have [Mpc/h]^3,
   #   so divide volume by (h/h70)^3 = (h/(h/0.7))^3 = 0.7^3
   # Mstar is in units of [M_o/h70^2]

       sfrtest = sfr**(hval/0.70)**2
       logstelltest = logstell + 2.*N.log10(hval/0.70) 
   if (addcolor==1):
     jj = N.nonzero(N.log10(sfrtest+1.e-16)<-0.49+(0.65+slopeval)*(logstelltest-10)+1.07*(redz-0.1)+shiftval)[0]
   if (addcolor==2):  
        jj = N.nonzero(N.log10(sfrtest+1.e-16)>=-0.49+(0.65+slopeval)*(logstelltest-10)+1.07*(redz-0.1)+shiftval)[0]
   logstell = logstell[jj]     

   nbin = 50
   nhist,bins = N.histogram(logstell,nbin,range=(8.3,12.3))
   bins += (bins[1]-bins[0])/2.
   bins = N.delete(bins,nbin)
   ngalact = nhist*1./(vol*(bins[1]-bins[0]))
   
   galkind =("all","red","blue")
  
  
   return(bins,ngalact,nhist.sum())
###
###  different models
###
   
def getphim_bc03(zcen=0.25,addcolor=0):
    """
    only for z>0.2
    Behroozi,Wechsler,Conroy
    from papers
    http://arxiv.org/abs/1207.6105 , http://arxiv.org/abs/1209.3013
    The Average Star Formation Histories of Galaxies in Dark Matter Halos from z = 0-8, 2013,ApJ,770,57
     On the Lack of Evolution in Galaxy Star Formation Efficiency, 2013 ApJ, 762, L31

    data from publicly available  
    behroozi-2013-data-compilation at www.peterbehroozi.com
    Stellar Mass functions: smf_ms/moustakas*.smf
Columns: Log10(stellar mass) (Msun), Log10(ND) (1/Mpc^3/dex), Err+ (dex), Err- (dex)
*OR*
Columns: Log10(stellar mass) (Msun), ND (1/Mpc^3/dex), Err+ , Err-
In the latter case, the data files are marked with "#errors: linear".

Assumptions:BC03 SPS models, Chabrier (2003) IMF, Blanton & Roweis (kcorrect) dust modeling.


    """
    if ((addcolor !=0)|(zcen>1)):
      return(N.array([1,1]),N.array([1,1]),0.,0.,0.,0.)
    znamelist = ("0.105","0.25","0.35","0.45","0.575","0.725","0.9")
    zvals = N.array([0.01,0.2,0.3,0.4,0.5,0.65,0.8,1.0])
    jj= N.nonzero(zcen>=zvals)[0]
    if (jj.size==0):
       return(N.array([1,1]),N.array([1,1]),0.,0.,0.,0.)
    if (jj.size==1):
        zmin =0.01
        zmax = 0.2
    if (jj.size>1):
        jj = jj.max()
        zmin =zvals[jj]
        zmax = zvals[jj+1]
    print "behroozi compilation"    
    ff =open("moustakas_z%s.smf"%(znamelist[jj]))
    phivals = N.loadtxt(ff)
    ff.close()
    # logm phi errplus errmin      
    ctypelist="all"
    #log phi errors are both positive
    logm = phivals[:,0]
    phi  = phivals[:,1]
    phip = phivals[:,2]
    phim = phivals[:,3]
    return(logm,phi,phip,phim,zmin,zmax)  


   
   
def getphibwc(zcen):
    """
    for all galaxies only, points in fig 3 of
    Behroozi, Peter S., Wechsler, Risa H., Conroy, Charlie
    The Average Star Formation Histories of Galaxies in Dark Matter Halos from z = 0-8
     2013,ApJ,770,57
     arXiv:1207.6105
    bc03, Blanton-Roweis dust, chabrier imf
    """
    zmid = -1. # flag
    phis = N.loadtxt("fig3_bwc12.dat")
    if (zcen<0.2):
        phis = phis[N.arange(66),:]
        zmid = 0.
    if ((zcen>0.3)&(zcen<0.75)):
        phis = phis[N.arange(66,88),:]
        zmid = 0.5
    if ((zcen>0.75)&(zcen<1.25)):    
        phis = phis[N.arange(88,103),:]
        zmid = 1.
    if (zmid>=0):    
      return(N.log10(phis[:,0]),phis[:,1],zmid)
    if (zmid<0):    
      return(N.array([1,1]),N.array([1,1]),zmid)
         
def getphisg(zcen=0.1,addcolor=0):
    """
    Moustakas, John, et al,
    PRIMUS: Constraints on Star Formation Quenching and Galaxy Merging, and the Evolution of the Stellar Mass Function from z = 0-1
         http://arxiv.org/abs/1301.1688
    table 3 h70 units
    fsps (Conroy/White/Gunn, Conroy/Gunn/White, Conroy/Gunn 2009,2010)
    Charlot Fall(2000) dust
    chabrier imf
    """
    if (zcen>0.2):
       return(N.array([1,1]), N.array([1,1]),0,0,0,0)
    ctypelist=("all","quiescent","starforming")
    logm = 9.0+ N.arange(31)*1./10.    
    if (addcolor==0): #only have fsps
        logphi = N.array([-1.899,-1.923,-1.970,-2.031,-2.055,-2.106,-2.144,-2.179,-2.188,-2.216,-2.234,-2.235,-2.262,-2.252,-2.285,-2.317,-2.365,-2.419,-2.504,-2.607,-2.728,-2.888,-3.104,-3.332,-3.606,-3.953,-4.363,-4.778,-5.255,-5.87,-6.49])
        logphi_plus = N.array([0.017,0.017,0.015,0.015,0.014,0.012,0.012,0.012,0.010,0.0086,0.0080,0.0069,0.0063,0.0056,0.0051,0.0047,0.0044,0.0041,0.0040,0.0039,0.0040,0.0043,0.0049,0.0059,0.0080,0.012,0.020,0.033,0.060,0.010,0.030])
        logphi_minus = N.array([-0.017,-0.016,-0.015,-0.014,-0.013,-0.012,-0.011,-0.012,-0.010,-0.0084,-0.0078,-0.0068,-0.0062,-0.0056,-0.0051,-0.0046,-0.0044,-0.0041,-0.0040,-0.0039,-0.0040,-0.0043,-0.0048,-0.0059,-0.0079,-0.012,-0.019,-0.031,-0.053,-0.010,-0.020])
    if (addcolor==1): #only have fsps
        logphi = N.array([-2.495,-2.486,-2.485,-2.523,-2.576,-2.603,-2.634,-2.642,-2.652,-2.655,-2.649,-2.614,-2.607,-2.5640,-2.5640,-2.5800,-2.6050,-2.6450,-2.7050,-2.7860,-2.8840,-3.0190,-3.2090,-3.4130,-3.6670,-4.002,-4.401,-4.806,-5.296,-5.93,-6.16])
        logphi_plus = N.array([0.048,0.044,0.038,0.037,0.033,0.030,0.026,0.028,0.021,0.018,0.015,0.013,0.011,0.0089,0.0077,0.0069,0.0062,0.0057,0.0053,0.0050,0.0049,0.0050,0.0055,0.0065,0.0085,0.013,0.021,0.034,0.063,0.10,0.40])
        logphi_minus = N.array([-0.043,-0.041,-0.035,-0.034,-0.031,-0.028,-0.025,-0.026,-0.020,-0.017,-0.015,-0.012,-0.011,-0.0087,-0.0076,-0.0068,-0.0061,-0.0056,-0.0052,-0.0050,-0.0049,-0.0050,-0.0054,-0.0064,-0.0084,-0.012,-0.020,-0.032,-0.056,-0.10,-0.20])
    if (addcolor==2): #only have fsps
        logphi = N.array([-2.026,-2.062,-2.129,-2.201,-2.211,-2.272,-2.313,-2.362,-2.371,-2.4120,-2.4450,-2.4700,-2.5240,-2.5410,-2.6090,-2.6600,-2.7370,-2.8110,-2.9340,-3.0770,-3.2500,-3.4720,-3.769,-4.102,-4.487,-4.930,-5.437,-5.98,-6.30,-6.77,-7.09])
        logphi_plus = N.array([0.018,0.017,0.015,0.014,0.014,0.012,0.012,0.011,0.011,0.0092,0.0090,0.0079,0.0074,0.0071,0.0066,0.0063,0.0062,0.0059,0.0061,0.0064,0.0071,0.0085,0.011,0.016,0.024,0.042,0.079,0.20,0.30,0.60,1.00])
        logphi_minus = N.array([-0.017,-0.016,-0.015,-0.014,-0.013,-0.012,-0.012,-0.011,-0.011,-0.0090,-0.0088,-0.0078,-0.0072,-0.0070,-0.0065,-0.0062,-0.0061,-0.0059,-0.0060,-0.0063,-0.0070,-0.0084,-0.010,-0.015,-0.023,-0.038,-0.067,-0.10,-0.20,-0.30,-0.40])
        
    phi = 10**logphi
    phip = phi*(10**logphi_plus -1)
    phim =phi*(1-10**logphi_minus)
    return(logm,phi,phip,phim,0.01,0.2)

        
    

def getphim(zcen=0.25,addcolor=0,ismf=0):
    """
    Moustakas, John, et al,
    PRIMUS: Constraints on Star Formation Quenching and Galaxy Merging, and the Evolution of the Stellar Mass Function from z = 0-1
     http://arxiv.org/abs/1301.1688
    table 4
    fsps (Conroy/White/Gunn, Conroy/Gunn/White, Conroy/Gunn 2009,2010)
    Charlot Fall(2000) dust
    chabrier imf

    only for z>0.2
    """
    if ((zcen<0.2)|(zcen>1.0)):
        return(N.array([1,1]),N.array([1,1]),0,0,0,0)
    ctypelist=("all","quiescent","starforming")
    ff = open("smf_%s_supergrid01.txt"%(ctypelist[addcolor]))
    #zlow 0 zhi 1 ngal 2 logm* 3 limit4 logphi5 logphierrm6 logphierrp7 logphierrcv8
    #log phi errors are both positive
    gals = N.loadtxt(ff,usecols=(0,1,3,5,6,7))
    ff.close()
    jj = N.nonzero((gals[:,0]<zcen)&(gals[:,1]>=zcen))[0]
    gals = gals[jj]
    zmin = gals[0,0]
    zmax = gals[0,1]
    logm    = gals[:,2]
    phi  = 10**gals[:,3]
    phim = phi*(1-10**(-gals[:,4]))
    phip = phi*(10**gals[:,5]-1)
    return(logm,phi,phip,phim,zmin,zmax)  

    
def getphit(zcen=0.45,addcolor=0):
   """
    Tomczak et al, 2014,
    Galaxy Stellar Mass Functions from ZFOURGE/CANDELS: An Excess of Low-mass Galaxies since z = 2 and the Rapid Buildup of Quiescent Galaxies
    arXiv:1309.5972
	2014,ApJ,783,85
    Another reference with with Tomczak et al is the paper detailing
    the way that the data/catalogs were put
    together for the survey, to appear in Straatman et al submitted.
    
   read data of tomczak et al 2014 table 1
   surrounding region of zcen chosen

   units are: log M, log Phi
   M* is in units of M_o/h70^2
       -if your M* is in units of M_o, multiply your M* by h70^2
   Phi is in units of [h70/Mpc]^3, 
       -if your Phi is in units of [Mpc^-3], divide by h70^3
   stellar masses using FAST (Kriek et al 2009)
   Bruzual & Charlot (2003) following an exponentially declining starformation history assuming a Chabrier (2003) initial mass function. They assume solar metallicity and allow Av to vary between [0, 4].    
   """
   # now need to find right redshift and type
   # first redshift
   zrange = N.array([0.2,0.5,0.75,1.0,1.25, 1.5,2.0,2.5,3.0])
   jjz = N.nonzero(zcen>=zrange)[0]
   if ((jjz.size==0)|(zcen>3)):
      return(N.array([1,1]),N.array([1,1]),0.,0.,0.,0.)
   if (jjz.size>1):
       jjz = jjz.max()
   zmin = zrange[jjz]
   zmax = zrange[jjz+1]
   print "using ZFOURGE range %3.2f < z < %3.2f "%(zrange[jjz],zrange[jjz+1])
   print "Bruzual Charlot used to calculate M*, solar Z, Av in [0,4] "
   
   colornamelist =("all","quiescent","starforming")
   ff = open("smf_zfourge_%s_supergrid.txt"%(colornamelist[addcolor]))
   phitom = N.loadtxt(ff,usecols=(0,1,3,5,6,7))
   ff.close()
   #zlo0 zhi 1 logm2 logphi3 logphim4 logphip5 
   jj = N.nonzero((zcen> phitom[:,0])&(zcen<=phitom[:,1]))[0]
   phitom = phitom[jj,:] # now have right redshift and right color sel
   logm = phitom[:,2]
   logphi = phitom[:,3]
   logphim = phitom[:,4]
   logphip = phitom[:,5]
   
   phi   = 10**logphi
   phip  = phi*(10**logphip-1)
   phim  = phi*(1-10**(-logphim))
   return(logm,phi,phip,phim,zmin,zmax)

def getphiuv(zcen=0.25,addcolor=0):
    """
   
   The Evolution of the Stellar Mass Functions of Star-Forming and Quiescent Galaxies to z = 4 from the COSMOS/UltraVISTA Survey
Adam Muzzin, Danilo Marchesini, Mauro Stefanon, Marijn Franx, Henry J. McCracken, Bo Milvang-Jensen, James S. Dunlop, J. P. U. Fynbo, Gabriel Brammer, Ivo Labbe, Pieter van Dokkum, 2013,
http://arxiv.org/abs/1303.4409
downloads at    cosmos2.phy.tufts.edu/~danilo/Downloads.html
 bc03, calzetti dust, kroupa imf
    """
    if (zcen<0.2):
        return(N.array([1,1]),N.array([1,1]),0,0,0,0)
    if (zcen>4.):
        return(N.array([1,1]),N.array([1,1]),0,0,0,0)
    
    ctypelist=("all","quiescent","starforming")
    zlist = N.array([0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0])
    jj = N.nonzero(zcen>zlist)[0]
    if (jj.size>1):
        jj = jj.max() #largest value of z less than zcen
    zmin = zlist[jj]
    zmax = zlist[jj+1]
    print "using COSMOS/Ultravista range %3.2f < z < %3.2f "%(zmin,zmax)

    
    ff = open("Vmax_%2.1fz%2.1f.dat"%(zmin,zmax))
    #logMs EMstar(1) logphi(2), eu(phi3) el(phi4) logphiq(5) ueq(6) leq(7) phisf(8) uesf(9) lesf(10)
    #log phi errors are both positive
    
    gals = N.loadtxt(ff,usecols=(0,2+3*addcolor,3+3*addcolor,4+3*addcolor))
    ff.close()
    jj = N.nonzero(gals[:,1]>-99)[0] #only where have measurements.
    gals = gals[jj,:]
    logm    = gals[:,0]
    # shift from kroupa to chabrier using 0.61/0.66, Madau and Dickinson, page 14
    logm += N.log10(0.61/0.66) 
    phi  = 10**gals[:,1]
    phim = phi*(1-10**(-gals[:,3])) 
    phip = phi*(10**gals[:,2]-1)
    return(logm,phi,phip,phim,zmin,zmax)  
   
      
def getphiuv_sch(zcen=0.25,addcolor=0):
    """
   THE EVOLUTION OF THE STELLAR MASS FUNCTIONS OF STAR-FORMING AND QUIESCENT GALAXIES TO z = 4 FROM THE COSMOS/UltraVISTA SURVEY,

Adam Muzzin, Danilo Marchesini, Mauro Stefanon, Marijn Franx, Henry J. McCracken, Bo Milvang-Jensen, James S. Dunlop, J. P. U. Fynbo, Gabriel Brammer, Ivo Labbe, PG van Dokkum, 2013, ApJ, 777, 1
http://arxiv.org/abs/1303.4409
schechter function fit 
downloads at   cosmos2.phy.tufts.edu/~danilo/Downloads.html

bc03,calzetti dust,kroupa imf.
    """
    if (zcen<0.2):
        return(N.array([1,1]),N.array([1,1]),0,0)
    if (zcen>4.):
        return(N.array([1,1]),N.array([1,1]),0,0)
    
    galkind = ("ALL","QUIESCENT","STARFORMING")
    ff = open("MaxLik_Schechter_%s.dat"%(galkind[addcolor]))
    sparams = N.loadtxt(ff)
    ff.close()
    #zlow,zhigh,nobj(2),mlim(3),m*(4),m*1su(5),m*1sl(6),m*1sutot(7),m*1sltot(8),phi*(9),phi*_1su(10)
    #phi*_1sl(11),phi*1_sutot(12),phi*1_sltot(13),alpha(14),alpha_1su(15),alpha_1sl(16),alpha1sutot(17),
    #alpha1sltot(18),
    #phi2*(19),phi2*_su(20)
    #phi2*_1sl(21),phi2*1_sutot(22),phi2*1_sltot(23),alpha2(24),alpha2_1su(25),alpha2_1sl(26),alpha2sutot(27),
    #alpha2sltot(28)
    doublefit=0
    jj = N.nonzero((zcen>sparams[:,0])&(zcen<=sparams[:,1])&(sparams[:,15]!=0))[0]#redshift and floating alpha
    sparams = sparams[jj,:] #now at right redshift and floating alpha
    if (sparams[:,0].size>1): #double schechter fit
        doublefit=1
        jj = N.nonzero(sparams[:,20]>-99)[0] #get double schechter fit
        sparams=sparams[jj,:]
    #now just one row
    zmin = sparams[0,0]
    zmax = sparams[0,1]
    logm = N.linspace(sparams[:,3],12,100) #log M, units M_o
    mstar     = sparams[0,4]
    phistar   = sparams[0,9]*1.e-4
    alpha     = sparams[0,14]
    phi       = N.log(10)*phistar*N.power(10.,(logm-mstar)*(1+alpha))*N.exp(-N.power(10,logm-mstar))
    if (doublefit==1): #double schechter fit
        phistar2 = sparams[0,19]*1.e-4
        alpha2   = sparams[0,24]
        phi     += N.log(10)*phistar2*N.power(10.,(logm-mstar)*(1+alpha2))*N.exp(-N.power(10,logm-mstar))
    logm += N.log10(0.61/0.66) #convert kroupa to chabrier
    return(logm,phi,zmin,zmax)
    
   
def getphihen(zcen=0.25,addcolor=0):
    """ 
  Henriques, Bruno M. B.; White, Simon D. M.; Thomas, Peter A.; Angulo, Raul; Guo, Qi; Lemson, Gerard; Springel, Volker; Overzier, Roderik
Galaxy formation in the Planck cosmology - I. Matching the observed evolution of star formation rates, colours and stellar masses
 2015, MNRAS,451,2663 
1410.0365: compilation of several measurements, details in appendix A2.
sigma8 = 0.829, H0 = 67.3 km/s/mpc, OmegaL = 0.685, Omegam = 0.315, Omegab = 0.0487 (fb = 0.155) and n = 0.96.
downloadable tables:
http://galformod.mpa-garching.mpg.de/public/LGalaxies/figures_and_data.php
points from figures 2 (all), figure 7 (quiescent and starforming)

shift from maraston to BC03 by adding 0.14 to logM*
    """
    #they have hval = 0.673 but have taken it out, units (mpc/h)^-3 and M*/h^2
    ctypelist=("all","quiescent","starforming")
    ff = open("henriques_%s.dat"%(ctypelist[addcolor]))
    gals = N.loadtxt(ff)
    ff.close()
    #vol is Mpc/h^3 so need to multiply by h^{-3}
    #M* is 
    startpoints = N.zeros(6*3,dtype='int')
    startpoints.shape=(6,3)
    startpoints[:,0]=N.array([0,17,17,30,41,49]) #all
    startpoints[:,1]=N.array([0,13,26,38,47,49]) #red
    startpoints[:,2]=N.array([0,13,26,38,48,55]) #blue

    zlist = N.array([0.1,0.4,1.0,2.0,3.0])
    jj = N.nonzero(abs(zcen-zlist)<0.2)[0]
    if (jj.size==0):
        return(N.array([1,1]),N.array([1,1]),0,0,0)
    if (jj.size>1):
        distval = abs(zcen-zlist)
        jj = N.nonzero(distval <=distval.min())[0]
        if (jj.size>1):
            jj = jj[0]    
    zmidh=zlist[jj]    
    jstart = startpoints[jj,addcolor]
    jend = startpoints[jj+1,addcolor]
    if (jend<=jstart):
        return(N.array([1,1]),N.array([1,1]),0,0,0) #just junk :), e.g.all for z=3. 
    mass = (10**gals[N.arange(jstart,jend),0]+10**gals[N.arange(jstart,jend),1])/2.
    phi  = gals[N.arange(jstart,jend),2]
    phip = phim = gals[N.arange(jstart,jend),3]
     #log will be later
    logm = N.log10(mass) 
    print "using Henriques z = %3.2f "%(zlist[jj])
    print "Maraston used to calculate M*, convert +0.14 dex"
    # shift from maraston to BC using 0.14, Henriques++1410.0365
    logm += 0.14 
    return(logm,phi,phip,phim,zmidh)  

def getphiv(zcen=0.45,addcolor=0):
    """
    Moutard et al, 2016
    The VIPERS Multi-Lambda Survey. II
Diving with massive galaxies in 22 square degrees since z = 1.5
    arxiv: 1602.05917 v3
    Use Schechter functions from table II (double)
    also have single schechter functions in table I 
    
   stellar masses using LePhare, metallicities 0.008 and 0.02,
   Bruzual & Charlot (2003),
   assuming a Chabrier (2003) initial mass function.
    exponentially declining star formation history, 9 possible decay rates between 0.1 and 30 Gyr.  See paper sec. 4.1 for dust
    prescription-three considered including Chabrier (2000).
    use equation 6
    phi(M*) dM* = exp(-M*/Mref) (phi1* (M*/Mref)^alpha1 + phi2*(M*/Mref)^alpha2) dM*/Mref
    
    """
    # now need to find right redshift and type
    # first redshift
    zrange = N.array([0.2,0.5,0.8,1.1,1.5])
    if ((zcen<zrange.min())|(zcen>zrange.max())):
      return(N.array([1,1]),N.array([1,1]),0.,0.)
    galkind =("all","quiescent","starforming")
    ff = open("viper_sche_%s.dat"%(galkind[addcolor]))
    sparams = N.loadtxt(ff)
    ff.close()
    #zlow,zhigh,nobj(2),mlim(3),m*(4),m*1su(5),m*1sl(6),phi*(7),phi*_1su(8)
    #phi*_1sl(9),alpha(10),alpha_1su(11),alpha_1sl(12),
    #phi2*(13),phi2*_su(14)
    #phi2*_1sl(15),
    doublefit=0
    jj = N.nonzero((zcen>sparams[:,0])&(zcen<=sparams[:,1]))[0]#redshift and floating alpha
    sparams = sparams[jj,:].flatten() #now at right redshift
    if (sparams[13] != -99): #double schechter fit
        doublefit=1
    #now just one row
    zmin = sparams[0]
    zmax = sparams[1]
    print "using VIPERS range %3.2f < z < %3.2f "%(zmin,zmax)
    print "Bruzual Charlot used to calculate M*, Chabrier IMF "
    
    logm = N.linspace(sparams[3],12,100) #log M, units M_o
    mstar     = sparams[4]
    phistar   = N.power(10,sparams[7])
    alpha     = sparams[10]
    phi = phistar*N.power(10.,(logm-mstar)*(alpha+1))*N.exp(-N.power(10,logm-mstar))
    print "doublefit=",doublefit
    if (doublefit==1):
       phistar2 = N.power(10,sparams[13])
       alpha2   = sparams[16]
       phi     += phistar2*N.power(10.,(logm-mstar)*(alpha2+1))*N.exp(-N.power(10,logm-mstar))
    phi = phi*N.log(10)
    return(logm,phi,zmin,zmax)
   

def plot4tog(zcen=0.45,fname="galshort.dat",hval=0.7,omm=0.31,slopeval=0.,shiftval=0.,boxside=-1,runname="runname",delz=0.02,ramin=16.98,ramax=20.17,decmin=13.23,decmax=16.33):   
    """
    three colors, four models, all together, with or without obs scatter
    """
    f,ax = plt.subplots(2,2,sharex=True, sharey=True)
    collist = ('k','r','b')
    galtype=("all","quiescent","starforming")
    smftype=("primus_bc03","bwc_comp","sdss_gal","primus_fsps","zfourge","cos/uv","hen15","vipers16")
    zminlist = N.zeros(len(smftype),dtype='float')
    zmaxlist = N.zeros(len(smftype),dtype='float')
    smfflag=N.zeros(len(smftype),dtype='int') #flag for what appears
    smarkerlist=('s','^','x','*','o','+','v')
    scollist=('c','y','g','m','darkgreen','thistle','pink','sandybrown')

    for i in range(3): #color
      #first get data  
      bin_centers,ngalact,ngal=getsimstellar(zcen,i,fname,hval,omm,slopeval,shiftval,boxside,delz,ramin,ramax,decmin,decmax,0)
      #now with scatter Behroozi/Wechsler/Conroy 2013
      ax[i%2,i/2].step(bin_centers, ngalact,collist[i],linestyle=':',label='simulation')
      bin_centers,ngalact,ngal=getsimstellar(zcen,i,fname,hval,omm,slopeval,shiftval,boxside,delz,ramin,ramax,decmin,decmax,1)      
      ax[i%2,i/2].step(bin_centers, ngalact,collist[i],label="scattered sim")
      #run through smf's
      hrat = hval/0.7 #h70, most people use h=0.7 in their analysis
      if (i==0):
          ismf = 0 
          logm,phi,phip,phim,bpmin,bpmax =getphim_bc03(zcen,i)
          if (logm[0]>1):
             phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
             phip *=(hrat*hrat*hrat)
             phim *=(hrat*hrat*hrat)
             logm -=  2.*N.log10(hval/0.70) 
 
             ax[i%2,i/2].errorbar(logm,phi,yerr=[phim,phip],fmt=' ',color=scollist[ismf],marker=smarkerlist[ismf],label="%s %3.2f<z<%3.2f"%(smftype[ismf],bpmin,bpmax))
             zminlist[ismf] = bpmin
             zmaxlist[ismf] = bpmax
             smfflag[ismf]=1
          ismf=1   
          logm,phi,zmid = getphibwc(zcen)
          if (logm[0]>1):
            phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
            logm -=  2.*N.log10(hval/0.70) 

            ax[i%2,i/2].plot(logm,phi,color=scollist[ismf],marker=smarkerlist[ismf],label="%s z=%3.2f"%(smftype[ismf],zmid))
            smfflag[ismf]=1
            zminlist[ismf] = zmid
      ismf=2      
      logm,phi,phip,phim,szmin,szmax =getphisg(zcen,i)
      if (logm[0]>1):
             phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
             phip *=(hrat*hrat*hrat)
             phim *=(hrat*hrat*hrat)
             logm -=  2.*N.log10(hval/0.70) 
             ax[i%2,i/2].errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="SDSS-FSPS %3.2f<z<%3.2f"%(szmin,szmax))
             smfflag[ismf]=1
             zminlist[ismf] = szmin
             zmaxlist[ismf] = szmax
      ismf=3       
      logm,phi,phip,phim,pzmin,pzmax=getphim(zcen,i)
      if (logm[0]>1):
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         phip *=(hrat*hrat*hrat)
         phim *=(hrat*hrat*hrat)
         logm -=  2.*N.log10(hval/0.70) 
          
         ax[i%2,i/2].errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="%s %3.2f<z<%3.2f"%(smftype[ismf],pzmin,pzmax))
         smfflag[ismf]=1
         zminlist[ismf] = pzmin
         zmaxlist[ismf] = pzmax
      ismf = 4
      #add zfourge
      logm,phi,phip,phim,zfmin,zfmax =getphit(zcen,i)
      if (logm[0]>1):
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         phip *=(hrat*hrat*hrat)
         phim *=(hrat*hrat*hrat)
         logm -=  2.*N.log10(hval/0.70) 
          
         ax[i%2,i/2].errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="%s %3.2f<z<%3.2f"%(smftype[ismf],zfmin,zfmax))
         smfflag[ismf]=1
         zminlist[ismf] = zfmin
         zmaxlist[ismf] = zfmax
      ismf = 5 #ultravista
      logm,phi,phip,phim,uzmin,uzmax =getphiuv(zcen,i)
      if (logm[0]>1):
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         phip *=(hrat*hrat*hrat)
         phim *=(hrat*hrat*hrat)
         logm -=  2.*N.log10(hval/0.70) 
  
         ax[i%2,i/2].errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="%s %3.2f<z<%3.2f "%(smftype[ismf],uzmin,uzmax))
         smfflag[ismf]=1
         zminlist[ismf] = uzmin
         zmaxlist[ismf] = uzmax
       #schechter function to ultravista
      logm,phi,uzmin,uzmax = getphiuv_sch(zcen,i)
      if (logm[0]>1):
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         phip *=(hrat*hrat*hrat)
         phim *=(hrat*hrat*hrat)
         logm -=  2.*N.log10(hval/0.70) 
      
         ax[i%2,i/2].plot(logm,phi,color=scollist[ismf],label="%s %3.2f<z<%3.2f "%(smftype[ismf],zfmin,zfmax))
         smfflag[ismf]=1   
      ismf = 6 #henriques #center of mass bin taken
      logm,phi,phip,phim,zhen = getphihen(zcen,i)
      if (logm[0]>1):
         phi    *= hval*hval*hval #units [h/Mpc]^3
         phip   *= hval*hval*hval
         phim   *= hval*hval*hval
         logm   -= 2*N.log10(hval) #units [M*/h^2]
         ax[i%2,i/2].errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="%s z=%3.2f "%(smftype[ismf],zhen))
         smfflag[ismf]=1
         zminlist[ismf] = zhen
      ismf = 7 #vipers 
      logm,phi,vzmin,vzmax = getphiv(zcen,i)
      if (logm[0]>1): #they assume h=0.70
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         logm -=  2.*N.log10(hval/0.70) 
      
         ax[i%2,i/2].plot(logm,phi,color=scollist[ismf],label="%s %3.2f<z<%3.2f "%(smftype[ismf],vzmin,vzmax))
         smfflag[ismf]=1
         zminlist[ismf]=vzmin
         zmaxlist[ismf] = vzmax 
         
      ax[i%2,i/2].set_xlim(8.2,12.)
      ax[i%2,i/2].set_ylim(1.e-5,0.02)
      ax[i%2,i/2].set_yscale("log")
      ax[i%2,i/2].text(8.5,1.e-4,'%s'%(galtype[i]))
      ax[i%2,i/2].text(8.5,5.e-5,r'$\bar{z}_{\rm sim}$=%3.2f'%(zcen))

      
    #trick it into putting legend in empty box  
    logm = N.array([6.0,6.01])
    phi = N.array([1.e-7,1.2e-7])
    phim = N.array([1.e-8,1.2e-8])
    phip = N.array([1.e-8,1.2e-8])    
    bin_centers = N.array([7,7.2])
    ngalact = N.array([1.e-8,1.e-8])
    
    ax[1,1].set_ylim(1.e-5,0.04)
    ax[1,1].set_yscale("log")
    
    ax[1,1].step(bin_centers, ngalact,'k',label="simulation")
    ax[1,1].step(bin_centers, ngalact,'k',linestyle=':',label="sim w/out obs scatter")
    ax[1,1].get_xaxis().set_visible(False)
    ax[1,1].get_yaxis().set_visible(False)    
    
    logm=N.array([5.,5.1])
    phi = 1.e-7*N.ones(2,dtype="float")
    phim = phip = 1.e-8
    i = 0
    for ismf in range(len(smftype)):
        if (smfflag[ismf]>0):
          if ((ismf==1)|(ismf==6)):
              ax[1,1].plot(logm,phi,marker=smarkerlist[ismf],color=scollist[ismf],linestyle='None',label="%s z=%3.2f"%(smftype[ismf],zminlist[ismf]))
          else:
            if (ismf !=7):       
              ax[1,1].errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="%s %3.2f<z<%3.2f"%(smftype[ismf],zminlist[ismf],zmaxlist[ismf]))
          if ((ismf==5)|(ismf==7)):
            ax[1,1].plot(logm,phi,color=scollist[ismf],label="%s %3.2f<z<%3.2f "%(smftype[ismf],zminlist[ismf],zmaxlist[ismf]))      

    
    ax[0,0].set_ylabel(" $\Phi$ [Mpc${}^{-3}$/dex]")      
    ax[1,1].legend(loc=3,fontsize='10',frameon=False)
    ax[1,1].axis('off')
    ax[1,0].set_xlabel(r'M* [log $M_\odot$]')
    f.subplots_adjust(hspace=0.001)
    f.subplots_adjust(wspace=0.001)
    plt.tight_layout()
    plt.savefig("smf4sims_%d_%s.pdf"%(zcen*100.5,runname))
      
    plt.close("all")  

def plot4sep(zcen=0.45,fname="galshort.dat",hval=0.7,omm=0.31,slopeval=0.,shiftval=0.,boxside=-1,runname="runname",delz=0.02,ramin=16.98,ramax=20.17,decmin=13.23,decmax=16.33):   
    """
    one color four models, all together, with or without obs scatter
    """
    hrat = hval/0.7
    collist = ('k','r','b')
    galtype=("all","quiescent","starforming")
    smftype=("primus_bc03","bwc_comp","sdss_gal","primus_fsps","zfourge","cos/uv","hen15","vipers16")
    smarkerlist=('s','^','x','*','o','+','v')
    scollist=('c','y','g','m','darkgreen','thistle','pink','sandybrown')
    for i in range(3): #color
      #set flags for each color separately  
      zminlist = N.zeros(len(smftype),dtype='float')
      zmaxlist = N.zeros(len(smftype),dtype='float')
      smfflag=N.zeros(len(smftype),dtype='int') #flag for what appears
      f,ax = plt.subplots(1,1)
      bin_centers,ngalact,ngal=getsimstellar(zcen,i,fname,hval,omm,slopeval,shiftval,boxside,delz,ramin,ramax,decmin,decmax,1)
      #with scatter
      ax.step(bin_centers, ngalact,collist[i],label="simulation")
      bin_centers,ngalact,ngal=getsimstellar(zcen,i,fname,hval,omm,slopeval,shiftval,boxside,delz,ramin,ramax,decmin,decmax,0)
      #no scatter
      ax.step(bin_centers, ngalact,collist[i],linestyle=':',label="sim, no obs scatter")
      #run through smf's
      if (i==0):
          ismf = 0 
          logm,phi,phip,phim,bpmin,bpmax =getphim_bc03(zcen,i)
          if (logm[0]>1):
             phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
             phip *=(hrat*hrat*hrat)
             phim *=(hrat*hrat*hrat)
             logm -=  2.*N.log10(hval/0.70) 
    
             ax.errorbar(logm,phi,yerr=[phim,phip],fmt=' ',color=scollist[ismf],marker=smarkerlist[ismf],label="%s %3.2f<z<%3.2f"%(smftype[ismf],bpmin,bpmax))
             zminlist[ismf] = bpmin
             zmaxlist[ismf] = bpmax
             smfflag[ismf]=1
          ismf=1   
          logm,phi,zmid = getphibwc(zcen)
          if (logm[0]>1):
            phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
            logm -=  2.*N.log10(hval/0.70) 
            ax.plot(logm,phi,color=scollist[ismf],marker=smarkerlist[ismf],linestyle='None',label="%s z=%3.2f"%(smftype[ismf],zmid))
            smfflag[ismf]=1
            zminlist[ismf] = zmid
      ismf=2      
      logm,phi,phip,phim,szmin,szmax =getphisg(zcen,i)
      if (logm[0]>1):
             phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
             phip *=(hrat*hrat*hrat)
             phim *=(hrat*hrat*hrat)
             logm -=  2.*N.log10(hval/0.70) 
          
             ax.errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="SDSS-FSPS %3.2f<z<%3.2f"%(szmin,szmax))
             smfflag[ismf]=1
             zminlist[ismf] = szmin
             zmaxlist[ismf] = szmax
      ismf=3       
      logm,phi,phip,phim,pzmin,pzmax=getphim(zcen,i)
      if (logm[0]>1):
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         phip *=(hrat*hrat*hrat)
         phim *=(hrat*hrat*hrat)
         logm -=  2.*N.log10(hval/0.70) 
          
         ax.errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="%s %3.2f<z<%3.2f"%(smftype[ismf],pzmin,pzmax))
         smfflag[ismf]=1
         zminlist[ismf] = pzmin
         zmaxlist[ismf] = pzmax
      ismf = 4
      #add zfourge
      logm,phi,phip,phim,zfmin,zfmax =getphit(zcen,i)
      if (logm[0]>1):
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         phip *=(hrat*hrat*hrat)
         phim *=(hrat*hrat*hrat)
         logm -=  2.*N.log10(hval/0.70) 
          
         ax.errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="%s %3.2f<z<%3.2f"%(smftype[ismf],zfmin,zfmax))
         smfflag[ismf]=1
         zminlist[ismf] = zfmin
         zmaxlist[ismf] = zfmax
      ismf = 5 #ultravista
      logm,phi,phip,phim,uzmin,uzmax =getphiuv(zcen,i)
      if (logm[0]>1):
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         phip *=(hrat*hrat*hrat)
         phim *=(hrat*hrat*hrat)
         logm -=  2.*N.log10(hval/0.70) 
          
         ax.errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="%s %3.2f<z<%3.2f "%(smftype[ismf],uzmin,uzmax))
         smfflag[ismf]=1
         zminlist[ismf] = uzmin
         zmaxlist[ismf] = uzmax
       #schechter function to ultravista
      logm,phi,uzmin,uzmax = getphiuv_sch(zcen,i)
      if (logm[0]>1):
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         logm -=  2.*N.log10(hval/0.70) 
          
         ax.plot(logm,phi,color=scollist[ismf],label="%s Schechter %3.2f<z<%3.2f "%(smftype[ismf],uzmin,uzmax))
         smfflag[ismf]=1   
      ismf = 6 #henriques
      logm,phi,phip,phim,zhen = getphihen(zcen,i)
      if (logm[0]>1):
      	 phi    *= hval*hval*hval #units [h/Mpc]^3
         phip   *= hval*hval*hval
         phim   *= hval*hval*hval
         logm   -= 2*N.log10(hval) #units [M*/h^2]
         ax.errorbar(logm,phi,yerr=[phim,phip],xerr=0.0,fmt=' ',marker=smarkerlist[ismf],color=scollist[ismf],label="%s z=%3.2f "%(smftype[ismf],zhen))
         smfflag[ismf]=1
         zminlist[ismf] = zhen
      ismf = 7 #vipers 
      logm,phi,vzmin,vzmax = getphiv(zcen,i)
      if (logm[0]>1): #they assume h=0.70
         phi *=(hrat*hrat*hrat) #they use (h70^-1 mpc)^3 for vol, h70^-2 Mo
         logm -=  2.*N.log10(hval/0.70) 
      
         ax.plot(logm,phi,color=scollist[ismf],label="%s Schechter %3.2f<z<%3.2f "%(smftype[ismf],vzmin,vzmax))
         smfflag[ismf]=1
         zminlist[ismf]=vzmin
         zmaxlist[ismf] = vzmax 
         
      ax.set_xlim(8.2,12.)
      ax.set_ylim(1.e-5,0.02)
      ax.set_yscale("log")
   
      ax.set_xlim(8.2,12.)
      ax.set_ylim(1.e-5,0.04)
      ax.set_yscale("log")
      ax.text(10.6,0.02,r'$\bar{z}_{\rm sim}$ = %3.2f'%(zcen))
      ax.text(10.6,0.013,'%s'%(galtype[i]),color=collist[i])
      ax.text(10.6,0.008,'%s'%(runname),color=collist[i])      
      ax.set_ylabel(" $\Phi$ [Mpc${}^{-3}$/dex]")      
      ax.legend(loc=3,fontsize='10',frameon=False)
      ax.set_xlabel(r'M* [log $M_\odot$]')

      plt.tight_layout()
      plt.savefig("smf4sims_%s_%d_%s.pdf"%(galtype[i],zcen*100.5,runname))
        
      plt.close("all")  


    
def getmsmh(fname="inputfile.dat",ratflag=1):
    """
     M*(Mh) for centrals
     ratflag=1 :M*/Mh as fn of Mh
     ratflag= 0: M*   as fn of Mh
    """        
   
    ff = open(fname)
    gals = N.loadtxt(ff,usecols=(0,5,6))
    ff.close()
    jj = N.nonzero(gals[:,1]==0)[0] #get centrals
    logmstar    = gals[jj,0]
    logmh    = gals[jj,2]
    mhbin    = N.linspace(logmh.min(),logmh.max(),40)
    mstarave = N.zeros(40,dtype="float")
    mstarlist =[]
    ngaltot = 0
    for i in range(39):
        jj = N.nonzero((logmh> mhbin[i])&(logmh<=mhbin[i+1]))[0]
        mstarave[i] = (10**logmstar[jj]).sum()/(jj.size +1.e-10)
        mstarlist.append(10**(logmstar[jj]-ratflag*logmh[jj]))
        ngaltot += jj.size
    mhbin +=(mhbin[1]-mhbin[0])/2.
    mhbin = N.delete(mhbin,39)        
    return(mhbin,mstarave,mstarlist,ngaltot)
        
    
def rhox(xval):
    """
    from M. White
    nfw profile, x=r/rs, rho_0 = 1
    """
    return(1./(xval*(1+xval)*(1+xval)))

def rhobar(xval):
    """ mean interior density
        from M. White
    """
    return(3*(N.log(1+xval)/(xval*xval*xval) -1./((1+xval)*xval*xval)))

def mmean(rho,rho0):
   """
   from M. White 
   Returns the mass enclosed within the radius at which the mean interior
   density falls to rho (times the critical density). 
   for defns like delta = 200c or 500omegab*/
   """    
   xlo = 1.e-10
   rlo = rho0*rhobar(xlo)-rho
   xhi = 100.
   rhi = rho0*rhobar(xhi)-rho
   xmid = (xlo+xhi)/2.0
   rmid = rho0*rhobar(xmid)-rho
   if (rmid*rhi < 0):
      xlo=xmid
      rlo=rmid
   else:
      xhi = xmid
      rhi = rmid
   while (xhi-xlo>1.e-3*xmid):
    xmid = (xlo+xhi)/2.0
    rmid = rho0*rhobar(xmid)-rho
    if (rmid*rhi < 0):
      xlo=xmid
      rlo=rmid
    else:
      xhi = xmid
      rhi = rmid
   tmp = 4*N.pi*rho0*(N.log(1+xmid)-xmid/(1+xmid))
   return(tmp)

    
def mvirtom200(logmval,zcen=0.45):
    """
    from M. White
    given Mvir
    convert to m200 to feed to moster mf
    from martin white, uses nwf profile to convert as in 2001 paper mass of halo
    http://arxiv.org/abs/astro-ph/0011495
    assume mvir and m200 not so different that concentration for
    one can be used for other.
    """
    omm = 0.31
    #use mvir to get c, difference with m200 too small to care
    mvirstart = 10**logmval
    c = 10*N.power(mvirstart/3.42e12,-0.2/(1+zcen))
    omz = omm/(omm + (1-omm)/(1+zcen)**3)
    DelC = 18*N.pi*N.pi+82*(omz-1.)-39*(omz-1)*(omz-1)
    rho0 = (200./3.)*N.power(c,3.0)/(N.log(1+c)-c/(1+c))
    mvir = 4*N.pi*(200./3.)*c*c*c
    mvirdivm200 = mmean(DelC*1.,rho0)/mvir
    return(1/mvirdivm200)
    
    
def getmoster(logmval,zcen=0.25,convflag=1):
    """
    units are in M* for everything, including Mh (no h)
    Moster, Naab & White, 2013
Galactic star formation and accretion histories from matching galaxies to dark matter haloes, MNRAS, 428, 3121
http://arxiv.org/abs/1205.5807
    Mh = M200c
    Bruzual-Charlot
    tuned to perez-gonzalez 08 
    http://arxiv.org/abs/0709.1354
    664 arcmin^2
    and santini(2011)
    http://arxiv.org/abs/1111.5728
    33 arcmin^2
    table 1
    """
    mvirconvert = mvirtom200(logmval,zcen) #m200 = mvirconvert*mvir
    if (convflag==0):
        mvirconvert=1.
    m10 = 11.590
    m11 = 1.195
    n10 = 0.0351
    n11 = -0.0247
    beta10 = 1.376
    beta11 = -0.826
    gamma10 = 0.608
    gamma11 = 0.329

    #errors
    em10 = 0.236
    em11 = 0.353
    en10 = 0.0058
    en11 = 0.0069
    ebeta10 = 0.153
    ebeta11 = 0.225
    egamma10 = 0.059
    egamma11 = 0.173

    
    zrat = zcen/(zcen+1)
    norm  = n10+n11*zrat
    M1    = m10 +m11*zrat  #this is log
    beta  = beta10+beta11*zrat
    gamma = gamma10+gamma11*zrat

    rat = 2*norm/(10**(beta*(M1-logmval)) + 10**(gamma*(logmval-M1)))
    #sign of beta is correct, logmval and m1 switched
    rat *=mvirconvert
    #scatter sigma_m(logm) = 0.15
    return(rat)

def fpb(xval,zcen=0.25):
    """
    x is log(Mh/M1)
    equations 3,4 and and start of section 5
    of behroozi, wechsler, conroy
    The Average Star Formation Histories of Galaxies in Dark Matter Halos from z=0-8, http://arxiv.org/abs/1207.6105
    From Peter B: alpha should be -alpha in Eq. 3.  Note also that the exponent of gamma is intended to be applied after the logarithm is taken, not before.
  #parameters in section 5 of behroozi et al  
    """
    a = 1/(1+zcen)
    nu = N.exp(-4*a*a)
    #alpha      = -1.474 + 1.339*(a-1)*nu #older version of BWC
    #delta      =  3.529 + (4.152*(a-1) +1.122*zcen)*nu
    #gamma      =  0.395 + (0.766*(a-1) +0.435*zcen)*nu
    alpha      = -1.412 + 0.731*(a-1)*nu
    delta      =  3.508 + (2.608*(a-1) -0.043*zcen)*nu
    gamma      =  0.316 + (1.319*(a-1) +0.279*zcen)*nu

    tmp = -N.log10(10**(alpha*xval) +1) + delta* (N.log10(1+N.exp(xval)))**gamma /(1+N.exp(10**(-xval)))

    return(tmp)    

    
def getms_pb(mhbin,zcen=0.25):
     """
     Behroozi, Wechsler, Conroy fitting function
        The Average Star Formation Histories of Galaxies in Dark Matter Halos from z=0-8, http://arxiv.org/abs/1207.6105
    From Peter B: alpha should be -alpha in Eq. 3.  Note also that the exponent of gamma is intended to be applied after the logarithm is taken, not before.
 
     Mh is Mvir (M_delta)    
     """
     #parameters in section 5 of behroozi et al
     a= 1/(1+zcen)
     nu = N.exp(-4*a*a)
     
     #logm1      = 11.539 + (-1.751*(a-1)  -0.329*zcen)*nu older version of BWC
     #logepsilon = -1.785 + (-0.074*(a-1) + -0.048*zcen)*nu -0.179*(a-1)
     logm1      = 11.514 + (-1.793*(a-1)  -0.251*zcen)*nu
     logepsilon = -1.777 + (-0.006*(a-1) + -0.0*zcen)*nu -0.119*(a-1)

     mstarmh   =logm1 + logepsilon + fpb(mhbin-logm1,zcen) - fpb(0,zcen)       
     return(10**(mstarmh-mhbin))

def plotboxwhisker_mult(outid="thisrun",fname="inputfile.dat",zcen=0.45,runname="runname"):
    # instead use M* vs Mh not M*/Mh
    #box is lower and upper quartiles
    # line is median
    # whisker is by choice 10,90 percentiles
    # outliers not plotted ("showfliers='false')
    from matplotlib.ticker import MultipleLocator,FormatStrFormatter
    fig,ax = plt.subplots(1,1)
    ax.set_xlim(10.7,15.)
    ax.set_ylim(1.e7,5.e13)
    majorLocator= MultipleLocator(1)
    minorLocator=MultipleLocator(0.2)
    majorFormatter=FormatStrFormatter('%d')
    mhbin,mstarave,mstarlist,ngaltot = getmsmh(fname,0) #just M* not M*/Mh
    ax.boxplot(mstarlist,whis=[10,90],showfliers=False,positions=mhbin,widths=0.1)
    #mhbin is log M halo
    msrat_moster = N.zeros(mhbin.size,dtype="float")
    for im in range(mhbin.size):
        msrat_moster[im] = getmoster(mhbin[im],zcen)
    ax.plot(mhbin,msrat_moster*(10**mhbin),'m-.',label='MNW 2013')
    for im in range(mhbin.size):
        msrat_moster[im] = getmoster(mhbin[im],zcen,0)
    ax.plot(mhbin,msrat_moster*(10**mhbin),color='cornflowerblue',linestyle=':',label=r'MNW13, no $M_{200} \rightarrow M_{vir}$')
    msrat_pb = getms_pb(mhbin,zcen)
    ax.plot(mhbin,msrat_pb*(10**mhbin),'m',label='BWC 2013')
    #plot 15 in Behroozi, Wechsler,Conroy 2013, 0-8 paper
    ax.set_yscale("log")    
    ax.set_xlabel(r'log$M_{\rm vir}$ $[M_\odot]$ ')
    ax.set_ylabel(r' $M_* $')
    ax.text(mhbin[25],1.e8,r'$\bar{z}$= %3.2f'%(zcen))
    ax.text(mhbin[25],5.e8,r'$%d$ galaxies'%(ngaltot))    
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.legend(loc=1,fontsize=12,frameon=False)
    plt.tight_layout()
    plt.savefig("stellar_halo_noratio_%s_%s.pdf"%(outid,runname))
    plt.close("all")
     
      
def plotboxwhisker(outid="thisrun",fname="inputfile.dat",zcen=0.45,runname="runname"):
    #box is lower and upper quartiles
    # line is median
    # whisker is by choice 10,90 percentiles
    # outliers not plotted ("showfliers='false')
    from matplotlib.ticker import MultipleLocator,FormatStrFormatter
    fig,ax = plt.subplots(1,1)
    ax.set_xlim(10.7,15.)
    ax.set_ylim(0.001,0.05)
    majorLocator= MultipleLocator(1)
    minorLocator=MultipleLocator(0.2)
    majorFormatter=FormatStrFormatter('%d')
    mhbin,mstarave,mstarlist,ngaltot = getmsmh(fname)
    ax.boxplot(mstarlist,whis=[10,90],showfliers=False,positions=mhbin,widths=0.1)
    #mhbin is log M halo
    msrat_moster = N.zeros(mhbin.size,dtype="float")
    for im in range(mhbin.size):
        msrat_moster[im] = getmoster(mhbin[im],zcen)
    ax.plot(mhbin,msrat_moster,'m-.',label='MNW 2013')
    for im in range(mhbin.size):
        msrat_moster[im] = getmoster(mhbin[im],zcen,0)
    ax.plot(mhbin,msrat_moster,color='cornflowerblue',linestyle=':',label=r'MNW13, no $M_{200} \rightarrow M_{vir}$')
    msrat_pb = getms_pb(mhbin,zcen)
    ax.plot(mhbin,msrat_pb,'m',label='BWC 2013')
    #plot 15 in Behroozi, Wechsler,Conroy 2013, 0-8 paper
    ax.set_yscale("log")    
    ax.set_xlabel(r'log$M_{\rm vir}$ $[M_\odot]$ ')
    ax.set_ylabel(r' $M_*/M_{\rm vir} $')
    ax.text(14,0.03,r'$\bar{z}$= %3.2f'%(zcen))
    ax.text(14,0.04,r'$%d$ galaxies'%(ngaltot))    
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.legend(loc=1,fontsize=12,frameon=False)
    plt.tight_layout()
    plt.savefig("stellar_halo_ratio_%s_%s.pdf"%(outid,runname))
    plt.close("all")

def plotcon(binx,biny,nhist2,ax,delchoice=0.5):
    """
    stolen from
    http://matplotlib.org/examples/pylab_examples/contour_demo.html
    """
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

    delta = delchoice
    nhist2 = N.rot90(nhist2) #setorigin=lower axes bottom to top
    nhist2 = N.flipud(nhist2) #nhist2.T() transpose for physics
    #http://oceanpython.org/2013/02/25/2d-histogram/
    cs = ax.contour(binx,biny,nhist2)
    return(cs)
    

def mstar_sfr(fname="inputfile.dat",hval=0.67,slopeval=0.,shiftval=0.,runname="runname"):
        """
        from M. White
        """
        ff = open(fname)
        gals = N.loadtxt(ff, usecols=(0,1,4))
        ff.close()
        #logm, sfr, zval, logMh
        zmin = gals[:,2].min()
        zmax = gals[:,2].max()
        zcen = gals[:,2].sum()/gals[:,2].size
        #convert to h70 stellar mass
        logm = gals[:,0]+2*N.log10(hval/0.70)
        logsfr = N.log10(gals[:,1]+1.e-8) + 2*N.log10(hval/0.70)
        fig,ax=plt.subplots(1,1)
        maxms = 12.
        minms = 8.5
        maxsf = 3.
        minsf = -5.
        nbin = 40
    # seems reversed
        nhist2,binx,biny = N.histogram2d(logm,logsfr,nbin,range=([minms,maxms],[minsf,maxsf]))
        nhist2 = N.arcsinh(nhist2)
        #nhist2 = nhist2*1./nhist2.sum()
        binx +=(binx[1]-binx[0])/2.
        biny +=(biny[1]-biny[0])/2.
        binx = N.delete(binx,nbin)
        biny = N.delete(biny,nbin)
        cs=plotcon(binx,biny,nhist2,ax)
        fig.colorbar(cs)        
        logms = N.arange(8.5,12.5,0.1)
        logsf = -0.49+ 0.65*(logms - 10)+1.07*(zcen-0.1)
        ax.plot(logms,logsf,'--',lw=2,label=r'log(sfr) =-0.49+ 0.65(log $M^*$ - 10)+1.07($\bar{z}$-0.1)')
        logsf = -0.49+ (0.65+slopeval)*(logms - 10)+1.07*(zcen-0.1)+shiftval
        ax.plot(logms,logsf,lw=2,label=r'log(sfr) =-0.49+ %3.2f(log $M^*$ - 10)+1.07($\bar{z}$-0.1)+%3.2f'%(slopeval+0.65,shiftval))
        
        ax.legend(loc=4)
        ax.set_xlabel('log M* $[M_\odot]$')
        ax.set_ylabel('log SFR ')
        ax.set_title(r'$\bar{z}_{\rm sim}$ = %3.2f'%(zcen))                         
        plt.tight_layout()
        plt.savefig("mstarsfr_%d_%s.pdf"%(zcen*100.5,runname))
          
        plt.close("all")

def ssfr_slice(fname="inputfile.dat",hval=0.67,runname="runname"):
        """
        slice of fixed mstar (log 10 bin width 0.5)
        """
        ff = open(fname)
        gals = N.loadtxt(ff, usecols=(0,1,4,6))
        ff.close()
        #logm, sfr, zval, logMh
        zmin = gals[:,2].min()
        zmax = gals[:,2].max()
        zcen = gals[:,2].sum()/gals[:,2].size
        #convert to h70 stellar mass
        logm  =  gals[:,0]+2*N.log10(hval/0.70)
        logsfr  = N.log10( gals[:,1]+1.e-8) + 2*N.log10(hval/0.70)
        fig,ax=plt.subplots(2,2)
        minssfr = -13
        maxssfr = -8
        for i in range(4):
          mstarlow = 9.5+i*0.5      
          jj = N.nonzero((logm>mstarlow)&(logm<=mstarlow +0.5))[0]
          ssfr = logsfr[jj]-logm[jj]
          nbin = 50
          #maxssfr = ssfr.max()
          #minssfr = ssfr.min()
          nhist,bins = N.histogram(ssfr,nbin,range=(minssfr,maxssfr))
          nhist = nhist*1./jj.size #fraction of all galaxies in M*
          #range, including those outside for minssfr,maxssfr range
          bins +=(bins[1]-bins[0])/2.
          bins = N.delete(bins,nbin)
          ax[i/2,i%2].set_xlim(minssfr,maxssfr)
          ax[i/2,i%2].step(bins, nhist,'k',where='mid')
          if (i!=1):
             ax[i/2,i%2].set_title("%3.1f<logM*<%3.1f "%(mstarlow,mstarlow+0.5),fontsize="10")
          if (i==1):
             ax[i/2,i%2].set_title(r'%3.1f<logM*<%3.1f,  $\bar{z}=$%3.2f '%(mstarlow,mstarlow+0.5,zcen),fontsize="10")

          ax[i/2,i%2].set_ylabel('frac of %d gals'%(ssfr.size))

          ax[i/2,i%2].set_xlabel('$log_{10}$ SSFR [yr${}^{-1}$]')
        plt.tight_layout()
        plt.savefig("ssfr_z%d_%s.pdf"%(zcen*100.5,runname))
          
        plt.close("all")
        

def runsuite(zcen=0.45, fname="inputfile.dat",hval=0.7,omm=0.31,slopeval=0.,shiftval=0, boxside=256,runname="runname",delz=0.1,ramin=-2,ramax=2,decmin=-2,decmax=2):
    """
    inputfile:
    logM*(0) [M_o], sfr(1) [M_o/yr],ra(2),dec(3),zred(4),ifsat(5),logmh(6)[M_o]
    """
    #4 smf's
    plot4sep(zcen,fname,hval,omm,slopeval,shiftval, boxside,runname,delz,ramin,ramax,decmin,decmax)
    plot4tog(zcen,fname,hval,omm, slopeval,shiftval,boxside,runname,delz,ramin,ramax,decmin,decmax)          
    plotboxwhisker("z%d"%(zcen*100.5),fname,zcen,runname)    #M*/Mh of Mh
    plotboxwhisker_mult("z%d"%(zcen*100.5),fname,zcen,runname)    #M* of Mh   
    mstar_sfr(fname,hval,slopeval,shiftval,runname)
    ssfr_slice(fname,hval,runname)
    #plotting M*/Mh all different redshifts has to be tailored by hand
    print "bwc_comp from Behroozi/Wechlsler/Conroy compilation Fig.3"
    print "    arxiv   http://arxiv.org/abs/1207.6105"
    print " "
    print "primus_bc03 from www.behroozi.com/data.html"
    print "    behroozi-2013-data-compilation calculated in "
    print "    http://arxiv.org/abs/1207.6105 , http://arxiv.org/abs/1209.3013"
    print " "
    print "sdss_gal from Moustakas++ 2013, table 3,http://arxiv.org/abs/1301.1688"
    print " "
    print "primus_fsps from  Moustakas++ 2013, table 4,http://arxiv.org/abs/1301.1688"
    print " "
    print "zfourge from Tomczak++2014, table 1, http://arxiv.org/abs/1309.5972"
    print "survey details, Straatman et al, submitted"
    print " "
    print "cosmos/ultravista from Muzzin++2013, paper http://arxiv.org/abs/1303.4409"
    print " data http://cosmos2.phy.tufts.edu/~danilo/Downloads.html"
    print " "
    print " hen15 Henriques compilation http://arxiv.org/abs/1410.0365"
    print " data http://galformod.mpa-garching.mpg.de/public/LGalaxies/figures_and_data.php "
    print " "
    print "vipers16 from    Moutard et al, 2016, table 2, http://arxiv.org/abs/1602.05917"
    print " "
    print "default color cut from Moustakas et al 2013, modified to"
    print "log(sfr) =-0.49+ (0.65+%3.2f)*(log M* - 10)+1.07(zbar-0.1)+%3.2f"%(slopeval,shiftval)
