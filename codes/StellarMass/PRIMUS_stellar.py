__author__ = "joanne cohn"
__email__  = "jcohn@berkeley.edu"
__version__= "1.1"

def overview():
   """
plots a simulated data set stellar mass function (histogram) along with that in
Moustakas et al, 2013 (points plus error bars)
for z<0.2, this is SDSS-GALEX data, for z>0.2 it is PRIMUS data

FILES NEEDED (5 including this one):
 1. import this file to run commands
 have files in this directory (ascii files):
 2. Mous_13_table4.txt 
 3. Mous_13_errplus.txt (new to 1.1) 
 4. Mous_13_errminus.txt  (new to 1.1)
 5. have your simulation data file fname (explained below), e.g. 'galaxies.dat'

USAGE:

 plot3sep(zcen,fname,hval,boxside,ramin,ramax,decmin,decmax,delz)

OUTPUT:
  stellar mass function plots for red, blue, all at redshift zcen
  colors separated by Moustakas et al eqn 2,
  N.log10(sfr) = -0.49+0.65*(log10 M*-10)+1.07*(z-0.1)
  where sfr is in M_o/yr, M_* is in M_o
  (more precisely used factor of h70^{-2} for both)

ARGUMENTS of plot3sep():
  zcen = central mass of simulation sample
  fname = filename of data, in ascii (e.g. 'galaxies.dat'), more below
  hval  =hubble constant (e.g. 0.67)
  boxside: two choices depending upon data in fname
    periodic box, boxside = box side length in Mpc/h
    light cone,  any number < 0, to flag that will use ra and dec and delta z

  other entries to plot3sep
  for periodic box:
      every entry after boxside (which is positive) is ignored
      can be left out
  for light cone:
     ramin, ramax, decmin, decmax : minimum and maximum ra and dec
     delz : galaxies are taken in a region zcen +/- delz for light cone sample
     omega_m is assumed to be 0.31, if you use a much different value
    you might need to change the lightcone volume.

INPUT SIMULATION DATA FILE:
  (e.g. 'galaxies.dat')
 for periodic box, this file is list of galaxies in ascii, one line per galaxy
  log10 M_* [M_o]   sfr [M_o/Gyr]

 for light cone this file is list of galaxies in ascii, one line per galaxy
  log10 M_* [M_o]   sfr [M_o/Gyr]    ra     dec   redshift

INPUT OBSERVATIONAL DATA FILE:
source for data in file Mous_13_table4.txt and z<0.2 data used is
arXiv:1301.1688, ApJ 2013, 767, 50
Moustakas, J.,Coil, A. L., Aird, J., Blanton, M. R., Cool, R. J., Eisenstein, D. J., Mendez, A. J., Wong, K. C., Zhu, G., Arnouts, S.
PRIMUS: Constraints on Star Formation Quenching and Galaxy Merging, and the Evolution of the Stellar Mass Function from z = 0-1
--there are also errors in there, not included here.
 
   """
import numpy as N
import os
import time
import sys
import matplotlib
import matplotlib.pyplot as plt



def chiofz(zval=0.45,omm=0.31):
  """
  comoving distance to redshift zval
  omega matter = omm
  use for volume if specify region with ra/dec
  """
  Nint = 300000
  zp1int = N.linspace(1,zval+1,Nint)
  ez = N.sqrt(omm*zp1int*zp1int*zp1int + (1-omm))
  tmp = 2997.925* N.trapz(1/ez, dx = zp1int[1]-zp1int[0])
  return(tmp)

    


def phi_moutab(zcen=0.45,addcolor=0):
   """
   read data of moustakas et al, 2013 take redshift range
   surrounding to zcen chosen
    

   addcolor=0 all
   addcolor = 1 red
   addcolor = 2 blue
   changed order in moustakas data file too
   units are: log M, log Phi
   M is in units of M_o/h70^2
       -if your M* is in units of M_o, multiply your M* by (0.7/h)^2
   Phi is in units of [h70/Mpc]^3, 
       -if your volume is in [Mpc/h]^3, multiply your volume by (1/h)^3 /(h/0.70)^3
       -*or* divide your Phi by same factor 
   """
   # now need to find right redshift and type
   # first redshift
   print "using FSPS stellar masses"
   zrange = N.array([0.01,0.2,0.3,0.4,0.5,0.65,0.8,1.0])
   jjz = N.nonzero(zcen>=zrange)[0]
   if ((zcen>1.0)|(jjz.size==0)):
     print "z = %3.2f not in range"%(zcen)
     return()
   if (jjz.size==1):
     # use GALEX-SDSS
     print "using GALEX-SDSS"
     logm_mou = 9.0+ N.arange(31)*1./10.  #Primus table 3, Moustakas ++13
    
     if (addcolor==0): #all galaxies
         logphi = N.array([-1.899,-1.923,-1.970,-2.031,-2.055,-2.106,-2.144,-2.179,-2.188,-2.216,-2.234,-2.235,-2.262,-2.252,-2.285,-2.317,-2.365,-2.419,-2.504,-2.607,-2.728,-2.888,-3.104,-3.332,-3.606,-3.953,-4.363,-4.778,-5.255,-5.87,-6.49])

         logphi_plus = N.array([0.017,0.017,0.015,0.015,0.014,0.012,0.012,0.012,0.010,0.0086,0.0080,0.0069,0.0063,0.0056,0.0051,0.0047,0.0044,0.0041,0.0040,0.0039,0.0040,0.0043,0.0049,0.0059,0.0080,0.012,0.020,0.033,0.060,0.010,0.030])
         logphi_minus = N.array([-0.017,-0.016,-0.015,-0.014,-0.013,-0.012,-0.011,-0.012,-0.010,-0.0084,-0.0078,-0.0068,-0.0062,-0.0056,-0.0051,-0.0046,-0.0044,-0.0041,-0.0040,-0.0039,-0.0040,-0.0043,-0.0048,-0.0059,-0.0079,-0.012,-0.019,-0.031,-0.053,-0.010,-0.020])         
     if (addcolor==1): # red
         logphi = N.array([-2.495,-2.486,-2.485,-2.523,-2.576,-2.603,-2.634,-2.642,-2.652,-2.655,-2.649,-2.614,-2.607,-2.5640,-2.5640,-2.5800,-2.6050,-2.6450,-2.7050,-2.7860,-2.8840,-3.0190,-3.2090,-3.4130,-3.6670,-4.002,-4.401,-4.806,-5.296,-5.93,-6.61])
         logphi_plus = N.array([0.048,0.044,0.038,0.037,0.033,0.030,0.026,0.028,0.021,0.018,0.015,0.013,0.011,0.0089,0.0077,0.0069,0.0062,0.0057,0.0053,0.0050,0.0049,0.0050,0.0055,0.0065,0.0085,0.013,0.021,0.034,0.063,0.10,0.40])       
         logphi_minus = N.array([-0.043,-0.041,-0.035,-0.034,-0.031,-0.028,-0.025,-0.026,-0.020,-0.017,-0.015,-0.012,-0.011,-0.0087,-0.0076,-0.0068,0-.0061,-0.0056,-0.0052,-0.0050,-0.0049,-0.0050,-0.0054,-0.0064,-0.0084,-0.012,-0.020,-0.032,-0.056,-0.10,-0.20])                   
              
     if (addcolor==2): #blue
          logphi = N.array([-2.026,-2.062,-2.129,-2.201,-2.211,-2.272,-2.313,-2.362,-2.371,-2.4120,-2.4450,-2.4700,-2.5240,-2.5410,-2.6090,-2.6600,-2.7370,-2.8110,-2.9340,-3.0770,-3.2500,-3.4720,-3.769,-4.102,-4.487,-4.930,-5.437,-5.98,-6.30,-6.77,-7.09])
          logphi_plus = N.array([0.018,0.017,0.015,0.014,0.014,0.012,0.012,0.011,0.011,0.0092,0.0090,0.0079,0.0074,0.0071,0.0066,0.0063,0.0062,0.0059,0.0061,0.0064,0.0071,0.0085,0.011,0.016,0.024,0.042,0.079,0.20,0.30,0.60,1.000])      
          logphi_minus = N.array([-0.017,-0.016,-0.015,-0.014,-0.013,-0.012,-0.012,0-.011,-0.011,-0.0090,-0.0088,-0.0078,-0.0072,-0.0070,-0.0065,-0.0062,-0.0061,-0.0059,-0.0060,-0.0063,-0.0070,-0.0084,-0.010,-0.015,-0.023,-0.038,-0.067,-0.10,-0.20,-0.30,-0.40])             

     return(logm_mou,logphi,logphi_plus,logphi_minus,jjz)
   if (jjz.size > 1):
      jjz = jjz.max()
   print "using PRIMUS range %3.2f < z < %3.2f "%(zrange[jjz],zrange[jjz+1])
   
   ff = open("Mous_13_table4.txt")
   phimou = N.loadtxt(ff, usecols=(0,jjz))
   ff.close()
   ff = open("Mous_13_errminus.txt")
   phimou_minus = N.loadtxt(ff, usecols=(0,jjz))
   ff.close()
   ff = open("Mous_13_errplus.txt")
   phimou_plus = N.loadtxt(ff, usecols=(0,jjz))
   ff.close()
   # Now color, just need to start at right place
   if (addcolor==0):
     jjkeep = N.arange(33) #all
   if (addcolor==1):
     jjkeep = N.arange(33,61)  #red
   if (addcolor==2):
     jjkeep = N.arange(61,94)   #blue
   logm_mou = phimou[jjkeep,0]
   logphi   = phimou[jjkeep,1]
   logphi_plus = phimou_plus[jjkeep,1]
   logphi_minus = phimou_minus[jjkeep,1]
   jj = N.nonzero(logphi<0)[0]
   logm_mou = logm_mou[jj]
   logphi   = logphi[jj]
   logphi_plus = logphi_plus[jj]
   logphi_minus = logphi_minus[jj]
   return(logm_mou,logphi,logphi_plus,logphi_minus,jjz)
 
def teststellar(zcen=0.45,addcolor=0,fname="galshort.dat",hval=0.67,boxside=100,ramin=-2,ramax=-2,decmin=2,decmax=-2,delz=0.02):
   """
   usage:
   teststellar(zval,addcolor,inputfil,hval,boxside,ra_min,ra_max,dec_min,dec_max,delz)

   zcen: central redshift for simulation
   fname:  galaxy data file, more below
   addcolor=0  all galaxies
   addcolor=1  red
   addcolor=2  blue
   boxside: positive for periodic box, negative for light cone
   boxside = side of box when periodic box
   **or**
   boxside<0 (e.g. -1) if light cone
   ra_min,ra_max,dec_min,dec_max :min/max ra and dec for light cone [ignored for per]
   delz = use z range zcen+delz > z < zcen-delz for light cone [ignored for per]
   this is n(M) not N(>M)
   gal units M_o
   primus table 4
   color log sfrmin -0.49 + 0.65 (logM* - 10) +1.07 *(z-0.1)
   input file either has entries

   galaxy data file fname entries on each line, one per galaxy

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
   sfr /= 1.e9 # to yr^-1
   if (boxside < 0):
       print "using light cone"
       ra = gals[:,2]
       dec = gals[:,3]
       redz = gals[:,4]
       #need ra, dec, redshift,delz
       chimax = chiofz(zcen+delz)  #[Mpc/h]
       chimin = chiofz(zcen-delz)  #[Mpc/h]
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
   #Moustakas box is in units of [Mpc/h70]^3, we have [Mpc/h]^3,
   #    so divide volume by (h/h70)^3 = (h/(h/0.7))^3 = 0.7^3
   vol = vol/(0.7*0.7*0.7)
   # Mstar is in units of [M_o/h70^2], multiply by 
   logstell = logstell + 2.*N.log10(hval/0.70) 
   sfr      = sfr*(hval/0.70)**2

   jj = N.arange(logstell.size)
   #note color cut is assuming h70's in units
   if (addcolor==1):
     jj = N.nonzero(N.log10(sfr+1.e-16)<-0.49+0.65*(logstell-10)+1.07*(redz-0.1))[0]
   if (addcolor==2):  
        jj = N.nonzero(N.log10(sfr+1.e-16)>=-0.49+0.65*(logstell-10)+1.07*(redz-0.1))[0]
   logstell = logstell[jj]     

   nbin = 50
   nhist,bins = N.histogram(logstell,nbin,range=(8.3,12.3))
   bins += (bins[1]-bins[0])/2.
   bins = N.delete(bins,nbin)
   ngalact = nhist*1./(vol*(bins[1]-bins[0]))
   
   galkind =("all","red","blue")
   logm_mou,logphi,logphip,logphim,jjz = phi_moutab(zcen,addcolor)
  
   return(bins,ngalact,logm_mou,logphi,logphip,logphim,jjz,nhist.sum())


def plot3sep(zcen=0.45,fname="galshort.dat",hval=0.67,boxside=100,ramin=-2,ramax=-2,decmin=2,decmax=-2,delz=0.02):
    """
    plots separately, see comments for teststellar to get inputs.
    if running for periodic box, only need zcen (central sim redshift), fname(galaxy file), hval (hubble constant, i.e. 0.67 for planck cosmology), boxside (e.g. 100 if
    box is 100 Mpc/h on a side)
    galshort format is above, in teststellar comments
    """
    coltype=("all","red","blue")
    collfull=("all galaxies","red galaxies","blue galaxies")
    collist=('k','r','b')
    cshapelist=('ks','rs','bs')
    zrange = N.array([0.01,0.2,0.3,0.4,0.5,0.65,0.8,1.0])    
    for i in range(3):
      f,ax = plt.subplots(1,1)
      bin_centers,ngalact,logm_mou,logphi,logphip,logphim,jjz,ngal=teststellar(zcen,i,fname,hval,boxside,ramin,ramax,decmin,decmax,delz)
      ax.step(bin_centers, ngalact,collist[i],where='mid',label=r'simulation $z_{sim}$= %3.2f'%(zcen))
      
      phiplus = 10**(logphi+logphip)-10**(logphi)
      phiminus = 10**(logphi+logphim)-10**(logphi)        
      if (jjz==0):
        #ax.plot(logm_mou,10**logphi,cshapelist[i],label=r'SDSS-GALEX %3.2f<z<%3.2f'%(zrange[jjz],zrange[jjz+1]))
        ax.errorbar(logm_mou,10**logphi,yerr=[-phiminus,phiplus],xerr=0.0,fmt=' ',marker='s',color=collist[i],label=r'SDSS-GALEX %3.2f<z<%3.2f'%(zrange[jjz],zrange[jjz+1]) )          
      if (jjz>0):  
          #ax.plot(logm_mou,10**logphi,cshapelist[i],label=r'PRIMUS %3.2f<z<%3.2f'%(zrange[jjz],zrange[jjz+1]),)
        ax.errorbar(logm_mou,10**logphi,yerr=[-phiminus,phiplus],xerr=0.0,fmt=' ',marker='s',color=collist[i],label=r'PRIMUS %3.2f<z<%3.2f'%(zrange[jjz],zrange[jjz+1]) )
      ax.text(9.2,2.e-4,coltype[i],color=collist[i])
      ax.text(9.2,1.e-4,'%d galaxies'%(ngal),color=collist[i])
      ax.set_xlim(8.8,12.0)
      ax.set_ylim(1.e-5,0.1)
      ax.legend(loc=3)
      ax.set_yscale("log")
      ax.set_xlabel("$M_* [M_o/h_{70}^2]$")
      ax.set_ylabel("$\Phi$[$h_{70}^3$/Mpc${}^3$]")
      plt.tight_layout()
      plt.savefig("mstar_z%d_%s.pdf"%(zcen*101,coltype[i]))
      plt.close("all")




