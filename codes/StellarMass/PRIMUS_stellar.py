__author__ = "joanne cohn"
__email__  = "jcohn@berkeley.edu"
__version__= "1.0"

def overview():
   """
plots a simulated data set stellar mass function (histogram) along with that in
Moustakas et al, 2013 (points)

FILES NEEDED (3 including this one):
 import this file to run commands
 have file "Mous_13_table4.txt" (in same directory)
 have your simulation data file fname (explained below), e.g. 'galaxies.dat'

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
source for data in file Mous_13_table4.txt is
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
   zrange = N.array([0.2,0.3,0.4,0.5,0.65,0.8,1.0])
   jjz = N.nonzero(zcen>=zrange)[0]
   if ((jjz.size ==0)|(zcen>1.0)):
     print "z = %3.2f not in range"%(zcen)
   if (jjz.size > 1):
      jjz = jjz.max()
   print "using PRIMUS range %3.2f < z < %3.2f "%(zrange[jjz],zrange[jjz+1])
   
   ff = open("Mous_13_table4.txt")
   phimou = N.loadtxt(ff, usecols=(0,jjz+1))
   ff.close()
   # Now color, just need to start at right place
   if (addcolor==0):
     jjkeep = N.arange(33) #all
   if (addcolor==1):
     jjkeep = N.arange(33,61)  #red
   if (addcolor==2):
     jjkeep = N.arange(61,93)   #blue
   logm_mou = phimou[jjkeep,0]
   logphi   = phimou[jjkeep,1]
   jj = N.nonzero(logphi<0)[0]
   logm_mou = logm_mou[jj]
   logphi   = logphi[jj]
   return(logm_mou,logphi,jjz)
 
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
   if (boxside>0):
       print "using periodic box, side %8.2f Mpc/h"%(boxside)
       vol = boxside*boxside*boxside

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
     jj = N.nonzero(N.log10(sfr+1.e-16)<-0.49+0.65*(logstell-10)+1.07*(zcen-0.1))[0]
   if (addcolor==2):  
        jj = N.nonzero(N.log10(sfr+1.e-16)>=-0.49+0.65*(logstell-10)+1.07*(zcen-0.1))[0]
   logstell = logstell[jj]     

   nbin = 50
   nhist,bins = N.histogram(logstell,nbin,range=(8.3,12.3))
   bins += (bins[1]-bins[0])/2.
   bins = N.delete(bins,nbin)
   ngalact = nhist*1./(vol*(bins[1]-bins[0]))
   
   galkind =("all","red","blue")
   logm_mou,logphi,jjz = phi_moutab(zcen,addcolor)
   return(bins,ngalact,logm_mou,logphi,jjz)


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
    for i in range(3):
      f,ax = plt.subplots(1,1)
      ax.set_xlim(8.8,12.0)
      ax.set_ylim(1.e-5,0.1)
      ax.set_yscale("log")        
      bin_centers,ngalact,logm_mou,logphi,jjz=teststellar(zcen,i,fname,hval,boxside,ramin,ramax,decmin,decmax,delz)
      ax.step(bin_centers, ngalact,collist[i],label=r'simulation $z_{sim}$= %3.2f'%(zcen))

      zrange = N.array([0.2,0.3,0.4,0.5,0.65,0.8,1.0])      
      ax.plot(logm_mou,10**logphi,cshapelist[i],label=r'PRIMUS %3.2f<z<%3.2f'%(zrange[jjz],zrange[jjz+1]))
      ax.text(9.2,2.e-4,collfull[i],color=collist[i])

      ax.legend(loc=3)
      ax.set_yscale("log")
      ax.set_xlabel("$M_* [M_o/h_{70}^2]$")
      ax.set_ylabel("$\Phi$[$h_{70}^3$/Mpc${}^3$]")
      plt.tight_layout()
      plt.savefig("mstar_z%d_%s.pdf"%(zcen*101,coltype[i]))
      plt.close("all")




