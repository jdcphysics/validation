__author__ = "joanne cohn"
__email__  = "jcohn@berkeley.edu"
__version__= "1.0"

def overview():
   """
plots a simulated data set B-band luminsoity function (histogram) along with that in
Beare et al (points)

FILES NEEDED (3 including this one):
 import this file to run commands
 have file "Beare_15.txt" (in same directory)
 have your simulation data file fname (explained below), e.g. 'galaxies.dat'

USAGE:

 plot3sep(zcen,fname,hval,boxside,ramin,ramax,decmin,decmax,delz)

OUTPUT:
  B-band Luminosity function plots for red, blue, all at redshift zcen
  colors separated by 
  M_U - M_B > 1.074 - 0.18z -0.03 (M_B + 19.4) [red]
  

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
 if have B band:
 
 for periodic box, this file is list of galaxies in ascii, one line per galaxy
  M_B  M_U
  units: [AB], absolute magnitudes

 for light cone this file is list of galaxies in ascii, one line per galaxy
  M_B  M_U    ra     dec   redshift
  units:  [AB], absolute magnitudes for M_B,M_U
  if don't have B band:
     just set M_U to something random and throw out color plots

INPUT OBSERVATIONAL DATA FILE:
source for data in file Beare_tab8910.txt
http://arxiv.org/abs/1511.01580
R.A.Beare, M.J.I.Brown, K.A.Pimbblet, F.Bian, Y-T Lin
The z<1.2 optical luminosity function from a sample of ~410000 galaxies in Bootes
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

    


def phi_beare(zcen=0.45,addcolor=0):
   """
   read data of Beare et al, 2015 take redshift range
   surrounding to zcen chosen
   tables 8,9,10 for all, red, blue 

   addcolor=0 all
   addcolor = 1 red
   addcolor = 2 blue
   units:
   Magnitudes:  MB - 5 log h70
   Phi:  1.e-3 (h70/Mpc)^3 /dex

   """
   # now need to find right redshift and type
   # first redshift
   zrange = N.array([0.2,0.4,0.6,0.8,1.0,1.2])
   jjz = N.nonzero(zcen>=zrange)[0]
   if ((jjz.size ==0)|(zcen>=1.2)):
     print "z = %3.2f not in range"%(zcen)
   if (jjz.size > 1):
      jjz = jjz.max()
   print "using BOOTES range %3.2f < z < %3.2f "%(zrange[jjz],zrange[jjz+1])
   
   ff = open("Beare_tab8910.txt")
   phi = N.loadtxt(ff, usecols=(0,2*(jjz+1),2*(jjz+1)+1))
   ff.close()
   # Now color, just need to start at right place
   if (addcolor==0):
     jjkeep = N.arange(22) #all
   if (addcolor==1):
     jjkeep = N.arange(22,46)  #red
   if (addcolor==2):
     jjkeep = N.arange(46,66)   #blue
   Bmid    = phi[jjkeep,0]+0.25/2.  #Bmax = Bmin+0.25
   phival = phi[jjkeep,1] *1.e-3
   phierr  = phi[jjkeep,2] *1.e-3
   jj = N.nonzero(phival>0)[0]
   Bmid =  Bmid[jj]
   phival   = phival[jj]
   phierr   = phierr[jj]
   return(Bmid,phival,phierr,jjz)
 
def testlum(zcen=0.45,addcolor=0,fname="galshort.dat",hval=0.67,boxside=100,ramin=-2,ramax=-2,decmin=2,decmax=-2,delz=0.02):
   """
   usage:
   testlum(zval,addcolor,inputfile,hval,boxside,ra_min,ra_max,dec_min,dec_max,delz)

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
   
   BOOTES tables 8,9,10
     colors separated by 
  M_U - M_B > 1.074 - 0.18z -0.03 (M_B + 19.4) [red]

   input file either has entries

   galaxy data file fname entries on each line, one per galaxy

   [boxside > 0, i.e. fixed time, give boxside in units of Mpc/h]
    M_B   M_U    [ABS]

   [boxside < 0,use any boxside value < 0, lightcone]
   M_B   M_U    [ABS], ra, dec, redshift 
   """
   ff = open(fname)
   gals = N.loadtxt(ff)
   ff.close()
   magB      = gals[:,0]
   magU      = gals[:,1]
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
       magU = magU[jj]
       magB = magB[jj]
   if (boxside>0):
       print "using periodic box, side %8.2f Mpc/h"%(boxside)
       vol = boxside*boxside*boxside



   #units:
   #volume for Phi is in units of [Mpc/h70]^3, we have [Mpc/h]^3,
   #    so divide volume by (h/h70)^3 = (h/(h/0.7))^3 = 0.7^3
   vol = vol/(0.7*0.7*0.7)
   # Magnitudes are in units of MB- 5 log h70
   magB = magB - 5*N.log10(hval/0.70) 

   #note color cut is assuming h70's in units, but for h = 0.65 change is -0.005 
   jj = N.arange(magB.size)
   if (addcolor==1):
     jj = N.nonzero(magU - magB > 1.074 - 0.18*zcen -0.03* (magB + 19.4) )[0]
   if (addcolor==2):  
     jj = N.nonzero(magU - magB <= 1.074 - 0.18*zcen -0.03* (magB + 19.4) )[0]
   magB = magB[jj]     


   nbin = 50
   nhist,bins = N.histogram(magB,nbin,range=(-24.00,-17.75))
   bins += (bins[1]-bins[0])/2.
   bins = N.delete(bins,nbin)
   ngalact = nhist*1./(vol*(bins[1]-bins[0]))
   
   galkind =("all","red","blue")
   Bmid,phi,phierr,jjz = phi_beare(zcen,addcolor)
   phi += 1.e-10
   return(bins,ngalact,Bmid,phi,phierr,jjz)


def plot3sep(zcen=0.45,fname="galshort.dat",hval=0.67,boxside=100,ramin=-2,ramax=-2,decmin=2,decmax=-2,delz=0.02):
    """
    plots separately, see comments for testlum to get inputs.
    if running for periodic box, only need zcen (central sim redshift), fname(galaxy file), hval (hubble constant, i.e. 0.67 for planck cosmology), boxside (e.g. 100 if
    box is 100 Mpc/h on a side)
    galshort format is above, in testlum comments
    """
    coltype=("all","red","blue")
    collfull=("all galaxies","red galaxies","blue galaxies")
    collist=('k','r','b')
    cshapelist=('ks','rs','bs')
    for i in range(3):
      f,ax = plt.subplots(1,1)
      ax.set_xlim(-17.75,-24.00)
      ax.set_ylim(1.e-6,0.01)
      ax.set_yscale("log")        
      bin_centers,ngalact,Bmid,phi,phierr,jjz=testlum(zcen,i,fname,hval,boxside,ramin,ramax,decmin,decmax,delz)
      ax.step(bin_centers, ngalact,collist[i],label=r'simulation $z_{sim}$= %3.2f'%(zcen))

      zrange = N.array([0.2,0.3,0.4,0.5,0.65,0.8,1.0])      
      ax.plot(Bmid,phi,cshapelist[i],label=r'BOOTES %3.2f<z<%3.2f'%(zrange[jjz],zrange[jjz+1]))
      ax.text(9.2,2.e-4,collfull[i],color=collist[i])

      ax.legend(loc=3)
      ax.set_yscale("log")
      ax.set_xlabel("$M_B-5 log $h_{70}$")
      ax.set_ylabel("$\Phi$[$h_{70}^3$/Mpc${}^3$]")
      plt.tight_layout()
      plt.savefig("Bband_z%d_%s.pdf"%(zcen*101,coltype[i]))
      plt.close("all")





