Please do email me if you have questions!
jcohn@berkeley.edu
   
runsuite in valid_suite.py 
generates 7 plots using runsuite(): below.
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
slopeval = in sfr-M* bimodal diagram, **change in** slope of line to separate star-forming and quiescent from PRIMUS
shiftval = change in shift of line between star forming and quiescent from PRIMUS

PRIMUS starforming and quiescent split by: 
log SFR = log sfrmin -0.49 + (0.65+slopeval) (logM* - 10) +1.07 *(z-0.1) + shiftval 

boxside = in Mpc/h for fixed time, any negative number if light cone
runname = string, such as "run0"
if lightcone, delz,ramin,ramax,decmin,decmax listed next.
if fixed time these arguments (delz, ramin, ramax, decmin, decmax)
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
http://arxiv.org/abs/1410.0365
 2015, MNRAS,451,2663 
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
