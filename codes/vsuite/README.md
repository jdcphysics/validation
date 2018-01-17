Please do email me if you have questions/corrections/suggestions.

jcohn@berkeley.edu

More documentation is in  http://arxiv.org/abs/1609.03956 (by J Cohn).  Please reference that paper <b>and datasets below</b>, at %%%%,  if you use this work.
   
how to get plots:

runsuite in valid_suite.py 

generates 7 plots using runsuite(): below.

1. 4 are stellar mass functions:
all, quiescent, star forming, and all on one page, compared to several observations described below

(quiescent/star forming division at log sfr = -0.49 + (0.65+slopeval) (logM* - 10) +1.07 *(z-0.1)+shiftval, for slopeval=shiftval=0,
Moustakas et al eq 2), although many of the observational papers listed use UVJ.

2. 1 is stellar mass-sfr diagram [can be compared with e.g., Moustakas et al 2013, but not overplotted with it]

3. 1 is ssfr in 4 stellar mass bins* (no cut on ra, dec for this)

4. 1 is stellar mass to halo mass diagram 

compared to Behroozi, Wechsler, Conroy 2013 (fit using Mvir)

Moster,Naab, White 2013 (fit using M200)

<b>To test use of the code</b><br>
in the directory vsuite/example are the outputs from running

runsuite(1/0.9947-1.,"inputstats_bolshoi_P_0.9947.dat",0.678,0.31,0.15,-0.8,250,"example")

with the file "inputstats_bolshoi_P_0.9947.dat", found at <a href="http://mwhite.berkeley.edu/vsuite_data">http://mwhite.berkeley.edu/vsuite_data/</a> .  
(You must uncompress it first, it is 22M compressed.)

If you run 
runsuite(1/0.9947-1.,"inputstats_bolshoi_P_0.9947.dat",0.678,0.31,0.15,-0.8,250,"tests")
you can compare your outputs (just name it something besides "example", here I've chosen "tests").

To get just stellar mass functions at a given redshift zchoose from observations (all, quiescent and star forming)

runsuite(zchoose,"inputstats_short.dat",0.678,0.31,0,0,250,"justobs") 
will produce 4 files called smf4sims[stuff]justobs.pdf
--giving observational stellar mass function data at redshift zchoose for
all, quiescent and starforming
with the hubble constant 0.678 and omega_m 0.31
(other figures will be produced but won't have much useful stuff with only 4 galaxies)
stellar mass to halo mass reads off redshift from file, so will have to change column 4 (5th column) to get M*(Mh) for other redshifts
besides the default redshift 0.


Again, if you use this program, please reference the papers and people who measured
all of these data!!
They are listed below at "%%%"


<b>USAGE:</b><br>
runsuite(zcen, "inputfile.dat",hval,omm,slopeval,shiftval, boxside,runname,delz,ramin,ramax,decmin,decmax):

zcen is central redshift

fname = "inputfile.dat" described below, can call it something else
 if you want.  ascii text.
 
hval  = hubble constant

omm = omega_matter (e.g. 0.31)

slopeval = in sfr-M* bimodal diagram, **change in** slope of line to separate star-forming and quiescent from PRIMUS - for simplicity can set to 0.

shiftval = change in shift of line between star forming and quiescent from PRIMUS -for simplicity can set to 0.

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
provided at [www.peterbehroozi.com/data.html,observational-data.tar.gz](www.peterbehroozi.com/data.html,observational-data.tar.gz)

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
[http://cosmos.phy.tufts.edu/~danilo/MuzzinEtal2013/Muzzin_et_al._(2013).html](http://cosmos.phy.tufts.edu/~danilo/MuzzinEtal2013/Muzzin_et_al._(2013).html)
paper is Muzzin et al, below, [http://arxiv.org/abs/1303.4409](http://arxiv.org/abs/1303.4409)

6.  henriques_all.dat, henriques_quiescent.dat,henriques_starforming.dat:
 points from figs 2 and 7 of Henriques et al
[http://arxiv.org/abs/1410.0365,](http://arxiv.org/abs/1410.0365) data provided at
 [http://galformod.mpa-garching.mpg.de/public/LGalaxies/figures_and_data.php](http://galformod.mpa-garching.mpg.de/public/LGalaxies/figures_and_data.php)

7. viper_sch_*.dat from Moutard et al,    Moutard et al, 2016
    The VIPERS Multi-Lambda Survey. II
Diving with massive galaxies in 22 square degrees since z = 1.5
    [http://arxiv.org/abs/1602.05917 v3](http://arxiv.org/abs/1602.05917)
    table 2.
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Full references for papers:

     Behroozi, P., Wechsler, R.H., Conroy, C.,
     The Average Star Formation Histories of Galaxies in Dark Matter Halos from z = 0-8,
     [http://arxiv.org/abs/1207.6105](http://arxiv.org/abs/1207.6105)
     2013,ApJ,770,57

     Behroozi,P., Wechsler,R.H., Conroy, C.,
     On the Lack of Evolution in Galaxy Star Formation Efficiency,
     [http://arxiv.org/abs/1209.3013](http://arxiv.org/abs/1209.3013)
     2013 ApJ, 762, L31

     Henriques, Bruno M. B.; White, Simon D. M.; Thomas, Peter A.; Angulo, Raul; Guo, Qi; Lemson, Gerard; Springel, Volker; Overzier, Roderik,
     Galaxy formation in the Planck cosmology - I. Matching the observed evolution of star formation rates, colours and stellar masses
     [http://arxiv.org/abs/1410.0365](http://arxiv.org/abs/1410.0365)
     2015, MNRAS,451,2663 
     data tables: [http://galformod.mpa-garching.mpg.de/public/LGalaxies/figures_and_data.php](http://galformod.mpa-garching.mpg.de/public/LGalaxies/figures_and_data.php)

    Moster, Naab & White, 
    Galactic star formation and accretion histories from matching galaxies to dark matter haloes
    [http://arxiv.org/abs/1205.5807] (http://arxiv.org/abs/1205.5807)
    MNRAS, 2013, 428, 3121
     
    Moustakas, John, et al,
    PRIMUS: Constraints on Star Formation Quenching and Galaxy Merging, and the Evolution of the Stellar Mass Function from z = 0-1
    [http://arxiv.org/abs/1301.1688](http://arxiv.org/abs/1301.1688)
    ApJ, 2013, 767, 50

    Moutard et al, 2016
    The VIPERS Multi-Lambda Survey. II
    Diving with massive galaxies in 22 square degrees since z = 1.5
    [http://arxiv.org/abs/1602.05917 v3](http://arxiv.org/abs/1602.05917) 
    
    Muzzin, A., Marchesini, D., Stefanon, M., Franx, M., McCracken, H.J., Milvang-Jensen, B., Dunlop, J.S., Fynbo,J.P.U.,  Brammer, G., Labbe, I., van Dokkum, P., 
    The Evolution of the Stellar Mass Functions of Star-Forming and Quiescent Galaxies to z = 4 from the COSMOS/UltraVISTA Survey
    [http://arxiv.org/abs/1303.4409](http://arxiv.org/abs/1303.4409)
    2013, ApJ, 777, 18

    Tomczak et al, 2014,
    Galaxy Stellar Mass Functions from ZFOURGE/CANDELS: An Excess of Low-mass Galaxies since z = 2 and the Rapid Buildup of Quiescent Galaxies
    [http://arxiv.org/abs/1309.5972](http://arxiv.org/abs/1309.5972)
    2014,ApJ,783,85

    A companion reference to Tomczak et al is the description of how the data/catalogs were put
    together for the survey, to appear in Straatman et al, submitted.
