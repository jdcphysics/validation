 validation
 
 
 in codes/StellarMass/
 
   code to compare simulated data stellar mass function with observed
   stellar mass function from PRIMUS, in several redshift bins from 0.2-1.0
   makes plot with both on it
  
   2 files plus your simulated data are needed:
   
     PRIMUS_stellar.py 
     
         histograms an input simulated data set (M*, sfr) with
         Moustakas et al 2013 PRIMUS data, redshifts 0.2-1.0.
         Moustakas et al 2013 SDSS-GALEX data, redshifts 0.01-0.2
         
     Mous_13_table4.txt 
     
         is 2013 PRIMUS data, in same directory 
         (typed in table 4 of Moustakas et al 2013, arXiv:1301.1688v1)

 in codes/Bband/
 
   code to compare simulated data B band luminosity function with observed
   stellar mass function from BOOTS, in several redshift bins from 0.2-1.1
   makes plots of both together.
   
   **only tested on 'fake' simulated data so far**
   
   2 files plus your simulated data are needed:
   
     BOOTES_lf.py 
     
         histograms an input simulated data set (M_B, M_U) with
         Beare et al 2015 BOOTES data, redshifts 0.2-1.1.
         
     Beare_tab8910.txt
     
         is 2015 BOOTES data, in same directory 
         (typed in table 8,9,10 of Beare et al 2015, arXiv:1511.01580)
