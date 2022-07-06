# TNS_ToolBox

 ========================================
 AUTHOR: Sylvain de la Torre         2022        
 ========================================

This C code implements the TNS model (Taruya, Nishimishi & Saito 2010) with
renormalised 1-loop bias model as used in Bautista, Paviot et al. 2020. It
offers the possibility of providing real-space P_dd, P_dt, P_tt as input or
using Halofit with Bel et al. 2019 prescriptions for P_dt and P_tt. The final
products are the monopole, quadrupole, and hexadecapole moments of the
redshift-space two-point correlation function.

*** Requirements ***

It needs gsl and fftw3 libraries to be installed.

*** Compilation *** 

In the code folder, there are bash scripts to compile the codes:
   - comp_RP: it compiles the RealPower code to generate the inputs power spectra and TNS corrections for the TNS model
   - comp_TNS: it compiles the TNS code that computes the TNS model
   - comp_TNSlib: it compiles the TNS library to be used in python wrappers for instance to run TNS prediction
	                     	       
RealPower and TNS codes can be ran from the command line.

** Contact
For any question, please send an email to sylvain.delatorre@lam.fr.