# TNS_ToolBox

AUTHOR: Sylvain de la Torre, 2022        

This C code implements the TNS model (Taruya, Nishimishi & Saito 2010) with
renormalised 1-loop bias model as used in Bautista, Paviot et al. 2020. It
offers the possibility of providing real-space P_dd, P_dt, P_tt as input or
using Halofit with Bel et al. 2019 prescriptions for P_dt and P_tt. The final
products are the monopole, quadrupole, and hexadecapole moments of the
redshift-space two-point correlation function.

*** Requirements ***

It needs gsl, fftw3, openmp libraries to be installed.

*** Compilation *** 

Compilation is done with
> make all

Individual codes can be compile as:
> make rp
> make tns
> make libtns

RealPower and TNS executable can be ran from the command line.
'RealPower' computes the ingredients of the TNS model. It has several options.
'TNS' computes the prediction of TNS model.
'libTNS.so' is a library that can be used externally in other codes or wrappers.  

** Contact
For any question, please send an email to sylvain.delatorre@lam.fr.
