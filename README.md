# TNS_ToolBox

This C code implements the TNS model (Taruya, Nishimishi & Saito 2010) with
renormalised 1-loop bias model as used in Bautista, Paviot et al. 2020. It
offers the possibility of providing real-space P_dd, P_dt, P_tt as input or
using Halofit with Bel et al. 2019 prescriptions for P_dt and P_tt. The final
products are the monopole, quadrupole, and hexadecapole moments of the
redshift-space two-point correlation function. A python wrapper is also provided.


## Requirements

The `gsl`, `fftw3`, `openmp` libraries should to be installed. For the python 
wrapper, `cffi` and `numpy` modules should be installed.


## Compilation

Compilation is done with
> make all

Individual codes are compiled as:
> make rp

> make tns

> make libtns

The `RealPower` and `TNS` executables can be ran from the command line.
`RealPower` computes the ingredients of the TNS model (it has several options),
`TNS` computes the prediction of the TNS model for the monopole, quadrupole, 
and hexadecapole of the correlation function, and `libTNS.so` is a library that 
can be used externally in other C codes or wrappers.  


## Python wrapper

The model can be used as a python module calling the `libTNS` library. The TNS 
python module needs first to be built with `cffi`. For this:
> cd pymodule

> python build.py

The wrapper class is in `wrapperTNS.py` and can be tested by running:

> python wrapperTNS.py

For the python wrapper to work in other python scripts, `libTNS.so` and 
`pyTNS.cpython-*.so` files need to be put in the same directory.


## Contact

For any question, please send an email to sylvain.delatorre@lam.fr.
