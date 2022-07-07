from pyTNS import ffi, lib
import numpy as np

class TNS:
    def __init__(self, input, smin, smax, nbins):
        self.input = input.encode("utf-8")
        self.smin = smin
        self.smax = smax
        self.nbins = nbins
        self.load()
    
    def load(self):
        lib.init_prediction(self.input)

    def predict(self,f,b1,b2,sv,apar,aper):
        d = lib.get_xil_LL(self.smin, self.smax, self.nbins, f, b1, b2, sv, apar, aper)
        return np.array(ffi.unpack(d,3*self.nbins))


pkfile = '/home/sdelatorre/EUCLID/IST/RSD/nonlinearpk/pk_nonlinear_0.9.dat'
model = TNS(pkfile, 0., 200., 40)
p = model.predict(0.5, 1.1, -0.5, 5, 1.0, 1.0)
print(p)
