import numpy as np
import h5py 
from astropy import constants as const
from astropy import units as u
kpc = const.kpc.cgs.value
Msun = const.M_sun.cgs.value
m_p = const.m_p.cgs.value
Myr = ((1.*u.Myr).to(u.s).value)
G = const.G.cgs.value
k_B = const.k_B.to("(cm2 g) / (s2 K)").value

class SnapQuantities:
    def __init__(self,f):
        self.boxsize = f["Parameters"].attrs["BoxSize"]
        self.pos = np.array(f["PartType0/Coordinates"])

        self.eth_spec = np.array(f["PartType0/InternalEnergy"])
        self.chem = np.array(f["PartType0/ChemicalAbundances"])
    
    def temp(self,ABHE=0.1):
        a = (5/3-1)*self.eth_spec/k_B
        b = (1+4*ABHE)*m_p/(1+ABHE+self.chem[:,1]-self.chem[:,0])
        return a.to_value()*b.to_value()
