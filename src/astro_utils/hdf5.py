import numpy as np
import h5py 
from astropy import constants as const
from astropy import units as u
from .units import _read_units_from_file
kpc = const.kpc.cgs.value
Msun = const.M_sun.cgs.value
m_p = const.m_p.cgs.value
Myr = ((1.*u.Myr).to(u.s).value)
G = const.G.cgs.value
k_B = const.k_B.to("(cm2 g) / (s2 K)").value

class SnapQuantities:
    def __init__(self,f):
        self.U = _read_units_from_file(f)
        self.boxsize = f["Parameters"].attrs["BoxSize"]
        self.pos = np.array(f["PartType0/Coordinates"])*self.U["UnitLength"]
        
        self.xctr = (np.array(f["PartType0/Coordinates"][:,0])-self.boxsize/2)*self.U["UnitLength"]
        self.yctr = (np.array(f["PartType0/Coordinates"][:,1])-self.boxsize/2)*self.U["UnitLength"]
        self.zctr = (np.array(f["PartType0/Coordinates"][:,2])-self.boxsize/2)*self.U["UnitLength"]
        
        self.radxy = np.sqrt(self.xctr**2+self.yctr**2)
        self.mass = np.array(f["PartType0/Masses"])*self.U["UnitMass"]
        self.dens = np.array(f['PartType0/Density'])*self.U["UnitDensity"]
        self.volume = self.mass/self.dens

        self.eth_spec = np.array(f["PartType0/InternalEnergy"])*self.U["UnitEspecific"]
        self.chem = np.array(f["PartType0/ChemicalAbundances"])

        if "nu_100TeV" in f["PartType0"].keys():
            self.nu_100TeV = np.array(f["PartType0/nu_100TeV"])*self.volume
        
        if "CosmicRaySpecificEnergy" in f["PartType0"].keys():
            self.ecr = np.array(f["PartType0/CosmicRaySpecificEnergy"]) * self.U["UnitEspecific"] * self.dens
    
    def temp(self,ABHE=0.1):
        a = (5/3-1)*self.eth_spec/k_B
        b = (1+4*ABHE)*m_p/(1+ABHE+self.chem[:,1]-self.chem[:,0])
        return a*b

