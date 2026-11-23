import h5py
import numpy as np
from astropy import constants as const
from astropy import units as u
pc = const.pc.cgs.value
kpc = const.kpc.cgs.value
Msol = const.M_sun.cgs.value
Myr = ((1.*u.Myr).to(u.s).value)

def _read_units_from_file(fname):
    print("reading file: ", fname)
    f = h5py.File(fname, "r")
    # get momentum boundaries from file
    try:
        UnitMass = f[u'Parameters'].attrs.get("UnitMass_in_g")
        UnitLen  = f[u'Parameters'].attrs.get("UnitLength_in_cm")
        UnitVelo = f[u'Parameters'].attrs.get("UnitVelocity_in_cm_per_s")
        #BoxSize  = f[u'Parameters'].attrs.get("BoxSize")
    except:
        print("Unit parameters not found!!!")
    f.close()
    return _compute_arepo_units(UnitMass, UnitLen, UnitVelo)

def _compute_arepo_units(UnitMass, UnitLength, UnitVelocity, print_units=False):
    UnitTime = UnitLength/UnitVelocity
    d = { "UnitLength"   : UnitLength,
          "UnitMass"     : UnitMass,
          "UnitVelocity" : UnitVelocity,
          "UnitTime"     : UnitLength/UnitVelocity,
          "UnitDensity"  : UnitMass/(UnitLength**3),
          "UnitColdens"  : UnitMass/(UnitLength**2),
          "UnitEnergy"   : UnitMass*UnitVelocity*UnitVelocity,
          "UnitEspecific": UnitVelocity*UnitVelocity,
          "UnitEdensity" : UnitMass / UnitLength / UnitTime**2,
          "UnitPressure" : UnitMass / UnitLength / UnitTime**2,
          "UnitBfield"   : np.sqrt(UnitMass) / np.sqrt(UnitLength) / UnitTime,
          "UnitGrav"     : UnitLength**3 / (UnitMass * UnitTime**2),
          "UnitPotential": UnitLength**2 / (UnitTime**2),
          # non cgs
          "UnitColdens_Msol_pc2" : UnitMass/(UnitLength**2) / Msol * pc*pc

    }
    if print_units:
        print()
        print("Arepo units")
        # find longest key name
        maxlen = len(max(d, key=len))
        for key in d:
            print(key.ljust(maxlen), ":", d[key])
        print()
    return d
