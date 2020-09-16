import spectraldata as nasa
from pprint import pprint
import numpy as np
import molecules as m
from itertools import chain
import setcover as sc
import os
import joblib

# molecules = m.construct_mega_query

n_spikes = 0
m_dict = dict()

PATH = "Titan/Titan/"
files = [f"{PATH}{f}"  for f in os.listdir(PATH)]
sf = nasa.SpectralFile(files, temp = 94)
sd = nasa.read(sf)
#for f in os.listdir("Titan/Titan/"):
#    sf = nasa.SpectralFile(f"Titan/Titan/{f}")
#    sd = nasa.read(sf)
#    spikes = nasa.identify_spikes(sd)
#    nasa.plot_spikes(spikes)
#    m_dict = dict(m_dict, **m.get_molecules_from_spikes(spikes))
###    # use these spikes + splatalogue to query JPL
###    # given that that works, figure out optimal choice
## print(len(m_dict))
## print(n_spikes)
#
#
## import pdb; pdb.set_trace()  # XXX BREAKPOINT
#
#sf = nasa.SpectralFile(f"Titan/Titan/Win0.clean1.contsub_Jy.rest.scom.c.txt")
#sd = nasa.read(sf)
# ghz_data = nasa.convert_units(sd, to='GHz')
spikes = nasa.identify_spikes(sd)
#nasa.plot_spikes(spikes)
# import pdb; pdb.set_trace()  # XXX BREAKPOINT
# nasa.plot_spikes(spikes)
if "all_mols.sav" not in os.listdir():
    molecule_dict = m.get_molecules_from_spikes(spikes)
    m.save_molecules(molecule_dict, "all_mols.sav")
else:
    molecule_dict = m.load_molecules("all_mols.sav")

# import pdb; pdb.set_trace()  # XXX BREAKPOINT
#molecule_dict = joblib.load("mol_dict.pkl")
#molecule_dict = {k: m.SpectralQuery(**v) for k, v in molecule_dict.items()}

#molecule_dict = sc.convert_intensities(molecule_dict, sd.temp)
# move unit conversion!

covers = sc.SetCovering(spikes, molecule_dict)

pprint(covers.likeliest_molecules())

import pdb; pdb.set_trace()  # XXX BREAKPOINT

pprint(covers.likeliest_sets())
import pdb; pdb.set_trace()  # XXX BREAKPOINT

