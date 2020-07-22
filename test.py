import spectraldata as nasa
from pprint import pprint
import numpy as np
import molecules as m
import setcover as sc
import os
import joblib

# molecules = m.construct_mega_query

n_spikes = 0
m_dict = dict()
for f in os.listdir("Titan/Titan/"):
    sf = nasa.SpectralFile(f"Titan/Titan/{f}")
    sd = nasa.read(sf)
    spikes = nasa.identify_spikes(sd)
#    nasa.plot_spikes(spikes)
#    n_spikes += spikes.index.shape[0]
#    print(n_spikes)
#    m_dict = dict(m_dict, **m.get_molecules_from_spikes(spikes))
##    # use these spikes + splatalogue to query JPL
##    # given that that works, figure out optimal choice
#print(len(m_dict))
#print(n_spikes)



#import pdb; pdb.set_trace()  # XXX BREAKPOINT

sf = nasa.SpectralFile(f"Titan/Titan/Win0.clean1.contsub_Jy.rest.scom.c.txt")
sd = nasa.read(sf)
# ghz_data = nasa.convert_units(sd, to='GHz')
spikes = nasa.identify_spikes(sd)
#import pdb; pdb.set_trace()  # XXX BREAKPOINT
#nasa.plot_spikes(spikes)
#molecule_dict = m.get_molecules_from_spikes(spikes)
molecule_dict = joblib.load("mol_dict.pkl")
molecule_dict = {k: m.SpectralQuery(**v) for k, v in  molecule_dict.items()}
# move unit conversion!

import time
rn = time.time()
possible_sets, spike_dict = sc.set_generation(spikes, molecule_dict)
print(time.time() - rn)


minset = sc.minimal_set(possible_sets)
pprint(minset)
pprint(sc.spike_explanations(
    sc.minimal_set(possible_sets, unique = False), spike_dict))

sorted_sets = sc.minimal_set_iterator(possible_sets)
pprint(next(sorted_sets))

import pdb; pdb.set_trace()  # XXX BREAKPOINT

#set_dict = set_cover_runner(spikes, molecule_dict)
#
#from itertools import combinations, product
#big_guy = {}
#for s in set_dict.keys():
#    big_guy.setdefault(s, []).append([x for z in [list(combinations(set_dict[s],y+1)) for y in range(len(set_dict[s]))] for x in z])
## Spike Seen at postion
#
#l = list(map(list,list(product(*big_guy.values()))[0]))
#
#SingleSpike = namedtuple("SingleSpike", ("frequency", "intensity"))
#
#def frequency_cover(molecule, spike):
#    # factor in measurement error to fully function
#    print(molecule.molecule)
#    if (
#        sum(
#            (molecule.frequency <= spike.frequency + molecule.uncertainty)
#            & (molecule.frequency >= spike.frequency - molecule.uncertainty)
#        )
#        > 0
#    ):
#        return (spike.frequency, molecule.molecule)
#
#
#def set_cover_runner(spikes, moles, method=frequency_cover):
#    if list(moles.values())[0].units != "GHz":
#        ValueError("bad units! Please convert the molecules to GHz")
#    spike_dict = {f: [] for f in spikes.frequency}
#    spike_array = np.array([spikes.frequency, spikes.intensity])
#    for s in spike_array.T:
#        spike = SingleSpike(*s)
#        coverer = partial(method, spike=spike)
#        # convert to pool when we convert to real classes
#        out = map(coverer, moles.values())
#        spike_dict[spike.frequency] = [x for x in out if x]
#    return spike_dict
#
#def set_generation(spikes, moles, approach='all'):
#  if approach=='all':
#    possible_moles = []
#    for mole in moles.keys():
#      for entry in moles[mole]:
#        for spike in spikes:
#          if entry[0]==spike[0]:
#            possible_moles.append(mole)
#            # break appropriately
#
#    spike_dict = dict([(s[0],[]) for s in spike])
#    for mole in moles:
#      for entry in moles[mole]:
#        if entry[0] in spike_dict.keys():
#          spike_dict[entry[0]].append(mole)
#
#    big_guy = []
#    for s in spike_dict.keys():
#      big_guy.append([x for z in [list(combinations(spike_dict[s],y+1)) for y in range(len(spike_dict[s]))] for x in z])
#
#list(product(*big_guy))
#big_guy=[[(10,), (20,),(10,20)], [(1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]]
#
#big_guy = [[a],[b,c],[d,e]]
#a    b.    d
#     c     e
#     bc    de
#
#
#import pdb; pdb.set_trace()  # XXX BREAKPOINT
#
## dict with freq as keys
