import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from functools import partial
import itertools as it
from multiprocessing import Pool
from molecules import SpectralQuery

SpikeHelper = namedtuple("SpikeHelper", ("frequency_low", "frequency_high", "intensity", "frequency"))
MolHelper = namedtuple("MolHelper", ("frequency", "uncertainty"))

def frequency_cover(molecule, spike):
    return ((molecule.frequency <= spike.frequency_high + molecule.uncertainty) &
        (molecule.frequency >= spike.frequency_low - molecule.uncertainty))

def broadcast(l, a):
    return [np.array([i for _ in range(a.shape[0])]).T for i in l]

def set_generation(spikes, moles, method=frequency_cover):
    if list(moles.values())[0].units != "GHz":
        ValueError("bad units! Please convert the molecules to GHz")
    spike_dict = {f: [] for f in spikes.frequency}
    spacing = (spikes.data[1:, 0] - spikes.data[:-1, 0]).mean() / 2
    spike_list = [spikes.frequency - spacing, spikes.frequency + spacing, spikes.intensity, spikes.frequency]
    out_mols = {}
    for mol in moles.values():
        #spike_broadcast = [np.array([s for _ in range(mol.frequency.shape[0])]).T for s in spike_list]
        spike_broadcast = broadcast(spike_list, mol.frequency)
        spike_matches = method(mol, SpikeHelper(*spike_broadcast)).sum(1)
        for i in range(spike_matches.shape[0]):
            if spike_matches[i]: # santerre trick
                spike_dict[list(spike_dict.keys())[i]].append(mol.molecule)
    spike_dict = {k: v if len(v) > 0 else ["unknown"] for k, v in spike_dict.items()}
    result = {}
    for s in spike_dict.keys():
        result.setdefault(s, []).append([x for z in [list(it.combinations(spike_dict[s],y+1)) for y in range(len(spike_dict[s]))] for x in z])
    return list(map(list, zip(*list(map(list,list(it.product(*result.values()))[0]))))), spike_dict

def get_mol_set(tl):
    return set(sum(it.chain(tl), ()))

def spike_explanations(mol_list, spike_dict):
    return {list(spike_dict.keys())[i]: list(it.chain(mol_list[i])) for i in range(len(mol_list))}

def minimal_set(possible_sets, unique = True):
    set_complexity = [len(set(sum(it.chain(x), ()))) for x in possible_sets]
    min_index = set_complexity.index(min(set_complexity))
    return get_mol_set(possible_sets[min_index]) if unique else possible_sets[min_index]

def minimal_set_iterator(possible_sets, unique = True):
    set_complexity = [len(set(sum(it.chain(x), ()))) for x in possible_sets]
    set_indexes = [set_complexity.index(x) for x in sorted(set_complexity)]
    list_ = [possible_sets[i] for i in set_indexes]
    for item in list_:
        yield get_mol_set(item) if unique else item

