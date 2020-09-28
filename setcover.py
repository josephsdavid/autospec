import numpy as np
import tqdm
import math
import matplotlib.pyplot as plt
from collections import namedtuple, defaultdict
from functools import partial
import itertools as it
from multiprocessing import Pool
import molecules as m
from molecules import SpectralQuery
import spectraldata as nasa

def convert_intensities(molecules, temp):
    def intensity_helper(mol, temp):
        mul = mol['intensity']/2.99792458e18
        planck = 6.626069e-27
        bolt = 1.3806505e-16
        exponent = planck*2.99792458e10
        exponent *= mol['elo'][0]
        exponent /= bolt
        exponent *= (1/300 - 1/temp)
        mid = np.exp(exponent).astype(np.float128)
        mid *= (300 / temp)**(mol['deg']/2 + 1)
        return mul*mid
    out = {k: x._asdict() for k, x in molecules.items()}
    for k, v in out.items():
        out[k]['intensity'] = intensity_helper(v, temp)
        out[k]['data'][:,-1] = out[k]['intensity']
    return {k:SpectralQuery(**v) for k, v in out.items()}


SpikeHelper = namedtuple("SpikeHelper", ("frequency_low", "frequency_high", "intensity", "frequency"))
MolHelper = namedtuple("MolHelper", ("frequency", "uncertainty"))
def reverse_dict(d):
    reversed_dict = defaultdict(list)
    for key, values in d.items():
        for value in values:
            reversed_dict[value].append(key)
    return reversed_dict

def frequency_cover(molecule, spike):
    return ((molecule.frequency <= spike.frequency_high + molecule.uncertainty) &
        (molecule.frequency >= spike.frequency_low - molecule.uncertainty))

def broadcast(l, a):
    return [np.array([i for _ in range(a.shape[0])]).T for i in l]


class SetCovering(object):
    def __init__(self, spikes, moles, method=frequency_cover, scoring = "intensity_score", temperature=None):
        self.sets, self.spike_dict = set_generation(spikes, moles, method = method, scoring=scoring, temperature = temperature)
        self.spikes = spikes
        self.method = method
        self.scoring = scoring
        self.molecule_dict = moles
        self.associations = [spike_explanations(x, self.spike_dict) for x in self.sets]

    def __getitem__(self, index):
        return get_mol_set(self.sets[index])

    def __len__(self):
        return len(self.sets)

    def minimal_set(self, unique = True):
        return minimal_set(self.sets, unique)

    def likeliest_sets(self):
        intensity_results = []
        for ps in self.sets:
            ps = set(it.chain(*ps))
            scores = [x[-1] for x in ps if not np.isnan(x[-1])]
            average_score = sum(scores) / len(scores)
            intensity_results.append([ps, average_score])
        return sorted(intensity_results, key=lambda x: x[-1], reverse=True)

    def likeliest_molecules(self):
        return sorted( set(it.chain(*it.chain(*self.sets))), key = lambda x: x[-1], reverse = True)



def set_generation(spikes, moles, method=frequency_cover, scoring = "intensity_score", temperature = None):
    if list(moles.values())[0].units != "GHz":
        ValueError("bad units! Please convert the molecules to GHz")
    spike_dict = {f: [] for f in spikes.frequency}
    spacing = (spikes.data[1:, 0] - spikes.data[:-1, 0]).mean() / 2
    spike_list = [spikes.frequency - spacing, spikes.frequency + spacing, spikes.intensity, spikes.frequency]
    out_mols = {}
    if temperature is not None:
        print("converting intensities")
        moles = convert_intensities(moles, temperature)

    for mol in tqdm.tqdm(moles.values(), desc = "Scoring molecules..." ):
        #spike_broadcast = [np.array([s for _ in range(mol.frequency.shape[0])]).T for s in spike_list]
        spike_broadcast = broadcast(spike_list, mol.frequency)
        spike_matches = method(mol, SpikeHelper(*spike_broadcast))
        for i in range(spike_matches.sum(1).shape[0]):
            if spike_matches.sum(1)[i]: # santerre trick
                # score now, make this to classes
                winner = mol.molecule
                spike_intensities = mol.intensity[np.where(spike_matches.sum(0) != 0)].astype(np.float128)
                I = np.power(np.array([10.]).astype(np.float128), mol.intensity.astype(np.float128))
                I = np.nansum(I)
                score = ((np.power(np.array([10.]).astype(np.float128), spike_intensities))/I)
                score = np.nansum(score)
                if np.isnan(score):
                    score=-1
                spike_dict[list(spike_dict.keys())[i]].append((mol.molecule, score))
    spike_dict = {k: v if len(v) > 0 else ["unknown"] for k, v in spike_dict.items()}
    result = {}
    for s in tqdm.tqdm(spike_dict.keys(), desc = "Generating combinations..."):
        result.setdefault(s, []).append([x for z in [list(it.combinations(spike_dict[s],y+1)) for y in range(len(spike_dict[s]))] for x in z])
    print("Restructuring data...")
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


def intensity_score(single_cover, spike_dict, molecule_dict, spikes, temp_adjust = None):
    # invert the spike_dict, so we have molecule: spikes
    inverse_spikes = reverse_dict(spike_dict)
    molecules = set(it.chain(single_cover)) if type(single_cover) is set else set(it.chain(*single_cover))
    res = []
    if temp_adjust:
        molecule_dict = convert_intensities(molecule_dict, temp_adjust)
    for m in molecules:
        if m not in inverse_spikes.keys():
            res.append((m, np.nan))
        else:
            spike_intensities = spikes.intensity[np.where(spikes.frequency == np.array(inverse_spikes[m]))]
            i = np.array(10 ** molecule_dict[m].intensity).sum()
            score = (10**spike_intensities/I).sum()
            res.append((m,score))
    return res

def batch_process(files, plot=False, temperature = None):
    cover_list = []
    for f in files:
        sf = nasa.SpectralFile(f)
        sd = nasa.read(sf)
        spikes = nasa.identify_spikes(sd)
        if plot:
            nasa.plot_spikes(spikes)
        molecules = m.get_molecules_from_spikes(spikes)
        cover_list.append(SetCovering(spikes, molecules, temperature = temperature))
    return cover_list

def summation_score(cover_list):
    mols = [x.likeliest_molecules() for x in cover_list]
    mols = sum(mols, [])
    result = {key: 0. for key in set([x[0] for x in mols])}
    for mol in tqdm.tqdm(mols, desc="Combining results..."):
        result[mol[0]] += mol[-1]
    out = [(key, value) for key, value in result.items()]
    return sorted(out, key = lambda x: x[-1], reverse = True)
