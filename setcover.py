import numpy as np
from collections import namedtuple
from functools import partial
import itertools as it
from multiprocessing import Pool

SingleSpike = namedtuple("SingleSpike", ("frequency", "intensity"))

def frequency_cover(molecule, spike):
    # factor in measurement error to fully function
    # the idea is for this to be extensible, so if we want other things in our
    # output we can do it, such as a confidence score. This will allow us to
    # rank sets of molecules or whatever
    if (
        sum(
            (molecule.frequency <= spike.frequency + molecule.uncertainty)
            & (molecule.frequency >= spike.frequency - molecule.uncertainty)
        )
        > 0
    ):
        return (spike.frequency, molecule.molecule)

def set_covering(spikes, moles, method=frequency_cover):
    if list(moles.values())[0].units != "GHz":
        ValueError("bad units! Please convert the molecules to GHz")
    spike_dict = {f: [] for f in spikes.frequency}
    spike_array = np.array([spikes.frequency, spikes.intensity])
    for s in spike_array.T:
        spike = SingleSpike(*s)
        coverer = partial(method, spike=spike)
        # convert to pool when we convert to real classes
        out = map(coverer, moles.values())
        # do something about if itsempty
        spike_dict[spike.frequency] = [x for x in out if x]
    out = {}
    for s in spike_dict.keys():
        out.setdefault(s, []).append([x for z in [list(it.combinations(spike_dict[s],y+1)) for y in range(len(spike_dict[s]))] for x in z])
    # returns candidate sets, remove list(map(list, zip(*l))) to getcandidates
    # by spike
    # a good test is set(sum([[len(y) for y in x] for x in p],[]))
    # fix this
    return list(map(list, zip(*list(map(list,list(it.product(*out.values()))[0])))))

