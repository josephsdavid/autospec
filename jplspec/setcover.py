import pandas as pd
from itertools import cycle
import plotly.express as px
import numpy as np
import math
import matplotlib.pyplot as plt
from collections import namedtuple, defaultdict
from functools import partial
import itertools as it
from multiprocessing import Pool
import plotly.graph_objects as go
from .molecules import SpectralQuery
import tqdm


def convert_intensities(molecules, temp):
    def intensity_helper(mol, temp):
        mul = mol["intensity"] / 2.99792458e18
        planck = 6.626069e-27
        bolt = 1.3806505e-16
        exponent = planck * 2.99792458e10
        exponent *= mol["elo"][0]
        exponent /= bolt
        exponent *= 1 / 300 - 1 / temp
        mid = np.exp(exponent).astype(np.float128)
        mid *= (300 / temp) ** (mol["deg"] / 2 + 1)
        return mul * mid

    out = {k: x._asdict() for k, x in molecules.items()}
    for k, v in out.items():
        out[k]["intensity"] = intensity_helper(v, temp)
        out[k]["data"][:, -1] = out[k]["intensity"]
    return {k: SpectralQuery(**v) for k, v in out.items()}


SpikeHelper = namedtuple(
    "SpikeHelper", ("frequency_low", "frequency_high", "intensity", "frequency")
)
MolHelper = namedtuple("MolHelper", ("frequency", "uncertainty"))


def reverse_dict(d):
    reversed_dict = defaultdict(list)
    for key, values in d.items():
        for value in values:
            reversed_dict[value].append(key)
    return reversed_dict


def frequency_cover(molecule, spike):
    return (molecule.frequency <= spike.frequency_high + molecule.uncertainty) & (
        molecule.frequency >= spike.frequency_low - molecule.uncertainty
    )

def entropy_score(molecule, spike):
    mol_intense = (molecule.intensity - molecule.intensity.min()) / (molecule.intensity.max() - molecule.intensity.min())
    spike_intense = (spike.intensity - spike.intensity.min()) / (spike.intensity.max() - spike.intensity.min())
    mol_ent = (10**mol_intense * np.divide(mol_intense, spike_intense,  where=spike_intense!=0)).sum()
    spike_ent = (10**spike_intense * np.divide(spike_intense, mol_intense,  where=mol_intense!=0)).sum()
    import pdb; pdb.set_trace()
    return mol_ent + spike_ent


def broadcast(l, a):
    return [np.array([i for _ in range(a.shape[0])]).T for i in l]


class SetCovering(object):
    def __init__(
        self,
        spikes,
        moles,
        method=frequency_cover,
        scoring="intensity_score",
        temperature=None,
    ):
        self.sets, self.spike_dict, self.plot_dict = set_generation(
            spikes, moles, method=method, scoring=scoring, temperature=temperature
        )
        self.spikes = spikes
        self.method = method
        self.scoring = scoring
        self.molecule_dict = moles
        self.associations = [spike_explanations(x, self.spike_dict) for x in self.sets]

    def __getitem__(self, index):
        return get_mol_set(self.sets[index])

    def __len__(self):
        return len(self.sets)

    def minimal_set(self, unique=True):
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
        return sorted(
            set(it.chain(*it.chain(*self.sets))), key=lambda x: x[-1], reverse=True
        )

    def visualize(self, lines=False):
        colorgen = cycle(px.colors.qualitative.Plotly)
        colors = {mol: next(colorgen) for mol in set(self.plot_dict.keys())}
        plot_dict = dict(
            sorted(
                self.plot_dict.items(),
                key=lambda x: x[1]["match_intensity"][0],
                reverse=False,
            )
        )
        fig = go.Figure(
            data=go.Scatter(
                x=self.spikes.data[:, 0],
                y=self.spikes.data[:, -1],
                legendgroup="Spectrum",
                name="Spectrum",
                opacity=0.2,
                mode="lines",
            )
        )
        if lines:
            for molecule, info in plot_dict.items():
                freq = info["match_freq"].tolist()
                intense = (info["match_intensity"] * lines).tolist()
            #    fig.add_trace(
            #        go.Scatter(
            #            x=[[f, f] for f in freq],
            #            y=[[-i,i] for i in intense] ,
            #            mode="lines+text",
            #            legendgroup=f"{molecule} spectrum",
            ##            name=molecule,
            #            textposition="bottom center",
            #            visible="legendonly",
            #            line=dict(color=colors[molecule])
            #        )
            #    )
                for f, i in zip(freq, intense):
                    fig.add_trace(
                        go.Scatter(
                            x=[f, f],
                            y=[-i, i] ,
                            mode="lines+text",
                            legendgroup=f"{molecule} spectrum",
                            name=molecule,
                            textposition="bottom center",
                            visible="legendonly",
                            line=dict(color=colors[molecule])
                        )
                    )
        else:
            mean = self.spikes.data[:, -1].mean()
            minim = self.spikes.data[:, -1].min()
            maxim = self.spikes.data[:, -1].max()
            for molecule, info in plot_dict.items():
                fig.add_trace(
                    go.Scatter(
                        x=list(info["match_freq"]),
                        y=list(info["match_intensity"] * (minim - maxim) ** 3),
                        mode="markers",
                        legendgroup=molecule,
                        name=f"{molecule}",
                        textposition="bottom center",
                        line=dict(color=colors[molecule])
                    )
                )
        return fig
    def visualize00(self, lines=False):


        plot_dict = dict(
            sorted(
                self.plot_dict.items(),
                key=lambda x: x[1]["match_intensity"][0],
                reverse=False,
            )
        )
        fig = go.Figure(
            data=go.Scatter(
                x=self.spikes.data[:, 0],
                y=self.spikes.data[:, -1],
                legendgroup="Spectrum",
                name="Spectrum",
                opacity=0.2,
                mode="lines",
            )
        )
        if lines:

            for molecule, info in plot_dict.items():
                freq = info["match_freq"].tolist()
                intense = (info["match_intensity"] * lines).tolist()
                fig.add_trace(
                    go.Scatter(
                        x=[[f, f] for f in freq],
                        y=[[-i,i] for i in intense] ,
                        mode="lines+text",
                        legendgroup=f"{molecule} spectrum",
                        name=molecule,
                        textposition="bottom center",
                        visible="legendonly"
                    )
                )

                for f, i in zip(freq, intense):
                    fig.add_trace(
                        go.Scatter(
                            x=[f, f],
                            y=[0, i] ,
                            mode="lines+text",
                            legendgroup=f"{molecule} spectrum",
                            name=molecule,
                            textposition="bottom center",
                            visible="legendonly"
                        )
                    )

        else:
            mean = self.spikes.data[:, -1].mean()
            minim = self.spikes.data[:, -1].min()
            maxim = self.spikes.data[:, -1].max()

            for molecule, info in plot_dict.items():
                fig.add_trace(
                    go.Scatter(
                        x=list(info["match_freq"]),
                        y=list(info["match_intensity"] * (minim - maxim) ** 3),
                        mode="markers",
                        legendgroup=molecule,
                        name=f"{molecule}",
                        textposition="bottom center",
                    )
                )

        return fig




def set_generation(
    spikes, moles, method=frequency_cover, scoring="intensity_score", temperature=None
):
    spike_dict = {f: [] for f in spikes.frequency}
    spacing = (spikes.data[1:, 0] - spikes.data[:-1, 0]).mean() / 2
    spike_list = [
        spikes.frequency - spacing,
        spikes.frequency + spacing,
        spikes.intensity,
        spikes.frequency,
    ]
    out_mols = {}
    plot_dict = {}
    if temperature is not None:
        print("converting intensities")
        moles = convert_intensities(moles, temperature)

    for mol in tqdm.tqdm(moles.values(), desc="identifying matches..."):
        # spike_broadcast = [np.array([s for _ in range(mol.frequency.shape[0])]).T for s in spike_list]
        spike_broadcast = SpikeHelper(*broadcast(spike_list, mol.frequency))
        spike_matches = method(mol, spike_broadcast)
        for i in range(spike_matches.sum(1).shape[0]):
            if spike_matches.sum(1)[i]:  # santerre trick
                # score now, make this to classes
                winner = mol.molecule
                molidxs = np.where(
                    (mol.frequency >= spikes.data[:, 0].min() - spacing)
                    & (mol.frequency <= spikes.data[:, 0].max() + spacing)
                )
                matchids = np.where(
                    (spike_matches.sum(0) != 0)
                    & (mol.frequency >= spikes.data[:, 0].min() - spacing)
                    & (mol.frequency <= spikes.data[:, 0].max() + spacing)
                )
                spike_intensities = mol.intensity[
                    np.where(
                        (spike_matches.sum(0) != 0)
                        & (mol.frequency >= spikes.data[:, 0].min() - spacing)
                        & (mol.frequency <= spikes.data[:, 0].max() + spacing)
                    )
                ]
                plot_dict[winner] = {
                    "mol_freq": mol.frequency[molidxs],
                    "mol_intensity": mol.frequency[molidxs],
                    "match_freq": mol.frequency[matchids],
                    "match_intensity": mol.intensity[matchids],
                }
                if scoring == "intensity_score":
                    I = np.power(10.0, mol.intensity[molidxs[0]]).sum()
                    score = ((np.power(10.0, spike_intensities))).sum() / I
                elif scoring == "entropy_score":
                    score = entropy_score(mol, spike_broadcast)
                else:
                    print("error!!")

                if score > 1.0:
                    import pdb

                    pdb.set_trace()  # XXX BREAKPOINT
                if np.isnan(score):
                    score = -1.0
                spike_dict[list(spike_dict.keys())[i]].append((mol.molecule, score))
    spike_dict = {k: v if len(v) > 0 else ["unknown"] for k, v in spike_dict.items()}
    result = {}
    for s in tqdm.tqdm(spike_dict.keys(), desc="identifying sets..."):
        result.setdefault(s, []).append(
            [
                x
                for z in [
                    list(it.combinations(spike_dict[s], y + 1))
                    for y in range(len(spike_dict[s]))
                ]
                for x in z
            ]
        )
    print("Restructuring data...")
    return (
        list(map(list, zip(*list(map(list, list(it.product(*result.values()))[0]))))),
        spike_dict,
        plot_dict,
    )


def get_mol_set(tl):
    return set(sum(it.chain(tl), ()))


def spike_explanations(mol_list, spike_dict):
    return {
        list(spike_dict.keys())[i]: list(it.chain(mol_list[i]))
        for i in range(len(mol_list))
    }


def minimal_set(possible_sets, unique=True):
    set_complexity = [len(set(sum(it.chain(x), ()))) for x in possible_sets]
    min_index = set_complexity.index(min(set_complexity))
    return get_mol_set(possible_sets[min_index]) if unique else possible_sets[min_index]


def minimal_set_iterator(possible_sets, unique=True):
    set_complexity = [len(set(sum(it.chain(x), ()))) for x in possible_sets]
    set_indexes = [set_complexity.index(x) for x in sorted(set_complexity)]
    list_ = [possible_sets[i] for i in set_indexes]
    for item in list_:
        yield get_mol_set(item) if unique else item


def intensity_score(single_cover, spike_dict, molecule_dict, spikes, temp_adjust=None):
    # invert the spike_dict, so we have molecule: spikes
    inverse_spikes = reverse_dict(spike_dict)
    molecules = (
        set(it.chain(single_cover))
        if type(single_cover) is set
        else set(it.chain(*single_cover))
    )
    res = []
    if temp_adjust:
        molecule_dict = convert_intensities(molecule_dict, temp_adjust)
    for m in molecules:
        if m not in inverse_spikes.keys():
            res.append((m, np.nan))
        else:
            spike_intensities = spikes.intensity[
                np.where(spikes.frequency == np.array(inverse_spikes[m]))
            ]
            I = np.array(10 ** molecule_dict[m].intensity).sum()
            score = (10 ** spike_intensities / I).sum()
            res.append((m, score))
    return res



def summation_score(cover_list):
    mols = [x.likeliest_molecules() for x in cover_list]
    mols = sum(mols, [])
    result = {key: 0.0 for key in set([x[0] for x in mols])}
    for mol in mols:
        result[mol[0]] += mol[-1]
    out = [(key, value) for key, value in result.items()]
    return sorted(out, key=lambda x: x[-1], reverse=True)
