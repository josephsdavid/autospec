import sys
import os
from itertools import cycle
import plotly.express as px
import numpy as np
from collections import namedtuple, defaultdict
import itertools as it
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


SpikeHelper = namedtuple("SpikeHelper", ("frequency_low", "frequency_high", "intensity", "frequency"))
MolHelper = namedtuple("MolHelper", ("frequency", "uncertainty"))


def reverse_dict(d):
    reversed_dict = defaultdict(list)
    for key, values in d.items():
        for value in values:
            reversed_dict[value].append(key)
    return reversed_dict


def frequency_cover(molecule, spike):
    return (molecule.frequency <= spike.frequency_high + molecule.uncertainty) & (
        molecule.frequency >= spike.frequency_low - molecule.uncertainty)


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
        threshold = 0.0):
        self.sets, self.spike_dict, self.plot_dict = set_generation(spikes, moles, method=method, scoring=scoring, temperature=temperature, threshold=threshold)
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
        set_complexity = [len(set(sum(it.chain(x), ()))) for x in self.sets]
        min_index = set_complexity.index(min(set_complexity))
        return get_mol_set(self.sets[min_index]) if unique else self.sets[min_index]


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

    def visualize(self, lines=True):
        colorgen = cycle(px.colors.qualitative.Plotly)
        colors = {mol: next(colorgen) for mol in set(self.plot_dict.keys())}
        plot_dict = dict(
            sorted(
                self.plot_dict.items(),
            )
        )
        fig = go.Figure()
        if lines:
            linelist = []
            spike_max = self.spikes.data[:,-1].max()
            spike_mean = self.spikes.data[:,-1].mean()
            for molecule, info in plot_dict.items():
                freq = info["match_freq"].tolist()
                intense = (info["match_intensity"] * lines).tolist()
                for f, i in zip(freq, intense):
                    fig.add_trace(
                        go.Scatter(
                            x=[f, f],
                            y=[spike_mean, spike_max] ,
                            mode="lines+text",
                            showlegend=True if molecule not in linelist else False,
                            legendgroup=f"{molecule} spectrum",
                            name=f"{molecule}, {info['score']:.3f}",
                            textposition="bottom center",
                            visible="legendonly",
                            line=dict(color=colors[molecule]),
                        )
                    )
                    linelist.append(molecule)
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
        fig.add_trace(
            go.Scatter(
                x=self.spikes.data[:, 0],
                y=self.spikes.data[:, -1],
                marker_color="black",
                legendgroup="Spectrum",
                name="Spectrum",
                opacity=0.4,
                mode="lines",
                showlegend=False
            )
        )
        return fig
    def save_results(self, path: str):
        fig = self.visualize(lines=True)
        if not os.path.exists(path):
            os.makedirs(path)
        fig.write_html(f"{path}/visualization.html")
        stdout = sys.stdout
        with open(f"{path}/results.txt", 'w') as outfile:
            sys.stdout = outfile
            print("*"*80, flush=True)
            print()
            mols = self.likeliest_molecules()
            for idx, mol in enumerate(mols):
                print(f"{idx}. {mol[0]} \t score: {mol[1]}")
        sys.stdout = stdout
        return




def set_generation(spikes, moles, method=frequency_cover, scoring="intensity_score", temperature=None, threshold = 0.05):
    spike_dict = {f: [] for f in spikes.frequency}
    spacing = (spikes.data[1:, 0] - spikes.data[:-1, 0]).mean() / 2 #what?
    spike_list = [
        spikes.frequency + spacing,
        spikes.frequency - spacing,
        spikes.intensity,
        spikes.frequency,
    ]
    out_mols = {}
    plot_dict = {}
    if temperature is not None:
        print("converting intensities")
        moles = convert_intensities(moles, temperature)
    #mol=moles["Ethyl Cyanide"]
    #spike_broadcast = SpikeHelper(*broadcast(spike_list, mol.frequency))
    #spike_matches = method(mol, spike_broadcast)
    #import pdb; pdb.set_trace()



    for mol in tqdm.tqdm(moles.values(), desc="identifying matches..."):
        #if mol.molecule == "Ethyl Cyanide":
        #  print('ahtath')
        #  import pdb; pdb.set_trace()

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
                if scoring == "intensity_score":
                    I = np.power(10.0, mol.intensity[molidxs[0]]).sum()
                    score = ((np.power(10.0, spike_intensities))).sum() / I
                else:
                    print("error!!")

                if np.isnan(score):
                    score = -1.0
                if score >= threshold:
                    plot_dict[winner] = {
                        "mol_freq": mol.frequency[molidxs],
                        "mol_intensity": mol.frequency[molidxs],
                        "match_freq": mol.frequency[matchids],
                        "match_intensity": mol.intensity[matchids],
                        "score": score
                    }
                    spike_dict[list(spike_dict.keys())[i]].append((mol.molecule, score))
    #import pdb; pdb.set_trace()
    spike_dict = {k: v if len(v) > 0 else ["unknown"] for k, v in spike_dict.items()}
    result = {}
    for s in tqdm.tqdm(spike_dict.keys(), desc="identifying sets..."):
        result.setdefault(s, []).append(
            [ x for z in [ list(it.combinations(spike_dict[s], y + 1))
                    for y in range(len(spike_dict[s])) ] for x in z ] )
    print("Restructuring data..233333")
    return (list(map(list, zip(*list(map(list, list(it.product(*result.values()))[0]))))), spike_dict, plot_dict, )


def get_mol_set(tl):
    return set(sum(it.chain(tl), ()))


def spike_explanations(mol_list, spike_dict):
    return {
        list(spike_dict.keys())[i]: list(it.chain(mol_list[i]))
        for i in range(len(mol_list))  }



