import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
import scipy.stats as ss
from pprint import pprint
from .utility import doctuple
import pdb

# this is probably not the way forward, but I think having minimalist classes
# makes a lot of sense, and is easier to manage, I also like immutability but
# that is just for messing around
# for unit conversions i think we can just use a helper dict



SpectralFile = doctuple(
    """
    A Spectral Data File
    Attributes:
        file_location: path to the file.
        name=None: the name of the planetary body.
        time=None: the the time of observation.
        date=None: the date of observation.
        temp=None: the temperature of the planetary body.
        units=GHz: the units of frequency used.
    functions:
        read(SpectralFile):
            returns a SpectralData object
        get_stats(read(SpectralFile)):
            returns a SpectralDataStats object
        identify_spikes(read(SpectralFile)):
            returns a Spike object
    """,
    "SpectralFile",
    ["file_location", "name", "time", "date", "temp", "units"],
    defaults=[None] * 4 + ["GHz"],
)


SpectralData = doctuple(
    """
    A Spectral Data Object
    Attributes:
        file_location: path to the file.
        name=None: the name of the planetary body.
        time=None: the the time of observation.
        date=None: the date of observation.
        temp=None: the temperature of the planetary body.
        units=GHz: the units of frequency used.
        data: the frequency and intensity data.
        frequency: frequency data.
        intensity: intensity data.
        num_entries: number of observations in the dataset.
    functions:
        read(SpectralFile):
            returns a SpectralData object
        get_stats(SpectralData):
            returns a SpectralDataStats object
        identify_spikes(SpectralData):
            returns a Spike object
        convert_units(SpectralData):
            converts the units of a SpectralData object
    """,
    "SpectralData",
    SpectralFile._fields + ("data", "frequency", "intensity", "num_entries"),
)


SpectralDataStats = doctuple(
    """
    A Spectral Data Object, with statistics
    Attributes:
        file_location: path to the file.
        name=None: the name of the planetary body.
        time=None: the the time of observation.
        date=None: the date of observation.
        temp=None: the temperature of the planetary body.
        units=GHz: the units of frequency used.
        data: the frequency and intensity data.
        frequency: frequency data.
        intensity: intensity data.
        num_entries: number of observations in the dataset.
        average_method: function used to calculate the mean
        std_method: function used to calculate the std deviation
        spectral_range: tuple of min and max frequencies
    functions:
        get_stats(SpectralData):
            returns a SpectralDataStats object
        identify_spikes(SpectralDataStats):
            returns a Spike object
        convert_units(SpectralDataStats):
            converts the units of a SpectralDataStats object

    """,
    "SpectralDataStats",
    SpectralData._fields
    + ("average_method", "std_method", "mean", "std", "spectral_range"),
)

# should spikes have *sd_stats?
Spike = doctuple(
    """doc holder""",
    "Spike",
    ["data", "frequency", "intensity", "index", "spike_method", "units"],
)

def three_sigma_spike(sd):
    return sd.intensity > (sd.mean + 3 * sd.std)


def convert_units(s, to='GHz'):
    # fix for SDS class, make pure
    unit_dict = {
        'Hz': 1,
        'kHz': 1e3,
        'MHz': 1e6,
        'GHz': 1e9
    }
    freq = s.frequency.copy()
    freq *= unit_dict[s.units]/unit_dict[to]
    out = s._asdict()
    out['frequency'] = freq
    out['units'] = to
    out['data'][:,0] = freq
    if 'uncertainty' in list(out.keys()):
        out['uncertainty'] *= unit_dict[s.units]/unit_dict[to]
    if type(s) is not SpectralDataStats:
        return type(s)(**out)
    else:
        return get_stats(SpectralData(*s[:-5]), s.average_method, s.std_method)

def trimmed_mean(x):
    return ss.trim_mean(x, proportiontocut=0.2)


def std(x):
    return ss.trimboth(x, 0.1).std()


def get_stats(sd: SpectralData, average_method=trimmed_mean, std_method=std):
    return SpectralDataStats(
        *sd,
        average_method,
        std_method,
        average_method(sd.intensity),
        std_method(sd.intensity),
        (sd.frequency.min(), sd.frequency.max()),
    )


def _read_file(p:str) -> np.ndarray:
    with open(p, 'r') as f:
        raw = [line.replace('\n', '') for line in f.readlines()]
    raw = [row for row in raw if row]
    metadata = ''.join([x for x in raw if '#' in x])
    data = np.array([[float(x) for x in y.split()] for y in raw if '#' not in y])
    if metadata:
        print("metadata detected")
        GHZ = 'MHz' not in metadata
        if not GHZ:
            print("cleaning units")
            data[:,0] /= 1000
    # fix this logic!
    data = data[1:,:]
    data = data[:-1,:]
    return data

def read(sf: SpectralFile):
    '''
    read: Extract data from a spectralFile
    '''
    if type(sf.file_location) is not list:
        data = _read_file(sf.file_location).astype(np.float32)
    else: #multiple file case
        data = np.vstack([np.genfromtxt(f) for f in sf.file_location])
    # *data.T is so cool
    return SpectralData(*sf, data, *data.T, data.shape[0])


def identify_spikes(sd, spike_method=three_sigma_spike, **kwargs):
    if type(sd) is not SpectralDataStats:
        sds = get_stats(sd, **kwargs)
    else:
        sds = sd
    ids = np.where(spike_method(sds))[0]
    return Spike(sds.data, sds.frequency[ids], sds.intensity[ids], ids, spike_method, sds.units)


def plot_spikes(sp):
    plt.plot(sp.data[:, 0], sp.data[:, 1], c="blue")
    for i in sp.frequency:
        plt.axvline(i, c="red", linestyle="--")
    plt.show()
    return

