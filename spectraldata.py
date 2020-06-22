import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
import scipy.stats as ss
from pprint import pprint
from utility import doctuple
import pdb

# this is probably not the way forward, but I think having minimalist classes
# makes a lot of sense, and is easier to manage, I also like immutability but
# that is just for messing around
# for unit conversions i think we can just use a helper dict


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
    # this doesnt really matter
    return ss.trim_mean(x, proportiontocut=0.4)


def std(x):
    return x.std()


def doctuple(doc, *args, **kwargs):
    # so we can use docstrings
    nt = namedtuple(*args, **kwargs)
    nt.__doc__ = doc
    return nt


SpectralFile = doctuple(
    """
    A Spectral Data File
    Attributes:
        file_location: path to the file.
        name=None: the name of the planetary body.
        time=None: the the time of observation.
        date=None: the date of observation.
        temp=None: the temperature of the planetary body.
        units=MHz: the units of frequency used.
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


def read(sf: SpectralFile):
    '''
    read: Extract data from a spectralFile
    '''
    data = np.genfromtxt(sf.file_location)
    # *data.T is so cool
    return SpectralData(*sf, data, *data.T, data.shape[0])


SpectralDataStats = doctuple(
    """
    A Spectral Data Object, with statistics
    Attributes:
        file_location: path to the file.
        name=None: the name of the planetary body.
        time=None: the the time of observation.
        date=None: the date of observation.
        temp=None: the temperature of the planetary body.
        units=MHz: the units of frequency used.
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


def get_stats(sd: SpectralData, average_method=trimmed_mean, std_method=std):
    return SpectralDataStats(
        *sd,
        average_method,
        std_method,
        average_method(sd.intensity),
        std_method(sd.intensity),
        (sd.frequency.min(), sd.frequency.max()),
    )


# should spikes have *sd_stats?
Spike = doctuple(
    """doc holder""",
    "Spike",
    ["data", "frequency", "intensity", "index", "spike_method"],
)


def three_sigma_spike(sd):
    return np.where(sd.intensity > (sd.mean + 3 * sd.std))[0]


def identify_spikes(sd, spike_method=three_sigma_spike, **kwargs):
    # the idea is that if the method uses the stats, we will have them, but if
    # the method does not use the stats, we run one useless computation and
    # preserve all the information
    if type(sd) is not SpectralDataStats:
       # warning_text = "Please either run get_stats(spectral data) or provide both average_method and std_method"
       # if "average_method" not in kwargs and "std_method" not in kwargs:
       #     message = f"using trimmed mean and untrimmed std for stats. {warning_text}"
       #     warnings.warn(message, SyntaxWarning)
       # elif "average_method" not in kwargs:
       #     message = f"using trimmed mean for stats. {warning_text}"
       #     warnings.warn(message, SyntaxWarning)
       # elif "std_method" not in kwargs:
       #     message = f"using untrimmed std for stats. {warning_text}"
       #     warnings.warn(message, SyntaxWarning)
        sds = get_stats(sd, **kwargs)
    else:
        sds = sd
    ids = spike_method(sds)
    return Spike(sds.data, sds.frequency[ids], sds.intensity[ids], ids, spike_method)


def plot_spikes(sp):
    plt.plot(sp.data[:, 0], sp.data[:, 1], c="blue")
    for i in sp.frequency:
        plt.axvline(i, c="red", linestyle="--")
    plt.show()
    return



# testing zone
sf = SpectralFile("Titan/Titan/Win0.clean1.contsub_Jy.rest.scom.c.txt")
#print(read(sf).data - convert_units(get_stats(read(sf))).data)
#import pdb; pdb.set_trace()  # XXX BREAKPOINT
##sf = read(sf)
##sd = get_stats(sf)
##pprint(sd)
##
##
### test out all our warnings
##spikes = identify_spikes(sf, three_sigma_spike, average_method=trimmed_mean, std_method = std)
##pdb.set_trace()  # XXX BREAKPOINT
##
##pprint(identify_spikes(sd, three_sigma_spike))
##pprint(spikes)
##
##plot_spikes(spikes)
##plot_spikes(identify_spikes(sd))
##pdb.set_trace()  # XXX BREAKPOINT
#
#
##out = []
##for f in os.listdir("Titan/Titan/"):
##    sf = SpectralFile(f"Titan/Titan/{f}")
##    sd = read(sf)
##    out += [sd]
##    sd_stats = get_stats(sd)
##    plot_spikes(identify_spikes(sd))
##
##import pdb; pdb.set_trace()  # XXX BREAKPOINT
##
