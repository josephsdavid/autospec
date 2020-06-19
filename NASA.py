import numpy as np
import scipy.stats as ss
from collections import namedtuple


Spike = namedtuple("spike", ['frequency','intensity','index'])
# this is bad i know
trimmed_mean = lambda x: ss.trim_mean(x, proportiontocut = 0.2)

class SpectralFile(object):
    def __init__(self, floc, name, time, temp, date, units):
    # Name of planet
    # file_loc
    # time
    # temp for connor algo
    # date
    # data (ordered dict?) List of tuples? numpy ext
    # num_entries
    # range_of_spectra
    # is_spectra_continuous
    # is it a sloping signal
    # units of messurment
        # throw error if no file location
        self.floc = floc
        # better name
        self.name = name
        # these default to none
        self.time = time
        self.date = date
        self.temp = temp
        # till here
        # default unit
        self.units = units
        #self.data = None
        self.num_entries = None
        self.spikes = None
        self.std = None
        self.average = None
        self.spectral_range = None
        self.is_signal_sloping = None
        self.is_spectra_continuous = None
        self.load_file()
        # put ornot
        self.get_spectral_range()
    ## Functions
    # open File
    # build object
    # Calcuate noise (average?) remove 5 sigma events and calculate avg
    def load_file(self):
        self.data = np.genfromtxt(self.floc)

    # add a get_data method!
    def get_file_location(self):
        return self.floc

    def get_num_entries(self):
        if not self.num_entries:
            self.num_entries = self.data.shape[-1]
        return self.num_entries

    def get_spectral_range(self):
        if not self.spectral_range:
            self.sprectral_range = (self.data[:,0].min(), self.data[:,0].max())
        return self.sprectral_range

    def get_is_signal_sloping(self):
        # ????
        # do notimplemented thing
        pass

    def get_units(self):
        return self.units

    def get_stats(self, average_method = trimmed_mean):
        # make better, make tests?
        self.average = average_method(self.data[:,1])
        self.std = self.data[:,1].std()
        return self.average, self.std

    def identify_spikes(self, spike_method):
        ids = spike_method(self)[0]
        self.spikes = Spike(self.data[ids,0], self.data[ids,1], ids[0])
        return self.spikes

def three_sigma_spike(sf , average_method = trimmed_mean):
    # argument for all these to be outside of the class
    avg, std = sf.get_stats(average_method)
    intensity = sf.data[:,1]
    return np.where(intensity > (avg + 3*std))

test = SpectralFile("Titan/Titan/Win0.clean1.contsub_Jy.rest.scom.c.txt", "Titan", time=None, temp=None, date=None, units="MhZ")

print(test.identify_spikes(lambda x: three_sigma_spike(x, average_method = lambda y: y.mean())))

import matplotlib.pyplot as plt
plt.plot(test.data[:,0],test.data[:,1], c='blue')
for i in test.spikes.frequency:
    plt.axvline(i,c='red', linestyle='--')
plt.show()

print(test.spikes.frequency)
print(test.spikes.intensity)
print(test.spikes.index)

import pdb; pdb.set_trace()  # XXX BREAKPOINT

print(test.data)



