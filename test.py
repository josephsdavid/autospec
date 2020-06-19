import scratch as nasa
import molecules as m
import os

molecules = construct_mega_query()

for f in os.listdir("Titan/Titan/"):
    sf = nasa.SpectralFile(f"Titan/Titan/{f}")
    sd = nasa.read(sf)
    sd = nasa.change_units(sd, 'GHz')
    sd_stats = nasa.get_stats(sd)
    print(sd_stats)
    nasa.plot_spikes(nasa.identify_spikes(sd))

import pdb; pdb.set_trace()  # XXX BREAKPOINT

