import spectraldata as nasa
import molecules as m
import os

#molecules = m.construct_mega_query

#for f in os.listdir("Titan/Titan/"):
#    sf = nasa.SpectralFile(f"Titan/Titan/{f}")
#    sd = nasa.read(sf)
#    sd = nasa.convert_units(sd, 'GHz')
#    sd_stats = nasa.get_stats(sd)
#    print(sd_stats)
#    import pdb; pdb.set_trace()  # XXX BREAKPOINT
#    nasa.plot_spikes(nasa.identify_spikes(sd))
#    # use these spikes + splatalogue to query JPL
#    # given that that works, figure out optimal choice
#
#import pdb; pdb.set_trace()  # XXX BREAKPOINT

sf = nasa.SpectralFile(f"Titan/Titan/Win0.clean1.contsub_Jy.rest.scom.c.txt")
sd = nasa.read(sf)
#ghz_data = nasa.convert_units(sd, to='GHz')
spikes = nasa.identify_spikes(sd)
mms = m.query_splatalogue(spikes)
out = []
for mm in mms:
    try:
        x = m.get_molecule_data(*mm)
        print(x)
        out.append(x)
    except:
        import pdb; pdb.set_trace()  # XXX BREAKPOINT
        pass


import pdb; pdb.set_trace()  # XXX BREAKPOINT


