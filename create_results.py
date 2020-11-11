import jplspec as jpl
from pprint import pprint
import numpy as np
import sys
from itertools import chain
import os

stdout = sys.stdout

# molecules = m.construct_mega_query
PATH = "data2/"
files = [f"{PATH}{f}" for f in os.listdir(PATH)]
#import pdb; pdb.set_trace()  # XXX BREAKPOINT
#
#f = "Titan/Titan/Win0.clean1.contsub_Jy.rest.scom.c.txt"
#data = jpl.read(jpl.SpectralFile(f))
#spikes = jpl.identify_spikes(data)
#molecules = jpl.get_molecules_from_spikes(spikes, save_path="win0.dmp")
#cover = jpl.SetCovering(spikes, molecules)
#fig = cover.visualize()
#fig2 = cover.visualize(lines=0.02)
#pprint(cover.likeliest_molecules())
import pdb; pdb.set_trace()  # XXX BREAKPOINT
#
#


for f in files:
    print(f)
    data = jpl.read(jpl.SpectralFile(f))
    import pdb; pdb.set_trace()  # XXX BREAKPOINT
    spikes = jpl.identify_spikes(data)
    import pdb; pdb.set_trace()  # XXX BREAKPOINT
    molecules = jpl.get_molecules_from_spikes(spikes)
    import pdb; pdb.set_trace()  # XXX BREAKPOINT
    cover = jpl.SetCovering(spikes, molecules)
    print("writing results")
    with open("results.txt", 'a') as fi:
        sys.stdout = fi
        print("*"*80, flush=True)
        print(f"Processing: {f}")
        print()
        mols = cover.likeliest_molecules()
        for idx, mol in enumerate(mols[:30]):
            print(f"{idx}. {mol[0]} \t score: {mol[1]}")
    sys.stdout = stdout


