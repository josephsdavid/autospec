import jplspec as jpl
from scipy import stats as ss
from pprint import pprint
import numpy as np
import sys
from itertools import chain
import os

stdout = sys.stdout


PATH = "data2/"
files = [f"{PATH}{f}" for f in os.listdir(PATH) if f=="all_win0_fluxJy.txt"]

for f in files:
    data = jpl.read(jpl.SpectralFile(f))
    spikes = jpl.identify_spikes(data)
    save_path = f.split("/")[-1].replace(".", "_").replace("_txt","") + "working_example"
    molecules = jpl.get_molecules_from_spikes(spikes, f"{save_path}.obj")
    cover = jpl.SetCovering(spikes, molecules, threshold=0.01)
    cover.save_results(save_path)


