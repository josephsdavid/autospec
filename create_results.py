import autospec
from scipy import stats as ss
import numpy as np
import sys
import os

stdout = sys.stdout


PATH = "present/"
fe="Thelen_2020_SPW37_spectrum.txt"
#"bonus/2013_00227_X1942_spectrum_f.txt
files = [f"{PATH}{f}" for f in os.listdir(PATH) if f == fe]
#files += [f"data/{f}" for f in os.listdir("data/") ]


for f in files:
    os.system(f"notify-send 'started!' {f}")
    data = autospec.read(autospec.SpectralFile(f))
    data = data._asdict()
    data['frequency'] = data['frequency'][np.where((data['frequency'] > 255.8) & (data['frequency']< 256.08)) ]
    data['intensity'] = data['intensity'][np.where((data['frequency'] > 255.8) & (data['frequency']< 256.08)) ]
    data2 = autospec.SpectralData(**data)
    spikes = autospec.identify_spikes(data2, std_method=lambda x: ss.trimboth(x, 0.1).std())
    autospec.plot_spikes(spikes)
#    import pdb; pdb.set_trace()
    save_path = f.split("/")[-1].replace(".", "_").replace("_txt","") + "working_example_zoom2"
    molecules = autospec.get_molecules_from_spikes(spikes, f"{save_path}.obj")
    cover = autospec.SetCovering(spikes, molecules, threshold=0.1)
    cover.save_results(save_path)
    os.system("notify-send 'done!'")


#
