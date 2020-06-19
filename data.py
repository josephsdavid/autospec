import numpy as np
from typing import *
from os import listdir


# subclass numpy array
class SpectralFile(np.ndarray):
    def __new__(cls, path: str) -> None:
        obj = np.genfromtxt(path).view(cls)
        # note the lack of np.copy!
        obj.frequency = obj[:,0]
        obj.flux = obj[:,1]
        return obj
    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.frequency = getattr(obj, 'frequency', None)
        self.flux = getattr(obj, 'flux', None)


class SpectralData(list):
    def __init__(self, directory: str) -> None:
        super().__init__((SpectralFile(f"{directory}/{x}") for x in listdir(directory)))
    def flux(self) -> np.ndarray:
        return np.concatenate([self[i].flux for i in range(self.__len__())])
    def frequency(self) -> np.ndarray:
        return np.concatenate([self[i].frequency for i in range(self.__len__())])
    def array(self) -> np.ndarray:
        return np.concatenate([self[i] for i in range(self.__len__())])



Spectra = Union[SpectralFile, SpectralData]


# this has side effects! and is overcomplcicated! copy doesnt work!
def _constructor(fn: Callable[...,SpectralFile]) -> Callable[..., Spectra]:
    def output(s: Spectra, *args, **kwargs):
        if isinstance(s, SpectralFile):
            return(fn(s, *args, **kwargs))
        else:
            for i in range(len(s)):
                s[i] = output(s[i], *args, **kwargs)
            return s
    return output


def _change_units(s: Spectra, mul: float = 1e-3) -> SpectralFile:
    s.frequency *= mul
    return s



def _temp_adjust(s: Spectra, fn: Callable[..., SpectralFile]) -> SpectralFile:
    s.flux = fn(s.flux)
    return s

change_units = _constructor(_change_units)
temp_adjust = _constructor(_temp_adjust)



sd = SpectralData("Titan/Titan")

print(sd[0])
print(sd[0].frequency)
print(sd[0].flux)
print(sd)
print(sd.frequency())
print(sd.flux())
print(sd.array())
sd = change_units(sd)
print(sd[0])
print(sd[0].frequency)
print(sd[0].flux)
print(sd)
print(sd.frequency())
print(sd.flux())
print(sd.array())
import pdb; pdb.set_trace()  # XXX BREAKPOINT

# if we assign, it still modifies the original:
x = change_units(sd)
print(x[0] - sd[0])

