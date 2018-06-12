#!/usr/bin/env python
#
import pandas as pd

from matplotlib import pylab as plt
plt.style.use('ggplot')

from numpy import log10

path = '/Users/mancinelli/Desktop/SS_Precursors/attenuation'


class Profile:
    def __init__(self, source_file):
        path2file = path + '/' + source_file
        self.df = pd.read_csv(path2file, delim_whitespace=True, names=["Depth","dum1","dum2","dum3","Qmu"])

    def setup_interpolation_functions(self):
        from scipy.interpolate import interp1d
        self.fun = interp1d(self.df.Depth.tolist(),self.df.Qmu.tolist())

        return self

    def plot(self):
        plt.plot(log10(self.df.Qmu), self.df.Depth)
        plt.show()


if __name__ == "__main__":
    prof = Profile('seismic_JF10_gs10mm_Temp_AD0p40_TP1350_ZG155').setup_interpolation_functions()
    print(prof.fun(299))
