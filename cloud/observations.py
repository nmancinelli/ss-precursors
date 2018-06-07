#!/usr/bin
#
import pandas as pd

path = '/Users/mancinelli/Desktop/SS_Precursors/cloud'

class PhaseVelocities:
    def __init__(self, lat, lon):

        self.measurements = {}

        for ii, period in enumerate([25,40,50,60,80,100,120,140, 180]):

            if period >= 120:
                filt=500
            else:
                filt=200

            if period == 180:
                tag = '6.3'
            else:
                tag = ''

            fname = '%s/chelm%s_selmedian5_filt%d_%ds_stripped' % (path, tag,filt,period)
            df0 = pd.read_csv(fname, names = ['Lat', 'Lon', 'C'], delim_whitespace=True)
            df0 = df0.query("Lat == %d and Lon == %d" %(lat, lon))

            self.measurements[period] = df0.C.tolist()[0]

        self.lat = lat
        self.lon = lon

if __name__ == "__main__":
    pv = PhaseVelocities(36,-110)
    print(pv.measurements)
