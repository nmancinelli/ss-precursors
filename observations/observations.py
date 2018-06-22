#!/usr/bin
#
import pandas as pd
from numpy import mean, std

path = '/Users/mancinelli/Desktop/SS_Precursors/observations/'

class PhaseVelocities:
    def __init__(self, lat, lon, source='Cloud'):

        self.measurements = {}
        self.lat = lat
        self.lon = lon

        if source == 'Cloud':
            self = self.get_velocities_cloud()
        else:
            self = self.get_velocities_ma()

    def get_velocities_cloud(self):

        for ii, period in enumerate([25,40,50,60,80,100,120,140, 180]):

            if period >= 120:
                filt=500
            else:
                filt=200

            if period == 180:
                tag = '6.3'
            else:
                tag = ''

            fname = '%s/cloud/chelm%s_selmedian5_filt%d_%ds_stripped' % (path, tag,filt,period)
            df0 = pd.read_csv(fname, names = ['Lat', 'Lon', 'C'], delim_whitespace=True)
            df0 = df0.query("Lat == %d and Lon == %d" %(self.lat, self.lon))

            self.measurements[period] = df0.C.tolist()[0]
        return self

    def get_velocities_ma(self):
        from subprocess import call
        from os import chdir

        chdir('%s/ma/rayleigh' % path)
        call('bash get.bash %f %f' % (self.lat, self.lon), shell=True)
        chdir('%s/..' % path)

        df = pd.read_csv('%s/ma/rayleigh/point.1' % path, names = ['freq_mHz','v1','v2'], delim_whitespace=True)

        df["Period"] = (df.freq_mHz / 1000.)**-1

        for _, row in df.iterrows():
            if row.v1 < 0.0:
                continue
            self.measurements[row.Period] = row.v1*1000.

        return self

class Region:
    def __init__(self,lat1,lat2,lon1,lon2,source):
        self.lats = range(lat1,lat2+1)
        self.lons = range(lon1,lon2+1)
        self.pvs  = {}

        for lat in self.lats:
            for lon in self.lons:
                pv = PhaseVelocities(lat, lon, source=source)
                self.pvs[lat,lon] = pv

        periods = self.pvs[lat1,lon1].measurements.keys()
        self.mean, self.std = {},{}

        for period in periods:
            pvlist = []
            for _, pv in self.pvs.items():
                pvlist.append(pv.measurements[period])

            self.mean[period] = mean(pvlist)
            self.std[period]  = std(pvlist)


if __name__ == "__main__":
    pv = PhaseVelocities(36,-110, "Ma")
    print(pv.measurements)
