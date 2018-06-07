#!/usr/bin/env python
#
alpha = 2.9e-5
Tp = 1350.0
g = 9.8
cp = 1250.0
crustal_thickness = 41.

class Geotherm:
    def __init__(self,
                 depths=None, temps=None, q0=None,
                 rhoH=None, rhoH_se=None,
                 xeno_depths=None, xeno_temps=None):

        self.q0 = q0
        self.rhoH = rhoH
        self.rhoH_se = rhoH_se
        self.alpha = alpha
        self.Tp = Tp
        self.g = g
        self.cp = cp
        self.crustal_thickness = crustal_thickness
        self.depths = depths
        self.temps = temps
        self.xeno_depths = xeno_depths
        self.xeno_temps = xeno_temps

    def save(self, filename='default_geotherm.pickle'):
        import pickle
        fout = open(filename, 'wb')
        pickle.dump(self, fout)
        fout.close()
