#!/usr/bin/env python
#
import pandas as pd
from numpy import arange

path = '/Users/mancinelli/Desktop/SS_Precursors/shen'

class CrustalModel:
    def __init__(self, filename):
        path2file = path + '/' + filename
        self.model = pd.read_csv(path2file, delim_whitespace=True, names=["Vs","Vp","Rho"])

        self.model["Depth"] = arange(len(self.model.Vp.tolist())) * 0.5

        tmp = self.model.query("Vp > 8.0")
        self.zmoho = min(tmp.Depth.tolist())

    def splice_to_card_file(self, cardfile, zmoho):
        pass

if __name__ == "__main__":
    model = CrustalModel('coloplat')
    print(model.zmoho)