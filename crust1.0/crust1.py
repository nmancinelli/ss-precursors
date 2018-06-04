#!/usr/bin/env python
#
class Crust1:
    def __init__(self):
        pass

    def get_crustal_model(self, lat, lon):
        import subprocess
        import pandas as pd
        subprocess.call("bash get_crustal_model.bash %f %f > /dev/null" % (lat, lon), shell=True)

        self.model=pd.read_csv('tmp2', skiprows=1, delim_whitespace=True,
                       names=['Density','Vp','Vs','Bottom','Label'])

        self.model.Label = ["Water", "Ice", "Seds1", "Seds2", "Seds3", "Crust1", "Crust2", "Crust3"]
        self.lat = lat
        self.lon = lon

        #cleanup
        subprocess.call("rm tmp tmp2", shell=True)

        return self

    def print(self):
        print("Crust 1.0 at lat, lon = %.2f %.2f" % (self.lat, self.lon))
        print(self.model)

if __name__ == "__main__":
    Crust1().get_crustal_model(lat=60,lon=-105).print()