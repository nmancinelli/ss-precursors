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
                       names=['Vp','Vs','Density','Bottom','Label'])

        self.model.Label = ["Water", "Ice", "Seds1", "Seds2", "Seds3", "Crust1", "Crust2", "Crust3"]
        self.lat = lat
        self.lon = lon

        #cleanup
        subprocess.call("rm tmp tmp2", shell=True)

        return self

    def print(self):
        print("Crust 1.0 at lat, lon = %.2f %.2f" % (self.lat, self.lon))
        print(self.model)
        return(self)

    def save(self,filename = 'crustal_model.pickle'):
        import pickle
        fout = open(filename,'wb')
        pickle.dump(self,fout)
        fout.close()

    def load(self, filename):
        import pickle
        fin = open(filename, 'rb')
        return pickle.load(fin)


if __name__ == "__main__":
    Crust1().get_crustal_model(lat=36, lon=-110).print().save("crustal_model_ColoradoPlateau.pickle")
    Crust1().get_crustal_model(lat=64.5, lon=-110).print().save("crustal_model_SlaveCraton.pickle")
    Crust1().get_crustal_model(lat=62, lon=96).print().save("crustal_model_Siberia.pickle")
    Crust1().get_crustal_model(lat=-27, lon=27).print().save("crustal_model_SouthAfrica.pickle")