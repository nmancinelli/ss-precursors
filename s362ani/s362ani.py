#!/usr/bin/env python
#
class Profile():
    def __init__(self, lat, lon, depths = range(1,410,4)):
        self._read_ref()
        self.lat=lat
        self.lon=lon
        self.df["Depth"] = (max(self.df.Radius) - self.df.Radius) /1000.0

    def get_model_perturbations(self):
        vpv,vph,vsv,vsh,eta = [],[],[],[],[]

        for irow, row in self.df.iterrows():
            depth = row.Depth

            df = self.get_model_perturbation(depth)
            vpv.append( float(df.vpv)/100. )
            vph.append( float(df.vph)/100. )
            vsv.append( float(df.vsv)/100. )
            vsh.append( float(df.vsh)/100. )
            eta.append( float(df.eta)/100. )

        self.df["Dlnvpv"] = vpv
        self.df["Dlnvph"] = vph
        self.df["Dlnvsv"] = vsv
        self.df["Dlnvsh"] = vsh
        self.df["Dlneta"] = eta

        return self

    def get_model_perturbation(self, depth):
        import subprocess
        import pandas as pd
        subprocess.call("bash get_model_values.bash %f %f %f" % (self.lat, self.lon, depth), shell=True)
        df=pd.read_csv('tmp2', delim_whitespace=True, names=['vsh','vsv','vph','vpv','eta','rho'])
        subprocess.call("rm tmp tmp2", shell=True)
        return df

    def _read_ref(self, filename="PROGRAMS/REF"):
        import pandas as pd
        df=pd.read_csv(filename, skiprows=3, delim_whitespace=True,
                       names=['Radius','Density','Vpv','Vsv','Qkappa','Qmu','Vph','Vsh','Eta'])
        self.df = df
        return self

if __name__ == "__main__":
    from matplotlib import pylab as plt
    plt.style.use('ggplot')
    profile = Profile(60,-110).get_model_perturbations()
    plt.plot(profile.df.Vsv, profile.df.Depth)
    plt.plot(profile.df.Vsv * (1. + profile.df.Dlnvsv), profile.df.Depth)
    plt.show()


