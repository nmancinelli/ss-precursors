#!/usr/bin/env python
#
from numpy import sqrt
VpVsRatioDefault=sqrt(3)

class Minos:
    def __init__(self,id="default"):
        self.id = "default"
        pass

    def run(self):
        import subprocess
        subprocess.call('bash run_minos.bash %s' % self.id, shell=True)
        return self

    def plot_output(self):
        from matplotlib import pylab as plt
        plt.style.use('ggplot')
        from numpy import array

        def main():
            plot_curves('%s.out' % self.id, '-')
            plt.show()

        def plot_curves(infile, ls='-'):
            fin = open(infile)
            a, b, c = [], [], []
            for line in fin.readlines():
                nfo = line.strip('\n').split()
                a.append(float(nfo[0]))
                b.append(float(nfo[1]))
                c.append(float(nfo[2]))

            period = 1000. / array(b)

            # plt.plot(period,a, ls, label='Phase Velocity')
            plt.plot(period, c, ls, label=infile)

            plt.ylabel('Rayleigh Group Velocity (km/s)')
            plt.xlabel('Period')

            plt.legend(fancybox=True)

            plt.xlim(0, 100)
            plt.ylim(0.5, 5)

        main()

class CardFile:
    def __init__(self, name='default'):
        pass

    def read(self, filename="STW105.txt"):
        import pandas as pd
        df=pd.read_csv(filename, skiprows=3, delim_whitespace=True,
                       names=['Radius','Density','Vpv','Vsv','Qkappa','Qmu','Vph','Vsh','Eta'])
        self.df = df
        return self

    def set_crustal_region(self, rmin, rmax, vs, vpvs_ratio):

        self.df.loc[(self.df.Radius >= rmin * 1000.) & (self.df.Radius <= rmax * 1000.), 'Vsv'] = vs
        self.df.loc[(self.df.Radius >= rmin * 1000.) & (self.df.Radius <= rmax * 1000.), 'Vsh'] = vs
        self.df.loc[(self.df.Radius >= rmin * 1000.) & (self.df.Radius <= rmax * 1000.), 'Vpv'] = vs*vpvs_ratio
        self.df.loc[(self.df.Radius >= rmin * 1000.) & (self.df.Radius <= rmax * 1000.), 'Vph'] = vs*vpvs_ratio

        return self

    def set_mantle_region(self, rmin, rmax, vs, vpvs_ratio):
        pass

    def plot(self, rmin = 0):
        from matplotlib import pylab as plt
        plt.style.use('ggplot')
        qry = "Radius >= %d " % (rmin*1000)
        df = self.df.query(qry)
        plt.plot(df.Vpv/1000, df.Radius/1000,label='Vpv (km/s)')
        plt.plot(df.Vsv/1000, df.Radius/1000,label='Vsv (km/s)')
        plt.ylabel('Radius (km)')
        plt.legend(loc=3)
        plt.show()

    def write(self,filename='default.txt'):
        fout = open(filename,'w')

        fout.write("""# radius [m] density [kg/m^3] vpv [m/s] vsv [m/s] Q kappa Q miu vph [m/s] vsh [m/s] eta [m/s] REF\n""")
        fout.write("""1 1. 1 1\n750  180  358  717  739\n""")

        for irow, row in self.df.iterrows():
            line = "%7.0f. %8.2f %8.2f %8.2f %8.1f %8.1f %8.2f %8.2f %8.5f\n" %\
                   (row.Radius, row.Density, row.Vpv, row.Vsv, row.Qkappa, row.Qmu, row.Vph, row.Vsh, row.Eta)
            fout.write(line)
        fout.close()

        return self


if __name__ == "__main__":
    cf= CardFile().read()
    cf.set_crustal_region(6340., 6371., 2500., VpVsRatioDefault)
    cf.set_crustal_region(6320., 6340., 3500., VpVsRatioDefault)
    cf.set_crustal_region(6300., 6320., 4000., VpVsRatioDefault).write()

    Minos().run().plot_output()
