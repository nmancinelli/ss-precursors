#!/usr/bin/env python
#
from numpy import sqrt
import subprocess
from os import chdir, getcwd
import pickle as pickle
import pandas as pd
from subprocess import call

workingdir = '/Users/mancinelli/Desktop/SS_Precursors/minos'
modeldir = '/Users/mancinelli/Desktop/SS_Precursors/s362ani'

VpVsRatioDefault=sqrt(3)

def calculate_anis_from_moduli(pvh_ratio, svh_ratio, bulk_modulus, shear_modulus, eta):
    # invert for C, N
    from numpy import zeros, matmul
    from numpy.linalg import inv

    M = zeros(4).reshape(2, 2)
    M[0, 0] = 1. + 4. * pvh_ratio ** -2 + 4. * eta * pvh_ratio ** -2
    M[0, 1] = -(4. + 8 * eta * svh_ratio ** 2)
    M[1, 0] = 1. + pvh_ratio ** -2 - 2. * eta * pvh_ratio ** -2
    M[1, 1] = 6. * svh_ratio ** 2 + 5. + 4. * eta * svh_ratio ** -2

    tmp = matmul(inv(M), [9 * bulk_modulus, 15 * shear_modulus])
    C = tmp[0]
    N = tmp[1]
    A = pvh_ratio ** -2 * C
    L = svh_ratio ** 2 * N
    F = eta * A - 2. * eta * L

    return C, N, A, L, F

    # string = 'C, N, A, L, F = %f %f %f %f %f' % (C, N, A, L, F)
    # print(string)
    # print('kappa + 4/3*mu, mu = %f %f' % ( bulk_modulus + 4./3. * shear_modulus, shear_modulus) )
    # print(M)

class Minos:

    def __init__(self,id="default", cardfile = None):
        self.id  = id
        self.cwd = getcwd()
        self.workingdir = workingdir
        self.cardfile = cardfile

    def run(self):
        chdir(self.workingdir)
        subprocess.call('bash run_minos.bash %s' % self.id, shell=True)
        chdir(self.cwd)

        return self

    def read_output(self):
        from matplotlib import pylab as plt
        plt.style.use('ggplot')
        from numpy import array

        def main():
            return read_curves('%s.out' % self.id, '-')

        def read_curves(infile, ls='-'):
            fin = open(infile)
            a, b, c = [], [], []
            for line in fin.readlines():
                nfo = line.strip('\n').split()
                a.append(float(nfo[0])) #phase
                b.append(float(nfo[1])) #freq
                c.append(float(nfo[2])) #group

            period = 1000. / array(b)

            return period, a

        chdir(self.workingdir)
        self.period, self.c = main()
        chdir(self.cwd)

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

        chdir(self.workingdir)
        main()
        chdir(self.cwd)

    def save(self,filename='minos.pickle'):
        fout=open(filename,'wb')
        pickle.dump(self, fout)
        fout.close()
        return

    def load(self,filename='minos.pickle'):
        fin=open(filename,'rb')
        self = pickle.load(fin)
        fin.close()
        return self

class CardFile:
    def __init__(self, name='default'):
        pass

    def read(self, filename="STW105.txt"):
        df=pd.read_csv(filename, skiprows=3, delim_whitespace=True,
                       names=['Radius','Density','Vpv','Vsv','Qkappa','Qmu','Vph','Vsh','Eta'])
        self.df = df

        self.__calculate_anisotropic_scaling()

        return self

    def regional_profile_from_kustow(self, lat=0, lon=0):
        print('Pulling regional profile...')
        curdir = getcwd()
        chdir(modeldir)

        names = ['dlnvsh','dlnvsv','dlnvph','dlnvpv','dlneta','dlnrho']  #vshout,vsvout,vphout,vpvout,etaout,rhoout
        for ii,row in self.df.iterrows():
            depth_in_km = (max(self.df.Radius) - row["Radius"]) / 1000.0

            if depth_in_km > 3000:
                continue

            call('bash get_model_values.bash %f %f %d' % (lat, lon, depth_in_km), shell=True)

            tmp = pd.read_csv('tmp2', delim_whitespace=True, names = names)

            self.df.Vpv[ii] *= (1.+tmp.dlnvpv/100.)
            self.df.Vph[ii] *= (1.+tmp.dlnvph/100.)
            self.df.Vsh[ii] *= (1.+tmp.dlnvsh/100.)
            self.df.Vsv[ii] *= (1.+tmp.dlnvsv/100.)

        chdir(curdir)
        print('done.')

        self.__calculate_anisotropic_scaling()

        return self


    def __calculate_anisotropic_scaling(self):
        self.df["Pvh_Ratio"] = self.df.Vpv / self.df.Vph
        self.df["Svh_Ratio"] = self.df.Vsv / self.df.Vsh

    def __calculate_anisotropic_constants(self):
        #p_scaling = self.df.Vph / self.df.Vpv
        #s_scaling = self.df.Vsh / self.df.Vsv

        C = self.df.Vpv**2 * self.df.Density
        A = self.df.Vph**2 * self.df.Density
        L = self.df.Vsv**2 * self.df.Density
        N = self.df.Vsh**2 * self.df.Density
        F = self.df.Eta*(A-2.*L)

        self.df["C"] = C
        self.df["A"] = A
        self.df["L"] = L
        self.df["N"] = N
        self.df["F"] = F

        kappa_voigt = 1./9. * (C + 4.*A - 4.*N + 4.*F)
        mu_voigt = 1./15. * (C + A + 6.*L + 5.*N - 2.*F)
        lambda_voigt = kappa_voigt - 2./3. * mu_voigt

        Vp_Voigt =  ( (lambda_voigt + 2. * mu_voigt) / self.df.Density )**0.5
        Vs_Voigt =  (mu_voigt / self.df.Density)**0.5

        self.df["VpVoigt"] = Vp_Voigt
        self.df["VsVoigt"] = Vs_Voigt

        return self

    def set_crustal_region(self, crustal_model):

        rearth0 = 6371000.

        rearth = rearth0 + crustal_model.loc[0,"Bottom"] * 1000.0
        rmoho = rearth0 - abs(crustal_model.loc[7, "Bottom"]) * 1000.
        cf1 = self.df.query("Radius < %d" % (rmoho))
        ii = len(cf1) - 1

        crustal_model = crustal_model.query("Density > 0.01")

        inds = crustal_model.index.tolist()[::-1]

        for kk,iclay in enumerate(inds[:-2]):
            rmoho = rearth0 - abs(crustal_model.loc[iclay, "Bottom"]) * 1000.
            if crustal_model.loc[iclay, "Density"] == 0:
                continue
            for iknot in [1, 0]:
                ii += 1

                if iknot == 0:
                    jj = inds[kk]
                else:
                    jj = inds[kk-1]

                if jj > 0:
                    rho = crustal_model.loc[jj, "Density"] * 1000.
                    vp = crustal_model.loc[jj, "Vp"] * 1000.
                    vs = crustal_model.loc[jj, "Vs"] * 1000.
                    cf1.loc[ii] = [rmoho, rho, vp, vs, 999., 999., vp, vs, 1., 1., 1.]
                else:
                    cf1.loc[ii] = cf1.loc[ii - 1]
                    cf1.loc[ii, 'Radius'] = rmoho
        ii += 1
        cf1.loc[ii] = [rearth, rho, vp, vs, 999., 999., vp, vs, 1., 1., 1.]

        self.df = cf1

        return self

    def set_crustal_region_shen(self, crustal_model, zmoho_in_km):

        rearth0 = 6371000.
        rmoho = rearth0 - zmoho_in_km*1000.

        cf1    = self.df.query("Radius < %d" % (rmoho))
        ii     = len(cf1) - 1

        crustal_model = crustal_model.query("Depth < %f" % zmoho_in_km).sort_values('Depth', ascending=False)

        for jj, row in crustal_model.iterrows():
            ii += 1
            radius = rearth0 - row.Depth*1000.0
            cf1.loc[ii] = [radius, row.Rho*1000., row.Vp*1000., row.Vs*1000., 99999., 300., row.Vp*1000., row.Vs*1000., 1., 1., 1.]

        self.df = cf1

        return self

    def set_mantle_region(self, vp_fun, vs_fun, rho_fun, zmax=350, zmoho=35.):
        for ii,row in self.df.iterrows():
            z = (max(self.df.Radius) - row["Radius"]) / 1000.0

            if z > zmax:
                continue

            vp, vs, rho = vp_fun(z), vs_fun(z), rho_fun(z)

            shear_modulus    = vs**2 * rho
            bulk_modulus     = vp**2 * rho - 4./3. * shear_modulus

            pvh_ratio = row["Pvh_Ratio"]
            svh_ratio = row["Svh_Ratio"]
            eta       = row["Eta"]

            if z <= zmoho:
                pvh_ratio = 1.
                svh_ratio = 1.
                eta = 1.


            C, N, A, L, F = calculate_anis_from_moduli(pvh_ratio, svh_ratio, bulk_modulus, shear_modulus, eta)

            self.df.loc[ii, 'Vpv'] = (C / rho) ** 0.5 * 1000.
            self.df.loc[ii, 'Vph'] = (A / rho) ** 0.5 * 1000.
            self.df.loc[ii, 'Vsv'] = (L / rho) ** 0.5 * 1000.
            self.df.loc[ii, 'Vsh'] = (N / rho) ** 0.5 * 1000.

            self.df.loc[ii, 'Density'] = rho * 1000.

            #print(z, self.df.loc[ii, 'Vpv'], self.df.loc[ii, 'Vph'], pvh_ratio, svh_ratio, eta)

        return self

    def set_mantle_attenuation(self, qmu_fun, zmax=350, zmoho=35.):
        for ii,row in self.df.iterrows():
            z = 6371. - row["Radius"]/1000.0
            if z > zmax or z < zmoho:
                continue

            qmu = qmu_fun(z)
            if qmu > 999999.:
                qmu = 999999.

            self.df.loc[ii, 'Qmu']    = qmu
            self.df.loc[ii, 'Qkappa'] = 999999.

        return self

    def plot(self, rmin = 0):
        from matplotlib import pylab as plt
        plt.style.use('ggplot')
        qry = "Radius >= %d " % (rmin*1000)
        df = self.df.query(qry)
        plt.subplot(1,2,1)
        plt.plot(df.Vsh/1000, df.Radius/1000,label='Vsh (km/s)')
        plt.plot(df.Vsv/1000, df.Radius/1000,label='Vsv (km/s)')
        plt.plot(df.VsVoigt / 1000, df.Radius / 1000, label='Vs Voigt (km/s)')
        plt.legend(loc=3)
        plt.subplot(1,2,2)
        plt.plot(df.Vph/1000, df.Radius/1000,label='Vph (km/s)')
        plt.plot(df.Vpv/1000, df.Radius/1000,label='Vpv (km/s)')
        plt.plot(df.VpVoigt / 1000, df.Radius / 1000, label='Vp Voigt (km/s)')
        plt.ylabel('Radius (km)')
        plt.legend(loc=3)
        plt.show()

    def write(self, ref_period=-1.0, filename=workingdir + '/default.txt', ):
        fout = open(filename,'w')

        nknots = len(self.df.Radius.tolist())

        fout.write("""# radius [m] density [kg/m^3] vpv [m/s] vsv [m/s] Q kappa Q miu vph [m/s] vsh [m/s] eta [m/s] REF\n""")
        fout.write("""1 %.1f 1 1\n%d  180  358  717  739\n""" % (ref_period, nknots) )

        for irow, row in self.df.iterrows():
            line = "%7.0f. %8.2f %8.2f %8.2f %8.1f %8.1f %8.2f %8.2f %8.5f\n" %\
                   (row.Radius, row.Density, row.Vpv, row.Vsv, row.Qkappa, row.Qmu, row.Vph, row.Vsh, row.Eta)
            fout.write(line)
        fout.close()

        return self

if __name__ == "__main__":
    cf= CardFile().read()
    print(cf.df.query("Radius > 6000000")[["Radius","Vpv","Vph","Pvh_Ratio"]])
    #cf.calculate_ani(1.01,1.02,5.0,3.0,1.01)
