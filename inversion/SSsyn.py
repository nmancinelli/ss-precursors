import matplotlib.pyplot as plt
from numpy import zeros, array, hanning
plt.style.use('ggplot')

class LayerCakeModel:

    def get_nodes_and_values(lcm):
        zs = lcm.zs


        return zs, vsh, rho

    def __init__(self, zs=None, vsh=None, rho=None):
        """

        :param zs: Nodes
        :param vsh: Specify vsh at each node
        :param rho: Specify rho at each node
        """
        from numpy import arange

        if zs is None:
            self.zs = arange(0, 350, 1)
        else:
            self.zs = zs.copy()

        if vsh is None:
            self.vsh = self.zs*0.0 + 5.
        else:
            self.vsh = vsh.copy()

        if rho is None:
            self.rho = self.zs*0.0 + 4.5
        else:
            self.rho = rho.copy()

        #Upsample here

        if True:
            from scipy.interpolate import interp1d

            fvsh = interp1d(zs, vsh)
            frho = interp1d(zs, rho)

            self.zs = arange(0, 410, 2.0)
            self.rho = frho(self.zs)
            self.vsh = fvsh(self.zs)

    def layerize(self):
        self.h = []
        self.vslay = []
        self.rholay = []
        self.ztop = []
        self.zbot = []
        self.zmid = []


        for ii, each in enumerate(self.zs):
            if ii == 0:
                continue

            h = self.zs[ii] - self.zs[ii - 1]
            if h == 0:
                continue

            self.ztop.append(self.zs[ii-1])
            self.zbot.append(self.zs[ii])
            self.zmid.append((self.zs[ii-1]+self.zs[ii]) / 2.)

            self.h.append(h)
            self.vslay.append((self.vsh[ii] + self.vsh[ii - 1]) / 2.)
            self.rholay.append((self.rho[ii] + self.rho[ii - 1]) / 2.)

        self.vslay = array(self.vslay)

        return self

    def compute_perturbations(self,lcm0):
        self.dvs_vs = []
        for ii, _  in enumerate(self.vslay):
            self.dvs_vs.append( (self.vslay[ii] - lcm0.vslay[ii]) / lcm0.vslay[ii]  )

        return self

    def compute_layers_from_model_parameters(self,dm):
        for ilay in range(len(self.vslay)):
            pert_pct = dm[ilay]
            if ilay == None:
                pass
            else:
                self.vslay[ilay] = self.vslay[ilay] * (1.+pert_pct)

        return self


    def plot(self):
        plt.plot(self.vsh, self.zs)
        plt.ylim(max(self.zs), 0)
        plt.show()
        return self


def rtcoefSH(vs1, den1, vs2, den2, hslow):
    from numpy import arcsin, cos
    j1 = arcsin(hslow * vs1)
    j2 = arcsin(hslow * vs2)
    a = vs1 * den1 * cos(j1) - vs2 * den2 * cos(j2)
    b = vs1 * den1 * cos(j1) + vs2 * den2 * cos(j2)

    c = 2 * den1 * vs1 * cos(j1)

    return a / b, c / b


class SSsyn:
    def __init__(self, layercake_model, ray_param_s_km=0.113):
        dt = 0.
        p = ray_param_s_km

        dts, impedance_contrasts = [], []

        for ii, h in enumerate(layercake_model.h):
            if ii + 1 == len(layercake_model.h):
                continue
            u = 1. / layercake_model.vslay[ii]
            dt = dt + 2. * u ** 2 * h / (u ** 2 - p ** 2) ** 0.5

            vs1 = layercake_model.vslay[ii + 1]
            vs2 = layercake_model.vslay[ii]
            den1 = layercake_model.rholay[ii+1]
            den2 = layercake_model.rholay[ii+1]

            rcoef, tcoef = rtcoefSH(vs1, den1, vs2, den2, p)

            dts.append(-dt)
            impedance_contrasts.append(rcoef)

        self.dts = dts.copy()
        self.impedance_contrasts = impedance_contrasts.copy()

    def construct_spike_train(self, dt=0.2):
        from numpy import zeros, arange
        t1 = -350
        t2 = 25
        npts = int((t2 - t1) / dt)
        s = zeros(npts)
        for jj, t in enumerate(self.dts):
            it = int((t - t1) / dt)
            s[it] += self.impedance_contrasts[jj]

        it0 = int(-t1 / dt)
        s[it0] = 1.0

        self.seis = s
        self.ts = arange(len(s)) * dt + t1
        self.dt = dt

        return self

    def convolve_with_wavelet(self):
        from scipy.signal import ricker, convolve, hilbert

        self.gwin = ricker(500, 21)
        self.seis = convolve(self.seis, self.gwin, 'same')

        return self

    def normalize(self):
        namp = max(abs(self.seis))
        self.seis = self.seis / namp
        return self

    def interpolate_to_time_points(self,time_points):
        from scipy.interpolate import interp1d
        f = interp1d(self.ts, self.seis)
        self.ts = time_points
        self.seis  = f(time_points)
        return self

    def plot(self):
        plt.plot(self.ts, self.seis)
        plt.show()


class FrechetDerivatives:
    def __init__(self, lcm0, time_points):
        from numpy import array

        #plt.figure(995)
        #plt.plot(lcm0.vslay, color='black', lw=5, alpha=0.4)
        #plt.show()

        sss0 = SSsyn(lcm0).construct_spike_train().convolve_with_wavelet().normalize().interpolate_to_time_points(time_points)

        nt = len(sss0.seis)
        nz = len(lcm0.zmid)

        self.nt = nt
        self.nz = nz

        FD = zeros(nt*nz).reshape(nt,nz)


        for iz, _ in enumerate(lcm0.zmid):
            lcm = LayerCakeModel(zs=lcm0.zs, vsh=lcm0.vsh, rho=lcm0.rho).layerize()
            lcm.vslay = lcm0.vslay.copy()
            lcm.vslay[iz:] = array(lcm0.vslay[iz:]) * (1. + 0.01)
            sss = SSsyn(lcm).construct_spike_train().convolve_with_wavelet().normalize().interpolate_to_time_points(time_points)
            damp = +sss.seis - sss0.seis

            #if iz == 150:
                #plt.figure(999)
                #plt.plot(lcm0.zmid,lcm.vslay,'blue')
                #plt.plot(lcm0.zmid,lcm0.vslay,'red')
                #plt.plot(sss.seis)
                #plt.plot(sss0.seis)
                #plt.show()

            if lcm.ztop[iz] > 0:
                FD[:,iz] = damp #* hanning(len(damp))

        self.FD = FD

    def plot(self):
        plt.imshow(self.FD, aspect = 'auto')
        plt.colorbar()
        plt.show()


def make_example_trace():
    import sys

    sys.path.append('..')
    from SSsyn import LayerCakeModel
    from minos.minos import Minos

    path = '/Users/mancinelli/Desktop/SS_Precursors/'
    starting_model = 'minos.min.anelas.dunite.ColoradoPlateau.pickle'

    # Load cardfile and construct layer_cake_model
    #
    cf = Minos().load(path + starting_model).cardfile
    cf.df['Z'] = 6371. - cf.df.Radius / 1000.
    df = cf.df.query('Z <= 500').sort_values('Z')
    df.Vsh = df.Vsh / 1000.
    df.Density = df.Density / 1000.

    df.Vsh[(df.Z > 160)] = df.Vsh[(df.Z > 160)] * 0.96

    df.Vsh[(df.Z > 180)] = df.Vsh[(df.Z > 180)] * 1.04

    lcm = LayerCakeModel(df.Z.tolist(), df.Vsh.tolist(), df.Density.tolist()).layerize()
    sss = SSsyn(lcm).construct_spike_train().convolve_with_wavelet().normalize()
    return sss

if __name__ == "__main__":
    sss = make_example_trace()
    sss.plot()
