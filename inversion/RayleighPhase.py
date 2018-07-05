
import pandas as pd

from scipy.interpolate import interp1d

from numpy import arange, zeros

from matplotlib import pylab as plt

def _read_derivative_file_and_interpolate(period):

    def __fun():
        return pd.read_csv(fname, delim_whitespace=True, names=['Radius', 'Vsv', 'Vpv', 'Vsh', 'Vph', 'eta', 'rho'])

    try:
        fname = '/Users/mancinelli/Desktop/SS_Precursors/inversion/frechet/cvfrechet_%ds' % int(period)
        df = __fun()
    except FileNotFoundError:
        fname = '/Users/mancinelli/Desktop/SS_Precursors/inversion/frechet/cvfrechet_%ds' % int((period + 1.))
        df = __fun()

    df['Z'] = 6371. - df.Radius / 1000.
    df.sort_values('Z')

    zs  = df.Z.tolist()
    vsv = df.Vsv.tolist()
    vsh = df.Vsh.tolist()
    vpv = df.Vpv.tolist()
    vph = df.Vph.tolist()

    fun_vsv = interp1d(zs, vsv)
    fun_vsh = interp1d(zs, vsh)
    fun_vpv = interp1d(zs, vpv)
    fun_vph = interp1d(zs, vph)

    return fun_vsv, fun_vsh, fun_vpv, fun_vph


class FrechetDerivatives:
    def __init__(self, periods=None, zs=None, hs=None):
        from numpy import sum

        if zs == None:
            zs = arange(0,410,1)
        if periods == None:
            periods = [25, 30, 40, 50, 60, 75, 80, 100, 125, 133, 140, 150, 180, 200]

        np = len(periods)
        nz = len(zs)

        FDvsv = zeros(np*nz).reshape(np,nz)
        FDvsh = zeros(np*nz).reshape(np,nz)
        FDvpv = zeros(np*nz).reshape(np,nz)
        FDvph = zeros(np*nz).reshape(np,nz)

        for ip, period in enumerate(periods):
            fun_vsv, fun_vsh, fun_vpv, fun_vph = _read_derivative_file_and_interpolate(period)
            FDvsv[ip,:] = fun_vsv(zs) * 1000.
            FDvsh[ip,:] = fun_vsh(zs) * 1000.
            FDvpv[ip,:] = fun_vpv(zs) * 1000.
            FDvph[ip,:] = fun_vph(zs) * 1000.

        #multipy by thickness of each layer
        for iz, h in enumerate(hs):
            assert(h>0.0)

            FDvsv[:,iz] =  sum(h*FDvsv[:,iz:],axis=1)
            FDvsh[:,iz] =  sum(h*FDvsh[:,iz:],axis=1)
            FDvpv[:,iz] =  sum(h*FDvpv[:,iz:],axis=1)
            FDvph[:,iz] =  sum(h*FDvph[:,iz:],axis=1)

            if zs[iz] < 0:
                FDvsv[:, iz] *= 0.
                FDvsh[:, iz] *= 0.
                FDvpv[:, iz] *= 0.
                FDvph[:, iz] *= 0.


        self.FDvsv = FDvsv
        self.FDvsh = FDvsh
        self.FDvpv = FDvpv
        self.FDvph = FDvph

    def plot(self):
        plt.imshow(self.FDvsv, aspect = 'auto', interpolation = 'nearest')
        plt.colorbar()
        plt.show()


if __name__ == "__main__":
    FrechetDerivatives().plot()