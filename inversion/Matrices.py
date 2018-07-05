import matplotlib.pyplot as plt
plt.style.use('ggplot')
from numpy import hstack, identity, vstack, zeros, matmul, array, arange, outer
import SSsyn
from numpy.linalg import inv

class ScalingMatrix:
    def __init__(self,nz):

        #sv, sh, pv, ph

        scalings = [1.0, 1.0, 0.55, 0.55]

        def mat(nz,scaling):
            return identity(nz)*scaling

        mat0 = mat(nz, scalings[0])
        mat1 = mat(nz, scalings[1])
        mat2 = mat(nz, scalings[2])
        mat3 = mat(nz, scalings[3])

        self.M = vstack((mat0,mat1,mat2,mat3))

class KernelMatrix:
    def __init__(self, periods, lcm, time_points):
        #Build SS part

        fd = SSsyn.FrechetDerivatives(lcm, time_points)
        B = fd.FD

        A = zeros(fd.nt * fd.nz).reshape(fd.nt, fd.nz)
        C=A.copy()
        D=A.copy()

        TOP = hstack((A,B,C,D))

        #Build RW part

        import RayleighPhase
        fd = RayleighPhase.FrechetDerivatives(periods=periods, zs=lcm.zmid, hs=lcm.h)

        E=fd.FDvsv
        F=fd.FDvsh
        G=fd.FDvpv
        H=fd.FDvph

        BOT = hstack((E, F, G, H))

        self.M = vstack((TOP,BOT))

    def plot(self):
        plt.imshow(self.M, aspect="auto", interpolation = "nearest")
        plt.show()

class ModelVector:



    def __init__(self, zs, dvs_vs, sigma=0.01):
        self.zs = zs.copy()
        self.dvs_vs = dvs_vs.copy()
        self.m = array(len(dvs_vs))
        self.m = dvs_vs.copy()
        self.C = identity(len(dvs_vs))

        for ii in range(len(zs)):
            for jj in range(len(zs)):
                    if abs(ii - jj) == 0:
                        self.C[ii, jj] = (1.0*sigma) ** 2
                    if abs(ii - jj) == 1:
                        self.C[ii,jj] = (2.0*sigma)**2
                    if abs(ii - jj) == 2:
                        self.C[ii,jj] = (0.000*sigma)**2
                    #if abs(ii - jj) == 3:
                    #    self.C[ii, jj] = (0.99 * sigma) ** 2

        self.Cinv = inv(self.C)


class DataVector:
    def __init__(self, timeseries, phasevelocity_perturbations, sigma_pv=0.025, sigma_ss=0.005):
        nts= len(timeseries)
        nd = len(timeseries) + len(phasevelocity_perturbations)
        self.d = zeros(nd)

        self.d[:nts] = timeseries
        self.d[nts:] = phasevelocity_perturbations

        self.C = identity(nd) * (sigma_pv)**2

        for ii in range(nts):
            self.C[ii,ii] = (sigma_ss)**2

        self.Cinv = inv(self.C)

class GMatrix:
    def __init__(self, periods, lcm, time_points):
        A = KernelMatrix(periods=periods, lcm=lcm, time_points=time_points).M
        B = ScalingMatrix(len(lcm.zmid)).M
        self.G = matmul(A, B)

    def plot(self):
        plt.imshow(self.G, aspect="auto", interpolation="nearest")
        plt.show()
        return self

class Inversion:

    def __init__(self):
        import sys
        sys.path.append('..')
        from SSsyn import LayerCakeModel
        from minos.minos import Minos
        from observations import observations as o

        plt.figure(1,(15,15))


        def _plot_iter(cmv, lcm, gm, dd, label):

            plt.subplot(2, 2, 1)
            plt.plot(array(reference_model.vslay) * (1. + array(cmv.dvs_vs)), lcm.zmid, lw=1)
            plt.ylim(max(lcm.zmid), 0)
            plt.ylabel('Depth (km)')
            plt.xlabel('Vs (km/s)')

            plt.subplot(2, 2, 2)
            a = list(periods)
            b = dd[-len(periods):]
            plt.plot(a, b, 'o', label=label)
            plt.ylabel('Fractional Residual')
            plt.xlabel('Period (s)')

            plt.subplot(2, 2, 3)
            plt.plot(times, dd[:-len(periods)], '-', label=label)
            plt.ylabel('Fractional Residual')
            plt.xlabel('Time from SS (s)')

            #plt.subplot(2,2,4)
            #plt.plot(array(cmv.m) - array(smv.m))
            #plt.ylabel('array(cmv.m) - array(smv.m)')

            sse = matmul( matmul(dd.T, dv.Cinv), dd )

            print('iteration, SSE = ', iteration, sse)

        #Set up path and starting model
        #
        path = '/Users/mancinelli/Desktop/SS_Precursors/'
        starting_model = 'minos.min.anelas.dunite.SlaveCraton.pickle'


        #Load cardfile and construct layer_cake_model
        #
        cf = Minos().load(path + starting_model).cardfile
        cf.df['Z'] = 6371. - cf.df.Radius/1000.
        df = cf.df.query('Z <= 500').sort_values('Z')
        df.Vsh = df.Vsh/1000.
        df.Density = df.Density/1000.
        lcm = LayerCakeModel(df.Z.tolist(), df.Vsh.tolist(), df.Density.tolist()).layerize()


        #Import phase velocity observations
        #
        pv = o.PhaseVelocities(64.5, -110, "Ma").compute_perturbations(starting_model)

        periods = pv.perturbations.keys()
        dc_cs = list(pv.perturbations.values())


        #Construct reference SS precursor waveform
        #
        sss = SSsyn.SSsyn(lcm).construct_spike_train().convolve_with_wavelet().normalize()

        times = arange(-250, 0, 1.0)


        tmp = SSsyn.make_example_trace()
        seis = tmp.interpolate_to_time_points(times).seis

        sss = sss.interpolate_to_time_points(times)


        #Caluclate optimal scaling
        d1 = sss.seis.copy()
        d2 = seis.copy()

        tmp1 = outer(d1, d1) + identity(len(d1))*0.001
        tmp2 = inv(tmp1)
        tmp3 = matmul(tmp2, d1)
        optscl = matmul(tmp3, d2)
        dseis = seis - sss.seis * optscl

        #Set up G matrix and data vector
        gm = GMatrix(periods, lcm, time_points=times)
        dv = DataVector(dseis,dc_cs)

        #Set up starting model vector
        reference_model = lcm
        lcm = lcm.compute_perturbations(reference_model)
        smv = ModelVector(lcm.h,lcm.dvs_vs)

        #calculate current model vector
        cmv = smv
        cmvlast = cmv

        dpred = matmul(gm.G, cmv.m)
        dd = dv.d - dpred

        iteration = -1
        _plot_iter(cmv, lcm, gm, dd, '-1')

        #vslay0 = lcm.vslay.copy()

        #iterate

        niter = 10

        for iteration in range(niter):

            #if iteration%2 == 0:
            #    dv = DataVector(dseis, dc_cs, sigma_ss=10.0)
            #else:
            #    dv = DataVector(dseis, dc_cs, sigma_pv=10.0)


            TMP = matmul( matmul( gm.G.T, dv.Cinv ), gm.G ) + smv.Cinv
            TERM1 = inv(TMP)

            TERM2 =  matmul( matmul( gm.G.T, dv.Cinv ), dd )
            TERM3 =  0.5*matmul( smv.Cinv, array(cmv.m) - array(cmvlast.m))
            TERM4 =  matmul( smv.Cinv, array(cmv.m) - array(smv.m))

            dm = matmul(TERM1, TERM2-TERM4)

            #if iteration % 2 == 1:
            #    from numpy import diff, argmax
            #    from numpy import hanning
            #    ii = argmax(abs(diff(dm) * hanning(len(diff(dm))) ))
            #    djump = dm[ii] * 1.0
            #    dm2 = dm.copy()*0.0
            #    dm2[ii:] = dm2[ii:] - djump/2.
            #    dm2[:ii] = dm2[:ii] + djump / 2.
            #    #plt.figure(1717)
            #    #plt.plot(dm)
            #    #plt.plot(dm2,'black')
            #    #plt.show()
            #    dm = dm2

            #mnew = cmv.dvs_vs + dm
            #mnew = 0.0 + dm

            #from sklearn.linear_model import Lasso
            #clf = Lasso(alpha=1.0)
            #clf.fit(gm.G,dv.d)
            #dm = clf.coef_
            #print(dm)


            cmvlast = cmv

            cmv = ModelVector(lcm.h, cmv.m + dm)

            dpred = matmul(gm.G, cmv.m)
            dd =  dv.d - dpred

            lcm=lcm.compute_layers_from_model_parameters(dm)

            gm = GMatrix(periods, lcm, time_points=times)

            if iteration%1 == 5 or iteration == niter - 1 :
                _plot_iter(cmv, lcm, gm, dd, iteration)



        plt.legend(loc='best', ncol=5)
        plt.savefig('mypost2.eps')



if __name__ == "__main__":
    Inversion()
