import matplotlib.pyplot as plt
plt.style.use('ggplot')
from numpy import hstack, identity, vstack, zeros, matmul, array, arange, outer, shape, argmax, exp
import SSsyn
from numpy.linalg import inv

from scipy.interpolate import interp1d

sigma_pv=0.003
sigma_ss=0.04

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
        self.Mss = TOP
        self.Mrp = BOT

    def plot(self):
        plt.imshow(self.M, aspect="auto", interpolation = "nearest")
        plt.show()

class ModelVector:
    def __init__(self, zs, dvs_vs, reference_vs_model=None, sigma_mantle=0.20):
        self.zs = zs.copy()
        self.dvs_vs = dvs_vs.copy()
        self.m = array(len(dvs_vs))
        self.m = dvs_vs.copy()
        self.C = identity(len(dvs_vs))
        self.reference_vs_model=reference_vs_model

        sigma_crust = 0.02

        z1 = 50
        z2 = 55
        z3 = 280
        z4 = 290

        s1 = sigma_crust
        s2 = sigma_crust
        s3 = sigma_mantle
        s4 = sigma_mantle

        fsigma = interp1d([z1,z2,z3,z4],[s1,s2,s3,s4],'linear',fill_value='extrapolate')

        for ii in range(len(zs)):
            for jj in range(len(zs)):
                    if abs(ii - jj) == 0:
                        if zs[ii]<50:
                            self.C[ii, jj] = sigma_crust**2
                        else:
                            sigma = fsigma(zs[ii])
                            self.C[ii, jj] = (1.0*sigma) ** 2
                            if abs(ii - jj) == 1:
                                self.C[ii, jj] = (0.50 * sigma) ** 2

        self.Cinv = inv(self.C)

        print(self.Cinv)

    def calculate_difference_between_layers(self):
        model = array(self.reference_vs_model) * (1. + array(self.m))

        self.mdiff = zeros(len(self.m))

        for ii, each in enumerate(model):
            if ii == 0:
                self.mdiff[ii] = 0.0
            else:
                a = (model[ii] - model[ii-1])
                b = (self.reference_vs_model[ii] - self.reference_vs_model[ii-1])
                self.mdiff[ii] = (a - b)  /  self.reference_vs_model[ii]

        self.mdiff *= 100.0  # convert to percent perturbations

        return self


    def plot(self):
        plt.figure(8345)
        #plt.plot(model)
        #plt.plot(self.reference_vs_model)
        plt.plot(self.mdiff*10.0, label = 'mdiff')
        plt.plot(self.m,'k',label='m vector')
        plt.legend(loc='best')
        plt.show()

        return self

class DataVector:
    def __init__(self, timeseries, phasevelocity_perturbations, number_of_model_parameters, type):
        nts= len(timeseries)
        nd = len(timeseries) + len(phasevelocity_perturbations)
        self.d = zeros(nd)

        self.d[:nts] = timeseries
        self.d[nts:nts+len(phasevelocity_perturbations)] = phasevelocity_perturbations

        self.C = identity(nd) * (sigma_pv)**2

        for ii in range(nts):
            self.C[ii,ii] = (sigma_ss)**2

        self.Cinv = inv(self.C)

        if type == 'rayleigh':
            self.Cinv = self.Cinv[nts:nts+len(phasevelocity_perturbations):,nts:nts+len(phasevelocity_perturbations)]
            self.d = self.d[nts:nts+len(phasevelocity_perturbations)]
            print('shape(self.Cinv) = ', shape(self.Cinv))
        elif type == 'ss':
            self.Cinv = self.Cinv[:nts,:nts].copy()
            self.d = self.d[:nts]

class GMatrix:
    def __init__(self, periods, lcm, time_points, number_of_model_parameters, type='full'):
        if type == 'full':
            A = KernelMatrix(periods=periods, lcm=lcm, time_points=time_points).M
        elif type == 'rayleigh':
            A = KernelMatrix(periods=periods, lcm=lcm, time_points=time_points).Mrp
        elif type == 'ss':
            A = KernelMatrix(periods=periods, lcm=lcm, time_points=time_points).Mss

        B = ScalingMatrix(len(lcm.zmid)).M.copy()
        G = matmul(A, B)

        self.G = G

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

        fig=plt.figure(1,(9,9))


        def _plot_iter(label):

            plt.figure(1)

            ax = fig.add_axes([0.10, 0.10, 0.30, 0.8])

            plt.plot(array(reference_model.vslay) * (1. + array(cmv.dvs_vs)), lcm.zmid, lw=3, label = label)
            plt.ylim(max(lcm.zmid), 0)
            plt.ylabel('Depth (km)')
            plt.xlabel('Vs (km/s)')

            #plt.subplot(2, 2, 2)
            #a = list(periods)
            #b = dd1
            #plt.plot(a, b, 'o', label=label)
            #plt.ylabel('Fractional Residual')
            #plt.xlabel('Period (s)')

            #plt.subplot(2, 2, 3)
            #plt.plot(times, dd2, '-', label=label)
            #plt.ylabel('Fractional Residual')
            #plt.xlabel('Time from SS (s)')

            sse = matmul( matmul(dd1.T, dv1.Cinv), dd1 )

            print('iteration, SSE = ', iteration, sse)

        #Set up path and starting model
        #
        path = '/Users/mancinelli/Desktop/SS_Precursors/'
        starting_model = 'minos.min.anelas.dunite.ColoradoPlateau.pickle'


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
        pv = o.PhaseVelocities(37, -110, "Cloud").compute_perturbations(starting_model)

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

        #Set up starting model vector
        reference_model = lcm
        lcm = lcm.compute_perturbations(reference_model)
        smv = ModelVector(lcm.zmid, lcm.dvs_vs, reference_model.vslay).calculate_difference_between_layers()


        #Set up G matrix and data vector
        gm1 = GMatrix(periods, lcm, times, len(smv.dvs_vs),'rayleigh')
        dv1 = DataVector(dseis,dc_cs, len(smv.dvs_vs),'rayleigh')
        gm2 = GMatrix(periods, lcm, times, len(smv.dvs_vs), 'ss')
        dv2 = DataVector(dseis, dc_cs, len(smv.dvs_vs), 'ss')


        #calculate current model vector
        cmv = smv
        cmvlast = cmv

        dpred1 = matmul(gm1.G, cmv.m)
        dd1 = dv1.d - dpred1

        dpred2 = matmul(gm2.G, cmv.mdiff)
        dd2 = dv2.d - dpred2


        iteration = -1
        _plot_iter('Starting Model')

        #vslay0 = lcm.vslay.copy()

        #iterate

        niter = 1

        for iteration in range(niter):
            # Set up starting model vector
            #reference_model = lcm
            #lcm = lcm.compute_perturbations(reference_model)
            #smv = ModelVector(lcm.zmid, lcm.dvs_vs, reference_model.vslay).calculate_difference_between_layers()
            dm = array(cmv.m)

            for case in ['b','a']:

                if case == 'a':
                    gm1 = GMatrix(periods, lcm, times, len(cmv.dvs_vs), 'rayleigh')
                    dv1 = DataVector(dseis, dc_cs, len(smv.dvs_vs), 'rayleigh')
                    dpred1 = matmul(gm1.G, cmv.m)
                    dd1 = dv1.d - dpred1



                    TMP = matmul( matmul( gm1.G.T, dv1.Cinv ), gm1.G ) + smv.Cinv
                    TERM1 = inv(TMP)

                    TERM2 =  matmul( matmul( gm1.G.T, dv1.Cinv ), dd1 )
                    TERM3 =  1.0*matmul( smv.Cinv, array(cmv.m) - array(cmvlast.m))
                    TERM4 =  1.0*matmul( smv.Cinv, array(cmv.m) - array(smv.m))

                    dm = matmul(TERM1, TERM2-TERM4)

                    cmv = ModelVector(lcm.zmid, cmv.m + dm, reference_model.vslay)

                    cmvlast = cmv.calculate_difference_between_layers()

                    #plt.figure(12345)
                    #plt.plot(dv.d)
                    #plt.plot(dpred,'k')
                    #plt.show()

                    lcm=lcm.compute_layers_from_model_parameters(dm)

                else:
                    if iteration > 1:
                        continue


                    gm2 = GMatrix(periods, lcm, times, len(cmv.dvs_vs), 'ss')
                    dv2 = DataVector(dseis,dc_cs, len(smv.dvs_vs), 'ss')

                    cmv = cmv.calculate_difference_between_layers()


                    dpred2 = matmul(gm2.G, cmv.mdiff)
                    dd2 = dv2.d - dpred2

                    from sklearn.linear_model import Lasso
                    for alpha in [5e-8]:
                        clf = Lasso(alpha=alpha)
                        clf.fit(gm2.G,dd2)
                        dmdiff = clf.coef_

                    for ii in range(len(dmdiff)):
                    #ii = argmax(abs(dmdiff))


                    #lcm.vslay[ii:] = lcm.vslay[ii:] * (1. + dmdiff[ii]/100.0)

                        cmv.m[ii:]  = cmv.m[ii:] + dmdiff[ii] / 100.0
                        dm[ii:] = dm[ii:] + dmdiff[ii] / 100.0


                    cmv = ModelVector(lcm.zmid, cmv.m, reference_model.vslay).calculate_difference_between_layers()
                    #cmv.mdiff *=0.0
                    #cmv.mdiff[ii] = dmdiff[ii]

                    lcm = lcm.compute_layers_from_model_parameters(dm)

                    #plt.figure(98765)
                    #plt.plot(cmv.m)
                    #
                    # plt.show()


            dpred2 = matmul(gm2.G, cmv.mdiff)
            dd2 = dv2.d - dpred2

            dpred1 = matmul(gm1.G, cmv.m)
            dd1 =  dv1.d - dpred1

            _plot_iter('Final Model')

        plt.legend(loc='best', ncol=1)

        sss_final = SSsyn.SSsyn(lcm).construct_spike_train().convolve_with_wavelet().normalize().interpolate_to_time_points(times)

        plt.subplot(2, 2, 2)
        plt.plot(list(periods), dv1.d, 'ok')
        yerr = dv1.d*0.0 + sigma_pv
        plt.errorbar(list(periods), dv1.d, yerr=yerr,
                         elinewidth=3, label='Observations (1$\sigma$)', color='black', capsize=5, linestyle=' ', capthick=3,
                         zorder=1)
        plt.ylabel('dc/c')
        plt.plot(list(periods), dpred1, 'or', label='Model Prediction')
        plt.xlabel('Period (s)')
        plt.xlim([20,160])
        plt.legend(loc='best', ncol=1)

        plt.subplot(2, 2, 4)
        plt.plot(times, dv2.d, '-k', label='Synthetic Observation')
        plt.ylabel('SS Precursor Amplitude')
        plt.xlabel('Time from SS (s)')
        plt.plot(times, sss_final.seis, '--r',label='Model Prediction')
        plt.xlim([-120,-45])
        plt.ylim([-0.05,0.05])
        plt.legend(loc='best', ncol=1)
        plt.savefig('img/jtinv_coloplat.eps')


if __name__ == "__main__":
    Inversion()
