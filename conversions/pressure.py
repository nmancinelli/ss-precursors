#!/usr/bin/env python
#
from numpy import log, exp
g = 9.8
beta = 8.7e-12
rho0 = 3350.

def PtoZ(Pgpa):
    """

    :param Pgpa: Pressure (GPa)
    :return: depth (km)
    """
    P = Pgpa * 1e9
    num = 1. - exp(-P*beta)
    den = g*beta*rho0
    return num/den / 1000.0

def ZtoP(zkm):
    """

    :param zkm: depth (km)
    :return: Pressure (GPa)
    """

    z=zkm*1000.0
    num = -log(1.-z*g*beta*rho0)
    den = beta
    return num/den * 1e-9

if __name__ == "__main__":
    from numpy import arange
    zs = arange(410)
    for z in zs:
        print(z, PtoZ(ZtoP(z)), ZtoP(z))