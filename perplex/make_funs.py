#!/usr/bin/env python
#
import pandas as pd
import subprocess
from numpy import arange, array, nan
from conversions import pressure
from scipy.interpolate import interp2d
import pickle
from glob import glob

table_file = glob('*_1.tab')[0]
print('table_file = ', table_file)

df = pd.read_csv(table_file, skiprows=13, header=None, delim_whitespace=True, names = ['rho','vp','vs'])

def cut_value_from_line(lineno=5):
    output = subprocess.check_output("sed -n -e %sp *.tab" % lineno, shell=True)
    return float(str(output).strip("\\n'").strip("b'").strip())


p1=cut_value_from_line(5)
dp=cut_value_from_line(6)
np=int(cut_value_from_line(7))

t1=cut_value_from_line(9)
dt=cut_value_from_line(10)
nt=int(cut_value_from_line(11))

ps = arange(np)*dp + p1
ts = arange(nt)*dt + t1 - 273.15

rhos = array(df.rho.tolist()).reshape(nt,np)
vps = array(df.vp.tolist()).reshape(nt,np)
vss = array(df.vs.tolist()).reshape(nt,np)

rhos[rhos<10] = nan
vps[vps>20]   = nan
vss[vss>20]   = nan

subprocess.call("ln -sf ../conversions .", shell=True)


zs = pressure.PtoZ(ps*10e-5)

frho=interp2d(ts,zs,rhos.T / 1000.0)
fvp =interp2d(ts,zs,vps.T)
fvs =interp2d(ts,zs,vss.T)


funs = [frho,fvp,fvs]
fout = open('funs.pickle','wb')
pickle.dump(funs, fout)
fout.close()

print('funs written to funs.pickle')
