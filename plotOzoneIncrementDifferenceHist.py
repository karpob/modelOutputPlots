import matplotlib
matplotlib.use('Agg')
import os, h5py, argparse, sys
import numpy as np
# add path so script will work outside pwd.
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from lib.graphics.maps import plotMapHist 
from netCDF4 import Dataset
def readProfileH5():
    rttovPath = '/discover/nobackup/bkarpowi/rt/rttov12_gcc7.2_openmp/'
    h5 = h5py.File(os.path.join(rttovPath,'rttov_test/profile-datasets-hdf/standard101lev_allgas.H5'))
    nprofiles = 6
    groups = []
    for i in range(1,nprofiles+1):
        groups.append("{:04d}".format(i))
    P = np.zeros([nprofiles,101])
    T = np.zeros([nprofiles,101])
    Q = np.zeros([nprofiles,101])
    CO2 = np.zeros([nprofiles,101])
    O3 = np.zeros([nprofiles,101])
    for i,g in enumerate(groups):
        P[i,:] = np.asarray(h5['PROFILES'][g]['P'])
        Q[i,:] = np.asarray(h5['PROFILES'][g]['Q'])
        T[i,:] = np.asarray(h5['PROFILES'][g]['T'])
        CO2[i,:] = np.asarray(h5['PROFILES'][g]['CO2'])
        O3[i,:] = np.asarray(h5['PROFILES'][g]['O3'])
        GasUnits = int(np.asarray(h5['PROFILES'][g]['GAS_UNITS']))

    return P, T, Q, CO2, O3, GasUnits

def dobson( p, o3 ):
    
    """
    -----------------------------------------------------------

    Purpose:
        Compute column ozone ammount in dobson units (DU), from a
        profile of O3 mixing ratio vs. pressure.
        1 DU = 2.69e16 molecules cm-2

    Input:
        p.......pressure, mb
        o3......ozone mixing ratio,  ppmv

    Output:
        dobs....column o3 abundance, DU
 
    Source: Mark Hervig IDL code from GATS website http://gwest.gats-inc.com/software/dobson.pro
    converted to python by Bryan Karpowicz.
    -----------------------------------------------------------
    """
    # some constants

    g        = 9.8         # gravity, m/s2
    avagadro = 6.02252e23  # molec/mol
    mdry     = 0.028964    # molec. wt. of dry air, kg/mol

    const    = 0.01 * avagadro / (g * mdry)

    dobs     = 0.0

    # sum o3 over height

    for j in range( 0, p.shape[0] -1 ):
        term = 0.5*( o3[j] + o3[j-1] ) *1.0e-6 * abs( p[j-1] - p[j] ) * const/2.69e16
        dobs = dobs + term 

    return dobs
def getLonLatLev(filename):
    d = Dataset(filename)

    lon = d.variables['lon']
    lat = d.variables['lat']
    lev = d.variables['lev']
    lon2d, lat2d = np.meshgrid(lon,lat)   
    return lon2d, lat2d, lev

def getVar(filename, vname):
    
    d = Dataset(filename)
    if (vname not in list(d.variables.keys())):
        print("Can't find {} try one of these instead:".format(vname))
        for dd in list(d.variables.keys()):
            print(dd)
    varVal = d.variables[vname]

    return varVal 
def getT(filename):
    
    d = Dataset(filename)
    for dd in list(d.variables.keys()):
        print(dd)

    lon = d.variables['lon']
    lat = d.variables['lat']
    Tv = d.variables['tv']
    w = np.asarray(d.variables['sphu'])*1000.0
    T = Tv - (w/6.0)
    return T

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'plot differences between two increments.')
    parser.add_argument('--control', help = 'control data file', required = True, dest = 'control')
    parser.add_argument('--experiment',help = 'experiment data file', required = True, dest='experiment')
    parser.add_argument('--output',help = 'output path', required = False, dest='output', default ='')
    a = parser.parse_args()


    # use given path and grab all nc diag files with instrument name in them.
    
    lon2d, lat2d, lev = getLonLatLev(a.control) 
    cntlOzone = getVar(a.control,'ozone')
    expOzone = getVar(a.experiment,'ozone') 

    cntlTv = getVar(a.control,'tv')
    expTv = getVar(a.experiment,'tv') 

    cntlTs = getVar(a.control,'ts')
    expTs = getVar(a.experiment,'ts') 

    cntlH = getVar(a.control,'sphu')
    expH = getVar(a.experiment,'sphu') 

    diffString = '{} - {}'.format(os.path.basename(a.experiment), os.path.basename(a.control))
    oo = a.output
    if( not os.path.exists(oo) ): os.makedirs(oo)

    plotMapHist(lat2d.flatten(), lon2d.flatten(), expTs[0,:,:].flatten()-cntlTs[0,:,:].flatten(),'diff-Ts: '+diffString, os.path.join(oo,'diff-Ts'), units ='Kelvin',plotRange=[-4,4])
    for l in range(0, lev.shape[0]): 
        plotMapHist(lat2d.flatten(), lon2d.flatten(), expOzone[0,l,:,:].flatten()-cntlOzone[0,l,:,:].flatten(), 'diff-PPMV Level {:07.3f}'.format(lev[l]), os.path.join(oo,'diff-PPMV{:07.3f}'.format(lev[l])), units ='ppmv')
        plotMapHist(lat2d.flatten(), lon2d.flatten(), expTv[0,l,:,:].flatten()-cntlTv[0,l,:,:].flatten(), 'diff-T Level {:07.3f}'.format(lev[l]), os.path.join(oo,'diff-T-Level_{:07.3f}'.format(lev[l])), units ='Kelvin')
        plotMapHist(lat2d.flatten(), lon2d.flatten(), expH[0,l,:,:].flatten()-cntlH[0,l,:,:].flatten(), 'diff-SpHu Level {:07.3f}'.format(lev[l]), os.path.join(oo,'diff-SPHU-Level_{:07.3f}'.format(lev[l])), units ='kg/kg')
