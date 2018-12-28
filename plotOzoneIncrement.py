import matplotlib
matplotlib.use('Agg')
import os,h5py
import numpy as np
from lib.maps import plotMap 
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

if __name__ == "__main__":
    """
    # Test for dobson unit...seems a bit high, but ballpark - what we're looking for.
    P, T, Q, CO2, O3, GasUnits = readProfileH5()
    for i in range(0,6):
        print(P[i,0:98].max())
        dob = dobson(P[i,0:97], O3[i,0:97])
        print(dob)
    """
    testFile = os.path.join('/discover','nobackup','bkarpowi', 'controlIncrements','x35_control.xinc.eta.20180709_23z.nc4')
    d = Dataset(testFile)
    for dd in list(d.variables.keys()):
        print(dd)

    lon = d.variables['lon']
    lat = d.variables['lat']
    ppmvOzone = d.variables['ozone']
    lon2d, lat2d = np.meshgrid(lon,lat)   
    #(1, 72, 361, 576) 
    """
    dobsonPixels = np.zeros([ppmvOzone.shape[2], ppmvOzone.shape[3]])
    for i in range(0,ppmvOzone.shape[2]):
        print("{} of {}".format(i,ppmvOzone.shape[2]))
        for j in range(0,ppmvOzone.shape[3]):
            dobsonPixels[i,j] = dobson(d.variables['lev'],ppmvOzone[0,:,i,j])
    print(ppmvOzone.shape) 
    print(lon.shape)
    print(lat.shape)
    print(lon2d.shape)
    print(lat2d.shape)
    
    plotMap(lat2d.flatten(), lon2d.flatten(), dobsonPixels.flatten(), 'Dobson Unit', 'whir', units ='Dobson Units')
    """
    for l in range(0,d.variables['lev'].shape[0]): 
        plotMap(lat2d.flatten(), lon2d.flatten(), ppmvOzone[0,l,:,:].flatten(), 'PPMV{:07.3f}'.format(d.variables['lev'][l]), 'PPMV{:07.3f}'.format(d.variables['lev'][l]), units ='ppmv')
