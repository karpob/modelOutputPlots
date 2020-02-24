import matplotlib
matplotlib.use('Agg')
import cartopy.crs as ccrs
import os, h5py, argparse, sys
import numpy as np
# add path so script will work outside pwd.
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from lib.maps import plotMap, plotMapHist 
from netCDF4 import Dataset
def getLonLat(filename):
    d = Dataset(filename)

    lon = d.variables['lon']
    lat = d.variables['lat']
    lon2d, lat2d = np.meshgrid(lon,lat)   
    return lon2d, lat2d

def getVar(filename, vname):
    
    d = Dataset(filename)
    if (vname not in list(d.variables.keys())):
        print("Can't find {} try one of these instead:".format(vname))
        for dd in list(d.variables.keys()):
            print(dd)
    varVal = d.variables[vname]
    varVal = np.asarray(varVal)
    return varVal 

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'plot the val')
    parser.add_argument('--experiment', help = 'val file', required = True, dest = 'experiment')
    parser.add_argument('--control', help = 'val file', required = True, dest = 'control')
    parser.add_argument('--output',help = 'output path', required = False, dest='output', default ='')
    a = parser.parse_args()

 
    lon2d, lat2d = getLonLat(a.experiment) 
    valOzoneExperiment = getVar(a.experiment,'TO3')
    valOzoneControl = getVar(a.control,'TO3')
    idxValid = np.where( valOzoneExperiment < 1e15 )
    diff = valOzoneExperiment-valOzoneControl 
    oo = a.output
    
    experimentName = os.path.basename(a.experiment).split('.')[0]
    
    experiments = experimentName

    if( not os.path.exists(oo) and oo != '' ): os.makedirs(oo)
         

    plotMap( lat2d.flatten(),\
             lon2d.flatten(),\
             diff.flatten(),\
             os.path.split(a.experiment)[-1], os.path.join(oo, os.path.split(a.experiment)[-1]),\
             units ='Dobson Units',\
             projectionSelected = ccrs.Mollweide(central_longitude=180), plotRange=[-20,20] )
