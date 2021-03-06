import matplotlib
matplotlib.use('Agg')
import os, h5py, argparse, glob, sys, socket
from datetime import timedelta, date, datetime
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from matplotlib import pyplot as plt
# add path so script will work outside pwd.
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from lib.graphics.maps import plotMapHist 

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
def dateRange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)
 
def getFiles(start, end, opsPath, experimentName):
    files = []
    # pull out years, etc from args.
    startYear, startMonth, startDay, startHour = int(start[0:4]), int(start[4:6]), int(start[6:8]), int(start[8:10])
    endYear, endMonth, endDay, endHour = int(end[0:4]), int(end[4:6]), int(end[6:8]), int(end[8:10])

    # basepath
    pathInit =  os.path.join(opsPath, experimentName, 'ana')
    
    # get start/end date, put in datetime
    startDate = date(startYear, startMonth, startDay)
    endDate = date(endYear, endMonth, endDay)

    for today in dateRange(startDate, endDate):
        for hour in ['03','09','15','21']:
            path = os.path.join(pathInit, today.strftime("Y%Y/M%m"))

            if not os.path.exists(path): print(path +'does not exist.' )
            else:
                if( len( glob.glob(path+'/'+experimentName+'.xinc09.eta.'+today.strftime("%Y%m%d")+'_'+hour+'z.nc4' ) ) > 0): 
                    files.append( glob.glob(path+'/'+experimentName+'.xinc09.eta.'+today.strftime("%Y%m%d")+'_'+hour+'z.nc4' )[0]   )
    return files    
def processFilePair(controlFile, experimentFile, mapsOn, diffsOnly):

    experimentName = os.path.basename( experimentFile ).split('.')[0]
    tE = os.path.basename( experimentFile ).split('.')[3]

    controlName = os.path.basename( controlFile ).split('.')[0]
    tC = os.path.basename( controlFile ).split('.')[3]
    print(tE, tC)
    experimentTime = datetime(int(tE[0:4]), int(tE[4:6]), int(tE[6:8]), int(tE[9:11]))
    controlTime = datetime(int(tC[0:4]), int(tC[4:6]), int(tC[6:8]), int(tC[9:11]))
    
    if(experimentTime != controlTime):
        print("No Bueno! diferent times, I won't take the difference!")
        sys.exit()


    lon2d, lat2d, lev = getLonLatLev(controlFile) 
    cntlOzone = getVar(controlFile,'ozone')
    expOzone = getVar(experimentFile,'ozone') 

    cntlTv = getVar(controlFile,'tv')
    expTv = getVar(experimentFile,'tv') 

    cntlTs = getVar(controlFile,'ts')
    expTs = getVar(experimentFile,'ts') 

    cntlH = getVar(controlFile,'sphu')
    expH = getVar(experimentFile,'sphu') 

    dTs = np.asarray(expTs[0,:,:] - cntlTs[0,:,:])
    dOzone = np.asarray(expOzone[0,:,:] - cntlOzone[0,:,:])
    dTv = np.asarray(expTv[0,:,:] - cntlTv[0,:,:]) 
    dH = np.asarray(expH[0,:,:] - cntlH[0,:,:])
 
    if(mapsOn): 
        diffString = '{} - {}'.format(os.path.basename(experimentFile), os.path.basename(controlFile))
        oo = os.path.basename(experimentFile).split('.')[-2]
        if( not os.path.exists(oo) ): os.makedirs(oo)
        plotMapHist(lat2d.flatten(), lon2d.flatten(), dTs.flatten(),'diff-Ts: '+diffString, os.path.join(oo,'diff-Ts'), units ='Kelvin')
        plotMapHist(lat2d.flatten(), lon2d.flatten(), np.asarray(expTs).flatten(),'Ts ', os.path.join(oo,'Ts_exp'), units ='Kelvin')
        plotMapHist(lat2d.flatten(), lon2d.flatten(), np.asarray(cntlTs).flatten(),'Ts ', os.path.join(oo,'Ts_ctl'), units ='Kelvin')
        for l in range(0, lev.shape[0]): 
            plotMapHist(lat2d.flatten(), lon2d.flatten(), dOzone[l,:,:].flatten(), 'diff-PPMV Level {:07.3f}'.format(lev[l]), os.path.join(oo,'diff-PPMV{:07.3f}'.format(lev[l])), units ='ppmv')
            plotMapHist(lat2d.flatten(), lon2d.flatten(), dTv[l,:,:].flatten(), 'diff-T Level {:07.3f}'.format(lev[l]), os.path.join(oo,'diff-T-Level_{:07.3f}'.format(lev[l])), units ='Kelvin')
            #plotMapHist(lat2d.flatten(), lon2d.flatten(), dH[0,l,:,:].flatten(), 'diff-SpHu Level {:07.3f}'.format(lev[l]), os.path.join(oo,'diff-SPHU-Level_{:07.3f}'.format(lev[l])), units ='kg/kg')
            plotMapHist(lat2d.flatten(), lon2d.flatten(), np.asarray(expOzone[0,l,:,:]).flatten(), 'PPMV Level {:07.3f}'.format(lev[l]), os.path.join(oo,'exp-PPMV{:07.3f}'.format(lev[l])), units ='ppmv')
            plotMapHist(lat2d.flatten(), lon2d.flatten(), np.asarray(expTv[0,l,:,:]).flatten(), 'T Level {:07.3f}'.format(lev[l]), os.path.join(oo,'exp-T-Level_{:07.3f}'.format(lev[l])), units ='Kelvin')
            #plotMapHist(lat2d.flatten(), lon2d.flatten(), expH[0,l,:,:].flatten(), 'SpHu Level {:07.3f}'.format(lev[l]), os.path.join(oo,'exp-SPHU-Level_{:07.3f}'.format(lev[l])), units ='kg/kg')
            plotMapHist(lat2d.flatten(), lon2d.flatten(), np.asarray(cntlOzone[0,l,:,:]).flatten(), 'PPMV Level {:07.3f}'.format(lev[l]), os.path.join(oo,'cntl-PPMV{:07.3f}'.format(lev[l])), units ='ppmv')
            plotMapHist(lat2d.flatten(), lon2d.flatten(), np.asarray(cntlTv[0,l,:,:]).flatten(), 'T Level {:07.3f}'.format(lev[l]), os.path.join(oo,'cntl-T-Level_{:07.3f}'.format(lev[l])), units ='Kelvin')
            #plotMapHist(lat2d.flatten(), lon2d.flatten(), np.asarray(cntlH[0,l,:,:]).flatten(), 'SpHu Level {:07.3f}'.format(lev[l]), os.path.join(oo,'cntl-SPHU-Level_{:07.3f}'.format(lev[l])), units ='kg/kg')
 
    dTvMean = np.mean(np.mean(dTv,axis=1),axis=1)
    dTvStd = np.std(np.std(dTv,axis=1),axis=1)
    dTvMax = np.max(np.max(dTv,axis=1),axis=1)
    dTvMin = np.min(np.min(dTv,axis=1),axis=1)
    dOzoneMean = np.mean(np.mean(dOzone,axis=1),axis=1)
    dOzoneStd = np.std(np.std(dOzone,axis=1),axis=1)
    dOzoneMax = np.max(np.max(dOzone,axis=1),axis=1)
    dOzoneMin = np.min(np.min(dOzone,axis=1),axis=1)

    expTvMean = np.mean(np.mean(expTv[0,:,:],axis=1),axis=1)
    expTvStd = np.std(np.std(expTv[0,:,:],axis=1),axis=1)
    expTvMax = np.max(np.max(expTv[0,:,:],axis=1),axis=1)
    expTvMin = np.min(np.min(expTv[0,:,:],axis=1),axis=1)

    expOzoneMean = np.mean(np.mean(expOzone[0,:,:],axis=1),axis=1)
    expOzoneStd = np.std(np.std(expOzone[0,:,:],axis=1),axis=1)
    expOzoneMax = np.max(np.max(expOzone[0,:,:],axis=1),axis=1)
    expOzoneMin = np.min(np.min(expOzone[0,:,:],axis=1),axis=1)

    cntlTvMean = np.mean(np.mean(cntlTv[0,:,:],axis=1),axis=1)
    cntlTvStd = np.std(np.std(cntlTv[0,:,:],axis=1),axis=1)
    cntlTvMax = np.max(np.max(cntlTv[0,:,:],axis=1),axis=1)
    cntlTvMin = np.min(np.min(cntlTv[0,:,:],axis=1),axis=1)

    cntlOzoneMean = np.mean(np.mean(cntlOzone[0,:,:],axis=1),axis=1)
    cntlOzoneStd = np.std(np.std(cntlOzone[0,:,:],axis=1),axis=1)
    cntlOzoneMax = np.max(np.max(cntlOzone[0,:,:],axis=1),axis=1)
    cntlOzoneMin = np.min(np.min(cntlOzone[0,:,:],axis=1),axis=1)

    stats2dTv = []
    stats2dOzone = []
    stats2dTvExp = []
    stats2dOzoneExp = []
    stats2dTvCntl = []
    stats2dOzoneCntl = []
    for l in range(0, lev.shape[0]):
        stats2dTv.append( [experimentTime,l, lev[l], dTvMean[l], dTvStd[l], dTvMax[l], dTvMin[l]] ) 
        stats2dOzone.append( [experimentTime,l, lev[l], dOzoneMean[l], dOzoneStd[l], dOzoneMax[l], dOzoneMin[l]] )
        stats2dTvExp.append([experimentTime,l, lev[l], expTvMean[l], expTvStd[l], expTvMax[l], expTvMin[l]] ) 
        stats2dOzoneExp.append( [experimentTime,l, lev[l], expOzoneMean[l], expOzoneStd[l], expOzoneMax[l], expOzoneMin[l]] )
        stats2dTvCntl.append([experimentTime,l, lev[l], cntlTvMean[l], cntlTvStd[l], cntlTvMax[l], cntlTvMin[l]] ) 
        stats2dOzoneCntl.append( [experimentTime,l, lev[l], cntlOzoneMean[l], cntlOzoneStd[l], cntlOzoneMax[l], cntlOzoneMin[l]] )
    statsDiff1d = [[experimentTime, np.mean(dTs), np.std(dTs), dTs.min(), dTs.max()]]  
    statsExp1d = [[experimentTime, np.mean(np.asarray(expTs)), np.std(np.asarray(expTs)), np.min(np.asarray(expTs)), np.max(np.asarray(expTs))]]  
    statsCntl1d = [[experimentTime, np.mean(np.asarray(cntlTs)), np.std(np.asarray(cntlTs)), np.min(np.asarray(cntlTs)), np.max(np.asarray(cntlTs))]]
  
    return statsDiff1d, statsExp1d, statsCntl1d, stats2dTv, stats2dOzone, stats2dTvExp, stats2dOzoneExp, stats2dTvCntl, stats2dOzoneCntl


def plot1dTimeseries( stats1d, varName  ):
    statArray = np.asarray(stats1d)    
    timeIndex = statArray[:,0]
    print(varName,statArray)
    df = pd.DataFrame(statArray, index = timeIndex, columns = ['date', 'mean', 'std','min', 'max'])
    df.plot(y=['mean','std','min','max'])
    plt.savefig(varName+'_timeseries.png')
    print(varName,df)

def plot2dTimeseries( stats, varName, units):
    statArray = np.asarray(stats)
    timeIndex = statArray[:,0]
    df = pd.DataFrame(statArray, index = timeIndex, columns = ['date','levelIdx', 'level', 'mean','std', 'max', 'min'])
    print(df)
    # use pandas to group data by model level 
    grouped = df.groupby('levelIdx')
    
    # hard code levels, smarter way to count levels would be nice.     
    levelCount = 72
    f, ax = plt.subplots(nrows=levelCount, ncols=1, figsize=(14,3*levelCount) )

    # iterate through each pandas group and plot time series of statistics for that level
    i = 0
    for n, g in grouped:
        g.plot(ax = ax[i], y=['mean','max','min'],legend=False)
        dd = g['date'].values
        m = np.array(g['mean'].values, dtype = np.float64)
        std = np.array(g['std'].values, dtype = np.float64)
        ax[i].fill_between(dd, m+std, m-std, alpha=0.3, facecolor='grey')
        ax[i].set_title('Level {:4f} '.format(np.max(g['level'].values) ) )
        ax[i].set_ylabel( units )
        ax[i].legend(bbox_to_anchor= (1.1,1),loc='upper right')
        i+=1
    plt.tight_layout()
    plt.savefig(varName+'_stats_timeseries.png')



if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'plot differences between two increments.')
    parser.add_argument('--control', help = 'control experiemnt name', required = True, dest = 'control')
    parser.add_argument('--start', help = 'start dtg YYYYMMDDhh', required = True, dest = 'start')
    parser.add_argument('--end', help = 'end dtg YYYYMMDDhh', required = True, dest = 'end')
    parser.add_argument('--experiment',help = 'test experiment name', required = True, dest='experiment')
    parser.add_argument('--archive',help = 'archive path', required = False, dest='archive', default ='/archive/u/bkarpowi')
    parser.add_argument('--no-maps', help="turn off maps/histograms.", dest='maps', action='store_false' )
    parser.add_argument('--only-diffs', help="only differences", dest='diffs', action='store_false' )
    a = parser.parse_args()
    allTsStats = []
    allTvStats = []
    allOzoneStats = []

    expTsStats = []
    expTvStats = []
    expOzoneStats = []

    cntlTsStats = []
    cntlTvStats = []
    cntlOzoneStats = []

    controlFiles = getFiles(a.start, a.end, a.archive, a.control)

    experimentFiles = getFiles(a.start, a.end, a.archive, a.experiment)
    if(len(controlFiles)!=len(experimentFiles)):
        print("Do not Pass go. We don't have the right number of files.")
        sys.exit()
    if('discover' in socket.gethostname()):
        for f in experimentFiles:
            v="ssh dirac 'dmget "+f+"'"
            print(v)
            os.system(v)
        for f in controlFiles:
            v="ssh dirac 'dmget "+f+"'"
            print(v)
            os.system(v)
 
    for i,fc in enumerate(controlFiles):
        print("Processing {}".format(os.path.basename(fc).split('.')[-2]))
        surfaceT, surfaceTexp, surfaceTcntl, Tv, Ozone, expTv, expOzone, cntlTv, cntlOzone = processFilePair(fc, experimentFiles[i], a.maps, a.diffs)
        allTsStats.extend(surfaceT)
        allTvStats.extend(Tv)
        allOzoneStats.extend(Ozone)

        expTsStats.extend(surfaceTexp)
        expTvStats.extend(expTv)
        expOzoneStats.extend(expOzone)

        cntlTsStats.extend(surfaceTcntl)
        cntlTvStats.extend(cntlTv)
        cntlOzoneStats.extend(cntlOzone)

    plot1dTimeseries( allTsStats, 'diff-Ts' )
    plot2dTimeseries( allTvStats, 'diff-Tv', 'Kelvin' )
    plot2dTimeseries( allOzoneStats, 'diff-Ozone', 'PPMV' )

    plot1dTimeseries( expTsStats, 'exp-Ts' )
    plot2dTimeseries( expTvStats, 'exp-Tv', 'Kelvin' )
    plot2dTimeseries( expOzoneStats, 'exp-Ozone', 'PPMV' )

    plot1dTimeseries( cntlTsStats, 'cntl-Ts' )
    plot2dTimeseries( cntlTvStats, 'cntl-Tv', 'Kelvin' )
    plot2dTimeseries( cntlOzoneStats, 'cntl-Ozone', 'PPMV' )
 
    print("Done!") 
