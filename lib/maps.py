import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors
import cartopy.crs as ccrs

import cartopy.feature as cfeature


def plotMap(lat, lon, data, title, graphicName, plotRange = None, units='', colorScheme = 'viridis', projectionSelected = ccrs.PlateCarree()):
    fig, ax1 = plt.subplots(nrows=1, figsize=(16,9), subplot_kw={'projection': ccrs.PlateCarree()}) 
    ax1.coastlines(resolution='110m', color='black') # plot coastlines on map
    ax1.set_global() # show the whole earth
    ax1.gridlines() # put lat/lon markers

    # if a plot range isn't specified, make it the max and min of the data.
    if(plotRange): pass 
    else: plotRange = ( data.min(), data.max() )

    # come up with normalization for data.
    norm = matplotlib.colors.Normalize(vmin=plotRange[0], vmax=plotRange[1]) #Nomalize color range to match set data range
    # plot scatter plot with data according to colormap and normalization specified
    ax1.scatter(lon, lat, c=data, cmap = colorScheme, norm=norm, s=2, marker="o", edgecolors='none', transform=ccrs.PlateCarree())

    statString = "Mean={:10.4f} Standard Deviation={:10.4f} Min={:10.4f} Max={:10.4f}  ".format( data.mean(), data.std(ddof=1), data.min(), data.max() )
    plt.title(statString, fontsize=12)

    fig.suptitle(title, fontsize=18) 
    bar = plt.cm.ScalarMappable(cmap=colorScheme, norm=norm)
    bar._A = []
    cbar = plt.colorbar(bar, ax=ax1)
    cbar.set_label(units)
    #plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    fig.savefig(graphicName+'.png')
    plt.close(fig)
