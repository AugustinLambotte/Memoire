import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy import interpolate
from datetime import date, timedelta
import cartopy.crs as ccrs
import os
import matplotlib.path as mpath
import cmocean
from scipy import signal

"""
    This script compute the Eddie Kinetic Energy (EKE) as in the Reynold's decomposition. u = u_mean + u_fluct.
"""
# Extracting data
u = []
v = []
for file in os.listdir('Data/gos/ke'):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2010 and int(file[:4]) <= 2020:
        u.append(np.nan_to_num(np.loadtxt(f'Data/gos/u/{file}')))
        v.append(np.nan_to_num(np.loadtxt(f'Data/gos/v/{file}')))


u_fluct = signal.detrend(u, axis = 0)
v_fluct = signal.detrend(v, axis = 0)
eke = 1/2 * (u_fluct**2 + v_fluct**2)
for file,i in zip(os.listdir('Data/gos/ke'),range(len(u))):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2010 and int(file[:4]) <= 2020:
        np.savetxt(f'Data/gos/eke/{file}',eke[i])
# Save the mensual mean of eke

date_data = np.loadtxt('Data/date_data.txt')
lon = np.loadtxt('Data/longitude.txt')
lat = np.loadtxt('Data/latitude.txt')

current_month_eke = [np.loadtxt(f'Data/gos/eke/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
for i in range(1,len(date_data)):
    eke = np.loadtxt(f'Data/gos/eke/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')
    current_month_eke.append(eke)
    
    if i == len(date_data)-1:
        np.savetxt(f'Data/gos/eke/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_eke, axis = 0))
        break

    if date_data[i,1] != date_data[i+1,1] : #last day of the month
        np.savetxt(f'Data/gos/eke/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_eke, axis = 0))
        current_month_eke = []

#Plotting data

projection = ccrs.LambertConformal(central_longitude = -18)
mask = np.loadtxt('Data/sit/month_mean/2010-10.txt')

for file in os.listdir('Data/gos/eke/month_mean'):
    try:
        year = int(file[:4])
        month = int(file[5:7])
        print(f"### - Saving eke: {year}-{month} - ###\n")
        fig = plt.figure(figsize=(9,7))
        axs = plt.axes(projection = projection)            
        xlim = [-43, 16]
        ylim = [61, 81]
        lower_space = 3 
        rect = mpath.Path([[xlim[0], ylim[0]],
                        [xlim[1], ylim[0]],
                        [xlim[1], ylim[1]],
                        [xlim[0], ylim[1]],
                        [xlim[0], ylim[0]],
                        ]).interpolated(20)
        proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
        rect_in_target = proj_to_data.transform_path(rect)
        axs.set_boundary(rect_in_target)
        axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
        
        axs.coastlines()
        axs.gridlines()
        levels = np.linspace(0,0.01,20)
        mean_eke = np.loadtxt(f'Data/gos/eke/month_mean/{year}-{month:02d}.txt')
        mean_eke = np.where(mask >= 0, mean_eke, np.nan)
        
        cs = axs.contourf(lon, lat, mean_eke, extend = 'both',levels = levels, cmap = "cmo.amp", transform=ccrs.PlateCarree())
        axs.set_title(" {} - {}".format(year,month), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [0,0.01])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots\Maps\GOS\EKE/{year}-{month:02d}.png")
        plt.close()
    except:
        pass