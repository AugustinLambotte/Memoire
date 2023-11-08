from netCDF4 import Dataset
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import xarray as xr
import pandas as pd
from datetime import datetime, date, timedelta
import cmocean
import matplotlib.path as mpath
import os
from global_land_mask import globe


####### - Functions definitions - #####

def extracting_data_sit(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc", lat_range = [60,83], lon_range = [-40,20]):
    """ 
        Given the file path "C:/.../..." to the .nc file it returns the longitude, latitude, sea_ice_thickness and time 
        restricted on the area defined by lat_range and lon_range

        output are in DataArray.
    """

    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]
    ds = xr.open_dataset(file, decode_times = False)
    lon = ds['Longitude'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    lat = ds['Latitude'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    sit = ds['Sea_Ice_Thickness']#.where((ds.Sea_Ice_Thickness != 0)&(ds.Sea_Ice_Concentration != 0)) 
    sit = sit.where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    sic = ds['Sea_Ice_Concentration']#.where((ds.Sea_Ice_Concentration != 0)& (ds.Sea_Ice_Thickness != 0)) 
    sic = sic.where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    
    sit_uncertainty = ds['Sea_Ice_Thickness_Uncertainty'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    time =  ds['Time']
    ds.close

    land = np.zeros(np.shape(np.array(lat)))
    for i in range(np.shape(np.array(lat))[0]):
        for j in range(np.shape(np.array(lat))[1]):
            if not np.isnan(np.array(lat)[i,j]):
                if globe.is_land(np.array(lat)[i,j],np.array(lon)[i,j]):
                    land[i,j] = 1
    
    np.savetxt('Data/land.txt',land)

    date_data = []
    for time_,sit_,sic_ in zip(time,sit,sic):
        year = (date(2010,10,1) + timedelta(days = int(time_) - int(time[0]))).year
        month = (date(2010,10,1) + timedelta(days = int(time_) - int(time[0]))).month
        day = (date(2010,10,1) + timedelta(days = int(time_) - int(time[0]))).day
        date_data.append([year,month,day,int(time_) - int(time[0])])
        np.savetxt(f'Data/sit/{year}-{month:02d}-{day:02d}.txt',np.array(sit_))
        np.savetxt(f'Data/sic/{year}-{month:02d}-{day:02d}.txt',np.array(sic_))
    np.savetxt('Data/longitude.txt',np.array(lon))
    np.savetxt('Data/latitude.txt',np.array(lat))
    np.savetxt('Data/date_data.txt',date_data, fmt = '%i')
    



def plot_mensual_mean(figsize = (9,7)):
    projection = ccrs.LambertConformal(central_longitude = -18)
    date_data = np.loadtxt('Data/date_data.txt')
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')

    # Save the mensual mean of the sit and sic
    current_month_sit = [np.loadtxt(f'Data/sit/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
    current_month_sic = [np.loadtxt(f'Data/sic/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
    for i in range(1,len(date_data)):
        sit = np.loadtxt(f'Data/sit/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')
        sic = np.loadtxt(f'Data/sic/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')

        current_month_sit.append(sit)
        current_month_sic.append(sic)
        if i == len(date_data)-1:
            np.savetxt(f'Data/sit/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_sit, axis = 0))
            np.savetxt(f'Data/sic/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_sic, axis = 0))

            break
        if date_data[i,1] != date_data[i+1,1] : #last day of the month
            np.savetxt(f'Data/sit/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_sit, axis = 0))
            np.savetxt(f'Data/sic/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_sic, axis = 0))

            current_month_sit = []
            current_month_sic = []

    for file in os.listdir('Data/sic/month_mean'):
        year = int(file[:4])
        month = int(file[5:7])
        print(f"### - Saving SIT: {year}-{month} - ###\n")
        fig = plt.figure(figsize=figsize)
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
        levels = np.linspace(0,4,1000)
        mean_sit = np.loadtxt(f'Data/sit/month_mean/{year}-{month:02d}.txt')
        mean_sic = np.loadtxt(f'Data/sic/month_mean/{year}-{month:02d}.txt')
        cs = axs.contourf(lon, lat, mean_sit, levels = levels, extend = 'both',cmap = "cmo.ice", transform=ccrs.PlateCarree())
        cs_ = axs.contour(lon, lat, mean_sic,[0.15], colors = 'red', transform=ccrs.PlateCarree())
        axs.set_title(" {} - {}".format(year,month), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [0,1,2,3,4])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots\Maps\Sea_ice\Thickness/{year}-{month:02d}.png")
        plt.close()
    
extracting_data_sit()
plot_mensual_mean()