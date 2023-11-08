import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy import interpolate
from datetime import date, timedelta
import cartopy.crs as ccrs
import os
import matplotlib.path as mpath
import cmocean

earth_radius = 6370*1e3

def extracting_data(lat_range = [60,83], lon_range = [-40,20]):
    """ 
        Given the file path "C:/.../..." to the .nc file it returns the longitude, latitude, sea_ice_thickness and time 
        restricted on the area defined by lat_range and lon_range
    """
    print("\n##############################")
    print("##### - Extracting data - ####")
    print("##############################\n")

    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]

    v_file = "C:/Users/Augustin/Downloads/Northward_sea_water_velocity"
    u_file = "C:/Users/Augustin/Downloads/Eastward_sea_water_velocity"
    
    v_ds = xr.open_dataset(v_file, decode_times = False)
    u_ds = xr.open_dataset(u_file, decode_times = False)
    u_gos= u_ds['ugos'].where((u_ds.longitude > lon_min) & (u_ds.longitude < lon_max) & (u_ds.latitude > lat_min) & (u_ds.latitude > 65.4 + (76.5-65.4)/(9+17) * (u_ds.longitude + 17)) & (u_ds.latitude < lat_max))
    v_gos= v_ds['vgos'].where((v_ds.longitude > lon_min) & (v_ds.longitude < lon_max) & (v_ds.latitude > lat_min) & (v_ds.latitude > 65.4 + (76.5-65.4)/(9+17) * (v_ds.longitude + 17)) & (v_ds.latitude < lat_max))
    #print(u_gos.sel(latitude = 60.125,longitude=-39.875,method = 'nearest'))
    lon_gos = u_ds.coords['longitude']
    lat_gos = u_ds.coords['latitude']


    
    #Interpolation over Cryosat2 grid

    lon_sit = np.loadtxt('Data/longitude.txt')
    lat_sit = np.loadtxt('Data/latitude.txt')
    date_data = np.loadtxt('Data/date_data.txt')
    
    
    for i in range(len(date_data)):
        ### - find the two date between which the averaging has to be done - ###
        year_cs2 = int(date_data[i,0])
        month_cs2 = int(date_data[i,1])
        day_cs2 = int(date_data[i,2])
        date_cs2 = date(year_cs2,month_cs2,day_cs2)

        if i == len(date_data)-1:
            next_year_cs2 = 2020
            next_month_cs2 = 7
            next_day_cs2 = 31
        else:
            next_year_cs2 = int(date_data[i+1,0])
            next_month_cs2 = int(date_data[i+1,1])
            next_day_cs2 = int(date_data[i+1,2])
        next_date_cs2 = date(next_year_cs2,next_month_cs2,next_day_cs2)
        print(f'--- saving GOS {year_cs2}-{month_cs2:02d}-{day_cs2:02d} ---')

        ### - Extracting mean value of the gos current for the two weeks of interest - ###
        two_weeks_u_gos = []
        two_weeks_v_gos = []
        for time in u_gos.time:
            time = int(time)
            corresponding_date = date(2010,10,1) + timedelta(days = time - int(u_gos.time[0]))
            if corresponding_date >= date_cs2 and corresponding_date < next_date_cs2:
                two_weeks_u_gos.append(u_gos.sel(time = time))
                two_weeks_v_gos.append(v_gos.sel(time = time))
        curr_ugos = np.nanmean(np.array(two_weeks_u_gos), axis = 0)
        curr_vgos = np.nanmean(np.array(two_weeks_v_gos), axis = 0)
        
        # Turn from m/s to degree_lat and degree_lon/s.
        lon_gos_mg,lat_gos_mg = np.meshgrid(lon_gos,lat_gos)
        curr_dlat = curr_vgos/(2*np.pi*earth_radius) * 360
        curr_dlon = curr_ugos/(2*np.pi*earth_radius*np.cos(lat_gos_mg/360*2*np.pi)) * 360
        
        ### - Interpolate over the spatial cs2 grid - ###
        points = [] # list of length NxM containing all the coordinates [lat,lon] of all points from si drift map
        values_dlon = []
        values_dlat = []
        for lat in range(len(lat_gos)):
            for lon in range(len(lon_gos)):
                if curr_dlon[lat,lon] !=0 and curr_dlat[lat,lon] != 0 and (not np.isnan(curr_dlon[lat,lon])) and (not np.isnan(curr_dlat[lat,lon])):
                    points.append([float(lat_gos[lat]),float(lon_gos[lon])]) 
                    values_dlon.append(curr_dlon[lat,lon])
                    values_dlat.append(curr_dlat[lat,lon])
        points = np.array(points)
        values_dlon = np.array(values_dlon)
        values_dlat = np.array(values_dlat)

        dlon_interp = interpolate.griddata(points, values_dlon, (lat_sit, lon_sit), method='linear')
        dlat_interp = interpolate.griddata(points, values_dlat, (lat_sit, lon_sit), method='linear')

        # Inverse conversion, degree/s => m/s
        u_gos_interp = dlon_interp/360 *2*np.pi * earth_radius * np.cos(lat_sit/360 * 2 *np.pi)
        v_gos_interp = dlat_interp/360 *2*np.pi * earth_radius
        #Computin kinetic energy with the new interpolated gos current values
        ke = 1/2 * (u_gos_interp**2 + v_gos_interp**2)
        land = np.loadtxt('Data/land.txt')
        np.savetxt(f'./Data/gos/u/{year_cs2}-{month_cs2:02d}-{day_cs2:02d}.txt',np.where(land == 0,u_gos_interp,np.nan))
        np.savetxt(f'./Data/gos/v/{year_cs2}-{month_cs2:02d}-{day_cs2:02d}.txt',np.where(land == 0,v_gos_interp,np.nan))
        np.savetxt(f'./Data/gos/ke/{year_cs2}-{month_cs2:02d}-{day_cs2:02d}.txt',np.where(land == 0,ke,np.nan))
    v_ds.close
    u_ds.close

    print('##### - Data extracted - #####\n')

def plot_mensual_mean(figsize = (9,7)):
    """
        Save and plot the mensual mean value of the SID data
    """
    projection = ccrs.LambertConformal(central_longitude = -18)
    date_data = np.loadtxt('Data/date_data.txt')
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')

    # Save the mensual mean of the sit and sic
    current_month_u = [np.loadtxt(f'Data/gos/u/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
    current_month_v = [np.loadtxt(f'Data/gos/v/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
    current_month_ke = [np.loadtxt(f'Data/gos/ke/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]

    for i in range(1,len(date_data)):
        u = np.loadtxt(f'Data/gos/u/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')
        v = np.loadtxt(f'Data/gos/v/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')
        ke = np.loadtxt(f'Data/gos/ke/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')

        current_month_u.append(u)
        current_month_v.append(v)
        current_month_ke.append(ke)
        if i == len(date_data)-1:
            np.savetxt(f'Data/gos/u/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_u, axis = 0))
            np.savetxt(f'Data/gos/v/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_v, axis = 0))
            np.savetxt(f'Data/gos/ke/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_ke, axis = 0))

            break
        if date_data[i,1] != date_data[i+1,1] : #last day of the month
            np.savetxt(f'Data/gos/u/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_u, axis = 0))
            np.savetxt(f'Data/gos/v/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_v, axis = 0))
            np.savetxt(f'Data/gos/ke/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_ke, axis = 0))

            current_month_u = []
            current_month_v = []
            current_month_ke = []
    
    ##### - Plot current - #####
    for file in os.listdir('Data/gos/v/month_mean'):
        year = int(file[:4])
        month = int(file[5:7])
        print(f"### - Saving GOS: {year}-{month} - ###\n")

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
        levels = np.linspace(0,0.5,100)
        mean_u = np.loadtxt(f'Data/gos/u/month_mean/{year}-{month:02d}.txt') 
        mean_v = np.loadtxt(f'Data/gos/v/month_mean/{year}-{month:02d}.txt')
        
        #Convert from m/s to degree/s
        dlon = mean_u/(2*np.pi*earth_radius*np.cos(lon/360 * 2*np.pi)) * 360
        dlat = mean_v/(2*np.pi*earth_radius) * 360
        current_intensity = np.sqrt(mean_u**2 + mean_v**2)
        cs = axs.contourf(lon, lat, current_intensity, levels = levels, extend = 'both',cmap = "cmo.speed", transform=ccrs.PlateCarree())
        
        #Vector plot
        current_intensity_vector_plot = np.sqrt(dlon**2 + dlat**2)
        dlon = mean_u/(current_intensity_vector_plot*10)
        dlat = mean_v/(current_intensity_vector_plot*10)
        axs.quiver(np.array(lon),np.array(lat),np.array(dlon),np.array(dlat),transform = ccrs.PlateCarree())

        axs.set_title(" {} - {}".format(year,month), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [0,0.1,0.2,0.3,0.4,0.5])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots\Maps\GOS\Current/{year}-{month:02d}.png")
        plt.close()


        ##### - Plot KE - #####
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
        levels = np.linspace(0,0.05,100)
        mean_ke = np.loadtxt(f'Data/gos/ke/month_mean/{year}-{month:02d}.txt') 
        cs = axs.contourf(lon, lat, mean_ke, levels = levels, extend = 'both',cmap = "cmo.speed", transform=ccrs.PlateCarree())
    
        axs.set_title(" {} - {}".format(year,month), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [0,0.01])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots\Maps\GOS\KE/{year}-{month:02d}.png")
        plt.close()
extracting_data()
plot_mensual_mean()