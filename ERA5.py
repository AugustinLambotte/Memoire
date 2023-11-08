import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy import interpolate
from datetime import date, timedelta
import cartopy.crs as ccrs
import os
import matplotlib.path as mpath
import cmocean
from global_land_mask import globe

earth_radius = 6370*1e3

def interp(lat_range = [60,83], lon_range = [-40,20]):
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

    lat_ERA5 = np.loadtxt('Data/ERA5/raw/lat_ERA5.txt')
    lon_ERA5 = np.loadtxt('Data/ERA5/raw/lon_ERA5.txt')

    lat_sit = np.loadtxt('Data/latitude.txt')
    lon_sit = np.loadtxt('Data/longitude.txt')

    mg_lon_ERA5, mg_lat_ERA5 = np.meshgrid(lon_ERA5,lat_ERA5)

    for file in os.listdir('Data\ERA5/raw/u'):
        print(file[:-4])
        u = np.loadtxt('Data\ERA5/raw/u/' + file)
        v = np.loadtxt('Data\ERA5/raw/v/' + file)

        # m/s => degree/s
        curr_dlat = u/(2*np.pi*earth_radius) * 360
        curr_dlon = v/(2*np.pi*earth_radius*np.cos(mg_lat_ERA5/360*2*np.pi)) * 360
        
        ### - Interpolate over the spatial cs2 grid - ###
        points = [] # list of length NxM containing all the coordinates [lat,lon] of all points from si drift map
        values_dlon = []
        values_dlat = []
        for lat in range(len(lat_ERA5)):
            for lon in range(len(lon_ERA5)):
                if curr_dlon[lat,lon] !=0 and curr_dlat[lat,lon] != 0 and (not np.isnan(curr_dlon[lat,lon])) and (not np.isnan(curr_dlat[lat,lon])):
                    points.append([float(lat_ERA5[lat]),float(lon_ERA5[lon])]) 
                    values_dlon.append(curr_dlon[lat,lon])
                    values_dlat.append(curr_dlat[lat,lon])
        points = np.array(points)
        values_dlon = np.array(values_dlon)
        values_dlat = np.array(values_dlat)

        dlon_interp = interpolate.griddata(points, values_dlon, (lat_sit, lon_sit), method='linear')
        dlat_interp = interpolate.griddata(points, values_dlat, (lat_sit, lon_sit), method='linear')

        # Inverse conversion, degree/s => m/s
        u_interp = dlon_interp/360 *2*np.pi * earth_radius * np.cos(lat_sit/360 * 2 *np.pi)
        v_interp = dlat_interp/360 *2*np.pi * earth_radius

        ke = 1/2 * (u_interp**2 + v_interp**2)
        land = np.loadtxt('Data/land.txt')

        np.savetxt(f'./Data/ERA5/interpolated/u/{file}',np.where(land == 0,u_interp,np.nan))
        np.savetxt(f'./Data/ERA5/interpolated/v/{file}',np.where(land == 0,v_interp,np.nan))
        np.savetxt(f'./Data/ERA5/interpolated/ke/{file}',np.where(land == 0,ke,np.nan))

def plot_mensual_mean(figsize = (9,7)):
    """
        Save and plot the mensual mean value of the SID data
    """
    projection = ccrs.LambertConformal(central_longitude = -18)
    date_data = np.loadtxt('Data/date_data.txt')
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')

    # Save the mensual mean of the sit and sic
    current_month_u = [np.loadtxt(f'Data/ERA5/interpolated/u/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
    current_month_v = [np.loadtxt(f'Data/ERA5/interpolated/v/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
    current_month_ke = [np.loadtxt(f'Data/ERA5/interpolated/ke/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]

    for i in range(1,len(date_data)):
        u = np.loadtxt(f'Data/ERA5/interpolated/u/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')
        v = np.loadtxt(f'Data/ERA5/interpolated/v/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')
        ke = np.loadtxt(f'Data/ERA5/interpolated/ke/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')

        current_month_u.append(u)
        current_month_v.append(v)
        current_month_ke.append(ke)
        if i == len(date_data)-1:
            np.savetxt(f'Data/ERA5/interpolated/u/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_u, axis = 0))
            np.savetxt(f'Data/ERA5/interpolated/v/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_v, axis = 0))
            np.savetxt(f'Data/ERA5/interpolated/ke/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_ke, axis = 0))

            break
        if date_data[i,1] != date_data[i+1,1] : #last day of the month
            np.savetxt(f'Data/ERA5/interpolated/u/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_u, axis = 0))
            np.savetxt(f'Data/ERA5/interpolated/v/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_v, axis = 0))
            np.savetxt(f'Data/ERA5/interpolated/ke/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_ke, axis = 0))

            current_month_u = []
            current_month_v = []
            current_month_ke = []
    """ ### - Mensual plots - ###  
    for file in os.listdir('Data/ERA5/interpolated/u/month_mean'):
        year = int(file[:4])
        month = int(file[5:7])
        print(f"### - Saving GOS: {year}-{month} - ###\n")

        ##### - Plot current - #####
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
        levels = np.linspace(0,10,100)
        mean_u = np.loadtxt(f'Data/ERA5/interpolated/u/month_mean/{year}-{month:02d}.txt') 
        mean_v = np.loadtxt(f'Data/ERA5/interpolated/v/month_mean/{year}-{month:02d}.txt')
        
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
        
        cb = plt.colorbar(cs, cax = cax, ticks = [0,5,10])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots\Maps\ERA5/vector/{year}-{month:02d}.png")
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
        levels = np.linspace(0,30,100)
        mean_ke = np.loadtxt(f'Data//ERA5/interpolated/ke/month_mean/{year}-{month:02d}.txt') 
        cs = axs.contourf(lon, lat, mean_ke, levels = levels, extend = 'both',cmap = "cmo.speed", transform=ccrs.PlateCarree())
    
        axs.set_title(" {} - {}".format(year,month), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [0,5,10,15,20,25,30])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots\Maps\ERA5\KE/{year}-{month:02d}.png")
        plt.close()
    """
    ### - Annual plots - ###
    for year_ in range(2011,2020):
        print(f"### - Saving wind: {year_} - ###\n")
        annual_era5ke = []
        for file in os.listdir('Data/ERA5/interpolated/ke/month_mean'):
            
            year = int(file[:4])
            month = int(file[5:7])
            if year == year_:
                annual_era5ke.append(np.loadtxt(f'Data/ERA5/interpolated/ke/month_mean/{year}-{month:02d}.txt'))
        print(len(annual_era5ke))
        mean_era5ke = np.nanmean(annual_era5ke,axis = 0)
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
        levels = np.linspace(0,30,100)
        mean_era5ke[mean_era5ke == 0] = np.nan
        cs = axs.contourf(lon, lat, mean_era5ke, levels = levels, extend = 'both',cmap = "cmo.speed", transform=ccrs.PlateCarree())
        axs.set_title(" {} - mean = {}J/kg".format(year_,round(np.nanmean(mean_era5ke),3)), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [0,10,20,30])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots/Maps/ERA5/KE/annual/{year_}.png")
        plt.close()

    ### - plot Total mean wind - ###
    
    print(f"### - Saving wind: {year_} - ###\n")
    annual_era5ke = []
    for file in os.listdir('Data/ERA5/interpolated/ke/month_mean'):
        year = int(file[:4])
        month = int(file[5:7])
        if year <= 2019 and year >= 2011:
            annual_era5ke.append(np.loadtxt(f'Data/ERA5/interpolated/ke/month_mean/{year}-{month:02d}.txt'))
    print(len(annual_era5ke))
    mean_era5ke = np.nanmean(annual_era5ke,axis = 0)
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
    levels = np.linspace(0,30,100)
    mean_era5ke[mean_era5ke == 0] = np.nan
    cs = axs.contourf(lon, lat, mean_era5ke, levels = levels, extend = 'both',cmap = "cmo.speed", transform=ccrs.PlateCarree())
    axs.set_title(" Mean over 2011 to 2019 period - spatial mean = {}J/kg".format(round(np.nanmean(mean_era5ke),3)), fontsize = 17)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [0,10,20,30])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots/Maps/ERA5/KE/annual/total_mean.png")
    plt.close()

#interp()
plot_mensual_mean()