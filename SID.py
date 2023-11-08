import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from datetime import datetime, date, timedelta
import matplotlib.path as mpath
import cmocean as cmo
import cartopy.crs as ccrs
import scipy
import os 
import seawater as sw

"""
    This script interpol the sea ice drift data over the sit grid and with a daily resolution (the original resolution of this dataset) save it in ./Data
"""
earth_radius = 6370 * 1e3
def extracting_SID(lat_range = [60,83], lon_range = [-40,20]):
    
    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]
    daily = True
    
    date_data = np.loadtxt('Data/date_data.txt')
    str_date_data = [f'{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}' for i in range(len(date_data))]
    dir = "C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/"
    print("\n##### - Extracting Sea Ice drift data - #####\n")
    for year in range(year_,year_end):
        day_ignored = 0
        if year == 2010:
            month_ = 10
        else:
            month_ = 1
        for month in range(month_,13):
            
            if month ==1:
                nb_days = 31
                month_name = "jan"
            if month ==2:
                month_name = "feb"
                if year == 2010:
                    nb_days = 28
                if year == 2011:
                    nb_days = 28
                if year == 2012:
                    nb_days = 29
                if year == 2013:
                    nb_days = 28
                if year == 2014:
                    nb_days = 28
                if year == 2015:
                    nb_days = 28
                if year == 2016:
                    nb_days = 29
                if year == 2017:
                    nb_days = 28
                if year == 2018:
                    nb_days = 28
                if year == 2019:
                    nb_days = 28
                if year == 2020:
                    nb_days = 29
            if month ==3:
                nb_days = 31
                month_name = "mar"
            if month ==4:
                nb_days = 30
                month_name = "april"
            if month ==5:
                nb_days = 31
                month_name = "may"
            if month ==6:
                nb_days = 30
                month_name = "june"
            if month ==7:
                nb_days = 31
                month_name = "july"
            if month ==8:
                nb_days = 31
                month_name = "aug"
            if month ==9:
                nb_days = 30
                month_name = "sept"
            if month ==10:
                nb_days = 31
                month_name = "oct"
            if month ==11:
                nb_days = 30
                month_name = "nov"
            if month ==12:
                nb_days = 31
                month_name = "dec"
            for day in range(1,nb_days+1):
                print(f'{year}-{month}-{day}')
                if f'{year}-{month:02d}-{day:02d}' in str_date_data or daily:
                   
                    file =  dir+f"{year}" + "/" + f"{month:02d}" + f"/ice_drift_nh_ease2-750_cdr-v1p0_24h-{year}{month:02d}{day:02d}1200.nc"
                    ds = xr.open_dataset(file, decode_times = False)
                    
                    dlat = (ds['lat1'] - ds['lat']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)
                    dlon = (ds['lon1'] - ds['lon']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)
                    
                    dlat = dlat.sel(time = int(ds.time))
                    dlon = dlon.sel(time = int(ds.time))
                    
                    if year == year_:
                        if month == month_:
                            if day == 1:
                                lat = ds['lat'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)
                                lon = ds['lon'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)

                                
                    ds.close
                    #Interpolation over SIT grid
                    points = [] # list of length NxM containing all the coordinates [lat,lon] of all points from si drift map
                    values_dlat = []
                    values_dlon = []
                    for i in range(len(lat)):
                        for j in range(len(lon[0])):
                            if dlat[i,j] !=0 and dlon[i,j] != 0 and not np.isnan(dlat[i,j]) and not np.isnan(dlon[i,j]):
                                points.append([lat[i,j],lon[i,j]])
                                values_dlat.append(dlat[i,j])
                                values_dlon.append(dlon[i,j])
                    points = np.array(points)
                    values_dlat = np.array(values_dlat)
                    values_dlon = np.array(values_dlon)
                    if len(values_dlat) > 2:
                        dlat_interp = interpolate.griddata(points, values_dlat, (lat_sit, lon_sit), method='linear')
                        dlon_interp = interpolate.griddata(points, values_dlon, (lat_sit, lon_sit), method='linear')
                        
                    else:
                        print('Manual interpolation')
                        #In this case, there are not enough data points (less than 3) to perform a classical interpolation
                        day_ignored += 1
                        manualy_interpolated_dlat = np.empty(lat_sit.shape)
                        manualy_interpolated_dlon = np.empty(lat_sit.shape)
                        manualy_interpolated_dlat[:] = np.nan
                        manualy_interpolated_dlon[:] = np.nan
                        for coords, temp_dX, temp_dY in zip(points,values_dlat,values_dlon):
                            idx_lat = int((np.abs(lat_sit - coords[0])).argmin())
                            idx_lon = int((np.abs(lon_sit - coords[1])).argmin())
                            
                            
                            manualy_interpolated_dlat[np.unravel_index(idx_lon,lon_sit.shape)] = temp_dX
                            manualy_interpolated_dlon[np.unravel_index(idx_lon,lon_sit.shape)] = temp_dY

                        N_drift = manualy_interpolated_dlat
                        E_drift = manualy_interpolated_dlon
                    

                    N_drift = (dlat_interp/360 * 2*np.pi) * earth_radius
                    E_drift = (dlon_interp/360 * 2*np.pi) * np.cos(lat_sit/360 * 2*np.pi)* earth_radius
                    land = np.loadtxt('Data/land.txt')
                
                    if daily:
                        np.savetxt(f'Data/N_drift/daily/{year}-{month:02d}-{day:02d}.txt',np.where(land ==0,N_drift,np.nan))
                        np.savetxt(f'Data/E_drift/daily/{year}-{month:02d}-{day:02d}.txt',np.where(land ==0,E_drift,np.nan))
                    else:
                        pass
                        """ np.savetxt(f'Data/X_drift/{year}-{month:02d}-{day:02d}.txt',lat_drift[f'y{year}'][month_name][-1])
                        np.savetxt(f'Data/Y_drift/{year}-{month:02d}-{day:02d}.txt',lon_drift[f'y{year}'][month_name][-1]) """
            #print(f"{day_ignored} days man interp during year {year}")

def bw_mean():
    """
        This function load the daily SID data and compute the biweekly mean as in the cryosat2 dataset.
        The daily SID data are drift in N or E direction but here we make the conversion to X and Y drift (m/days)
    """
    date_data = np.loadtxt('Data/date_data.txt')
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')
    for i in range(len(date_data)):
        two_weeks_X_drift = []
        two_weeks_Y_drift = []

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
        
        print(f'--- saving SID {year_cs2}-{month_cs2:02d}-{day_cs2:02d} ---')
        for file in os.listdir('Data/N_drift/daily'):

            year_sid = int(file[:4])
            month_sid = int(file[5:7])
            day_sid = int(file[8:10])
            date_sid = date(year_sid,month_sid,day_sid)
            if date_sid >= date_cs2 and date_sid < next_date_cs2:
                #extracting sea ice drift
                N_drift = np.loadtxt('Data/N_drift/daily/' + file) # m/days
                E_drift = np.loadtxt('Data/E_drift/daily/' + file) # m/days

                # Passing in grid direction (X and Y)
                X_drift = E_drift * np.sin(lon/360 * 2*np.pi) + N_drift * np.cos(lon/360 * 2*np.pi) # m/days
                Y_drift = -E_drift * np.cos(lon/360 * 2*np.pi) + N_drift * np.sin(lon/360 * 2*np.pi) # m/days
                   
                two_weeks_X_drift.append(X_drift)
                two_weeks_Y_drift.append(Y_drift)
                
        X_drift = np.nanmean(np.array(two_weeks_X_drift),axis = 0)
        Y_drift = np.nanmean(np.array(two_weeks_Y_drift),axis = 0)

        np.savetxt(f'Data/X_drift/{year_cs2}-{month_cs2:02d}-{day_cs2:02d}.txt',X_drift)
        np.savetxt(f'Data/Y_drift/{year_cs2}-{month_cs2:02d}-{day_cs2:02d}.txt',Y_drift)

def plot_mensual_mean(figsize = (9,7)):
    """
        Save and plot the mensual mean value of the SID data
    """
    projection = ccrs.LambertConformal(central_longitude = -18)
    date_data = np.loadtxt('Data/date_data.txt')
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')

    # Save the mensual mean of the sit and sic
    current_month_X_drift = [np.loadtxt(f'Data/X_drift/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
    current_month_Y_drift = [np.loadtxt(f'Data/Y_drift/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
    for i in range(1,len(date_data)):
        X_drift = np.loadtxt(f'Data/X_drift/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')
        Y_drift = np.loadtxt(f'Data/Y_drift/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')

        current_month_X_drift.append(X_drift)
        current_month_Y_drift.append(Y_drift)
        if i == len(date_data)-1:
            np.savetxt(f'Data/X_drift/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_X_drift, axis = 0))
            np.savetxt(f'Data/Y_drift/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_Y_drift, axis = 0))

            break
        if date_data[i,1] != date_data[i+1,1] : #last day of the month
            np.savetxt(f'Data/X_drift/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_X_drift, axis = 0))
            np.savetxt(f'Data/Y_drift/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_Y_drift, axis = 0))

            current_month_X_drift = []
            current_month_Y_drift = []

    for file in os.listdir('Data/Y_drift/month_mean'):
        year = int(file[:4])
        month = int(file[5:7])
        print(f"### - Saving SID: {year}-{month} - ###\n")
        ### - Velocity field - ###
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
        levels = np.linspace(0,0.5,1000)
        mean_X_drift = np.loadtxt(f'Data/X_drift/month_mean/{year}-{month:02d}.txt') /(24*60*60) #Passing from m/days to m/s
        mean_Y_drift = np.loadtxt(f'Data/Y_drift/month_mean/{year}-{month:02d}.txt') /(24*60*60) #Passing from m/days to m/s
        drift_intensity = np.sqrt(mean_X_drift**2 + mean_Y_drift**2)
        cs = axs.contourf(lon, lat, drift_intensity, levels = levels, extend = 'both',cmap = "cmo.speed", transform=ccrs.PlateCarree())
        
        #Vector plot
        #Transform from X and Y displacement to lat and lon displacement

        dN = mean_X_drift*np.cos(lon/360 * 2*np.pi) + mean_Y_drift*np.sin(lon/360 * 2*np.pi)
        dE = mean_X_drift*np.sin(lon/360 * 2*np.pi) - mean_Y_drift*np.cos(lon/360 * 2*np.pi)

        dlat = dN/(2*np.pi * earth_radius) * 360
        dlon = dE/(2*np.pi * earth_radius * np.cos(lat/360 * 2*np.pi)) * 360
        
        drift_intensity = np.sqrt(dlat**2+dlon**2)

        dlat = dlat/drift_intensity * 10
        dlon = dlon/drift_intensity * 10
        axs.quiver(np.array(lon),np.array(lat),np.array(dlon),np.array(dlat),scale = 400,transform = ccrs.PlateCarree())

        axs.set_title(" {} - {}".format(year,month), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [0,0.5])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots\Maps\Sea_ice\Drift/{year}-{month:02d}.png")
        plt.close()

        ############################
        ### - divergence field - ###
        ############################

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
        levels = np.linspace(-1e-6,1e-6,100)
        mean_X_drift = np.loadtxt(f'Data/X_drift/month_mean/{year}-{month:02d}.txt') /(24*60*60) #Passing from m/days to m/s
        mean_Y_drift = np.loadtxt(f'Data/Y_drift/month_mean/{year}-{month:02d}.txt') /(24*60*60) #Passing from m/days to m/s
        
        #The two following arrays will be filled with grid cell distances in X and Y direction
        X_dist = np.zeros((np.shape(lat)[0],np.shape(lat)[1]-1))
        Y_dist = np.zeros((np.shape(lat)[0]-1,np.shape(lat)[1]))
        for i in range(np.shape(X_dist)[0]):
            for j in range(np.shape(X_dist)[1]):
                X_dist[i,j] = np.nan_to_num(sw.extras.dist([lat[i,j],lat[i,j+1]],[lon[i,j],lon[i,j+1]])[0])*1000 # [m]
        
        for i in range(np.shape(Y_dist)[0]):
            for j in range(np.shape(Y_dist)[1]):
                Y_dist[i,j] = np.nan_to_num(sw.extras.dist([lat[i,j],lat[i+1,j]],[lon[i,j],lon[i+1,j]])[0])*1000 # [m]

        X_dist_cumsum = np.zeros(np.shape(lat)) 
        Y_dist_cumsum = np.zeros(np.shape(lat)) 
        for i in range(np.shape(X_dist_cumsum)[0]):
            for j in range(np.shape(Y_dist_cumsum)[1]):
                if j != 0:
                    X_dist_cumsum[i,j] = np.sum(X_dist[i,:j])
                if i != 0:
                    Y_dist_cumsum[i,j] = np.sum(Y_dist[:i,j])

        
        X_div = np.array([np.gradient(np.nan_to_num(mean_X_drift[line,:]),X_dist_cumsum[line,:]) for line in range(np.shape(X_drift)[0])])
        Y_div = np.transpose([np.gradient(np.nan_to_num(mean_Y_drift[:,col]),Y_dist_cumsum[:,col]) for col in range(np.shape(Y_drift)[1])])
        
        divergence = X_div + Y_div #[m/days]
        cs = axs.contourf(lon, lat, divergence, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
        
        #Vector plot
        #Transform from X and Y displacement to lat and lon displacement

        dN = mean_X_drift*np.cos(lon/360 * 2*np.pi) + mean_Y_drift*np.sin(lon/360 * 2*np.pi)
        dE = mean_X_drift*np.sin(lon/360 * 2*np.pi) - mean_Y_drift*np.cos(lon/360 * 2*np.pi)

        dlat = dN/(2*np.pi * earth_radius) * 360
        dlon = dE/(2*np.pi * earth_radius * np.cos(lat/360 * 2*np.pi)) * 360
        
        drift_intensity = np.sqrt(dlat**2+dlon**2)

        dlat = dlat/drift_intensity * 10
        dlon = dlon/drift_intensity * 10
        axs.quiver(np.array(lon),np.array(lat),np.array(dlon),np.array(dlat),scale = 400,transform = ccrs.PlateCarree())

        axs.set_title(" {} - {}".format(year,month), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [-1e-6,1e-6])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots\Maps\Sea_ice\Drift/divergence/{year}-{month:02d}.png")
        plt.close()

    ###########################################
    ### - Total averaged divergence field - ###
    ###########################################
    mean_X = []
    mean_Y = []
    for file in os.listdir('Data/Y_drift'):
        if file == 'month_mean':
            continue
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2019:
            mean_X.append(np.loadtxt('Data/X_drift/'+file)/(24*60*60)) #Passing from m/days to m/s)
            mean_Y.append(np.loadtxt('Data/Y_drift/'+file)/(24*60*60)) #Passing from m/days to m/s)

    mean_X_drift = np.nanmean(np.array(mean_X),axis = 0)
    mean_Y_drift = np.nanmean(np.array(mean_Y),axis = 0)

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
    levels = np.linspace(-1e-6,1e-6,100)

    
    X_div = [np.gradient(np.nan_to_num(mean_X_drift[line,:]),X_dist_cumsum[line,:]) for line in range(np.shape(X_drift)[0])]
    Y_div = [np.gradient(np.nan_to_num(mean_Y_drift[:,col]),Y_dist_cumsum[:,col]) for col in range(np.shape(Y_drift)[1])]
    
    X_div = np.array(X_div)
    Y_div = np.array(np.transpose(Y_div))
    
    divergence = X_div + Y_div #[m/days]
    cs = axs.contourf(lon, lat, divergence, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    
    #Vector plot
    #Transform from X and Y displacement to lat and lon displacement

    dN = mean_X_drift*np.cos(lon/360 * 2*np.pi) + mean_Y_drift*np.sin(lon/360 * 2*np.pi)
    dE = mean_X_drift*np.sin(lon/360 * 2*np.pi) - mean_Y_drift*np.cos(lon/360 * 2*np.pi)

    dlat = dN/(2*np.pi * earth_radius) * 360
    dlon = dE/(2*np.pi * earth_radius * np.cos(lat/360 * 2*np.pi)) * 360
    
    drift_intensity = np.sqrt(dlat**2+dlon**2)

    dlat = dlat/drift_intensity * 10
    dlon = dlon/drift_intensity * 10
    axs.quiver(np.array(lon),np.array(lat),np.array(dlon),np.array(dlat),scale = 400,transform = ccrs.PlateCarree())

    axs.set_title(" Divergence averaged over 2011 to 2019 period", fontsize = 20)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-1e-6,1e-6])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Sea_ice\Drift/divergence/2011_2019_averaged.png")
    plt.close()


    
year_ = 2010
year_end = 2021
lat_sit = np.loadtxt('Data/latitude.txt')
lon_sit = np.loadtxt('Data/longitude.txt')
extracting_SID()
bw_mean()
plot_mensual_mean()