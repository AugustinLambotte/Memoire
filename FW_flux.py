import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os 
from datetime import timedelta, date
from scipy import interpolate, stats
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import math
import seawater as sw
"""
    This script compute the correlation coefficient for every pixels (80km^2) over 
    the EGC between the intensity of the EGC (KE) and the fresh water flux.
"""
earth_radius = 6370*1e3

def save_siv_transp_fw(method='personal'):
    # For each cell NOT on the edge of the grid, we compute the net transport (positive when net import of sea ice) by considering the four cells 
    # on the edges of the cell. We stock this in Cell_transport. In Cell_siv_var we stock the variation of sea ice volume for each cell. We divide
    # them by the number of days between the two date to have a 'daily' siv variation.

    date_data = np.loadtxt('Data/date_data.txt')
    for day in range(np.shape(date_data)[0]-1):
        print(f'{int(date_data[day][0])}-{int(date_data[day][1])}-{int(date_data[day][2])}')
        #time_gap = date_data[day,-1] - date_data[day-1,-1]

        lat = np.loadtxt('Data/latitude.txt')
        lon = np.loadtxt('Data/longitude.txt')

        #The two following arrays will be filled with grid cell distances in X and Y direction
        X_dist = np.zeros((np.shape(lat)[0],np.shape(lat)[1]-1))
        Y_dist = np.zeros((np.shape(lat)[0]-1,np.shape(lat)[1]))
        for i in range(np.shape(X_dist)[0]):
            for j in range(np.shape(X_dist)[1]):
                X_dist[i,j] = np.nan_to_num(sw.extras.dist([lat[i,j],lat[i,j+1]],[lon[i,j],lon[i,j+1]])[0])*1000 # [m]
        
        for i in range(np.shape(Y_dist)[0]):
            for j in range(np.shape(Y_dist)[1]):
                Y_dist[i,j] = np.nan_to_num(sw.extras.dist([lat[i,j],lat[i+1,j]],[lon[i,j],lon[i+1,j]])[0])*1000 # [m]
    
        """ print(np.shape(X_dist))
        print(np.shape(Y_dist))
        plt.imshow(X_dist)
        plt.show()
        plt.imshow(Y_dist)
        plt.show() """
        
        # Load Sea ice drift data (m/days)
        X_drift = np.loadtxt(f'Data/X_drift/{int(date_data[day][0])}-{int(date_data[day][1]):02d}-{int(date_data[day][2]):02d}.txt')
        Y_drift = np.loadtxt(f'Data/Y_drift/{int(date_data[day][0])}-{int(date_data[day][1]):02d}-{int(date_data[day][2]):02d}.txt')

        # Load Sea ice thickness and concentration data for the two weeks of interest
        sit = np.loadtxt(f'Data/sit/{int(date_data[day][0])}-{int(date_data[day][1]):02d}-{int(date_data[day][2]):02d}.txt')
        sic = np.loadtxt(f'Data/sic/{int(date_data[day][0])}-{int(date_data[day][1]):02d}-{int(date_data[day][2]):02d}.txt')

        # Load Sea ice thickness and concentration data for the two following weeks 
        sit_next = np.loadtxt(f'Data/sit/{int(date_data[day+1][0])}-{int(date_data[day+1][1]):02d}-{int(date_data[day+1][2]):02d}.txt')
        sic_next = np.loadtxt(f'Data/sic/{int(date_data[day+1][0])}-{int(date_data[day+1][1]):02d}-{int(date_data[day+1][2]):02d}.txt')
        
        transport = np.zeros((np.shape(X_drift)[0],np.shape(X_drift)[1]))
        siv_var = np.zeros((np.shape(X_drift)[0],np.shape(X_drift)[1]))
        if method == 'personal':
            for line in range(np.shape(X_drift)[0]):
                for col in range(np.shape(X_drift)[1]):
                    if line == 0 or line == np.shape(X_drift)[0] -1 or col == 0 or col == np.shape(X_drift)[1] - 1:
                        transport[line,col] = 0
                        siv_var[line,col] = 0

                    else:
                        from_above_si_drift = -Y_drift[line-1,col]
                        from_left_si_drift = X_drift[line,col-1] 
                        from_right_si_drift = -X_drift[line,col+1]
                        from_below_si_drift = Y_drift[line+1,col]

                        from_above_si_drift = np.nan_to_num(from_above_si_drift)
                        from_left_si_drift = np.nan_to_num(from_left_si_drift)
                        from_right_si_drift = np.nan_to_num(from_right_si_drift)
                        from_below_si_drift = np.nan_to_num(from_below_si_drift)

                        from_above_transport = from_above_si_drift * sit[line-1,col] * sic[line-1,col] * X_dist[line-1,col-1] #[m^3] 
                        from_below_transport = from_below_si_drift * sit[line+1,col] * sic[line+1,col] * X_dist[line-1,col-1] #[m^3]
                        from_right_transport = from_right_si_drift * sit[line,col+1] * sic[line,col+1] * Y_dist[line-1,col-1] #[m^3]
                        from_left_transport = from_left_si_drift * sit[line,col-1] * sic[line,col-1] * Y_dist[line-1,col-1] #[m^3]            
                        in_cell = 0 #m^3 of sea ice coming into the cell
                        if from_above_transport > 0:
                            in_cell += from_above_transport
                        if from_left_transport > 0:
                            in_cell += from_left_transport
                        if from_right_transport > 0:
                            in_cell += from_right_transport
                        if from_below_transport > 0:
                            in_cell += from_below_transport
                        
                        out_cell = (abs(np.nan_to_num(X_drift[line,col])) * Y_dist[line-1,col-1] + abs(np.nan_to_num(Y_drift[line,col])) * X_dist[line-1,col-1]) * sit[line,col] * sic[line,col] #m^3 of sea ice leaving the cell
                        
                        transport[line,col] = in_cell - out_cell
                        siv_var[line,col] = (sit_next[line,col] * sic_next[line,col] - sit[line,col] * sic[line,col])*80000**2 #Positive when siv increase over the cell
            
            net_mass_flux_of_sea_ice_melt = transport-siv_var        
        
        if method == 'divergence':
            X_dist_cumsum = np.zeros(np.shape(lat)) 
            Y_dist_cumsum = np.zeros(np.shape(lat)) 
            for i in range(np.shape(X_dist_cumsum)[0]):
                for j in range(np.shape(Y_dist_cumsum)[1]):
                    if j != 0:
                        X_dist_cumsum[i,j] = np.sum(X_dist[i,:j])
                    if i != 0:
                        Y_dist_cumsum[i,j] = np.sum(Y_dist[:i,j])

            
            X_div = np.array([np.gradient(np.nan_to_num(X_drift[line,:]),X_dist_cumsum[line,:]) for line in range(np.shape(X_drift)[0])])
            Y_div = np.transpose([np.gradient(np.nan_to_num(Y_drift[:,col]),Y_dist_cumsum[:,col]) for col in range(np.shape(Y_drift)[1])])
            
            divergence = X_div + Y_div #[m/days]
            plt.imshow(divergence)
            plt.show()
            transport = divergence*sit*sic*80e3 #[m^3/days]
            siv_var = (sit_next * sic_next - sit * sic) * 80e3**2 #[m^3/days]
            net_mass_flux_of_sea_ice_melt = transport - siv_var #[m^3/days]
            plt.subplot(221)
            plt.title('X_dist_cumsum')
            plt.imshow(X_dist_cumsum)

            plt.subplot(222)
            plt.title('Y_dist_cumsum')
            plt.imshow(Y_dist_cumsum)

            plt.subplot(223)
            plt.title('Y_drift')
            plt.imshow(Y_drift)

            plt.subplot(224)
            plt.title('X_drift')
            plt.imshow(X_drift)

            plt.show()
        
        net_mass_flux_of_sea_ice_melt /= 80000**2 #from m^3/80km^2 to m^3/m^2
        transport /= 80000**2 #from m^3/80km^2 to m^3/m^2
        siv_var /= 80000**2 #from m^3/80km^2 to m^3/m^2
        np.savetxt(f'Data/nmfsim/{int(date_data[day][0])}-{int(date_data[day][1]):02d}-{int(date_data[day][2]):02d}.txt',net_mass_flux_of_sea_ice_melt)
        np.savetxt(f'Data/siv_var/{int(date_data[day][0])}-{int(date_data[day][1]):02d}-{int(date_data[day][2]):02d}.txt',transport)
        np.savetxt(f'Data/transport/{int(date_data[day][0])}-{int(date_data[day][1]):02d}-{int(date_data[day][2]):02d}.txt',siv_var)

def plot_fw_flux(figsize = (9,7)):
    projection = ccrs.LambertConformal(central_longitude = -18)
    date_data = np.loadtxt('Data/date_data.txt')
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')

    ### - Save the mensual mean of the fw_flux - ###
    current_month_nmfsim = [np.loadtxt(f'Data/nmfsim/{int(date_data[0,0])}-{int(date_data[0,1]):02d}-{int(date_data[0,2]):02d}.txt')]
    for i in range(1,len(date_data)-1):
        nmfsim = np.loadtxt(f'Data/nmfsim/{int(date_data[i,0])}-{int(date_data[i,1]):02d}-{int(date_data[i,2]):02d}.txt')

        current_month_nmfsim.append(nmfsim)
        if i == len(date_data)-1:
            np.savetxt(f'Data/nmfsim/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_nmfsim, axis = 0))

            break
        if date_data[i,1] != date_data[i+1,1] or date_data[i+1,2] > 20: #last day of the month
            np.savetxt(f'Data/nmfsim/month_mean/{int(date_data[i,0])}-{int(date_data[i,1]):02}.txt',np.nanmean(current_month_nmfsim, axis = 0))

            current_month_nmfsim = []

    ### - plot mensual nmfsim - ###
    for file in os.listdir('Data/nmfsim/month_mean'):
        year = int(file[:4])
        month = int(file[5:7])
        print(f"### - Saving nmfsim: {year}-{month} - ###\n")
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
        levels = np.linspace(-1,1,100)
        mean_fw = np.loadtxt(f'Data/nmfsim/month_mean/{year}-{month:02d}.txt')
        mean_fw[mean_fw == 0] = np.nan
        cs = axs.contourf(lon, lat, mean_fw, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
        axs.set_title(" {} - {}".format(year,month), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [-1,0,1])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots/Maps/nmfsim/{year}-{month:02d}.png")
        plt.close()
   
    ### - plot annual FW flux - ###
    for year_ in range(2011,2020):
        print(f"### - Saving fw: {year_} - ###\n")
        annual_nmfsim = []
        for file in os.listdir('Data/nmfsim/month_mean'):
            """ if file == 'month_mean':
                continue """
            year = int(file[:4])
            month = int(file[5:7])
            if year == year_:
                annual_nmfsim.append(np.loadtxt(f'Data/nmfsim/month_mean/{year}-{month:02d}.txt'))
        print(len(annual_nmfsim))
        mean_nmfsim = np.nanmean(annual_nmfsim,axis = 0)
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
        levels = np.linspace(-0.3,0.3,100)
        mean_nmfsim[mean_nmfsim == 0] = np.nan
        cs = axs.contourf(lon, lat, mean_nmfsim, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
        axs.set_title(" {} - mean = {}m".format(year_,round(np.nanmean(mean_nmfsim),3)), fontsize = 30)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [-0.3,0,0.3])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots/Maps/nmfsim/annual/{year_}.png")
        plt.close()

    ### - plot Total mean FW flux - ###
    
    print(f"### - Saving fw: {year_} - ###\n")
    annual_nmfsim = []
    for file in os.listdir('Data/sic/month_mean'):
        year = int(file[:4])
        month = int(file[5:7])
        if year <= 2019 and year >= 2011:
            annual_nmfsim.append(np.loadtxt(f'Data/nmfsim/month_mean/{year}-{month:02d}.txt'))
    print(len(annual_nmfsim))
    mean_nmfsim = np.nanmean(annual_nmfsim,axis = 0)
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
    levels = np.linspace(-0.3,0.3,100)
    mean_nmfsim[mean_nmfsim == 0] = np.nan
    cs = axs.contourf(lon, lat, mean_nmfsim, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    axs.set_title(" Mean over 2011 to 2019 period - spatial mean = {}m".format(round(np.nanmean(mean_nmfsim),3)), fontsize = 17)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.3,0,0.3])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots/Maps/nmfsim/total_mean.png")
    plt.close()

#save_siv_transp_fw(method = 'personal')
plot_fw_flux()