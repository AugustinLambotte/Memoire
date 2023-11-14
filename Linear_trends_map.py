import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import xarray as xr
from scipy import stats
"""
    This script is used to plot on a map the value of the trend over all the time span for the EGC energy and the fresh water flux and for each pixel
"""
earth_radius = 6370*1e3

def trend(mode):
    """
        return maps of the trend of nmfsim and EGC energy
    """
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')
    
    nmfsim = []
    ke = []
    sit = []
    sic = []
    X_drift = []
    Y_drift = []
    transp = []
    siv_var = []
    eke = []
    era5ke = []
    date_data = np.loadtxt('Data/date_data.txt')

    date_data_2011_to_2019 = []
    for date in date_data:
        if date[0] >= 2011 and date[0] < 2020:
            date_data_2011_to_2019.append(date)
    date_data_2011_to_2019 = np.array(date_data_2011_to_2019)


    for file in os.listdir('Data/nmfsim'):
        if file == 'month_mean':
            continue
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2019:
            
            nmfsim.append(np.loadtxt(f'Data/nmfsim/{file}'))
            ke.append(np.loadtxt(f'Data/gos/ke/{file}'))
            sit.append(np.loadtxt(f'Data/sit/{file}'))
            sic.append(np.loadtxt(f'Data/sic/{file}'))
            X_drift.append(np.loadtxt(f'Data/X_drift/{file}'))
            Y_drift.append(np.loadtxt(f'Data/Y_drift/{file}'))
            transp.append(np.loadtxt(f'Data/Transport/{file}'))
            siv_var.append(np.loadtxt(f'Data/SIV_var/{file}'))
            eke.append(np.loadtxt(f'Data/gos/eke/{file}'))
            era5ke.append(np.loadtxt(f'Data/ERA5/interpolated/ke/{file}'))

    nmfsim = np.array(nmfsim)
    ke = np.nan_to_num(np.array(ke))
    mean_ke = np.nanmean(ke,axis = 0)
    sic= np.nan_to_num(np.array(sic))
    nmfsim= np.nan_to_num(np.array(nmfsim))
    sit= np.nan_to_num(np.array(sit))
    X_drift = np.nan_to_num(X_drift)
    Y_drift = np.nan_to_num(Y_drift)
    transp= np.nan_to_num(np.array(transp))
    siv_var= np.nan_to_num(np.array(siv_var))
    eke= np.nan_to_num(np.array(eke))
    era5ke= np.nan_to_num(np.array(era5ke))
    drift_magnitude = np.sqrt(X_drift**2+Y_drift**2)

    trend_nmfsim = np.zeros((np.shape(nmfsim)[1:]))
    trend_ke = np.zeros((np.shape(nmfsim)[1:]))
    trend_sit = np.zeros((np.shape(nmfsim)[1:]))
    trend_sic = np.zeros((np.shape(nmfsim)[1:]))
    trend_drift = np.zeros((np.shape(nmfsim)[1:]))
    trend_transp = np.zeros((np.shape(nmfsim)[1:]))
    trend_siv_var = np.zeros((np.shape(nmfsim)[1:]))
    trend_eke = np.zeros((np.shape(nmfsim)[1:]))
    trend_era5ke = np.zeros((np.shape(nmfsim)[1:]))

    p_nmfsim = np.zeros((np.shape(nmfsim)[1:]))
    p_ke = np.zeros((np.shape(nmfsim)[1:]))
    p_sit = np.zeros((np.shape(nmfsim)[1:]))
    p_sic = np.zeros((np.shape(nmfsim)[1:]))
    p_drift = np.zeros((np.shape(nmfsim)[1:]))
    p_transp = np.zeros((np.shape(nmfsim)[1:]))
    p_siv_var = np.zeros((np.shape(nmfsim)[1:]))
    p_eke = np.zeros((np.shape(nmfsim)[1:]))
    p_era5ke = np.zeros((np.shape(nmfsim)[1:]))
    for l in range(np.shape(trend_nmfsim)[0]):
        for c in range(np.shape(trend_nmfsim)[1]):
            slope = stats.linregress(date_data_2011_to_2019[:,-1],nmfsim[:,l,c]).slope
            p_nmfsim[l,c] = stats.linregress(date_data_2011_to_2019[:,-1],nmfsim[:,l,c]).pvalue
            trend_nmfsim[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

            slope = stats.linregress(date_data_2011_to_2019[:,-1],ke[:,l,c]).slope
            p_ke[l,c] = stats.linregress(date_data_2011_to_2019[:,-1],ke[:,l,c]).pvalue
            trend_ke[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])
            
            slope = stats.linregress(date_data_2011_to_2019[:,-1],sit[:,l,c]).slope
            p_sit[l,c] = stats.linregress(date_data_2011_to_2019[:,-1],sit[:,l,c]).pvalue
            trend_sit[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])
            
            slope = stats.linregress(date_data_2011_to_2019[:,-1],sic[:,l,c]).slope
            p_sic[l,c] = stats.linregress(date_data_2011_to_2019[:,-1],sic[:,l,c]).pvalue
            trend_sic[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

            slope = stats.linregress(date_data_2011_to_2019[:,-1],drift_magnitude[:,l,c]).slope
            p_drift[l,c] = stats.linregress(date_data_2011_to_2019[:,-1],drift_magnitude[:,l,c]).pvalue
            trend_drift[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

            slope = stats.linregress(date_data_2011_to_2019[:,-1],transp[:,l,c]).slope
            p_transp[l,c] = stats.linregress(date_data_2011_to_2019[:,-1],transp[:,l,c]).pvalue
            trend_transp[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

            slope = stats.linregress(date_data_2011_to_2019[:,-1],siv_var[:,l,c]).slope
            p_siv_var[l,c] = stats.linregress(date_data_2011_to_2019[:,-1],siv_var[:,l,c]).pvalue
            trend_siv_var[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

            slope = stats.linregress(date_data_2011_to_2019[:,-1],eke[:,l,c]).slope
            p_eke[l,c] = stats.linregress(date_data_2011_to_2019[:,-1],eke[:,l,c]).pvalue
            trend_eke[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

            slope = stats.linregress(date_data_2011_to_2019[:,-1],era5ke[:,l,c]).slope
            p_era5ke[l,c] = stats.linregress(date_data_2011_to_2019[:,-1],era5ke[:,l,c]).pvalue
            trend_era5ke[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])
    """ trend = np.polyfit(date_data_2011_to_2019[:,-1],ke[:,-5,-5], deg = 1)
    plt.plot(np.linspace(2011,2020,len(sit)),ke[:,-5,-5], color = 'blue',label = 'kinetic energy')
    plt.plot(np.linspace(2011,2020,len(sit)),[trend[1] + trend[0] * day for day in date_data_2011_to_2019[:,-1]], label = 'linear regression',color = 'red')
    plt.fill_between(np.linspace(2011,2020,len(sit)),ke[:,-5,-5], y2 = [trend[1] + trend[0] * day for day in date_data_2011_to_2019[:,-1]], label = 'eddy kinetic energy')
    plt.xlim(2011,2020)
    plt.xlabel('year')
    plt.ylabel('J/kg')
    plt.legend()
    plt.grid()
    plt.show()
    quit() """

    trend_ke[trend_ke == 0] = np.nan
    trend_nmfsim[trend_nmfsim == 0] = np.nan
    trend_sit[trend_sit == 0] = np.nan
    trend_sic[trend_sic == 0] = np.nan
    trend_drift[trend_drift == 0] = np.nan
    mean_ke[mean_ke == 0] = np.nan
    trend_transp[trend_transp == 0] = np.nan
    trend_siv_var[trend_siv_var == 0] = np.nan
    trend_eke[trend_eke == 0] = np.nan
    trend_era5ke[trend_era5ke == 0] = np.nan

    p_nmfsim[p_nmfsim == 1] = np.nan
    p_ke[p_ke == 1] = np.nan
    p_sit[p_sit == 1] = np.nan
    p_sic[p_sic == 1] = np.nan
    p_drift[p_drift == 1] = np.nan
    p_transp[p_transp == 1] = np.nan
    p_siv_var[p_siv_var == 1] = np.nan
    p_eke[p_eke == 1] = np.nan
    p_era5ke[p_era5ke == 1] = np.nan
    if mode == 'map':
        return trend_nmfsim, trend_ke, trend_sit, trend_sic, trend_drift, trend_transp, trend_siv_var, trend_eke, trend_era5ke, p_nmfsim, p_ke, p_sit, p_sic, p_drift, p_transp, p_siv_var, p_eke, p_era5ke, mean_ke, lon, lat 
    if mode == 'plot':
        return nmfsim, ke, sit, sic, drift_magnitude, transp, siv_var, eke, mean_ke,

def map():
    """
        plot the results of trend
    """
    trend_nmfsim, trend_ke, trend_sit, trend_sic, trend_drift, trend_transp, trend_siv_var, trend_eke,trend_era5ke,p_nmfsim, p_ke, p_sit, p_sic, p_drift, p_transp, p_siv_var, p_eke,p_era5ke, mean_ke, lon, lat = trend(mode = 'map')

    #### - net mass flux of sea ice melt flux trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.set_boundary(rect_in_target)
    axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
    
    axs.coastlines()
    axs.gridlines()
    levels = np.linspace(-0.1,0.1,100)
    cs = axs.contourf(lon, lat, trend_nmfsim, levels=levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p_nmfsim,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    
    cs_ = axs.contour(lon, lat, p_nmfsim,[0.05], colors = 'black',transform=ccrs.PlateCarree())  

    axs.set_title(f"Trend of net mass flux of sea ice melt from 2011 to 2019", fontsize = 22)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.1,0,0.1])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/Maps/Trends/nmfsim.png')
    plt.close()

    #### - kinetic energy trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.set_boundary(rect_in_target)
    axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
    
    axs.coastlines()
    axs.gridlines()

    #trend_ke *= 1e+3 #Passing from J/kg to 1e-3J/kg
    levels = np.linspace(-0.01,0.01,100)
    cs = axs.contourf(lon, lat, trend_ke, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p_ke,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p_ke,[0.05], colors = 'black',transform=ccrs.PlateCarree()) 
    axs.set_title(f"Trend of kinetic energy from 2011 to 2019", fontsize = 22)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.01,0,0.01])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/Maps/Trends/KE.png')
    plt.close()

    #### - SIT trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.set_boundary(rect_in_target)
    axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
    
    axs.coastlines()
    axs.gridlines()

    levels = np.linspace(-1.2,1.2,100)
    cs = axs.contourf(lon, lat, trend_sit, levels=levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p_sit,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p_sit,[0.05], colors = 'black',transform=ccrs.PlateCarree()) 
    axs.set_title(f"Trend of SIT from 2011 to 2019", fontsize = 22)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-1.2,-1,-0.5,0,0.5,1,1.2])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/Maps/Trends/SIT.png')
    plt.close()
    
    #### - SIC trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.set_boundary(rect_in_target)
    axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
    
    axs.coastlines()
    axs.gridlines()
    levels = np.linspace(-0.35,0.35,100)
    cs = axs.contourf(lon, lat, trend_sic, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p_sic,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p_sic,[0.05], colors = 'black',transform=ccrs.PlateCarree()) 

    axs.set_title(f"Trend of SIC from 2011 to 2019", fontsize = 22)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.3,-0.2,-0.1,0,0.1,0.2,0.3])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/Maps/Trends/SIC.png')
    plt.close()

    #### - SID trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.set_boundary(rect_in_target)
    axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
    
    axs.coastlines()
    axs.gridlines()
    levels = np.linspace(-6e3,6e3,200)
    cs = axs.contourf(lon, lat, trend_drift, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p_drift,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p_drift,[0.05], colors = 'black',transform=ccrs.PlateCarree()) 
    
    axs.set_title(f"Trend of SID from 2011 to 2019", fontsize = 22)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-6e3,0,6e3])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/Maps/Trends/SID.png')
    plt.close()
    
    #### - Transport trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.set_boundary(rect_in_target)
    axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    
    axs.coastlines()
    axs.gridlines()
    levels = np.linspace(-0.15,0.15,200)
    cs = axs.contourf(lon, lat, trend_transp, levels = levels,extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p_transp,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p_transp,[0.05], colors = 'black',transform=ccrs.PlateCarree()) 
    
    axs.set_title(f"Trend of net SI transport from 2011 to 2019", fontsize = 22)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.15,0,0.15])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/Maps/Trends/Transp.png')
    plt.close()
   

    #### - siv_var - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.set_boundary(rect_in_target)
    axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
    
    axs.coastlines()
    axs.gridlines()
    levels = np.linspace(-0.004,0.004,200)
    cs = axs.contourf(lon, lat, trend_siv_var,levels = levels,extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p_siv_var,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p_siv_var,[0.05], colors = 'black',transform=ccrs.PlateCarree()) 
    
    axs.set_title(f"Trend of net SI transport from 2011 to 2019", fontsize = 22)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.004,0,0.004])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/Maps/Trends/siv_var.png')
    plt.close()

    #### - EKE - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.set_boundary(rect_in_target)
    axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
    
    axs.coastlines()
    axs.gridlines()
    levels = np.linspace(-0.004,0.004,200)

    """ mask = np.loadtxt('Data/sit/month_mean/2010-10.txt')
    trend_eke = np.where(mask >= 0, trend_eke, np.nan) """
    cs = axs.contourf(lon, lat, trend_eke,levels = levels,extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    
    cs_non_signi = axs.contourf(lon, lat, p_eke,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p_eke,[0.05], colors = 'black',transform=ccrs.PlateCarree()) 
    
    axs.set_title(f"Trend of EKE from 2011 to 2019", fontsize = 22)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.004,0,0.004])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/Maps/Trends/eke.png')
    plt.close()

    #### - ERA5KE - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.set_boundary(rect_in_target)
    axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
    
    axs.coastlines()
    axs.gridlines()
    levels = np.linspace(-3,3,200)


    cs = axs.contourf(lon, lat, trend_era5ke,levels = levels,extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p_era5ke,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p_era5ke,[0.05], colors = 'black',transform=ccrs.PlateCarree()) 
    
    axs.set_title(f"Trend of wind KE from 2011 to 2019", fontsize = 22)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-3,0,3])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/Maps/Trends/era5ke.png')
    plt.close()

def plot():
    nmfsim, ke, sit, sic, drift_magnitude, transp, siv_var, eke, mean_ke = trend(mode='plot')
    lat = np.loadtxt('Data/latitude.txt')
    lon = np.loadtxt('Data/longitude.txt')
    area = np.where((mean_ke > np.nanmean(mean_ke)) & (lat>70) & (lat <80),1,0) # Area where the mean flow is located and between the good lat range
    
    nmfsim = np.where(area == 1,nmfsim,np.nan)
    ke = np.where(area == 1,ke,np.nan)
    sit= np.where(area == 1,sit,np.nan)
    sic= np.where(area == 1,sic,np.nan)
    drift_magnitude= np.where(area == 1,drift_magnitude,np.nan)
    transp= np.where(area == 1,transp,np.nan)
    siv_var= np.where(area == 1,siv_var,np.nan)
    eke= np.where(area == 1,eke,np.nan)
    date_data = np.loadtxt('Data/date_data.txt')
    print(np.shape(ke))
    print(np.shape(date_data))
    date_to_delete = []
    for i in range(len(date_data)):
        if date_data[i,0] == 2010 or date_data[i,0] == 2020:
            date_to_delete.append(i)
    date_data = np.delete(date_data,date_to_delete,axis = 0)
    print(np.shape(date_data))

    # Save the mensual mean of the sit and sic
    current_month_ke = [ke[0]]
    ke_annual = []
    current_month_nmfsim = [ke[0]]
    nmfsim_annual = []
    for i in range(1,len(date_data)):
        current_month_ke.append(ke[i])
        current_month_nmfsim.append(nmfsim[i])

        if i == len(date_data)-1:
            ke_annual.append(np.nanmean(current_month_ke))
            nmfsim_annual.append(np.nanmean(current_month_nmfsim))
            break

        if date_data[i,0] != date_data[i+1,0] : #last day of the month
            ke_annual.append(np.nanmean(current_month_ke))
            nmfsim_annual.append(np.nanmean(current_month_nmfsim))
            current_month_ke = []
            current_month_nmfsim = []

    print(np.shape(ke_annual))
    print(np.shape(nmfsim_annual))

    """r_ke_nmfsim = stats.pearsonr(ke_annual,nmfsim_annual)[0]
    p_ke_nmfsim = stats.pearsonr(ke_annual,nmfsim_annual)[1]

    print(r_ke_nmfsim,p_ke_nmfsim)
    plt.title(f'nmfsim and kinetic energy averaged over the area. r= {r_ke_nmfsim}, p = {p_ke_nmfsim}')
    plt.plot(np.linspace(2011,2020,len(nmfsim_annual)),nmfsim_annual,color = 'blue',label = 'nmfsim')
    plt.twinx()
    plt.plot(np.linspace(2011,2020,len(ke_annual)),ke_annual,color = 'red',label = 'Kinetic energy')
    plt.legend()
    plt.grid()
    plt.xlabel('year')
    plt.ylabel('m^3')
    plt.show()

    r_eke_nmfsim = stats.pearsonr(np.nanmean(eke,axis = (1,2)),np.nanmean(nmfsim,axis = (1,2)))[0]
    p_eke_nmfsim = stats.pearsonr(np.nanmean(eke,axis = (1,2)),np.nanmean(nmfsim,axis = (1,2)))[1]

    print(r_eke_nmfsim,p_eke_nmfsim)
    plt.title(f'nmfsim and eddy kinetic energy averaged over the area. r= {r_eke_nmfsim}, p = {p_eke_nmfsim}')
    plt.plot(np.linspace(2011,2020,len(nmfsim)),np.nanmean(nmfsim,axis = (1,2)),color = 'blue',label = 'nmfsim')
    plt.twinx()
    plt.plot(np.linspace(2011,2020,len(eke)),np.nanmean(eke,axis = (1,2)),color = 'red',label = 'eddy Kinetic energy')
    plt.legend()
    plt.grid()
    plt.xlabel('year')
    plt.ylabel('m^3')
    plt.show()

    plt.title('Kinetic energy averaged over the area')
    plt.plot(np.linspace(2011,2020,len(ke)),np.nanmean(ke,axis = (1,2)))
    plt.grid()
    plt.xlabel('year')
    plt.ylabel('J/kg')
    plt.show()

    plt.title('SIT averaged over the area')
    plt.plot(np.linspace(2011,2020,len(sit)),np.nanmean(sit,axis = (1,2)))
    plt.grid()
    plt.xlabel('year')
    plt.ylabel('m^3')
    plt.show()

    plt.title('Sic averaged over the area')
    plt.plot(np.linspace(2011,2020,len(sic)),np.nanmean(sic,axis = (1,2)))
    plt.grid()
    plt.xlabel('year')
    plt.ylabel('%')
    plt.show()"""
    drift_magnitude = np.where((lat>=75)&(lat<=80.5),drift_magnitude,0)
    plt.close()
    plt.imshow(np.nanmean(drift_magnitude,axis = 0))
    plt.show()
    plt.title('Sea ice drift magnitude averaged over the area')
    plt.plot(np.linspace(2011,2020,len(drift_magnitude)),np.nanmean(drift_magnitude,axis = (1,2)))
    plt.grid()
    plt.xlabel('year')
    plt.ylabel('m/days')
    plt.show()

    plt.title('Net transport averaged over the area')
    plt.plot(np.linspace(2011,2020,len(transp)),np.nanmean(transp,axis = (1,2)))
    plt.grid()
    plt.xlabel('year')
    plt.ylabel('m^3')
    plt.show()

    plt.title('SIV_var averaged over the area')
    plt.plot(np.linspace(2011,2020,len(siv_var)),np.nanmean(siv_var,axis = (1,2)))
    plt.grid()
    plt.xlabel('year')
    plt.ylabel('m^3')
    plt.show()

    plt.title('Eddy kinetic energy averaged over the area')
    plt.plot(np.linspace(2011,2020,len(eke)),np.nanmean(eke,axis = (1,2)))
    plt.grid()
    plt.xlabel('year')
    plt.ylabel('J/kg')
    plt.show()

plot()