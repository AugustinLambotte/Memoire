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

"""
    This script compute the correlation coefficient for every pixels (80km^2) over 
    the EGC between the intensity of the EGC (KE) and the fresh water flux.
"""
""" 
def correlation():
    
    #    This function compute the correlation coefficient for each cells based on data stored in Data/bw
    
    ke = []
    nmfsim = []
    siv_var = []
    transp = []
    eke = []
    era5ke = []
    for file in os.listdir('Data/nmfsim'):
        if file == 'month_mean':
            continue
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2019:
            ke.append(np.loadtxt(f'Data/gos/ke/{file}'))
            nmfsim.append(np.loadtxt(f'Data/nmfsim/{file}'))
            siv_var.append(np.loadtxt(f'Data/siv_var/{file}'))
            transp.append(np.loadtxt(f'Data/transport/{file}'))
            eke.append(np.loadtxt(f'Data/gos/eke/{file}'))
            era5ke.append(np.loadtxt(f'Data/ERA5/interpoalted/ke/{file}'))
    ke = np.nan_to_num(np.array(ke))
    nmfsim = np.nan_to_num(np.array(nmfsim))
    siv_var = np.nan_to_num(np.array(siv_var))
    transp = np.nan_to_num(np.array(transp))
    eke = np.nan_to_num(np.array(eke))
    era5ke = np.nan_to_num(np.array(era5ke))

    # WITH KE
    corr_ke_nmfsim = np.zeros(np.shape(ke[0]))
    p_value_ke_nmfsim = np.zeros(np.shape(ke[0]))

    corr_ke_siv_var = np.zeros(np.shape(ke[0]))
    p_value_ke_siv_var = np.zeros(np.shape(ke[0]))

    corr_ke_transp = np.zeros(np.shape(ke[0]))
    p_value_ke_transp = np.zeros(np.shape(ke[0]))

    for line in range(np.shape(corr_ke_nmfsim)[0]):
        for col in range(np.shape(corr_ke_nmfsim)[1]):
            corr_ke_nmfsim[line,col] = stats.pearsonr(ke[:,line,col],nmfsim[:,line,col])[0]
            p_value_ke_nmfsim[line,col] = stats.pearsonr(ke[:,line,col],nmfsim[:,line,col])[1]

            corr_ke_siv_var[line,col] = stats.pearsonr(ke[:,line,col],siv_var[:,line,col])[0]
            p_value_ke_siv_var[line,col] = stats.pearsonr(ke[:,line,col],siv_var[:,line,col])[1]

            corr_ke_transp[line,col] = stats.pearsonr(ke[:,line,col],transp[:,line,col])[0]
            p_value_ke_transp[line,col] = stats.pearsonr(ke[:,line,col],transp[:,line,col])[1]

    np.savetxt(f'Data/correlations/ke_nmfsim/r.txt',corr_ke_nmfsim)
    np.savetxt(f'Data/correlations/ke_nmfsim/p.txt',p_value_ke_nmfsim)

    np.savetxt(f'Data/correlations/ke_siv_var/r.txt',corr_ke_siv_var)
    np.savetxt(f'Data/correlations/ke_siv_var/p.txt',p_value_ke_siv_var)

    np.savetxt(f'Data/correlations/ke_transp/r.txt',corr_ke_transp)
    np.savetxt(f'Data/correlations/ke_transp/p.txt',p_value_ke_transp)

    # WITH EKE
    corr_ke_nmfsim = np.zeros(np.shape(ke[0]))
    p_value_ke_nmfsim = np.zeros(np.shape(ke[0]))

    corr_ke_siv_var = np.zeros(np.shape(ke[0]))
    p_value_ke_siv_var = np.zeros(np.shape(ke[0]))

    corr_ke_transp = np.zeros(np.shape(ke[0]))
    p_value_ke_transp = np.zeros(np.shape(ke[0]))
    
    for line in range(np.shape(corr_ke_nmfsim)[0]):
        for col in range(np.shape(corr_ke_nmfsim)[1]):
            corr_ke_nmfsim[line,col] = stats.pearsonr(eke[:,line,col],nmfsim[:,line,col])[0]
            p_value_ke_nmfsim[line,col] = stats.pearsonr(eke[:,line,col],nmfsim[:,line,col])[1]

            corr_ke_siv_var[line,col] = stats.pearsonr(eke[:,line,col],siv_var[:,line,col])[0]
            p_value_ke_siv_var[line,col] = stats.pearsonr(eke[:,line,col],siv_var[:,line,col])[1]

            corr_ke_transp[line,col] = stats.pearsonr(eke[:,line,col],transp[:,line,col])[0]
            p_value_ke_transp[line,col] = stats.pearsonr(eke[:,line,col],transp[:,line,col])[1]

    np.savetxt(f'Data/correlations/eke_nmfsim/r.txt',corr_ke_nmfsim)
    np.savetxt(f'Data/correlations/eke_nmfsim/p.txt',p_value_ke_nmfsim)

    np.savetxt(f'Data/correlations/eke_siv_var/r.txt',corr_ke_siv_var)
    np.savetxt(f'Data/correlations/eke_siv_var/p.txt',p_value_ke_siv_var)

    np.savetxt(f'Data/correlations/eke_transp/r.txt',corr_ke_transp)
    np.savetxt(f'Data/correlations/eke_transp/p.txt',p_value_ke_transp)
"""

def correlation_annual():
    """
        This function compute the correlation coefficient for each cells based on data stored in Data/bw.
        the difference with correlation() is that a running average over one year is apply on the data.
    """
    ke = []
    nmfsim = []
    siv_var = []
    transp = []
    eke = []
    era5ke = []
    for file in os.listdir('Data/nmfsim'):
        if file == 'month_mean':
            continue
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2019:
            ke.append(np.loadtxt(f'Data/gos/ke/{file}'))
            nmfsim.append(np.loadtxt(f'Data/nmfsim/{file}'))
            siv_var.append(np.loadtxt(f'Data/siv_var/{file}'))
            transp.append(np.loadtxt(f'Data/transport/{file}'))
            eke.append(np.loadtxt(f'Data/gos/eke/{file}'))
            era5ke.append(np.loadtxt(f'Data/ERA5/interpolated/ke/{file}'))
    ke = np.nan_to_num(np.array(ke))
    nmfsim = np.nan_to_num(np.array(nmfsim))
    siv_var = np.nan_to_num(np.array(siv_var))
    transp = np.nan_to_num(np.array(transp))
    eke = np.nan_to_num(np.array(eke))
    era5ke = np.nan_to_num(np.array(era5ke))
    
    
    #Extracting date of data from 2011 to 2019
    date_data = np.loadtxt('Data/date_data.txt')
    date_data_2011_to_2019 = []
    for date_ in date_data:
        if date_[0] >= 2011 and date_[0] <2020:
            date_data_2011_to_2019.append(date_)
    date_data_2011_to_2019 = np.array(date_data_2011_to_2019)
    
    ### - One year running average - ###

    # Storing the running mean in 'name'_averaged
    ke_averaged = []
    nmfsim_averaged = []
    eke_averaged = []
    era5ke_averaged = []

    for i in range(len(date_data_2011_to_2019)):
        if date_data_2011_to_2019[i,0] > 2011:
            year = int(date_data_2011_to_2019[i,0])
            month = int(date_data_2011_to_2019[i,1])
            day = int(date_data_2011_to_2019[i,2])
            
            ke_last_year = []
            nmfsim_last_year = []
            eke_last_year = []
            era5ke_last_year = []
            for j in range(i):
                year_ = int(date_data_2011_to_2019[j,0])
                month_ = int(date_data_2011_to_2019[j,1])
                day_ = int(date_data_2011_to_2019[j,2])
                time_dif = (date(year,month,day) - date(year_,month_,day_)).days

                if time_dif <= 365:
                    ke_last_year.append(ke[j])
                    nmfsim_last_year.append(nmfsim[j])
                    eke_last_year.append(eke[j])
                    era5ke_last_year.append(era5ke[j])

            ke_averaged.append(np.nanmean(ke_last_year,axis = 0))
            nmfsim_averaged.append(np.nanmean(nmfsim_last_year,axis = 0))
            eke_averaged.append(np.nanmean(eke_last_year,axis = 0))
            era5ke_averaged.append(np.nanmean(era5ke_last_year,axis = 0))
    ke_averaged = np.array(ke_averaged)
    nmfsim_averaged = np.array(nmfsim_averaged)
    eke_averaged = np.array(eke_averaged)
    era5ke_averaged = np.array(era5ke_averaged)
    
    # Freshwater flux vs KE
    corr_ke_nmfsim = np.zeros(np.shape(ke_averaged[0]))
    p_value_ke_nmfsim = np.zeros(np.shape(ke_averaged[0]))

    for line in range(np.shape(corr_ke_nmfsim)[0]):
        for col in range(np.shape(corr_ke_nmfsim)[1]):
            corr_ke_nmfsim[line,col] = stats.pearsonr(ke_averaged[:,line,col],nmfsim_averaged[:,line,col])[0]
            p_value_ke_nmfsim[line,col] = stats.pearsonr(ke_averaged[:,line,col],nmfsim_averaged[:,line,col])[1]

    np.savetxt(f'Data/correlations/averaged/ke_nmfsim/r.txt',corr_ke_nmfsim)
    np.savetxt(f'Data/correlations/averaged/ke_nmfsim/p.txt',p_value_ke_nmfsim)

    # Freshwater flux vs EKE
    corr_eke_nmfsim = np.zeros(np.shape(eke_averaged[0]))
    p_value_eke_nmfsim = np.zeros(np.shape(eke_averaged[0]))

    for line in range(np.shape(corr_eke_nmfsim)[0]):
        for col in range(np.shape(corr_eke_nmfsim)[1]):
            corr_eke_nmfsim[line,col] = stats.pearsonr(eke_averaged[:,line,col],nmfsim_averaged[:,line,col])[0]
            p_value_eke_nmfsim[line,col] = stats.pearsonr(eke_averaged[:,line,col],nmfsim_averaged[:,line,col])[1]

    np.savetxt(f'Data/correlations/averaged/eke_nmfsim/r.txt',corr_eke_nmfsim)
    np.savetxt(f'Data/correlations/averaged/eke_nmfsim/p.txt',p_value_eke_nmfsim)

    # Ke vs EKE
    corr_ke_eke = np.zeros(np.shape(eke_averaged[0]))
    p_value_ke_eke = np.zeros(np.shape(eke_averaged[0]))

    for line in range(np.shape(corr_eke_nmfsim)[0]):
        for col in range(np.shape(corr_eke_nmfsim)[1]):
            corr_ke_eke[line,col] = stats.pearsonr(eke_averaged[:,line,col],ke_averaged[:,line,col])[0]
            p_value_ke_eke[line,col] = stats.pearsonr(eke_averaged[:,line,col],ke_averaged[:,line,col])[1]

    np.savetxt(f'Data/correlations/averaged/ke_eke/r.txt',corr_ke_eke)
    np.savetxt(f'Data/correlations/averaged/ke_eke/p.txt',p_value_ke_eke)

    # WITH ERA5_KE
    corr_era5ke_ke = np.zeros(np.shape(era5ke_averaged[0]))
    p_value_era5ke_ke = np.zeros(np.shape(era5ke_averaged[0]))

    corr_era5ke_eke = np.zeros(np.shape(era5ke_averaged[0]))
    p_value_era5ke_eke = np.zeros(np.shape(era5ke_averaged[0]))

    for line in range(np.shape(corr_era5ke_ke)[0]):
        for col in range(np.shape(corr_era5ke_ke)[1]):
            corr_era5ke_ke[line,col] = stats.pearsonr(era5ke_averaged[:,line,col],ke_averaged[:,line,col])[0]
            p_value_era5ke_ke[line,col] = stats.pearsonr(era5ke_averaged[:,line,col],ke_averaged[:,line,col])[1]

            corr_era5ke_eke[line,col] = stats.pearsonr(era5ke_averaged[:,line,col],eke_averaged[:,line,col])[0]
            p_value_era5ke_eke[line,col] = stats.pearsonr(era5ke_averaged[:,line,col],eke_averaged[:,line,col])[1]

    np.savetxt(f'Data/correlations/averaged/era5ke_ke/r.txt',corr_era5ke_ke)
    np.savetxt(f'Data/correlations/averaged/era5ke_ke/p.txt',p_value_era5ke_ke)

    np.savetxt(f'Data/correlations/averaged/era5ke_eke/r.txt',corr_era5ke_eke)
    np.savetxt(f'Data/correlations/averaged/era5ke_eke/p.txt',p_value_era5ke_eke)

def correlation_annual_lag_time():
    """
        This function compute the correlation coefficient for each cells based on data stored in Data/bw.
        the difference with correlation() is that a running average over one year is apply on the data.
    """
    ke = []
    nmfsim = []
    siv_var = []
    transp = []
    eke = []
    for file in os.listdir('Data/nmfsim'):
        if file == 'month_mean':
            continue
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2019:
            ke.append(np.loadtxt(f'Data/gos/ke/{file}'))
            nmfsim.append(np.loadtxt(f'Data/nmfsim/{file}'))
            siv_var.append(np.loadtxt(f'Data/siv_var/{file}'))
            transp.append(np.loadtxt(f'Data/transport/{file}'))
            eke.append(np.loadtxt(f'Data/gos/eke/{file}'))
    
    ke = np.nan_to_num(np.array(ke))
    nmfsim = np.nan_to_num(np.array(nmfsim))
    siv_var = np.nan_to_num(np.array(siv_var))
    transp = np.nan_to_num(np.array(transp))
    eke = np.nan_to_num(np.array(eke))
    #Extracting date of data from 2011 to 2019
    date_data = np.loadtxt('Data/date_data.txt')
    date_data_2011_to_2019 = []
    for date_ in date_data:
        if date_[0] >= 2011 and date_[0] <2020:
            date_data_2011_to_2019.append(date_)
    date_data_2011_to_2019 = np.array(date_data_2011_to_2019)
    
    #One year running average.

    # Storing the running mean in 'name'_averaged
    ke_averaged = []
    nmfsim_averaged = []
    eke_averaged = []

    for i in range(len(date_data_2011_to_2019)):
        if date_data_2011_to_2019[i,0] > 2011:
            year = int(date_data_2011_to_2019[i,0])
            month = int(date_data_2011_to_2019[i,1])
            day = int(date_data_2011_to_2019[i,2])
            
            ke_last_year = []
            nmfsim_last_year = []
            eke_last_year = []
            for j in range(i):
                year_ = int(date_data_2011_to_2019[j,0])
                month_ = int(date_data_2011_to_2019[j,1])
                day_ = int(date_data_2011_to_2019[j,2])
                time_dif = (date(year,month,day) - date(year_,month_,day_)).days

                if time_dif <= 365:
                    ke_last_year.append(ke[j])
                    nmfsim_last_year.append(nmfsim[j])
                    eke_last_year.append(eke[j])

            ke_averaged.append(np.nanmean(ke_last_year,axis = 0))
            nmfsim_averaged.append(np.nanmean(nmfsim_last_year,axis = 0))
            eke_averaged.append(np.nanmean(eke_last_year,axis = 0))
    
    ke_averaged = np.array(ke_averaged)
    nmfsim_averaged = np.array(nmfsim_averaged)
    eke_averaged = np.array(eke_averaged)
    
    #Impose lag time
    ke_averaged = ke_averaged[24:]
    eke_averaged = eke_averaged[24:]
    nmfsim_averaged = nmfsim_averaged[:-24]
    # Freshwater flux vs KE
    corr_ke_nmfsim = np.zeros(np.shape(ke_averaged[0]))
    p_value_ke_nmfsim = np.zeros(np.shape(ke_averaged[0]))

    for line in range(np.shape(corr_ke_nmfsim)[0]):
        for col in range(np.shape(corr_ke_nmfsim)[1]):
            corr_ke_nmfsim[line,col] = stats.pearsonr(ke_averaged[:,line,col],nmfsim_averaged[:,line,col])[0]
            p_value_ke_nmfsim[line,col] = stats.pearsonr(ke_averaged[:,line,col],nmfsim_averaged[:,line,col])[1]

    np.savetxt(f'Data/correlations/lag/averaged/ke_nmfsim/r.txt',corr_ke_nmfsim)
    np.savetxt(f'Data/correlations/lag/averaged/ke_nmfsim/p.txt',p_value_ke_nmfsim)

    # Freshwater flux vs EKE
    corr_eke_nmfsim = np.zeros(np.shape(eke_averaged[0]))
    p_value_eke_nmfsim = np.zeros(np.shape(eke_averaged[0]))

    for line in range(np.shape(corr_eke_nmfsim)[0]):
        for col in range(np.shape(corr_eke_nmfsim)[1]):
            corr_eke_nmfsim[line,col] = stats.pearsonr(eke_averaged[:,line,col],nmfsim_averaged[:,line,col])[0]
            p_value_eke_nmfsim[line,col] = stats.pearsonr(eke_averaged[:,line,col],nmfsim_averaged[:,line,col])[1]

    np.savetxt(f'Data/correlations/lag/averaged/eke_nmfsim/r.txt',corr_eke_nmfsim)
    np.savetxt(f'Data/correlations/lag/averaged/eke_nmfsim/p.txt',p_value_eke_nmfsim)

def plot_r_p_annual():
    """
        This function plots and save correlation map based on data in Data/bw/r_... and Data/bw/p_...
    """
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')
    ke = []
    for file in os.listdir('Data/gos/ke/month_mean'):
        ke.append(np.loadtxt('Data/gos/ke/month_mean/' + file))
    ke = np.array(ke)
    mean_ke = np.nanmean(ke,axis = 0)

    ######## - KE vs nmfsim - #########
    r = np.loadtxt('Data\correlations/averaged/ke_nmfsim/r.txt')
    p = np.loadtxt('Data\correlations/averaged/ke_nmfsim\p.txt')


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
    levels = np.linspace(-0.6,0.6,1000)
    
    cs = axs.contourf(lon, lat, r, extend = 'both', levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'black',transform=ccrs.PlateCarree())

    axs.set_title("Correlation between kinetic energy and freshwater flux", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.6,0,0.6])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/averaged/ke_nmfsim.png")
    plt.close()

    ######## - EKE vs nmfsim - #########
    r = np.loadtxt('Data\correlations/averaged/eke_nmfsim/r.txt')
    p = np.loadtxt('Data\correlations/averaged/eke_nmfsim\p.txt')
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
    levels = np.linspace(-0.6,0.6,100)
    cs = axs.contourf(lon, lat, r, extend = 'both', levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'black',transform=ccrs.PlateCarree())    

    axs.set_title("Correlation between eddy kinetic energy and freshwater flux", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.6,0,0.6])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/averaged/eke_nmfsim.png")

    ######## - KE vs EKE - #########
    r = np.loadtxt('Data\correlations/averaged/ke_eke/r.txt')
    p = np.loadtxt('Data\correlations/averaged/ke_eke\p.txt')
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
    levels = np.linspace(-1,1,1000)
    cs = axs.contourf(lon, lat, r, extend = 'both', levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'black',transform=ccrs.PlateCarree())    

    axs.set_title("Correlation between eddy kinetic energy and kinetic energy", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-1,0,1])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/averaged/ke_eke.png")

    ######## - KE vs ERA5KE - #########
    r = np.loadtxt('Data\correlations/averaged/era5ke_ke/r.txt')
    p = np.loadtxt('Data\correlations/averaged/era5ke_ke\p.txt')
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
    levels = np.linspace(-1,1,1000)
    cs = axs.contourf(lon, lat, r, extend = 'both', levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'black',transform=ccrs.PlateCarree())    

    axs.set_title("Correlation between wind kinetic energy and gos kinetic energy", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-1,0,1])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/averaged/era5ke_ke.png")

    ######## - EKE vs ERA5KE - #########
    r = np.loadtxt('Data\correlations/averaged/era5ke_eke/r.txt')
    p = np.loadtxt('Data\correlations/averaged/era5ke_eke\p.txt')
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
    levels = np.linspace(-1,1,1000)
    cs = axs.contourf(lon, lat, r, extend = 'both', levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'black',transform=ccrs.PlateCarree())    

    axs.set_title("Correlation between wind kinetic energy and gos eddy kinetic energy", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-1,0,1])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/averaged/era5ke_eke.png")

def plot_r_p_annual_lag():
    """
        This function plots and save correlation map based on data in Data/bw/r_... and Data/bw/p_...
    """
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')
    ke = []
    for file in os.listdir('Data/gos/ke/month_mean'):
        ke.append(np.loadtxt('Data/gos/ke/month_mean/' + file))
    ke = np.array(ke)
    mean_ke = np.nanmean(ke,axis = 0)

    ######## - KE vs nmfsim - #########
    r = np.loadtxt('Data\correlations/lag/averaged/ke_nmfsim/r.txt')
    p = np.loadtxt('Data\correlations/lag/averaged/ke_nmfsim\p.txt')


    
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
    levels = np.linspace(-0.6,0.6,1000)
    
    cs = axs.contourf(lon, lat, r, extend = 'both', levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    
    cs_non_signi = axs.contourf(lon, lat, p,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'black',transform=ccrs.PlateCarree())

    axs.set_title("Correlation between kinetic energy and freshwater flux with two years lag time", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.6,0,0.6])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/lag/averaged/ke_nmfsim.png")
    
    ######## - EKE vs nmfsim - #########
    r = np.loadtxt('Data\correlations/lag/averaged/eke_nmfsim/r.txt')
    p = np.loadtxt('Data\correlations/lag/averaged/eke_nmfsim\p.txt')
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
    levels = np.linspace(-0.6,0.6,1000)
    cs = axs.contourf(lon, lat, r, extend = 'both', levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_non_signi = axs.contourf(lon, lat, p,[0,0.05,1],hatches = ['','//'],colors = 'none', transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'black',transform=ccrs.PlateCarree())    

    axs.set_title("Correlation between eddy kinetic energy and freshwater flux with two years lag time", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.6,0,0.6])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/lag/averaged/eke_nmfsim.png")

def plot_r_p():
    """
        This function plots and save correlation map based on data in Data/bw/r_... and Data/bw/p_...
    """
    lon = np.loadtxt('Data/longitude.txt')
    lat = np.loadtxt('Data/latitude.txt')
    ke = []
    for file in os.listdir('Data/gos/ke/month_mean'):
        ke.append(np.loadtxt('Data/gos/ke/month_mean/' + file))
    ke = np.array(ke)
    mean_ke = np.nanmean(ke,axis = 0)

    ######## - KE vs nmfsim - #########
    r = np.loadtxt('Data\correlations\ke_nmfsim/r.txt')
    p = np.loadtxt('Data\correlations\ke_nmfsim\p.txt')
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
    levels = np.linspace(-0.25,0.25,1000)
    cs = axs.contourf(lon, lat, r, extend = 'both', levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red',transform=ccrs.PlateCarree())
    #cs_ = axs.contour(lon, lat, p,[0.01], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    axs.set_title("Correlation between kinetic energy and freshwater flux", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.25,0,0.25])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/ke_nmfsim.png")

    ######## - KE vs siv_variation - #########
    r = np.loadtxt('Data\correlations\ke_siv_var/r.txt')
    p = np.loadtxt('Data\correlations\ke_siv_var\p.txt')
    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -18))
    #fig, axs = plt.plots(nrows = 1, ncols = 1, figsize=figsize, subplot_kw={'projection': projection})
    #axs.set_extent([-47, 16, 60, 85], crs = ccrs.PlateCarree())
    
    xlim = [-43, 16]
    ylim = [61, 81]
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
    levels = np.linspace(-0.5,0.5,1000)
    cs = axs.contourf(lon, lat, r, levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
    axs.set_title("Correlation between kinetic energy and siv variation", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.5,0,0.5])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/ke_siv_var.png")

    ######## - KE vs transport - #########
    r = np.loadtxt('Data\correlations\ke_transp/r.txt')
    p = np.loadtxt('Data\correlations\ke_transp\p.txt')
    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -18))
    
    xlim = [-43, 16]
    ylim = [61, 81]
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
    levels = np.linspace(-0.5,0.5,1000)
    cs = axs.contourf(lon, lat, r, levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
    axs.set_title("Correlation between kinetic energy and net sea ice transport", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.5,0,0.5])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/ke_transp.png")

    
    ######## - EKE vs nmfsim - #########
    r = np.loadtxt('Data\correlations\eke_nmfsim/r.txt')
    p = np.loadtxt('Data\correlations\eke_nmfsim\p.txt')
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
    levels = np.linspace(-0.25,0.25,1000)
    cs = axs.contourf(lon, lat, r, extend = 'both', levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red',transform=ccrs.PlateCarree())
    #cs_ = axs.contour(lon, lat, p,[0.01], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    axs.set_title("Correlation between eddies kinetic energy and freshwater flux", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.25,0,0.25])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/eke_nmfsim.png")

    ######## - EKE vs siv_variation - #########
    r = np.loadtxt('Data\correlations\eke_siv_var/r.txt')
    p = np.loadtxt('Data\correlations\eke_siv_var\p.txt')
    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -18))
    #fig, axs = plt.plots(nrows = 1, ncols = 1, figsize=figsize, subplot_kw={'projection': projection})
    #axs.set_extent([-47, 16, 60, 85], crs = ccrs.PlateCarree())
    
    xlim = [-43, 16]
    ylim = [61, 81]
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
    levels = np.linspace(-0.5,0.5,1000)
    cs = axs.contourf(lon, lat, r, levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
    axs.set_title("Correlation between eddies kinetic energy and siv variation", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.5,0,0.5])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/eke_siv_var.png")

    ######## - EKE vs transport - #########
    r = np.loadtxt('Data\correlations\eke_transp/r.txt')
    p = np.loadtxt('Data\correlations\eke_transp\p.txt')
    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -18))
    
    xlim = [-43, 16]
    ylim = [61, 81]
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
    levels = np.linspace(-0.5,0.5,1000)
    cs = axs.contourf(lon, lat, r,levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
    axs.set_title("Correlation between eddies kinetic energy and net sea ice transport", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.5,0,0.5])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots\Maps\Correlations/eke_transp.png")

correlation_annual()
plot_r_p_annual()
#plot_r_p_annual() 
#plot_r_p()