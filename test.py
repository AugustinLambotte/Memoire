import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import date, timedelta
import os
# - profile -



iteration = 0
graph, (plot1, plot2) = plt.subplots(1, 2)

for folder in os.listdir('Data/argo_download/coriolis'):
    iteration += 1
    if iteration == 10:
        break
    print(folder)
    for file_name in os.listdir(f'Data/argo_download/coriolis/{folder}/profiles'):
        file = xr.open_dataset(f'C:/Users/Augustin/Desktop/Physique/Master Climatologie/Mémoire/Code 2.0/Data/argo_download/coriolis/{folder}/profiles/{file_name}')
        latitude = np.array(file.LATITUDE.sel(N_PROF = 0))
        longitude = np.array(file.LONGITUDE.sel(N_PROF = 0))
        file.close()
        year = int(file.JULD_LOCATION.sel(N_PROF = 0).dt.year)
        month = int(file.JULD_LOCATION.sel(N_PROF = 0).dt.month)
        day = int(file.JULD_LOCATION.sel(N_PROF = 0).dt.day)
        hour = int(file.JULD_LOCATION.sel(N_PROF = 0).dt.hour)
        minute = int(file.JULD_LOCATION.sel(N_PROF = 0).dt.minute)
        print('#################################')
        print(f"------- Date : {year}-{month}-{day} --------")
        print(f"--- Coords : [{latitude}N ; {longitude} E] ---")
        psal = np.array(file.PSAL.sel(N_PROF = 0))
        temp = np.array(file.TEMP.sel(N_PROF = 0))
        pres = np.array(file.PRES.sel(N_PROF = 0))
        if np.min(psal) < 34.5 or np.max(psal) >35.6 or np.min(pres) < 0 or np.max(temp) > 15 or np.max(pres) > 2000: #Some dataset have incoherent values
            continue
        if latitude >= 66.5: 
            plot1.plot(psal,pres,color = 'blue',alpha = 0.1)
            plot2.plot(temp,pres,color = 'blue',alpha = 0.1)
        if latitude < 66.5: 
            plot1.plot(psal,pres,color = 'red',alpha = 0.1)
            plot2.plot(temp,pres,color = 'red',alpha = 0.1)
plot1.invert_yaxis()
plot2.invert_yaxis()

plot1.grid()
plot2.grid()

plot1.set_title("Practical salinity", fontsize = 10)
plot2.set_title("Temperature", fontsize = 10)

plt.show()
