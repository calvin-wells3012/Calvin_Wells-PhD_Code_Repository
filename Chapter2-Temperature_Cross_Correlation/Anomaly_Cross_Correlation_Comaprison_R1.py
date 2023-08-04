## This program plots the depth averaged (first 30m from surface) current speed spacially and creates a video

import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
from mpl_toolkits.basemap import Basemap, shiftgrid

## load data file names    

all_files = glob.glob('../Altimetry/Aviso/Raw/*.nc')
all_files.sort()

anomaly_dates = pd.read_csv('../Measured_Temperature/Processed/FFT/*.csv')
anomaly_dates = pd.to_datetime(anomaly_dates['Date'])
startdate = datetime.datetime(1950,1,1,0,0,0)

ssh_anomaly = np.empty([len(anomaly_dates),65,69])
u_vel_anomaly = np.empty([len(anomaly_dates),65,69])
v_vel_anomaly = np.empty([len(anomaly_dates),65,69])

#Loops through each netcdf file
for timestep, datafile in enumerate(all_files): 
    fh = Dataset(datafile,mode='r')
    time = np.array(fh.variables['time'][:])
    lons = fh.variables['longitude'][:]
    lats = fh.variables['latitude'][:]
    lats = np.flip(lats)

    date = [None] * len(time)
    i=0
    while i < len(time):
        time_day = int(np.array(fh.variables['time'][i]))
        day = startdate + datetime.timedelta(days=time_day)
        date[i] = day
        i+=1

    date_int = np.empty(len(date))

    j=0
    while j< len(date):
        time_step = date[j]
        time_step = int(time_step.strftime('%Y%m%d'))
        date_int[j] = time_step
        j+=1
        
    index_list = np.empty(len(anomaly_dates))

    k=0
    while k< len(anomaly_dates):
        anomaly_step = anomaly_dates[k]
        anomaly_step = int(anomaly_step.strftime('%Y%m%d'))
        index = np.abs(date_int - anomaly_step).argmin()
        index_list[k] = index
        k+=1
    
    cc_ssh = np.empty(len(anomaly_dates))
    n = 0
    while n < len(anomaly_dates):
        ssh = np.array(fh.variables['sla'][index_list[n],:,:])
        ssh[ssh < -100] = np.nan
        ssh = np.flip(ssh, 0)
        ssh_anomaly[n] = ssh
        ssh01 = ssh[28:46,8:29]

        
        ssh_norm = (ssh-np.nanmean(ssh))/np.nanstd(ssh)
        
        ssh_mean = pd.read_csv("SSH_mean.csv")
        ssh_mean = np.array(ssh_mean)
        ssh_mean01 = ssh_mean[28:46,8:29]
        xy = ssh_mean01 * ssh01
        xy_avg = np.nanmean(xy)
        x_avg = np.nanmean(ssh_mean01)
        y_avg = np.nanmean(ssh01)
        x_std = np.nanstd(ssh_mean01)
        y_std = np.nanstd(ssh01)
        cc01 = (xy_avg - (x_avg*y_avg))/(x_std*y_std)
        cc_ssh[n] = cc01
        
        n+=1

        fig02 = plt.figure(figsize=(11,7))
        ax = fig02.add_subplot(1,1,1)
        m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax, resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
        border_color = 'black'    
        m.drawcoastlines()    
        plt.imshow(ssh_norm, vmin = -4, vmax = 4, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
        m.fillcontinents(color='grey')
        cbar1 = plt.colorbar()
        cbar1.ax.set_ylabel(r"Normalised SSH [m]",  labelpad=20, rotation=270) 
        plt.plot(32.65, -27.53, 'ro', markersize=5, color = 'deeppink')
        plt.xlabel('Longitude [deg]')
        plt.ylabel('Latitude [deg]')
        plt.xticks(np.arange(30, 45, step=1))
        plt.yticks(np.arange(-34, -20, step=1))   
        plt.plot([31, 36, 36, 31, 31], [-30, -30, -26, -26, -30], color='black', linestyle = 'dashed')
        plt.tight_layout() 
        plt.savefig('SSH_Correlation_Anomaly'+str(n)+'.png')    


ssh_mean_norm = (ssh_mean-np.nanmean(ssh_mean))/np.nanstd(ssh_mean)

##Plotting    
fig01 = plt.figure(figsize=(11,7))
ax = fig01.add_subplot(1,1,1)
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax, resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
border_color = 'black'    
m.drawcoastlines()    
plt.imshow(ssh_mean_norm, vmin = -4, vmax = 4, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
m.fillcontinents(color='grey')
cbar1 = plt.colorbar()
cbar1.ax.set_ylabel(r"Normalised SSH (conditional average) [m]",  labelpad=20, rotation=270) 
plt.plot(32.65, -27.53, 'ro', markersize=5, color = 'deeppink')
plt.xlabel('Longitude [deg]')
plt.ylabel('Latitude [deg]')
plt.xticks(np.arange(30, 45, step=1))
plt.yticks(np.arange(-34, -20, step=1))
plt.plot([31, 36, 36, 31, 31], [-30, -30, -26, -26, -30], color='black', linestyle = 'dashed')
plt.tight_layout()
plt.savefig('SSH_CondAvg_Alt_2deg_peak_regional.png')    


    

    
    
    
    
    
    
    
    
    
    
    
    

        
