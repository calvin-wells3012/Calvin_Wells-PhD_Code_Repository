## This program plots the depth averaged (first 30m from surface) current speed spacially and creates a video
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import matplotlib.colors
from matplotlib.colors import LinearSegmentedColormap


## load data file names    

all_files = glob.glob('../Raw/Anomalies_R1/Local_Upwelling/*.nc')
all_files.sort()

anomaly_dates = pd.read_csv('Local_Upwelling_Dates.csv')
anomaly_dates = pd.to_datetime(anomaly_dates['Date'])
startdate = datetime.datetime(1950,1,1,0,0,0)

temp_anomaly = np.empty([len(anomaly_dates),133,121])
u_vel_anomaly = np.empty([len(anomaly_dates),133,121])
v_vel_anomaly = np.empty([len(anomaly_dates),133,121])

i = 0
while i < len(all_files):

    
    fh = Dataset(all_files[i],mode='r')
    time = np.array(fh.variables['time'][:])
    lons = fh.variables['longitude'][:]
    lats = fh.variables['latitude'][:]

    date = [None] * len(time)
    k=0
    while k < len(time):
        time_day = int(np.array(fh.variables['time'][k]))
        day = startdate + datetime.timedelta(hours=time_day)
        date[k] = day
        k+=1

    date_int = np.empty(len(date))

    j=0
    while j< len(date):
        time_step = date[j]
        time_step = int(time_step.strftime('%Y%m%d'))
        date_int[j] = time_step
        j+=1
        
    index_list = np.empty(len(anomaly_dates))

    anomaly_step = anomaly_dates[i]
    anomaly_step = int(anomaly_step.strftime('%Y%m%d'))
    index = np.abs(date_int - anomaly_step).argmin()
    index_list[i] = index
    
    temp = np.array(fh.variables['thetao'][index,4,:,:])
    u_vel = np.array(fh.variables['uo'][index,4,:,:])
    v_vel = np.array(fh.variables['vo'][index,4,:,:])
    temp[temp < -100] = np.nan
    
    temp_mean = np.nanmean(temp)
    
    temp = temp - temp_mean
    
    temp_anomaly[i] = temp
    u_vel[u_vel < -100] = np.nan
    u_vel_anomaly[i] = u_vel
    v_vel[v_vel < -100] = np.nan
    v_vel_anomaly[i] = v_vel
        
    i+=1

temp_sum = temp_anomaly.sum(axis=0)
temp_mean = temp_sum / len(anomaly_dates)

u_vel_sum = u_vel_anomaly.sum(axis=0)
u_vel_mean = u_vel_sum / len(anomaly_dates)

v_vel_sum = v_vel_anomaly.sum(axis=0)
v_vel_mean = v_vel_sum / len(anomaly_dates)

cc = np.zeros([len(anomaly_dates),23,25])
cc_max = np.zeros(len(anomaly_dates))
cc_mean = np.zeros(len(anomaly_dates))

i = 0
while i < len(anomaly_dates):
    
    ai = np.array(u_vel_mean[29:52,20:45])
    bj = np.array(v_vel_mean[29:52,20:45])

    ci = np.array(u_vel_anomaly[i,29:52,20:45])
    dj = np.array(v_vel_anomaly[i,29:52,20:45])

    cc_01 = ((ai * ci) + (bj * dj)) / (np.sqrt(np.square(ai)+np.square(bj)) * np.sqrt(np.square(ci)+np.square(dj)))
    cc[i] = cc_01

    cc_max_01 = np.nanmax(cc_01)
    cc_mean_01 = np.nanmean(cc_01)
    
    cc_max[i] = cc_max_01
    cc_mean[i] = cc_mean_01
    
    i+=1
    
    
    
    
    
    
    
    
    
    
    

        
