## This code correlates each SSH time to the average SSH at peak of the anomalies

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib.dates as mdates
import glob
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy import signal
import cv2
from skimage import measure

## load data file names    

all_files = glob.glob('../Raw/Anomalies_R1/*.nc')
all_files.sort()

anomaly_dates = pd.read_csv('Remote_Upwelling_Dates.csv')
#anomaly_dates = pd.to_datetime(anomaly_dates['Date'])
anomaly_dates = pd.to_datetime(anomaly_dates['Remote upwelling date'])
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
    lats = np.flip(lats)

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
temp_sum = np.flip(temp_sum, 0)
temp_mean = temp_sum / len(anomaly_dates)

u_vel_sum = u_vel_anomaly.sum(axis=0)
u_vel_sum = np.flip(u_vel_sum, 0)
u_vel_mean = u_vel_sum / len(anomaly_dates)

v_vel_sum = v_vel_anomaly.sum(axis=0)
v_vel_sum = np.flip(v_vel_sum, 0)
v_vel_mean = v_vel_sum / len(anomaly_dates)

current_spd_mean = np.sqrt(u_vel_mean**2 + v_vel_mean**2)

cc_u_vel = np.empty(len(u_vel_anomaly))
cc_v_vel = np.empty(len(v_vel_anomaly))

k = 0
while k < len(cc_u_vel):
#while k < 1:
    
    u_vel_upwell = u_vel_anomaly[k]
    u_vel_upwell = np.flip(u_vel_upwell, 0)
    u_vel_upwell = u_vel_upwell
    
    xy = u_vel_mean * u_vel_upwell
    xy_avg = np.nanmean(xy)
    
    x_avg = np.nanmean(u_vel_mean)
    y_avg = np.nanmean(u_vel_upwell)
    x_std = np.nanstd(u_vel_mean)
    y_std = np.nanstd(u_vel_upwell)
    
    cc01 = (xy_avg - (x_avg*y_avg))/(x_std*y_std)
    cc_u_vel[k] = cc01
    
    v_vel_upwell = v_vel_anomaly[k]
    v_vel_upwell = np.flip(v_vel_upwell, 0)
    v_vel_upwell = v_vel_upwell
    
    xy = v_vel_mean * v_vel_upwell
    xy_avg = np.nanmean(xy)
    
    x_avg = np.nanmean(v_vel_mean)
    y_avg = np.nanmean(v_vel_upwell)
    x_std = np.nanstd(v_vel_mean)
    y_std = np.nanstd(v_vel_upwell)
    
    cc02 = (xy_avg - (x_avg*y_avg))/(x_std*y_std)
    cc_v_vel[k] = cc02
        
    k+=1

# #anomaly_df = pd.read_csv("../Measured_Temperature/Processed/Anomaly_Dates_2deg.csv")
# anomaly_df = pd.read_csv("../Measured_Temperature/Processed/FFT/Anomaly_Dates_2deg.csv")
# anomaly_date = anomaly_df['Date']
# anomaly_date = [datetime.datetime.strptime(d,"%Y/%m/%d").date() for d in anomaly_date]
# index_list = np.empty(len(anomaly_date))

# j=0
# while j< len(anomaly_date):
#     anomaly_step = anomaly_date[j]
#     anomaly_step = int(anomaly_step.strftime('%Y%m%d'))
#     index = np.abs(date_alt_int - anomaly_step).argmin()
#     index_list[j] = index
#     j+=1
# temp_anomaly = np.empty(len(anomaly_date))
# j=0
# while j < len(anomaly_date):
#     x = int(index_list[j])
#     temp_x = cc_ssh01[x]
#     temp_anomaly[j] = temp_x
#     j+=1

# zero_line = np.zeros(len(date_alt))
# #
# fig01, ax = plt.subplots(figsize=(12,5))
# ax.plot(date_alt,cc_ssh01, 'b', linewidth=0.8)
# #ax.plot(date_alt,cc_ssh02, 'r')
# ax.plot(anomaly_date,temp_anomaly, 'r', linestyle = ' ', marker = 'o', label = 'Anomalies')
# ax.plot(date_alt,zero_line, 'k', linestyle = '--')
# plt.ylabel('Correlation Coefficient')
# formatter = mdates.DateFormatter("%Y")
# ax.xaxis.set_major_formatter(formatter)
# locator = mdates.YearLocator()
# ax.xaxis.set_major_locator(locator)
# plt.setp( ax.xaxis.get_majorticklabels(), rotation=90 )
# plt.ylim((-1,1))
# plt.savefig('SSH_Correlation_R0.PDF')

