## This code correlates each SSH time to the average SSH at peak of the anomalies

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib.dates as mdates
import glob
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, shiftgrid
from PyEMD import EMD, Visualisation
from scipy import signal
import cv2
from skimage import measure

fh = Dataset('../Altimetry/Aviso/Raw/*.nc',mode='r')
startdate = datetime.datetime(1950,1,1,0,0,0)
time_alt = np.array(fh.variables['time'][:])
lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
lats_flip = np.flip(lats)

date_alt = [None] * len(time_alt)
i=0
while i < len(time_alt):
    time_day = int(np.array(fh.variables['time'][i]))
    day = startdate + datetime.timedelta(days=time_day)
    date_alt[i] = day
    i+=1
#date_alt = date_alt[70:7669]

date_alt_int = np.empty(len(date_alt))

i=0
while i < len(date_alt):
    date_dly = int(date_alt[i].strftime('%Y%m%d'))
    date_alt_int[i] = date_dly
    i+=1

ssh_mean = pd.read_csv("SSH_mean.csv")
#ssh_mean = ssh_mean.fillna(0)
ssh_mean = np.array(ssh_mean)
ssh_mean = ssh_mean[28:46,8:29]
#ssh_mean[ssh_mean > 0] = 0

##Load in altimetry SSH


cc_ssh01 = np.empty(len(date_alt))
cc_ssh02 = np.empty(len(date_alt))

k = 0
while k < len(cc_ssh01):
    ssh01 = np.array(fh.variables['sla'][k,:,:])
    ssh01[ssh01 < -100] = np.nan
    ssh01 = np.flip(ssh01, 0)
    ssh01 = ssh01[28:46,8:29]
    
    xy = ssh_mean * ssh01
    xy_avg = np.nanmean(xy)
    
    x_avg = np.nanmean(ssh_mean)
    y_avg = np.nanmean(ssh01)
    x_std = np.nanstd(ssh_mean)
    y_std = np.nanstd(ssh01)
    
    cc01 = (xy_avg - (x_avg*y_avg))/(x_std*y_std)
    cc_ssh01[k] = cc01
    
    x_sum = np.nansum(ssh01)
    y_sum = np.nansum(ssh_mean)
    xy_sum = np.nansum(xy)
    x_2_sum = np.nansum(ssh01**2)
    y_2_sum = np.nansum(ssh_mean**2)
    n = 335
    
    cc02 = ((n*xy_sum)-(x_sum*y_sum))/np.sqrt((n*x_2_sum-x_sum**2)*(n*y_2_sum-y_sum**2))
    cc_ssh02[k] = cc02
        
    k+=1

anomaly_df = pd.read_csv("../Measured_Temperature/Processed/FFT/*.csv")
anomaly_date = anomaly_df['Date']
anomaly_date = [datetime.datetime.strptime(d,"%Y/%m/%d").date() for d in anomaly_date]
index_list = np.empty(len(anomaly_date))

j=0
while j< len(anomaly_date):
    anomaly_step = anomaly_date[j]
    anomaly_step = int(anomaly_step.strftime('%Y%m%d'))
    index = np.abs(date_alt_int - anomaly_step).argmin()
    index_list[j] = index
    j+=1
temp_anomaly = np.empty(len(anomaly_date))
j=0
while j < len(anomaly_date):
    x = int(index_list[j])
    temp_x = cc_ssh01[x]
    temp_anomaly[j] = temp_x
    j+=1

zero_line = np.zeros(len(date_alt))
#
fig01, ax = plt.subplots(figsize=(12,5))
ax.plot(date_alt,cc_ssh01, 'b', linewidth=0.8)
#ax.plot(date_alt,cc_ssh02, 'r')
ax.plot(anomaly_date,temp_anomaly, 'r', linestyle = ' ', marker = 'o', label = 'Anomalies')
ax.plot(date_alt,zero_line, 'k', linestyle = '--')
plt.ylabel('Correlation Coefficient')
formatter = mdates.DateFormatter("%Y")
ax.xaxis.set_major_formatter(formatter)
locator = mdates.YearLocator()
ax.xaxis.set_major_locator(locator)
plt.setp( ax.xaxis.get_majorticklabels(), rotation=90 )
plt.ylim((-1,1))
plt.savefig('SSH_Correlation_R0.PDF')

