## This program takes co-ordinates in LAT LON from an input file and generates a CSV that can be read in by Delft3D as a temperature boundary input file.

import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import xarray as xr
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy.interpolate import griddata
import csv
from numpy import savetxt

#text_file = open('../Raw/*.txt', 'r')
#lines = text_file.readlines()

## load data file names    
all_files = glob.glob('../Raw/*.nc')
all_files.sort()
datafile = all_files[0]

julian_dates = pd.read_csv('data.csv')
julian_dates = np.array(julian_dates)

##Loops through each netcdf file

fh = xr.open_dataset(datafile)
track_no = np.array(fh.variables['track'][:])
date_test = pd.DataFrame(fh.variables['j1'][:100])
#obs_no = np.array(fh.variables['n'][:100])
#amp = fh.variables['A'][:]
cyc = np.array(fh.variables['cyc'][:])
lat = np.array(fh.variables['lat'][:])
lon = np.array(fh.variables['lon'][:])-360

##Filter latitudes
lat_idx = np.where(np.logical_and(lat>-30, lat<-26))
lat_idx_np = np.array(lat_idx[0])
lat_filter01 = np.empty(len(lat_idx_np))
lon_filter01 = np.empty(len(lat_idx_np))
cyc_filter01 = np.empty(len(lat_idx_np))
track_filter01 = np.empty(len(lat_idx_np))
date_filter01 = np.empty(len(lat_idx_np))
i=0
while i < len(lat_idx_np):
    x = int(lat_idx_np[i])
    lat_filter01[i] = lat[x]
    lon_filter01[i] = lon[x]
    cyc_filter01[i] = cyc[x]
    track_filter01[i] = track_no[x]
    date_filter01[i] = julian_dates[x]
    i+=1

##Filter longitudes
lon_idx = np.where(np.logical_and(lon_filter01>32, lon_filter01<36))
lon_idx_np = np.array(lon_idx[0])
lat_filter02 = np.empty(len(lon_idx_np))
lon_filter02 = np.empty(len(lon_idx_np))
cyc_filter02 = np.empty(len(lon_idx_np))
track_filter02 = np.empty(len(lon_idx_np))
date_filter02 = np.empty(len(lon_idx_np))
i=0
while i < len(lon_idx_np):
    x = int(lon_idx_np[i])
    lat_filter02[i] = lat_filter01[x]
    lon_filter02[i] = lon_filter01[x]
    cyc_filter02[i] = cyc_filter01[x]
    track_filter02[i] = track_filter01[x]
    date_filter02[i] = date_filter01[x]
    i+=1

##Filter cyclonic eddies
cyc_idx = np.where(cyc_filter02==-1)
cyc_idx_np = np.array(cyc_idx[0])
lat_filter03 = np.empty(len(cyc_idx_np))
lon_filter03 = np.empty(len(cyc_idx_np))
cyc_filter03 = np.empty(len(cyc_idx_np))
track_filter03 = np.empty(len(cyc_idx_np))
date_filter03 = np.empty(len(cyc_idx_np))
i=0
while i < len(cyc_idx_np):
    x = int(cyc_idx_np[i])
    lat_filter03[i] = lat_filter02[x]
    lon_filter03[i] = lon_filter02[x]
    cyc_filter03[i] = cyc_filter02[x]
    track_filter03[i] = track_filter02[x]
    date_filter03[i] = date_filter02[x]
    i+=1

date_days = date_filter03 - 2448623
startdate = datetime.datetime(1992,1,1,0,0,0)
date_real = [None] * len(date_filter03)
i=0
while i < len(date_filter03):
    time_day = date_days[i]
    day = startdate + datetime.timedelta(days=time_day)
    date_real[i] = day
    i+=1

date_int = np.empty(len(date_real))
j=0
while j< len(date_real):
    time_step = date_real[j]
    time_step = int(time_step.strftime('%Y%m%d'))
    date_int[j] = time_step
    j+=1
    
anomaly_dates = pd.read_csv('Anomaly_Dates_2deg_Julian.csv')
anomaly_dates = pd.to_datetime(anomaly_dates['Date'])
startdate = datetime.datetime(1950,1,1,0,0,0)

index_list = np.empty(len(anomaly_dates))
k=0
while k < len(anomaly_dates):
    anomaly_step = anomaly_dates[k]
    anomaly_step = int(anomaly_step.strftime('%Y%m%d'))
    index = np.abs(date_int - anomaly_step).argmin()
    index_list[k] = index
    k+=1
    
anomaly_date_eddy = np.empty(len(anomaly_dates))
anomaly_track = np.empty(len(anomaly_dates))
m = 0
while m < len(anomaly_dates):
    z = int(index_list[m])
    track = track_filter03[z]
    anomaly_track[m] = track
    m+=1

