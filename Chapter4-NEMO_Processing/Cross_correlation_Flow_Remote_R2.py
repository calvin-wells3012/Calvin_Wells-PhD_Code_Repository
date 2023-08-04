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
from scipy.signal import correlate2d

## load data file names    

all_files = glob.glob('../Raw/Anomalies_R1/*.nc')
all_files.sort()

anomaly_dates = pd.read_csv('Remote_Upwelling_Dates.csv')
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
    #lats = np.flip(lats)

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
#temp_sum = np.flip(temp_sum, 0)
temp_mean = temp_sum / len(anomaly_dates)

u_vel_sum = u_vel_anomaly.sum(axis=0)
#_vel_sum = np.flip(u_vel_sum, 0)
u_vel_mean = u_vel_sum / len(anomaly_dates)

v_vel_sum = v_vel_anomaly.sum(axis=0)
#v_vel_sum = np.flip(v_vel_sum, 0)
v_vel_mean = v_vel_sum / len(anomaly_dates)

    
#%%

# extract the x- and y-components of the vectors
# ai = np.array(u_vel_mean[88:114,67:93])
# bj = np.array(v_vel_mean[88:114,67:93])

# ci = np.array(u_vel_anomaly[0,88:114,67:93])
# dj = np.array(v_vel_anomaly[0,88:114,67:93])

cc = np.zeros([len(anomaly_dates),26,25])
cc_max = np.zeros(len(anomaly_dates))
cc_mean = np.zeros(len(anomaly_dates))

i = 0
while i < len(anomaly_dates):
    
    ai = np.array(u_vel_mean[70:96,56:81])
    bj = np.array(v_vel_mean[70:96,56:81])

    ci = np.array(u_vel_anomaly[i,70:96,56:81])
    dj = np.array(v_vel_anomaly[i,70:96,56:81])

    cc_01 = ((ai * ci) + (bj * dj)) / (np.sqrt(np.square(ai)+np.square(bj)) * np.sqrt(np.square(ci)+np.square(dj)))
    cc[i] = cc_01

    cc_max_01 = np.nanmax(cc_01)
    cc_mean_01 = np.nanmean(cc_01)
    
    cc_max[i] = cc_max_01
    cc_mean[i] = cc_mean_01
    
    i+=1
    
Y, X = np.mgrid[np.min(lats):np.max(lats):133j, np.min(lons):np.max(lons):121j]

meridians = np.arange(10.,351.,2.)
parallels = np.arange(-351.,-10.,2.)
meridians = np.arange(10.,351.,2.)
parallels = np.arange(-351.,-10.,2.)

lat_min = -31
lat_max = -22
lon_min = 30
lon_max = 38    

fig, ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree()},figsize=(10,6))
im0 = ax.imshow(np.flip(cc[16],0), vmin = -1, vmax = 1, cmap='bwr', extent=[lons[44], lons[93], lats[58], lats[107]], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax, resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax.coastlines()
ax.add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)  
ax.plot(32.65, -27.53, 'ro', markersize=5, color = 'lime')
ax.set_xlim(lon_min, lon_max)
ax.set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
ax.set_title('At time of remote upwelling')    
    
    

        
