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

c_00 = (1.0, 1.0, 1.0)
c_01 = (204/255, 255/255, 255/255)
c_02 = (155/255, 255/255, 254/255)
c_03 = (148/255, 238/255, 247/255)
c_04 = (135/255, 220/255, 234/255)
c_05 = (141/255, 202/255, 227/255)
c_06 = (139/255, 181/255, 219/255)
c_07 = (136/255, 164/255, 210/255)
c_08 = (132/255, 146/255, 200/255)
c_09 = (127/255, 126/255, 192/255)
c_10 = (170/255, 101/255, 144/255)
c_11 = (212/255, 75/255, 100/255)
c_12 = (255/255, 50/255, 50/255)
c_13 = (255/255, 75/255, 40/255)
c_14 = (255/255, 101/255, 26/255)
c_15 = (255/255, 104/255, 23/255)
c_16 = (254/255, 153/255, 0/255)
c_17 = (255/255, 177/255, 0/255)                  
c_18 = (255/255, 203/255, 0/255)
c_19 = (254/255, 229/255, 0/255)
c_20 = (253/255, 255/255, 7/255)
        
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [c_00,c_01,c_02,c_03,c_04,c_05,c_06,c_07,c_08,c_09,c_10
                                                                    ,c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,c_19,c_20])  



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
u_vel_mean_01 = u_vel_sum / len(anomaly_dates)

v_vel_sum = v_vel_anomaly.sum(axis=0)
v_vel_sum = np.flip(v_vel_sum, 0)
v_vel_mean_01 = v_vel_sum / len(anomaly_dates)

current_spd_mean_01 = np.sqrt(u_vel_mean_01**2 + v_vel_mean_01**2)


########################################################################

all_files = glob.glob('../Raw/Anomalies_R1/*.nc')
all_files.sort()

anomaly_dates = pd.read_csv('Remote_Upwelling_Dates.csv')
anomaly_dates = pd.to_datetime(anomaly_dates['Date'])
#nomaly_dates = pd.to_datetime(anomaly_dates['Remote upwelling date'])
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
u_vel_mean_02 = u_vel_sum / len(anomaly_dates)

v_vel_sum = v_vel_anomaly.sum(axis=0)
v_vel_sum = np.flip(v_vel_sum, 0)
v_vel_mean_02 = v_vel_sum / len(anomaly_dates)

current_spd_mean_02 = np.sqrt(u_vel_mean_02**2 + v_vel_mean_02**2)

###########################################################################

Y, X = np.mgrid[np.min(lats):np.max(lats):133j, np.min(lons):np.max(lons):121j]

meridians = np.arange(10.,351.,2.)
parallels = np.arange(-351.,-10.,2.)
meridians = np.arange(10.,351.,2.)
parallels = np.arange(-351.,-10.,2.)

lat_min = -31
lat_max = -22
lon_min = 30
lon_max = 38    

fig, ax = plt.subplots(1,2, subplot_kw={'projection':ccrs.PlateCarree()},figsize=(10,6))
im0 = ax[0].imshow(current_spd_mean_01, vmin = 0, vmax = 0.7, cmap=cmap, extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0], resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax[0].coastlines()
ax[0].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)  
ax[0].plot(32.65, -27.53, 'ro', markersize=5, color = 'lime')
ax[0].set_xlim(lon_min, lon_max)
ax[0].set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
quiv = ax[0].quiver(X[::2,::2], Y[::2,::2], np.flip(u_vel_mean_01, 0)[::2,::2], np.flip(v_vel_mean_01, 0)[::2,::2], zorder=1, scale = 10)
qk = plt.quiverkey(quiv, 0.42, 0.37, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 8}, labelsep = 0.15)
t = qk.text.set_backgroundcolor('w')
ax[0].set_title('At time of remote upwelling')

im1 = ax[1].imshow(current_spd_mean_02, vmin = 0, vmax = 0.7, cmap=cmap, extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
n = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1], resolution = 'l', llcrnrlat=-29, urcrnrlat = -22, llcrnrlon = 32, urcrnrlon = 38)
ax[1].coastlines()
ax[1].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)
ax[1].plot(32.65, -27.53, 'ro', markersize=5, color = 'lime')
ax[1].set_xlim(lon_min, lon_max)
ax[1].set_ylim(lat_min, lat_max)
n.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
n.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
quiv = ax[1].quiver(X[::2,::2], Y[::2,::2], np.flip(u_vel_mean_02, 0)[::2,::2], np.flip(v_vel_mean_02, 0)[::2,::2], zorder=1, scale = 10)
qt = plt.quiverkey(quiv, 0.84, 0.37, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 8}, labelsep = 0.15)
t = qt.text.set_backgroundcolor('w')
ax[1].set_title('At time of anomaly peak')

cbar1 = fig.colorbar(im1, ax=ax[:], orientation='horizontal', pad = 0.095)
cbar1.ax.tick_params(labelsize=7)
cbar1.ax.set_xlabel('Average current speed [m/s]', fontsize = 9)

plt.savefig('Anomaly_Currents_CondAvg_R0.png')   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

        
