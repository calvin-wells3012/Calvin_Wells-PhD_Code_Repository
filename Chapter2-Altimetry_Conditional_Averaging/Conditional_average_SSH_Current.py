## This program plots the depth averaged (first 30m from surface) current speed spacially and creates a video

import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits.axes_grid1 import make_axes_locatable

## load data file names    

all_files = glob.glob('../../Raw/*.nc')
all_files.sort()

anomaly_dates = pd.read_csv('../../../../Measured_Temperature/Processed/FFT/*.csv')
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
    
    n = 0
    while n < len(anomaly_dates):
        ssh = np.array(fh.variables['sla'][index_list[n]+0,:,:])
        u_vel = np.array(fh.variables['ugosa'][index_list[n]+0,:,:])
        v_vel = np.array(fh.variables['vgosa'][index_list[n]+0,:,:])
        ssh[ssh < -100] = np.nan
        ssh_anomaly[n] = ssh
        u_vel[u_vel < -100] = np.nan
        u_vel_anomaly[n] = u_vel
        v_vel[v_vel < -100] = np.nan
        v_vel_anomaly[n] = v_vel
        n+=1

ssh_sum = ssh_anomaly.sum(axis=0)
ssh_sum = np.flip(ssh_sum, 0)
ssh_mean = ssh_sum / len(anomaly_dates)

u_vel_sum = u_vel_anomaly.sum(axis=0)
u_vel_sum = np.flip(u_vel_sum, 0)
u_vel_mean = u_vel_sum / len(anomaly_dates)

v_vel_sum = v_vel_anomaly.sum(axis=0)
v_vel_sum = np.flip(v_vel_sum, 0)
v_vel_mean = v_vel_sum / len(anomaly_dates)

current_spd_mean = np.sqrt(u_vel_mean**2 + v_vel_mean**2)
    
#fig01 = plt.figure(figsize=(11,7))
#ax = fig01.add_subplot(1,1,1)
#m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax, resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
#border_color = 'black'    
#m.drawcoastlines()    
#plt.imshow(ssh_mean, vmin = -0.2, vmax = 0.2, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
#m.fillcontinents(color='grey')
#cbar1 = plt.colorbar()
#cbar1.ax.set_ylabel(r"Average SSH [m]",  labelpad=20, rotation=270) 
#plt.plot(32.65, -27.53, 'ro', markersize=5, color = 'deeppink')
#plt.xlabel('Longitude [deg]', fontsize=14)
#plt.ylabel('Latitude [deg]', fontsize=14)
#plt.xticks(np.arange(30, 45, step=1))
#plt.yticks(np.arange(-34, -20, step=1))
#ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
Y, X = np.mgrid[np.min(lats):np.max(lats):65j, np.min(lons):np.max(lons):69j]

fig, ax = plt.subplots(1,2,figsize=(10,6))
im0 = ax[0].imshow(ssh_mean, vmin = -0.2, vmax = 0.2, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0], resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
m.drawcoastlines()
m.fillcontinents(color='grey')
border_color = 'black'   
ax[0].plot(32.65, -27.53, 'ro', markersize=3, color = 'deeppink')
ax[0].set_xlabel("Longitude [deg]", fontsize = 7)
ax[0].set_ylabel("Latitude [deg]", fontsize = 7)
ax[0].set_xticks(range(30, 45, 1))
ax[0].set_yticks(range(-34, -20, 1))
ax[0].set_xticklabels(range(30, 45, 1), fontsize=7)
ax[0].set_yticklabels(range(-34, -20, 1), fontsize=7)
quiv = ax[0].quiver(X, Y, np.flip(u_vel_mean, 0), np.flip(v_vel_mean, 0), zorder=1)
cbar0 = fig.colorbar(im0, ax=ax[0], orientation='horizontal', pad = 0.095)
cbar0.ax.tick_params(labelsize=7)
cbar0.ax.set_xlabel('Average SSH [m]', fontsize = 9)

qk = plt.quiverkey(quiv, 0.44, 0.37, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 8}, labelsep = 0.15)
t = qk.text.set_backgroundcolor('w')

im1 = ax[1].imshow(current_spd_mean, vmin = 0, vmax = 0.25, cmap='jet', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
n = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1], resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
n.drawcoastlines()
n.fillcontinents(color='grey')
ax[1].plot(32.65, -27.53, 'ro', markersize=3, color = 'deeppink')
ax[1].set_xlabel("Longitude [deg]", fontsize = 7)
ax[1].set_ylabel("Latitude [deg]", fontsize = 7)
ax[1].set_xticks(range(30, 45, 1))
ax[1].set_yticks(range(-34, -20, 1))
ax[1].set_xticklabels(range(30, 45, 1), fontsize=7)
ax[1].set_yticklabels(range(-34, -20, 1), fontsize=7)
cbar1 = fig.colorbar(im1, ax=ax[1], orientation='horizontal', pad = 0.095)
cbar1.ax.tick_params(labelsize=7)
cbar1.ax.set_xlabel('Average geostrophic current speed [m/s]', fontsize = 9)
quiv = ax[1].quiver(X, Y, np.flip(u_vel_mean, 0), np.flip(v_vel_mean, 0), zorder=1)

qt = plt.quiverkey(quiv, 0.865, 0.37, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 8}, labelsep = 0.15)
t = qt.text.set_backgroundcolor('w')

#plt.tight_layout()
#plt.savefig('SSH_Currents_CondAvg_2deg.PDF')   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

        
