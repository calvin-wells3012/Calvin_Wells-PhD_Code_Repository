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

ssh_anomaly__15 = np.empty([len(anomaly_dates),65,69])
u_vel_anomaly__15 = np.empty([len(anomaly_dates),65,69])
v_vel_anomaly__15 = np.empty([len(anomaly_dates),65,69])

ssh_anomaly__10 = np.empty([len(anomaly_dates),65,69])
u_vel_anomaly__10 = np.empty([len(anomaly_dates),65,69])
v_vel_anomaly__10 = np.empty([len(anomaly_dates),65,69])

ssh_anomaly__5 = np.empty([len(anomaly_dates),65,69])
u_vel_anomaly__5 = np.empty([len(anomaly_dates),65,69])
v_vel_anomaly__5 = np.empty([len(anomaly_dates),65,69])

ssh_anomaly = np.empty([len(anomaly_dates),65,69])
u_vel_anomaly = np.empty([len(anomaly_dates),65,69])
v_vel_anomaly = np.empty([len(anomaly_dates),65,69])

ssh_anomaly_10 = np.empty([len(anomaly_dates),65,69])
u_vel_anomaly_10 = np.empty([len(anomaly_dates),65,69])
v_vel_anomaly_10 = np.empty([len(anomaly_dates),65,69])

ssh_anomaly_5 = np.empty([len(anomaly_dates),65,69])
u_vel_anomaly_5 = np.empty([len(anomaly_dates),65,69])
v_vel_anomaly_5 = np.empty([len(anomaly_dates),65,69])

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
        ssh__15 = np.array(fh.variables['sla'][index_list[n]-15,:,:])
        u_vel__15 = np.array(fh.variables['ugosa'][index_list[n]-15,:,:])
        v_vel__15 = np.array(fh.variables['vgosa'][index_list[n]-15,:,:])
        ssh__15[ssh__15 < -100] = np.nan
        ssh_anomaly__15[n] = ssh__15
        u_vel__15[u_vel__15 < -100] = np.nan
        u_vel_anomaly__15[n] = u_vel__15
        v_vel__15[v_vel__15 < -100] = np.nan
        v_vel_anomaly__15[n] = v_vel__15
        
        ssh__10 = np.array(fh.variables['sla'][index_list[n]-10,:,:])
        u_vel__10 = np.array(fh.variables['ugosa'][index_list[n]-10,:,:])
        v_vel__10 = np.array(fh.variables['vgosa'][index_list[n]-10,:,:])
        ssh__10[ssh__10 < -100] = np.nan
        ssh_anomaly__10[n] = ssh__10
        u_vel__10[u_vel__10 < -100] = np.nan
        u_vel_anomaly__10[n] = u_vel__10
        v_vel__10[v_vel__10 < -100] = np.nan
        v_vel_anomaly__10[n] = v_vel__10
        
        ssh__5 = np.array(fh.variables['sla'][index_list[n]-5,:,:])
        u_vel__5 = np.array(fh.variables['ugosa'][index_list[n]-5,:,:])
        v_vel__5 = np.array(fh.variables['vgosa'][index_list[n]-5,:,:])
        ssh__5[ssh__5 < -100] = np.nan
        ssh_anomaly__5[n] = ssh__5
        u_vel__5[u_vel__5 < -100] = np.nan
        u_vel_anomaly__5[n] = u_vel__5
        v_vel__5[v_vel__5 < -100] = np.nan
        v_vel_anomaly__5[n] = v_vel__5
        
        ssh = np.array(fh.variables['sla'][index_list[n],:,:])
        u_vel = np.array(fh.variables['ugosa'][index_list[n],:,:])
        v_vel = np.array(fh.variables['vgosa'][index_list[n],:,:])
        ssh[ssh < -100] = np.nan
        ssh_anomaly[n] = ssh
        u_vel[u_vel < -100] = np.nan
        u_vel_anomaly[n] = u_vel
        v_vel[v_vel < -100] = np.nan
        v_vel_anomaly[n] = v_vel
        
        ssh_10 = np.array(fh.variables['sla'][index_list[n]+10,:,:])
        u_vel_10 = np.array(fh.variables['ugosa'][index_list[n]+10,:,:])
        v_vel_10 = np.array(fh.variables['vgosa'][index_list[n]+10,:,:])
        ssh_10[ssh_10 < -100] = np.nan
        ssh_anomaly_10[n] = ssh_10
        u_vel_10[u_vel_10 < -100] = np.nan
        u_vel_anomaly_10[n] = u_vel_10
        v_vel_10[v_vel_10 < -100] = np.nan
        v_vel_anomaly_10[n] = v_vel_10
        
        ssh_5 = np.array(fh.variables['sla'][index_list[n]+5,:,:])
        u_vel_5 = np.array(fh.variables['ugosa'][index_list[n]+5,:,:])
        v_vel_5 = np.array(fh.variables['vgosa'][index_list[n]+5,:,:])
        ssh_5[ssh_5 < -100] = np.nan
        ssh_anomaly_5[n] = ssh_5
        u_vel_5[u_vel_5 < -100] = np.nan
        u_vel_anomaly_5[n] = u_vel_5
        v_vel_5[v_vel__5 < -100] = np.nan
        v_vel_anomaly_5[n] = v_vel_5
        
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

ssh_sum__5 = ssh_anomaly__5.sum(axis=0)
ssh_sum__5 = np.flip(ssh_sum__5, 0)
ssh_mean__5 = ssh_sum__5 / len(anomaly_dates)
u_vel_sum__5 = u_vel_anomaly__5.sum(axis=0)
u_vel_sum__5 = np.flip(u_vel_sum__5, 0)
u_vel_mean__5 = u_vel_sum__5 / len(anomaly_dates)
v_vel_sum__5 = v_vel_anomaly__5.sum(axis=0)
v_vel_sum__5 = np.flip(v_vel_sum__5, 0)
v_vel_mean__5 = v_vel_sum__5 / len(anomaly_dates)

ssh_sum__10 = ssh_anomaly__10.sum(axis=0)
ssh_sum__10 = np.flip(ssh_sum__10, 0)
ssh_mean__10 = ssh_sum__10 / len(anomaly_dates)
u_vel_sum__10 = u_vel_anomaly__10.sum(axis=0)
u_vel_sum__10 = np.flip(u_vel_sum__10, 0)
u_vel_mean__10 = u_vel_sum__10 / len(anomaly_dates)
v_vel_sum__10 = v_vel_anomaly__10.sum(axis=0)
v_vel_sum__10 = np.flip(v_vel_sum__10, 0)
v_vel_mean__10 = v_vel_sum__10 / len(anomaly_dates)

ssh_sum__15 = ssh_anomaly__15.sum(axis=0)
ssh_sum__15 = np.flip(ssh_sum__15, 0)
ssh_mean__15 = ssh_sum__15 / len(anomaly_dates)
u_vel_sum__15 = u_vel_anomaly__15.sum(axis=0)
u_vel_sum__15 = np.flip(u_vel_sum__15, 0)
u_vel_mean__15 = u_vel_sum__15 / len(anomaly_dates)
v_vel_sum__15 = v_vel_anomaly__15.sum(axis=0)
v_vel_sum__15 = np.flip(v_vel_sum__15, 0)
v_vel_mean__15 = v_vel_sum__15 / len(anomaly_dates)

ssh_sum_5 = ssh_anomaly_5.sum(axis=0)
ssh_sum_5 = np.flip(ssh_sum_5, 0)
ssh_mean_5 = ssh_sum_5 / len(anomaly_dates)
u_vel_sum_5 = u_vel_anomaly_5.sum(axis=0)
u_vel_sum_5 = np.flip(u_vel_sum_5, 0)
u_vel_mean_5 = u_vel_sum_5 / len(anomaly_dates)
v_vel_sum_5 = v_vel_anomaly_5.sum(axis=0)
v_vel_sum_5 = np.flip(v_vel_sum_5, 0)
v_vel_mean_5 = v_vel_sum_5 / len(anomaly_dates)

ssh_sum_10 = ssh_anomaly_10.sum(axis=0)
ssh_sum_10 = np.flip(ssh_sum_10, 0)
ssh_mean_10 = ssh_sum_10 / len(anomaly_dates)
u_vel_sum_10 = u_vel_anomaly_10.sum(axis=0)
u_vel_sum_10 = np.flip(u_vel_sum_10, 0)
u_vel_mean_10 = u_vel_sum_10 / len(anomaly_dates)
v_vel_sum_10 = v_vel_anomaly_10.sum(axis=0)
v_vel_sum_10 = np.flip(v_vel_sum_10, 0)
v_vel_mean_10 = v_vel_sum_10 / len(anomaly_dates)


Y, X = np.mgrid[np.min(lats):np.max(lats):65j, np.min(lons):np.max(lons):69j]

fig, ax = plt.subplots(2,3,figsize=(9,7))

im0 = ax[0,0].imshow(ssh_mean__15, vmin = -0.2, vmax = 0.2, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0,0], resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
m.drawcoastlines()
m.fillcontinents(color='grey')
border_color = 'black'   
ax[0,0].plot(32.65, -27.53, 'ro', markersize=2, color = 'deeppink')
ax[0,0].set_xlabel("Longitude [deg]", fontsize = 6)
ax[0,0].set_ylabel("Latitude [deg]", fontsize = 6)
ax[0,0].set_xticks(range(30, 45, 1))
ax[0,0].set_yticks(range(-34, -20, 1))
ax[0,0].set_xticklabels(range(30, 45, 1), fontsize=5)
ax[0,0].set_yticklabels(range(-34, -20, 1), fontsize=5)
quiv0 = ax[0,0].quiver(X, Y, np.flip(u_vel_mean__15, 0), np.flip(v_vel_mean__15, 0), zorder=1)
qq = plt.quiverkey(quiv0, 0.322, 0.335, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 5}, labelsep = 0.1)
t = qq.text.set_backgroundcolor('w')

im1 = ax[0,1].imshow(ssh_mean__10, vmin = -0.2, vmax = 0.2, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
n = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0,1], resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
n.drawcoastlines()
n.fillcontinents(color='grey')
border_color = 'black'   
ax[0,1].plot(32.65, -27.53, 'ro', markersize=2, color = 'deeppink')
ax[0,1].set_xlabel("Longitude [deg]", fontsize = 6)
ax[0,1].set_ylabel("Latitude [deg]", fontsize = 6)
ax[0,1].set_xticks(range(30, 45, 1))
ax[0,1].set_yticks(range(-34, -20, 1))
ax[0,1].set_xticklabels(range(30, 45, 1), fontsize=5)
ax[0,1].set_yticklabels(range(-34, -20, 1), fontsize=5)
quiv1 = ax[0,1].quiver(X, Y, np.flip(u_vel_mean__10, 0), np.flip(v_vel_mean__10, 0), zorder=1)
qk = plt.quiverkey(quiv1, 0.595, 0.335, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 5}, labelsep = 0.1)
t = qk.text.set_backgroundcolor('w')

im2 = ax[1,0].imshow(ssh_mean, vmin = -0.2, vmax = 0.2, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
n = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1,0], resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
n.drawcoastlines()
n.fillcontinents(color='grey')
border_color = 'black'   
ax[1,0].plot(32.65, -27.53, 'ro', markersize=2, color = 'deeppink')
ax[1,0].set_xlabel("Longitude [deg]", fontsize = 6)
ax[1,0].set_ylabel("Latitude [deg]", fontsize = 6)
ax[1,0].set_xticks(range(30, 45, 1))
ax[1,0].set_yticks(range(-34, -20, 1))
ax[1,0].set_xticklabels(range(30, 45, 1), fontsize=5)
ax[1,0].set_yticklabels(range(-34, -20, 1), fontsize=5)
quiv2 = ax[1,0].quiver(X, Y, np.flip(u_vel_mean, 0), np.flip(v_vel_mean, 0), zorder=1)
ql = plt.quiverkey(quiv2, 0.87, 0.335, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 5}, labelsep = 0.1)
t = ql.text.set_backgroundcolor('w')

im3 = ax[1,1].imshow(ssh_mean_5, vmin = -0.2, vmax = 0.2, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
n = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1,1], resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
n.drawcoastlines()
n.fillcontinents(color='grey')
border_color = 'black'   
ax[1,1].plot(32.65, -27.53, 'ro', markersize=2, color = 'deeppink')
ax[1,1].set_xlabel("Longitude [deg]", fontsize = 6)
ax[1,1].set_ylabel("Latitude [deg]", fontsize = 6)
ax[1,1].set_xticks(range(30, 45, 1))
ax[1,1].set_yticks(range(-34, -20, 1))
ax[1,1].set_xticklabels(range(30, 45, 1), fontsize=5)
ax[1,1].set_yticklabels(range(-34, -20, 1), fontsize=5)
quiv3 = ax[1,1].quiver(X, Y, np.flip(u_vel_mean_5, 0), np.flip(v_vel_mean_5, 0), zorder=1)
qt = plt.quiverkey(quiv3, 0.322, 0.66, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 5}, labelsep = 0.1)
t = qt.text.set_backgroundcolor('w')

im4 = ax[0,2].imshow(ssh_mean__5, vmin = -0.2, vmax = 0.2, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
n = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0,2], resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
n.drawcoastlines()
n.fillcontinents(color='grey')
border_color = 'black'   
ax[0,2].plot(32.65, -27.53, 'ro', markersize=2, color = 'deeppink')
ax[0,2].set_xlabel("Longitude [deg]", fontsize = 6)
ax[0,2].set_ylabel("Latitude [deg]", fontsize = 6)
ax[0,2].set_xticks(range(30, 45, 1))
ax[0,2].set_yticks(range(-34, -20, 1))
ax[0,2].set_xticklabels(range(30, 45, 1), fontsize=5)
ax[0,2].set_yticklabels(range(-34, -20, 1), fontsize=5)
quiv4 = ax[0,2].quiver(X, Y, np.flip(u_vel_mean__5, 0), np.flip(v_vel_mean__5, 0), zorder=1)
qk = plt.quiverkey(quiv4, 0.595, 0.66, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 5}, labelsep = 0.1)
t = qk.text.set_backgroundcolor('w')

im5 = ax[1,2].imshow(ssh_mean_10, vmin = -0.2, vmax = 0.2, cmap='bwr', extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
n = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1,2], resolution = 'l', llcrnrlat=-34, urcrnrlat = -20, llcrnrlon = 30, urcrnrlon = 45)
n.drawcoastlines()
n.fillcontinents(color='grey')
border_color = 'black'   
ax[1,2].plot(32.65, -27.53, 'ro', markersize=2, color = 'deeppink')
ax[1,2].set_xlabel("Longitude [deg]", fontsize = 6)
ax[1,2].set_ylabel("Latitude [deg]", fontsize = 6)
ax[1,2].set_xticks(range(30, 45, 1))
ax[1,2].set_yticks(range(-34, -20, 1))
ax[1,2].set_xticklabels(range(30, 45, 1), fontsize=5)
ax[1,2].set_yticklabels(range(-34, -20, 1), fontsize=5)
quiv5 = ax[1,2].quiver(X, Y, np.flip(u_vel_mean_10, 0), np.flip(v_vel_mean_10, 0), zorder=1)
qk = plt.quiverkey(quiv5, 0.87, 0.66, 0.5, r'$0.5 m/s$', labelpos ='S',coordinates='figure',fontproperties={'weight': 'bold','size': 5}, labelsep = 0.1)
t = qk.text.set_backgroundcolor('w')

cbar0 = fig.colorbar(im0, ax=ax[:], orientation='horizontal', pad = 0.08)
cbar0.ax.tick_params(labelsize=6)
cbar0.ax.set_xlabel('Average SSH [m]', fontsize = 9)

#plt.tight_layout()
#plt.savefig('SSH_Currents_CondAvg_2deg.PDF')   
plt.savefig('SSH_CondAvg_Alt_2deg_time_lag.png')  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

        
