import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
import xarray as xr
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

all_files = glob.glob('../../../Raw/Anomalies/*.nc')
all_files.sort()

mean_traj = np.array(pd.read_csv('../../Clustering_Analysis_20m_R1/Clusters_07/Mean_Trajectories_06Clusters.csv', header = None))
mean_lons_00 = mean_traj[0]
mean_lons_01 = mean_traj[1]
mean_lons_02 = mean_traj[2]
mean_lons_03 = mean_traj[3]
mean_lons_04 = mean_traj[4]
mean_lons_05 = mean_traj[5]
mean_lats_00 = mean_traj[6]
mean_lats_01 = mean_traj[7]
mean_lats_02 = mean_traj[8]
mean_lats_03 = mean_traj[9]
mean_lats_04 = mean_traj[10]
mean_lats_05 = mean_traj[11]

######################################################

anomaly_data = pd.read_csv('Anomaly_Dates_2deg_Class00_00.csv')
anomaly_num = np.array(anomaly_data['Anomaly'])
u_vel = np.zeros([len(anomaly_data),133,121])
v_vel = np.zeros([len(anomaly_data),133,121])
x = 0
while x < 1:    
    lon_traj_pos = mean_lons_00[x]
    lat_traj_pos = mean_lats_00[x]    
    i = 0
    while i < 6:    
        fh = xr.open_dataset(all_files[anomaly_num[i]])    
        time =  fh.variables['time'][:]
        date = str(time[40])
        date = date[28:38]
        depth = np.array(fh.variables['depth'][:])
        depth_layer = 12
        lons = np.array(fh.variables['longitude'][:])
        lats = np.array(fh.variables['latitude'][:])
        lats = np.flip(lats)
        u_vel_anomaly = np.array(fh.variables['uo'][40-x,depth_layer,:,:])
        u_vel_anomaly = np.flip(u_vel_anomaly,0)
        u_vel[i] = u_vel_anomaly
        v_vel_anomaly = np.array(fh.variables['vo'][40-x,depth_layer,:,:])
        v_vel_anomaly = np.flip(v_vel_anomaly,0)
        v_vel[i] = v_vel_anomaly    
        i+=1

    u_vel_sum = u_vel.sum(axis=0)
    u_vel_mean_00 = u_vel_sum / len(anomaly_data)
    v_vel_sum = v_vel.sum(axis=0)
    v_vel_mean_00 = v_vel_sum / len(anomaly_data)
    current_spd_mean_00 = np.sqrt(u_vel_mean_00**2 + v_vel_mean_00**2)
    Y, X = np.mgrid[np.min(lats):np.max(lats):133j, np.min(lons):np.max(lons):121j]
    x+=1
    
######################################################

anomaly_data = pd.read_csv('Anomaly_Dates_2deg_Class02_01.csv')
anomaly_num = np.array(anomaly_data['Anomaly'])
u_vel = np.zeros([len(anomaly_data),133,121])
v_vel = np.zeros([len(anomaly_data),133,121])
x = 0
while x < 1:    
    lon_traj_pos = mean_lons_01[x]
    lat_traj_pos = mean_lats_01[x]    
    i = 0
    while i < 10:    
        fh = xr.open_dataset(all_files[anomaly_num[i]])    
        time =  fh.variables['time'][:]
        date = str(time[40])
        date = date[28:38]
        depth = np.array(fh.variables['depth'][:])
        depth_layer = 12
        lons = np.array(fh.variables['longitude'][:])
        lats = np.array(fh.variables['latitude'][:])
        lats = np.flip(lats)
        u_vel_anomaly = np.array(fh.variables['uo'][40-x,depth_layer,:,:])
        u_vel_anomaly = np.flip(u_vel_anomaly,0)
        u_vel[i] = u_vel_anomaly
        v_vel_anomaly = np.array(fh.variables['vo'][40-x,depth_layer,:,:])
        v_vel_anomaly = np.flip(v_vel_anomaly,0)
        v_vel[i] = v_vel_anomaly    
        i+=1

    u_vel_sum = u_vel.sum(axis=0)
    u_vel_mean_01 = u_vel_sum / len(anomaly_data)
    v_vel_sum = v_vel.sum(axis=0)
    v_vel_mean_01 = v_vel_sum / len(anomaly_data)
    current_spd_mean_01 = np.sqrt(u_vel_mean_01**2 + v_vel_mean_01**2)
    Y, X = np.mgrid[np.min(lats):np.max(lats):133j, np.min(lons):np.max(lons):121j]
    x+=1    
    
######################################################

anomaly_data = pd.read_csv('Anomaly_Dates_2deg_Class03_02.csv')
anomaly_num = np.array(anomaly_data['Anomaly'])
u_vel = np.zeros([len(anomaly_data),133,121])
v_vel = np.zeros([len(anomaly_data),133,121])
x = 0
while x < 1:    
    lon_traj_pos = mean_lons_02[x]
    lat_traj_pos = mean_lats_02[x]    
    i = 0
    while i < 9:    
        fh = xr.open_dataset(all_files[anomaly_num[i]])    
        time =  fh.variables['time'][:]
        date = str(time[40])
        date = date[28:38]
        depth = np.array(fh.variables['depth'][:])
        depth_layer = 12
        lons = np.array(fh.variables['longitude'][:])
        lats = np.array(fh.variables['latitude'][:])
        lats = np.flip(lats)
        u_vel_anomaly = np.array(fh.variables['uo'][40-x,depth_layer,:,:])
        u_vel_anomaly = np.flip(u_vel_anomaly,0)
        u_vel[i] = u_vel_anomaly
        v_vel_anomaly = np.array(fh.variables['vo'][40-x,depth_layer,:,:])
        v_vel_anomaly = np.flip(v_vel_anomaly,0)
        v_vel[i] = v_vel_anomaly    
        i+=1

    u_vel_sum = u_vel.sum(axis=0)
    u_vel_mean_02 = u_vel_sum / len(anomaly_data)
    v_vel_sum = v_vel.sum(axis=0)
    v_vel_mean_02 = v_vel_sum / len(anomaly_data)
    current_spd_mean_02 = np.sqrt(u_vel_mean_02**2 + v_vel_mean_02**2)
    Y, X = np.mgrid[np.min(lats):np.max(lats):133j, np.min(lons):np.max(lons):121j]
    x+=1    
    
######################################################

anomaly_data = pd.read_csv('Anomaly_Dates_2deg_Class04_03.csv')
anomaly_num = np.array(anomaly_data['Anomaly'])
u_vel = np.zeros([len(anomaly_data),133,121])
v_vel = np.zeros([len(anomaly_data),133,121])
x = 0
while x < 1:    
    lon_traj_pos = mean_lons_03[x]
    lat_traj_pos = mean_lats_03[x]    
    i = 0
    while i < 26:    
        fh = xr.open_dataset(all_files[anomaly_num[i]])    
        time =  fh.variables['time'][:]
        date = str(time[40])
        date = date[28:38]
        depth = np.array(fh.variables['depth'][:])
        depth_layer = 12
        lons = np.array(fh.variables['longitude'][:])
        lats = np.array(fh.variables['latitude'][:])
        lats = np.flip(lats)
        u_vel_anomaly = np.array(fh.variables['uo'][40-x,depth_layer,:,:])
        u_vel_anomaly = np.flip(u_vel_anomaly,0)
        u_vel[i] = u_vel_anomaly
        v_vel_anomaly = np.array(fh.variables['vo'][40-x,depth_layer,:,:])
        v_vel_anomaly = np.flip(v_vel_anomaly,0)
        v_vel[i] = v_vel_anomaly    
        i+=1

    u_vel_sum = u_vel.sum(axis=0)
    u_vel_mean_03 = u_vel_sum / len(anomaly_data)
    v_vel_sum = v_vel.sum(axis=0)
    v_vel_mean_03 = v_vel_sum / len(anomaly_data)
    current_spd_mean_03 = np.sqrt(u_vel_mean_03**2 + v_vel_mean_03**2)
    Y, X = np.mgrid[np.min(lats):np.max(lats):133j, np.min(lons):np.max(lons):121j]
    x+=1    
    
######################################################

anomaly_data = pd.read_csv('Anomaly_Dates_2deg_Class01_04.csv')
anomaly_num = np.array(anomaly_data['Anomaly'])
u_vel = np.zeros([len(anomaly_data),133,121])
v_vel = np.zeros([len(anomaly_data),133,121])
x = 0
while x < 1:    
    lon_traj_pos = mean_lons_04[x]
    lat_traj_pos = mean_lats_04[x]    
    i = 0
    while i < 8:    
        fh = xr.open_dataset(all_files[anomaly_num[i]])    
        time =  fh.variables['time'][:]
        date = str(time[40])
        date = date[28:38]
        depth = np.array(fh.variables['depth'][:])
        depth_layer = 12
        lons = np.array(fh.variables['longitude'][:])
        lats = np.array(fh.variables['latitude'][:])
        lats = np.flip(lats)
        u_vel_anomaly = np.array(fh.variables['uo'][40-x,depth_layer,:,:])
        u_vel_anomaly = np.flip(u_vel_anomaly,0)
        u_vel[i] = u_vel_anomaly
        v_vel_anomaly = np.array(fh.variables['vo'][40-x,depth_layer,:,:])
        v_vel_anomaly = np.flip(v_vel_anomaly,0)
        v_vel[i] = v_vel_anomaly    
        i+=1

    u_vel_sum = u_vel.sum(axis=0)
    u_vel_mean_04 = u_vel_sum / len(anomaly_data)
    v_vel_sum = v_vel.sum(axis=0)
    v_vel_mean_04 = v_vel_sum / len(anomaly_data)
    current_spd_mean_04 = np.sqrt(u_vel_mean_04**2 + v_vel_mean_04**2)
    Y, X = np.mgrid[np.min(lats):np.max(lats):133j, np.min(lons):np.max(lons):121j]
    x+=1    
    
######################################################

anomaly_data = pd.read_csv('Anomaly_Dates_2deg_Class05_05.csv')
anomaly_num = np.array(anomaly_data['Anomaly'])
u_vel = np.zeros([len(anomaly_data),133,121])
v_vel = np.zeros([len(anomaly_data),133,121])
x = 0
while x < 1:    
    lon_traj_pos = mean_lons_05[x]
    lat_traj_pos = mean_lats_05[x]    
    i = 0
    while i < 4:    
        fh = xr.open_dataset(all_files[anomaly_num[i]])    
        time =  fh.variables['time'][:]
        date = str(time[40])
        date = date[28:38]
        depth = np.array(fh.variables['depth'][:])
        depth_layer = 12
        lons = np.array(fh.variables['longitude'][:])
        lats = np.array(fh.variables['latitude'][:])
        lats = np.flip(lats)
        u_vel_anomaly = np.array(fh.variables['uo'][40-x,depth_layer,:,:])
        u_vel_anomaly = np.flip(u_vel_anomaly,0)
        u_vel[i] = u_vel_anomaly
        v_vel_anomaly = np.array(fh.variables['vo'][40-x,depth_layer,:,:])
        v_vel_anomaly = np.flip(v_vel_anomaly,0)
        v_vel[i] = v_vel_anomaly    
        i+=1

    u_vel_sum = u_vel.sum(axis=0)
    u_vel_mean_05 = u_vel_sum / len(anomaly_data)
    v_vel_sum = v_vel.sum(axis=0)
    v_vel_mean_05 = v_vel_sum / len(anomaly_data)
    current_spd_mean_05 = np.sqrt(u_vel_mean_05**2 + v_vel_mean_05**2)
    Y, X = np.mgrid[np.min(lats):np.max(lats):133j, np.min(lons):np.max(lons):121j]
    x+=1    
    
#%%    
im_ratio = current_spd_mean_00.shape[0]/current_spd_mean_00.shape[1]

current_spd_mean_00[current_spd_mean_00 > 0] = np.nan
current_spd_mean_00[current_spd_mean_00 < 0] = np.nan

current_spd_mean_04[current_spd_mean_04 > 0] = np.nan
current_spd_mean_04[current_spd_mean_04 < 0] = np.nan

meridians = np.arange(10.,351.,2.)
parallels = np.arange(-351.,-10.,2.)

fig, ax = plt.subplots(1,2,figsize=(12,7))
fig.subplots_adjust(wspace = 0.1, hspace = 0.3)

im0 = ax[0].imshow(current_spd_mean_00, cmap=cmap, vmin = 0, vmax = 1, extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear', zorder = -1);
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0], resolution = 'l', llcrnrlat=-30, urcrnrlat = -23, llcrnrlon = 32, urcrnrlon = 38)
m.drawcoastlines()
m.fillcontinents(color='lightgrey')
border_color = 'black'   
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2)
ax[0].plot(32.7, -27.416, 'ro', markersize=6, color = 'deeppink')
quiv = ax[0].quiver(X[::3, ::3], Y[::3, ::3], np.flip(u_vel_mean_00, 0)[::3, ::3], np.flip(v_vel_mean_00, 0)[::3, ::3], zorder=0, scale = 10)  
ax[0].set_title('Cluster 1', fontsize = 14)

im4 = ax[1].imshow(current_spd_mean_04, cmap=cmap, vmin = 0, vmax = 1, extent=[min(lons), max(lons), min(lats), max(lats)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1], resolution = 'l', llcrnrlat=-30, urcrnrlat = -23, llcrnrlon = 32, urcrnrlon = 38)
m.drawcoastlines()
m.fillcontinents(color='lightgrey')
border_color = 'black'   
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2)
ax[1].plot(32.7, -27.416, 'ro', markersize=6, color = 'deeppink')
quiv = ax[1].quiver(X[::3, ::3], Y[::3, ::3], np.flip(u_vel_mean_04, 0)[::3, ::3], np.flip(v_vel_mean_04, 0)[::3, ::3], zorder=0, scale = 10)  
ax[1].set_title('Cluster 5', fontsize = 14)

qt = plt.quiverkey(quiv, 0.88, 0.885, 0.5, r'$0.5 m/s$', labelpos ='W',coordinates='figure',fontproperties={'weight': 'bold','size': 12}, labelsep = 0.15)
t = qt.text.set_backgroundcolor('w')

qt = plt.quiverkey(quiv, 0.47, 0.887, 0.5, r'$0.5 m/s$', labelpos ='W',coordinates='figure',fontproperties={'weight': 'bold','size': 12}, labelsep = 0.15)
t = qt.text.set_backgroundcolor('w')

plt.savefig('CMEMS_Currents_CondAvg_R2.PDF')
plt.savefig('CMEMS_Currents_CondAvg_R2.png')
#plt.close()








        
