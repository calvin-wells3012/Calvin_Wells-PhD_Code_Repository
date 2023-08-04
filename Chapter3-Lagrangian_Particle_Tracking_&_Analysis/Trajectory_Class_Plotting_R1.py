## 2D Depth plotting of mean trajectory classes and individual trajectories

import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import xarray as xr
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.dates as mdates
import cv2
bm = Basemap()
from scipy import signal
import math

all_files = glob.glob('../Paticle_Tracking_Runs_20m_R3/*.nc')
all_files.sort()

lats_median = np.zeros([63,11])
lons_median = np.zeros([63,11])
depths_median = np.zeros([63,11])

lats_traj = np.zeros([63,500,11])
lons_traj = np.zeros([63,500,11])
depths_traj = np.zeros([63,500,11])

dist_list = np.zeros([63])

i = 0
##Loops through each netcdf file
for timestep, datafile in enumerate(all_files):
    fr = xr.open_dataset(datafile)
    time = pd.DataFrame(np.array(fr.variables['time'][:].values, dtype='datetime64'))
    lons = np.array(fr.variables['lon'][:,:11])
    lats = np.array(fr.variables['lat'][:,:11])
    depths = np.array(fr.variables['z'][:,:11])
    depths = -1*depths
    fr.close()

    lats_median_01 = np.nanmedian(lats, axis=0)
    lons_median_01 = np.nanmedian(lons, axis=0)
    depths_median_01 = np.nanmedian(depths, axis=0)        
    
    lats_median[i] = lats_median_01
    lons_median[i] = lons_median_01
    depths_median[i] = depths_median_01
    i+=1


##############################################################
clusters = 6

data = pd.read_csv('Clusters_07/Mean_Trajectories_06Clusters.csv', header = None)
lon_mean_traj = np.array(data[:clusters])
lat_mean_traj = np.array(data[clusters:])

euclid_dist_list = np.zeros([len(lons_median),clusters])
k = 0
while k < clusters:
    i = 0
    while i < len(lons_median):
    
        lat_mean_traj_inp = lat_mean_traj[k] - lat_mean_traj[0,0]
        lon_mean_traj_inp = lon_mean_traj[k] - lon_mean_traj[0,0]
        lons_inp = lons_median[i] - lons_median[0,0]
        lats_inp = lats_median[i] - lats_median[0,0]
        
        euclid_dist = np.sqrt((lons_inp-lon_mean_traj_inp)**2+(lats_inp-lat_mean_traj_inp)**2)
        euclid_dist_list[i,k] = np.sum(euclid_dist)
        i+=1
    k+=1

class_list = np.zeros(len(euclid_dist_list))
j = 0
while j < len(euclid_dist_list):
    class_var = np.argmin(euclid_dist_list[j])
    class_list[j] = class_var
    j+=1

lat_mean_traj = np.zeros([clusters,11])
lon_mean_traj = np.zeros([clusters,11])
depth_mean_traj = np.zeros([clusters,11])

RMSE_list = np.zeros(clusters)
class_list = pd.DataFrame(class_list)
i = 0
while i < clusters:    
    globals()["class_"+str(i)] = class_list[class_list[0] == i].index.values
    globals()["class_"+str(i)+"_lons"] = np.zeros([len(globals()["class_"+str(i)]),11])
    globals()["class_"+str(i)+"_lats"] = np.zeros([len(globals()["class_"+str(i)]),11])
    globals()["class_"+str(i)+"_depths"] = np.zeros([len(globals()["class_"+str(i)]),11])
    j = 0
    while j < len(globals()["class_"+str(i)+"_lons"]):
        globals()["class_"+str(i)+"_lons"][j] = lons_median[globals()["class_"+str(i)][j]]
        lon_mean_traj[i] = np.mean(globals()["class_"+str(i)+"_lons"],0)
        globals()["class_"+str(i)+"_lats"][j] = lats_median[globals()["class_"+str(i)][j]]
        lat_mean_traj[i] = np.mean(globals()["class_"+str(i)+"_lats"],0)
        globals()["class_"+str(i)+"_depths"][j] = depths_median[globals()["class_"+str(i)][j]]
        depth_mean_traj[i] = np.mean(globals()["class_"+str(i)+"_depths"],0)
        
        j+=1    
    i+=1

i = 0
while i < clusters:    
    fig, ax = plt.subplots(figsize=(9,7))
    ax.plot(lon_mean_traj[i],depth_mean_traj[i], c = 'k')
    ax.plot(globals()["class_"+str(i)+"_lons"].T,globals()["class_"+str(i)+"_depths"].T)
    ax.set_xlabel("Longitude [deg]", fontsize = 10)
    ax.set_ylabel("Depth [m]", fontsize = 10)
    ax.set_xticks(range(32, 36, 1))
    ax.set_yticks(range(-70, -30, 5))
    ax.set_xticklabels(range(32, 36, 1), fontsize=9)
    ax.set_yticklabels(range(-70, -30, 5), fontsize=9)
    ax.set_title("Class "+str(i))
    #plt.savefig('Cross_Section_Class_'+str(i)+'_trajectory.png')   
    i+=1

# fig, ax = plt.subplots(figsize=(9,7))
# ax.plot(lon_mean_traj.T,depth_mean_traj.T)
# ax.set_xlabel("Longitude [deg]", fontsize = 10)
# ax.set_ylabel("Latitude [deg]", fontsize = 10)
# ax.set_xticks(range(30, 40, 1))
# ax.set_yticks(range(-30, -20, 1))
# ax.set_xticklabels(range(30, 40, 1), fontsize=9)
# ax.set_yticklabels(range(-30, -20, 1), fontsize=9)
# plt.savefig('Cross_Section_Class_trajectories.png') 





