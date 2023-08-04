## Classification using trajectory origin and normalising variables

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

all_files = glob.glob('../../Paticle_Tracking_Runs_20m_R3/*.nc')
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

lat_seed_traj = lats_median[::6]
lon_seed_traj = lons_median[::6]

clusters = 11
euclid_dist_list = np.zeros([len(lons_median),int(clusters)])

k = 0
while k < clusters:
    i = 0
    while i < len(lons_median):
    
        lat_seed_traj_inp = lat_seed_traj[k] - lat_seed_traj[0,0]
        lon_seed_traj_inp = lon_seed_traj[k] - lon_seed_traj[0,0]
    
        lons_inp = lons_median[i] - lons_median[0,0]
        lats_inp = lats_median[i] - lats_median[0,0]
        euclid_dist = np.sqrt((lons_inp-lon_seed_traj_inp)**2+(lats_inp-lat_seed_traj_inp)**2)
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

class_list = pd.DataFrame(class_list)
i = 0
while i < clusters:    
    globals()["class_"+str(i)] = class_list[class_list[0] == i].index.values
    globals()["class_"+str(i)+"_lons"] = np.zeros([len(globals()["class_"+str(i)]),11])
    globals()["class_"+str(i)+"_lats"] = np.zeros([len(globals()["class_"+str(i)]),11])
    j = 0
    while j < len(globals()["class_"+str(i)+"_lons"]):
        globals()["class_"+str(i)+"_lons"][j] = lons_median[globals()["class_"+str(i)][j]]
        lon_mean_traj[i] = np.mean(globals()["class_"+str(i)+"_lons"],0)
        globals()["class_"+str(i)+"_lats"][j] = lats_median[globals()["class_"+str(i)][j]]
        lat_mean_traj[i] = np.mean(globals()["class_"+str(i)+"_lats"],0)
        j+=1    
    i+=1

##############################################################


euclid_dist_list = np.zeros([len(lons_median),int(clusters)])
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

RMSE_list = np.zeros(clusters)
class_list = pd.DataFrame(class_list)
i = 0
while i < clusters:    
    globals()["class_"+str(i)] = class_list[class_list[0] == i].index.values
    globals()["class_"+str(i)+"_lons"] = np.zeros([len(globals()["class_"+str(i)]),11])
    globals()["class_"+str(i)+"_lats"] = np.zeros([len(globals()["class_"+str(i)]),11])
    j = 0
    while j < len(globals()["class_"+str(i)+"_lons"]):
        globals()["class_"+str(i)+"_lons"][j] = lons_median[globals()["class_"+str(i)][j]]
        lon_mean_traj[i] = np.mean(globals()["class_"+str(i)+"_lons"],0)
        globals()["class_"+str(i)+"_lats"][j] = lats_median[globals()["class_"+str(i)][j]]
        lat_mean_traj[i] = np.mean(globals()["class_"+str(i)+"_lats"],0)
        j+=1    
    lon_diff = globals()["class_"+str(i)+"_lons"] - lon_mean_traj[i]
    lat_diff = globals()["class_"+str(i)+"_lats"] - lat_mean_traj[i]
    dist = np.sqrt(lon_diff**2 + lat_diff**2)
    SE = np.square(dist)
    MSE = np.mean(SE,1)
    RMSE = np.sqrt(MSE)
    RMSE_sum = np.sum(RMSE)
    RMSE_list[i] = RMSE_sum
    i+=1

RMSE_total = np.sum(RMSE_list)
RMSE_print = np.zeros([1,11])
RMSE_print[0,0] = RMSE_total

mean_traj = np.append(lon_mean_traj,lat_mean_traj, axis=0)

############ Write in code to save mean_traj
np.savetxt('Mean_Trajectories_11Clusters.csv', mean_traj, delimiter=',')

mean_dist_list = np.zeros([int(clusters),int(clusters)])
k = 0
while k < clusters:
    i = 0
    while i < clusters:    
        lat_mean_traj_inp = lat_mean_traj[k] - lat_mean_traj[0,0]
        lon_mean_traj_inp = lon_mean_traj[k] - lon_mean_traj[0,0]
        lons_inp = lon_mean_traj[i] - lon_mean_traj[0,0]
        lats_inp = lat_mean_traj[i] - lat_mean_traj[0,0]
        
        mean_dist = np.sqrt((lons_inp-lon_mean_traj_inp)**2+(lats_inp-lat_mean_traj_inp)**2)        
        mean_dist_list[i,k] = np.sum(mean_dist)
        mean_dist_list = np.tril(mean_dist_list)
        mean_dist_list[mean_dist_list == 0] = 9999
        closest = np.where(mean_dist_list == np.amin(mean_dist_list))
        i+=1
    k+=1

lon_mean_avg01 = lon_mean_traj[int(closest[0])]
lon_mean_avg02 = lon_mean_traj[int(closest[1])]
lat_mean_avg01 = lat_mean_traj[int(closest[0])]
lat_mean_avg02 = lat_mean_traj[int(closest[1])]
lon_merged = (lon_mean_avg01 + lon_mean_avg02) / 2
lat_merged = (lat_mean_avg01 + lat_mean_avg02) / 2

i = 0
while i < clusters:    
    fig, ax = plt.subplots(figsize=(9,7))
    ax.plot(lon_mean_traj[i],lat_mean_traj[i], c = 'k')
    ax.plot(globals()["class_"+str(i)+"_lons"].T,globals()["class_"+str(i)+"_lats"].T)
    m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax, resolution = 'l', llcrnrlat=-30, urcrnrlat = -20.5, llcrnrlon = 30, urcrnrlon = 40)
    m.drawcoastlines()
    m.fillcontinents(color='grey')
    border_color = 'black'  
    ax.plot(32.65, -27.53, 'ro', markersize=8, color = 'deeppink')
    ax.set_xlabel("Longitude [deg]", fontsize = 10)
    ax.set_ylabel("Latitude [deg]", fontsize = 10)
    ax.set_xticks(range(30, 40, 1))
    ax.set_yticks(range(-30, -20, 1))
    ax.set_xticklabels(range(30, 40, 1), fontsize=9)
    ax.set_yticklabels(range(-30, -20, 1), fontsize=9)
    i+=1

lon_mean_traj_new = np.array(lon_mean_traj)
lon_mean_traj_new[int(closest[0])] = lon_merged
lon_mean_traj_new = np.delete(lon_mean_traj_new,int(closest[1]),0)

lat_mean_traj_new = np.array(lat_mean_traj)
lat_mean_traj_new[int(closest[0])] = lat_merged
lat_mean_traj_new = np.delete(lat_mean_traj_new,int(closest[1]),0)

mean_traj_new = np.append(lon_mean_traj_new,lat_mean_traj_new, axis=0)
np.savetxt('Mean_Trajectories_10Clusters.csv', mean_traj_new, delimiter=',')
np.savetxt('Class_List_11Clusters.csv', class_list, delimiter=',')







