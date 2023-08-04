import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import scipy.stats as stats
import seaborn as sns
from mpl_toolkits.basemap import Basemap, shiftgrid
bm = Basemap()
from scipy import interpolate
from sklearn import preprocessing
from pyproj import Proj
from scipy.interpolate import interp1d

convert = Proj("+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

all_files = glob.glob('../Paticle_Tracking_Runs_20m_R5/*.nc')
all_files.sort()

mean = np.zeros([63,11])
median = np.zeros([63,11])
std = np.zeros([63,11])
temp_norm_10 = np.zeros([63,500])
temp_10_all = np.zeros([63,500])
dist = np.zeros([63,11])
temp = np.zeros([63,11])
depth = np.zeros([63,11])

i = 0
while i < 63:
    
    fr = xr.open_dataset(all_files[i])  
    time_part = pd.DataFrame(np.array(fr.variables['time'][:].values, dtype='datetime64'))
    lons_part = fr.variables['lon'][:]
    lons = np.nanmean(lons_part, axis = 0)
    lons = lons[:11]
    lats_part = fr.variables['lat'][:]
    lats = np.nanmean(lats_part, axis = 0)
    lats = lats[:11]
    coord = convert(lons,lats) 
    x =  coord[0]
    x = x - x[0]
    y =  coord[1]
    y = y - y[0]
    dist_01 = np.sqrt(x**2+y**2)
    dist[i] = dist_01
    temp_part = fr.variables['t'][:]
    depth_part = fr.variables['z'][:]
    depth_01 = np.nanmean(depth_part, axis = 0)
    depth_01 = depth_01[:11]
    depth[i] = depth_01
    fr.close()
    
    temp_actual_mean = np.nanmean(temp_part,0)
    temp_actual_mean = temp_actual_mean[:11]
    temp[i] = temp_actual_mean
    
    temp_00 = np.array(temp_part[:,0])
    temp_01 = np.array(temp_part[:,1])
    temp_02 = np.array(temp_part[:,2])
    temp_03 = np.array(temp_part[:,3])
    temp_04 = np.array(temp_part[:,4])
    temp_05 = np.array(temp_part[:,5])
    temp_06 = np.array(temp_part[:,6])
    temp_07 = np.array(temp_part[:,7])
    temp_08 = np.array(temp_part[:,8])
    temp_09 = np.array(temp_part[:,9])
    temp_10 = np.array(temp_part[:,10])

    temp_cent = np.zeros([11,500])
    temp_cent[0] = temp_00 - temp_00
    temp_cent[1] = temp_00 - temp_01
    temp_cent[2] = temp_00 - temp_02
    temp_cent[3] = temp_00 - temp_03
    temp_cent[4] = temp_00 - temp_04
    temp_cent[5] = temp_00 - temp_05
    temp_cent[6] = temp_00 - temp_06
    temp_cent[7] = temp_00 - temp_07
    temp_cent[8] = temp_00 - temp_08
    temp_cent[9] = temp_00 - temp_09
    temp_cent[10] = temp_00 - temp_10
    
    temp_mean = np.mean(temp_cent,1)
    mean[i] = temp_mean
    
    temp_median = np.median(temp_cent,1)
    median[i] = temp_median
    
    temp_std = np.std(temp_cent,1)
    std[i] = temp_std
    
    i+=1

temp_chng = np.zeros([63,11])
depth_chng = np.zeros([63,11])

n = 0
while n < len(mean):
    
    temp_01 = mean[n]
    tc = temp_01 - temp_01[-1]
    temp_chng[n] = tc*-1
    
    depth_01 = depth[n]
    dc = depth_01[-1] - depth_01
    depth_chng[n] = dc
    
    n+=1

dist = dist/1000

temp_50 = np.zeros([len(dist),46])
temp_int = np.zeros([len(dist),46])
dist_int = np.zeros([len(dist),46])
x_new = np.linspace(0, 450, num=46, endpoint=True)
upwelling_int = np.zeros([len(dist),46])

k = 0
while k < len(dist):
    dist_shelf = 50
    temp_line = temp[k]
    upwelling_line = depth_chng[k]
    dist_line = dist[k]
    int_t = interp1d(dist_line, temp_line, bounds_error=False, fill_value=0)
    int_u = interp1d(dist_line, upwelling_line, bounds_error=False, fill_value=0)
    int_d = interp1d(dist_line, dist_line, bounds_error=False, fill_value=0)
    temp_int_01 = int_t(x_new)
    upwelling_int_01 = int_u(x_new)
    temp_int[k] = temp_int_01
    upwelling_int[k] = upwelling_int_01
    dist_int[k] = x_new
    temp_50_norm = temp_int_01[5]
    if temp_50_norm == 0:
        temp_50_norm = temp_int_01[4]
    temp_50_01 = temp_int_01 - temp_50_norm
    temp_50[k] = temp_50_01
    k+=1

temp_50[temp_50 < -10] = np.nan

###########################################


fig, ax = plt.subplots(figsize=(12,5))
ax.plot(dist.T,temp.T, c = 'k', alpha=0.7)
ax.set_xlabel("Pathway distance from Sodwana [km]", fontsize = 12)
ax.set_ylabel("Mean particle temperature [deg C]", fontsize = 12)
plt.xlim(0,450)
plt.ylim(20,30)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.invert_xaxis()
plt.savefig('Particle_Temp.png')

fig, ax = plt.subplots(figsize=(12,5))
ax.plot(dist.T,temp_chng.T, c = 'k', alpha=0.3)
ax.set_xlabel("Pathway distance from Sodwana [km]", fontsize = 12)
ax.set_ylabel("Mean particle temperature change [deg C]", fontsize = 12)
plt.xlim(0,450)
#plt.ylim(-4,2)
ax.invert_xaxis()

fig, ax = plt.subplots(figsize=(12,5))
ax.plot(dist_int.T,temp_50.T, c = 'k', alpha=0.3)
ax.set_xlabel("Pathway distance from Sodwana [km]", fontsize = 12)
ax.set_ylabel("Mean particle temperature change [deg C]", fontsize = 12)
plt.xlim(0,50)
plt.ylim(-3,1)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.invert_xaxis()
plt.savefig('Particle_Temp_Change.png')

fig, ax = plt.subplots(figsize=(12,5))
ax.plot(dist.T,-1*mean.T, c = 'k', alpha=0.3)
ax.set_xlabel("Pathway distance from Sodwana [km]", fontsize = 12)
ax.set_ylabel("Particle temperature centered at Anomaly [deg C]", fontsize = 12)
plt.xlim(0,450)
#plt.ylim(-4,2)
ax.invert_xaxis()

fig, ax = plt.subplots(figsize=(12,5))
ax.plot(dist.T,depth_chng.T, c='k')
ax.set_xlabel("Pathway distance from Sodwana [km]", fontsize = 12)
ax.set_ylabel("Mean particle upwelled distance [m]", fontsize = 12)
plt.xlim(0,450)
plt.ylim(-5,25)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.invert_xaxis()
plt.savefig('Upwelling_Dist.png')

#%%
# fig, ax = plt.subplots(figsize=(10,10))
# ax.scatter(depth_chng[:,0],temp_chng[:,0], c ='k')
# #ax.scatter(depth_chng[:],temp_chng[:])
# ax.set_xlabel("Average upwelled distance over ten days [m]", fontsize = 12)
# ax.set_ylabel("Mean particle temperaure change [deg C]", fontsize = 12)
# plt.xlim(-5,40)
# plt.ylim(-4,2)
# plt.savefig('Scatter_Temp_Upwelling.png')


# fig, ax = plt.subplots(figsize=(12,5))
# ax.plot(depth_chng.T,temp_chng.T, c = 'k', alpha=0.3)
# ax.set_xlabel("Average upwelled distance over ten days [m]", fontsize = 12)
# ax.set_ylabel("Mean particle temperaure change [deg C]", fontsize = 12)
# ax.invert_xaxis()

#%%

from matplotlib.collections import LineCollection

#cs = np.array([-2,2])
anom = 2

x = np.array(dist[anom])
x_01 = np.full(2,-10)
x_01 = np.append(x_01,x)

y = np.array(depth_chng[anom])
y_01 = np.full(2,50)
y_01 = np.append(y_01,y)

cols = np.array(temp_chng[anom])
cols_01 = np.zeros(13)
cols_01[2:] = cols
cols_01[0] = 3
cols_01[1] = -3

points = np.array([x_01, y_01]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, ax = plt.subplots()
lc = LineCollection(segments, cmap='jet')
lc.set_array(cols_01)
lc.set_linewidth(2)
line = ax.add_collection(lc)
fig.colorbar(line,ax=ax)
plt.xlim(0,450)
plt.ylim(0,20)
ax.invert_xaxis()

#%%

fig, ax = plt.subplots(figsize=(12,5))

i = 0
while i < 63: 
    
    x = np.array(dist[i])
    x_01 = np.full(2,-10)
    x_01 = np.append(x_01,x)

    y = np.array(depth_chng[i])
    y_01 = np.full(2,50)
    y_01 = np.append(y_01,y)

    cols = np.array(temp_chng[i])
    cols_01 = np.zeros(13)
    cols_01[2:] = cols
    cols_01[0] = 3
    cols_01[1] = -3

    points = np.array([x_01, y_01]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap='jet')
    lc.set_array(cols_01)
    #lc.set_linewidth(2)
    line = ax.add_collection(lc)
    #fig.colorbar(line,ax=ax)
    plt.xlim(0,450)
    plt.ylim(0,20)
    
    i+=1



depth_median = np.median(upwelling_int, 0)
depth_95 = np.percentile(upwelling_int, 95, 0)
depth_05 = np.percentile(upwelling_int, 5, 0)


ax.plot(x_new,depth_median, c='k', linewidth = 2)
ax.plot(x_new,depth_05, c='k', linewidth = 2, linestyle = '--')
ax.plot(x_new,depth_95, c='k', linewidth = 2, linestyle = '--')
ax.set_xlim([2, 500])
ax.set_ylim([-5, 20])
ax.invert_xaxis()
ax.set_xlabel("Pathway distance from Sodwana [km]", fontsize = 14)
ax.set_ylabel("Mean veritcal displacement [m]", fontsize = 14)
ax.tick_params(labelsize=12 )

cb = fig.colorbar(line,ax=ax)
cb.set_label(label='Mean particle temperature change [$^\circ$C]', fontsize = 14)
plt.tight_layout()
plt.savefig('upwelling_temp.PDF')

#%%

temp_final = np.array(temp_chng[:,0])
x = np.zeros(21)
p = np.zeros(21)

i = 0
while i < len(p)*5:
    
    k = int(i/5)
    
    x[k] = i
    p[k] = np.percentile(temp_final, i)    
    i+=5
    
      
    

fig01, ax = plt.subplots(figsize=(8,8))
ax.plot(p, x, color='k')
plt.xlabel('Mean particle temperature change [$^\circ$C]', fontsize = 14)
plt.ylabel('% Exceedance', fontsize = 14)
ax.set_xlim([-4, 1])
ax.set_ylim([0, 100])
ax.tick_params(axis='both', which='major', labelsize=14)
plt.savefig('temp_change.PDF')
