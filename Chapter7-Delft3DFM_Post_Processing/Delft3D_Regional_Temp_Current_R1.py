## This program takes co-ordinates in LAT LON from an input file and generates a CSV that can be read in by Delft3D as a temperature boundary input file.
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import xarray as xr
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy.interpolate import griddata
import scipy
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
all_files = glob.glob('../2004_02/01b_cdb.dsproj_data/FlowFM/input/output/FlowFM_map.nc')
all_files.sort()
datafile = all_files[0]
    
###############################################

depth_layer_30 = 44

fh = xr.open_dataset(datafile)

time = pd.DataFrame(fh.variables['time'][:])
lon = np.array(fh.variables['mesh2d_face_x'][:])
lat = np.array(fh.variables['mesh2d_face_y'][:])
co_ords = np.array([lon,lat])
co_ords = np.transpose(co_ords)

grid_x, grid_y = np.mgrid[np.min(lon):np.max(lon):2000j, np.min(lat):np.max(lat):2000j]
grid_x_vel, grid_y_vel = np.mgrid[np.min(lon):np.max(lon):250j, np.min(lat):np.max(lat):250j]


## 9 Days before anomaly ##


ts = 19

delft_temp_05 = np.array(fh.variables['mesh2d_tem1'][ts,:,depth_layer_30])

delft_temp_grid_05 = griddata(co_ords, delft_temp_05, (grid_x, grid_y), method='linear')
delft_temp_grid_05 = np.flip(delft_temp_grid_05, 1)
delft_temp_grid_05 = np.transpose(delft_temp_grid_05)

delft_uvel_05 = np.array(fh.variables['mesh2d_ucx'][ts,:,depth_layer_30])
delft_uvel_grid_05 = griddata(co_ords, delft_uvel_05, (grid_x_vel, grid_y_vel), method='linear')
#delft_uvel_grid_05 = np.flip(delft_uvel_grid_05, 1)
#delft_uvel_grid_05 = np.transpose(delft_uvel_grid_05)

delft_vvel_05 = np.array(fh.variables['mesh2d_ucy'][ts,:,depth_layer_30])
delft_vvel_grid_05 = griddata(co_ords, delft_vvel_05, (grid_x_vel, grid_y_vel), method='linear')
#delft_vvel_grid_05 = np.flip(delft_vvel_grid_05, 1)
#delft_vvel_grid_05 = np.transpose(delft_vvel_grid_05)

delft_curr_grid_05 = np.sqrt(delft_uvel_grid_05**2 + delft_vvel_grid_05**2)
delft_curr_grid_05 = np.flip(delft_curr_grid_05, 1)
delft_curr_grid_05 = np.transpose(delft_curr_grid_05)

date = np.array(time)
date_str = str(date[ts])
date_str_05 = date_str[2:12]

## 15 Days before anomaly ##

ts = 14

delft_temp_06 = np.array(fh.variables['mesh2d_tem1'][ts,:,depth_layer_30])

delft_temp_grid_06 = griddata(co_ords, delft_temp_06, (grid_x, grid_y), method='linear')
delft_temp_grid_06 = np.flip(delft_temp_grid_06, 1)
delft_temp_grid_06 = np.transpose(delft_temp_grid_06)

delft_uvel_06 = np.array(fh.variables['mesh2d_ucx'][ts,:,depth_layer_30])
delft_uvel_grid_06 = griddata(co_ords, delft_uvel_06, (grid_x_vel, grid_y_vel), method='linear')

delft_vvel_06 = np.array(fh.variables['mesh2d_ucy'][ts,:,depth_layer_30])
delft_vvel_grid_06 = griddata(co_ords, delft_vvel_06, (grid_x_vel, grid_y_vel), method='linear')

delft_curr_grid_06 = np.sqrt(delft_uvel_grid_06**2 + delft_vvel_grid_06**2)
delft_curr_grid_06 = np.flip(delft_curr_grid_06, 1)
delft_curr_grid_06 = np.transpose(delft_curr_grid_06)

date = np.array(time)
date_str = str(date[ts])
date_str_06 = date_str[2:12]

## 18 Days before anomaly ##

ts = 9

delft_temp_07 = np.array(fh.variables['mesh2d_tem1'][ts,:,depth_layer_30])

delft_temp_grid_07 = griddata(co_ords, delft_temp_07, (grid_x, grid_y), method='linear')
delft_temp_grid_07 = np.flip(delft_temp_grid_07, 1)
delft_temp_grid_07 = np.transpose(delft_temp_grid_07)

delft_uvel_07 = np.array(fh.variables['mesh2d_ucx'][ts,:,depth_layer_30])
delft_uvel_grid_07 = griddata(co_ords, delft_uvel_07, (grid_x_vel, grid_y_vel), method='linear')

delft_vvel_07 = np.array(fh.variables['mesh2d_ucy'][ts,:,depth_layer_30])
delft_vvel_grid_07 = griddata(co_ords, delft_vvel_07, (grid_x_vel, grid_y_vel), method='linear')

delft_curr_grid_07 = np.sqrt(delft_uvel_grid_07**2 + delft_vvel_grid_07**2)
delft_curr_grid_07 = np.flip(delft_curr_grid_07, 1)
delft_curr_grid_07 = np.transpose(delft_curr_grid_07)

date = np.array(time)
date_str = str(date[ts])
date_str_07 = date_str[2:12]

## load data file names    
all_files = glob.glob('../2004_01/01b_cdb.dsproj_data/FlowFM/input/output/FlowFM_map.nc')
all_files.sort()
datafile = all_files[0]
    
###############################################

depth_layer_30 = 44

fh = xr.open_dataset(datafile)
time = pd.DataFrame(fh.variables['time'][:])

## 21 Days before anomaly ##

ts = 29

delft_temp_08 = np.array(fh.variables['mesh2d_tem1'][ts,:,depth_layer_30])

delft_temp_grid_08 = griddata(co_ords, delft_temp_08, (grid_x, grid_y), method='linear')
delft_temp_grid_08 = np.flip(delft_temp_grid_08, 1)
delft_temp_grid_08 = np.transpose(delft_temp_grid_08)

delft_uvel_08 = np.array(fh.variables['mesh2d_ucx'][ts,:,depth_layer_30])
delft_uvel_grid_08 = griddata(co_ords, delft_uvel_08, (grid_x_vel, grid_y_vel), method='linear')

delft_vvel_08 = np.array(fh.variables['mesh2d_ucy'][ts,:,depth_layer_30])
delft_vvel_grid_08 = griddata(co_ords, delft_vvel_08, (grid_x_vel, grid_y_vel), method='linear')

delft_curr_grid_08 = np.sqrt(delft_uvel_grid_08**2 + delft_vvel_grid_08**2)
delft_curr_grid_08 = np.flip(delft_curr_grid_08, 1)
delft_curr_grid_08 = np.transpose(delft_curr_grid_08)

date = np.array(time)
date_str = str(date[ts])
date_str_08 = date_str[2:12]

#%%

tempmin = 24
tempmax = 27

lat_min = -29
lat_max = -22
lon_min = 32
lon_max = 38    

meridians = np.arange(10.,351.,2.)
parallels = np.arange(-351.,-10.,2.)
meridians = np.arange(10.,351.,2.)
parallels = np.arange(-351.,-10.,2.)


fig, ax = plt.subplots(2,4, subplot_kw={'projection':ccrs.PlateCarree()},figsize=(9,8))

im0 = ax[0,0].imshow(delft_temp_grid_08, vmin = tempmin, vmax = tempmax, cmap='jet', extent=[np.min(lon), np.max(lon), np.min(lat), np.max(lat)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0,0], resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax[0,0].coastlines()
ax[0,0].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1) 
ax[0,0].plot(32.68, -27.49, 'ro', markersize=3, color = 'deeppink')
ax[0,0].set_xlim(lon_min, lon_max)
ax[0,0].set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
ax[0,0].set_title(date_str_08)
quiv = ax[0,0].quiver(grid_x_vel[::2,::2], grid_y_vel[::2,::2], delft_uvel_grid_08[::2,::2], delft_vvel_grid_08[::2,::2], scale=10, zorder=0, color='white')

im = ax[0,1].imshow(delft_temp_grid_07, vmin = tempmin, vmax = tempmax, cmap='jet', extent=[np.min(lon), np.max(lon), np.min(lat), np.max(lat)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0,1], resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax[0,1].coastlines()
ax[0,1].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)  
ax[0,1].plot(32.68, -27.49, 'ro', markersize=3, color = 'deeppink')
ax[0,1].set_xlim(lon_min, lon_max)
ax[0,1].set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
ax[0,1].set_title(date_str_07)
quiv = ax[0,1].quiver(grid_x_vel[::2,::2], grid_y_vel[::2,::2], delft_uvel_grid_07[::2,::2], delft_vvel_grid_07[::2,::2], scale=10, zorder=0, color='white')

im = ax[0,2].imshow(delft_temp_grid_06, vmin = tempmin, vmax = tempmax, cmap='jet', extent=[np.min(lon), np.max(lon), np.min(lat), np.max(lat)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0,2], resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax[0,2].coastlines()
ax[0,2].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1) 
ax[0,2].plot(32.68, -27.49, 'ro', markersize=3, color = 'deeppink')
ax[0,2].set_xlim(lon_min, lon_max)
ax[0,2].set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
ax[0,2].set_title(date_str_06)
quiv = ax[0,2].quiver(grid_x_vel[::2,::2], grid_y_vel[::2,::2], delft_uvel_grid_06[::2,::2], delft_vvel_grid_06[::2,::2], scale=10, zorder=0, color='white')

im = ax[0,3].imshow(delft_temp_grid_05, vmin = tempmin, vmax = tempmax, cmap='jet', extent=[np.min(lon), np.max(lon), np.min(lat), np.max(lat)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[0,3], resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax[0,3].coastlines()
ax[0,3].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1) 
ax[0,3].plot(32.68, -27.49, 'ro', markersize=3, color = 'deeppink')
ax[0,3].set_xlim(lon_min, lon_max)
ax[0,3].set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
ax[0,3].set_title(date_str_05)
quiv = ax[0,3].quiver(grid_x_vel[::2,::2], grid_y_vel[::2,::2], delft_uvel_grid_05[::2,::2], delft_vvel_grid_05[::2,::2], scale=10, zorder=0, color='white')

im0 = ax[1,0].imshow(delft_curr_grid_08, vmin = 0, vmax = 1.4, cmap=cmap, extent=[np.min(lon), np.max(lon), np.min(lat), np.max(lat)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1,0], resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax[1,0].coastlines()
ax[1,0].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)  
ax[1,0].plot(32.68, -27.49, 'ro', markersize=3, color = 'deeppink')
ax[1,0].set_xlim(lon_min, lon_max)
ax[1,0].set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
ax[1,0].set_title(date_str_08)
quiv = ax[1,0].quiver(grid_x_vel[::2,::2], grid_y_vel[::2,::2], delft_uvel_grid_08[::2,::2], delft_vvel_grid_08[::2,::2], scale=10, zorder=0, color='k')


im0 = ax[1,1].imshow(delft_curr_grid_07, vmin = 0, vmax = 1.4, cmap=cmap, extent=[np.min(lon), np.max(lon), np.min(lat), np.max(lat)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1,1], resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax[1,1].coastlines()
ax[1,1].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)   
ax[1,1].plot(32.68, -27.49, 'ro', markersize=3, color = 'deeppink')
ax[1,1].set_xlim(lon_min, lon_max)
ax[1,1].set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
ax[1,1].set_title(date_str_07)
quiv = ax[1,1].quiver(grid_x_vel[::2,::2], grid_y_vel[::2,::2], delft_uvel_grid_07[::2,::2], delft_vvel_grid_07[::2,::2], scale=10, zorder=0, color='k')

im0 = ax[1,2].imshow(delft_curr_grid_06, vmin = 0, vmax = 1.4, cmap=cmap, extent=[np.min(lon), np.max(lon), np.min(lat), np.max(lat)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1,2], resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax[1,2].coastlines()
ax[1,2].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)  
ax[1,2].plot(32.68, -27.49, 'ro', markersize=3, color = 'deeppink')
ax[1,2].set_xlim(lon_min, lon_max)
ax[1,2].set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
ax[1,2].set_title(date_str_06)
quiv = ax[1,2].quiver(grid_x_vel[::2,::2], grid_y_vel[::2,::2], delft_uvel_grid_06[::2,::2], delft_vvel_grid_06[::2,::2], scale=10, zorder=0, color='k')

im0 = ax[1,3].imshow(delft_curr_grid_05, vmin = 0, vmax = 1.4, cmap=cmap, extent=[np.min(lon), np.max(lon), np.min(lat), np.max(lat)], interpolation='bilinear');
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, ax = ax[1,3], resolution = 'l', llcrnrlat=lat_min, urcrnrlat = lat_max, llcrnrlon = lon_min, urcrnrlon = lon_max)
ax[1,3].coastlines()
ax[1,3].add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)  
ax[1,3].plot(32.68, -27.49, 'ro', markersize=3, color = 'deeppink')
ax[1,3].set_xlim(lon_min, lon_max)
ax[1,3].set_ylim(lat_min, lat_max)
m.drawmeridians(meridians,labels=[True,False,False,True],zorder=-2, fontsize=9)
m.drawparallels(parallels,labels=[False,True,True,False],zorder=-2, fontsize=9)
ax[1,3].set_title(date_str_05)
quiv = ax[1,3].quiver(grid_x_vel[::2,::2], grid_y_vel[::2,::2], delft_uvel_grid_05[::2,::2], delft_vvel_grid_05[::2,::2], scale=10, zorder=0, color='k')

cbar0 = fig.colorbar(im, ax=ax[0,:], orientation='horizontal', pad = 0.08)
cbar0.ax.tick_params(labelsize=9)
cbar0.ax.set_xlabel('Temperature [deg C]', fontsize = 11)

cbar0 = fig.colorbar(im0, ax=ax[1,:], orientation='horizontal', pad = 0.08)
cbar0.ax.tick_params(labelsize=9)
cbar0.ax.set_xlabel('Current Speed [m/s]', fontsize = 11)




