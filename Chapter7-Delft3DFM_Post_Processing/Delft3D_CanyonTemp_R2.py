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
import cv2

## load in Delft3D map file
all_files = glob.glob('../2004_02/01b_cdb.dsproj_data/FlowFM/input/output/FlowFM_map.nc')
all_files.sort()
datafile = all_files[0]

fh = xr.open_dataset(datafile)
#time = pd.DataFrame(fh.variables['time'][:])
time = np.array(fh.variables['time'][:])
x = np.array(fh.variables['mesh2d_face_x'][:])
y = np.array(fh.variables['mesh2d_face_y'][:])

## Reduces maop extent to speed up interpolation from points onto map grid
lon_min = 32.647
lon_max = 32.78

lat_min = -27.559
lat_max = -27.405

x[x < lon_min] = np.nan
x[x > lon_max] = np.nan
y[y < lat_min] = np.nan
y[y > lat_max] = np.nan    

co_ords = np.array([x,y])
co_ords = np.transpose(co_ords)
co_ords = pd.DataFrame(co_ords)
co_ords.loc[co_ords.isnull().any(axis=1), :] = np.nan
co_ords = np.array(co_ords)

#################################

## Specify model depth layer
#depthLayer = 44 ## 30m depth
depthLayer = 41 ## 60m depth
#depthLayer = 39 ## 90m depth

timestep = 18

co_ords_zoom = co_ords[np.logical_not(np.isnan(co_ords[:,0]))]

##Extracts variables from map file
x_vel = np.array(fh.variables['mesh2d_ucx'][timestep,:,depthLayer])
x_vel_zoom_60 = x_vel[np.logical_not(np.isnan(co_ords[:,0]))]
y_vel = np.array(fh.variables['mesh2d_ucy'][timestep,:,depthLayer])
y_vel_zoom_60 = y_vel[np.logical_not(np.isnan(co_ords[:,0]))]
z_vel = np.array(fh.variables['mesh2d_ucz'][timestep,:,depthLayer])
z_vel_zoom_60 = z_vel[np.logical_not(np.isnan(co_ords[:,0]))]
z_vel_zoom_60[4106] = z_vel_zoom_60[4105]


temp = np.array(fh.variables['mesh2d_tem1'][timestep,:,depthLayer])
temp_zoom_60 = temp[np.logical_not(np.isnan(co_ords[:,0]))]
temp_zoom_60[4106] = temp_zoom_60[4105]

ts = str(time[timestep])
ts = ts[:10]

## Interpolates model output onto user defined grid
grid_x, grid_y = np.mgrid[np.min(co_ords_zoom[:,0]):np.max(co_ords_zoom[:,0]):2000j, np.min(co_ords_zoom[:,1]):np.max(co_ords_zoom[:,1]):2000j]

temp_grid_60 = griddata(co_ords_zoom, temp_zoom_60, (grid_x, grid_y), method='linear')
temp_grid_60 = np.flip(temp_grid_60, 1)
temp_grid_60 = np.transpose(temp_grid_60)

# grid_y_vel, grid_x_vel = np.mgrid[np.min(y):np.max(y):100j, np.min(x):np.max(x):100j]

# ## Normalises current speed if required
# normalized_x_vel = x_vel / np.sqrt(x_vel**2 + y_vel**2)
# normalized_y_vel = y_vel / np.sqrt(x_vel**2 + y_vel**2)

# normalized_current_speed = np.sqrt(normalized_x_vel**2 + normalized_y_vel**2)

## Specify model depth layer
depthLayer = 44 ## 30m depth
#depthLayer = 41 ## 60m depth
#depthLayer = 39 ## 90m depth

timestep = 18

co_ords_zoom = co_ords[np.logical_not(np.isnan(co_ords[:,0]))]

##Extracts variables from map file
x_vel = np.array(fh.variables['mesh2d_ucx'][timestep,:,depthLayer])
x_vel_zoom_30 = x_vel[np.logical_not(np.isnan(co_ords[:,0]))]
y_vel = np.array(fh.variables['mesh2d_ucy'][timestep,:,depthLayer])
y_vel_zoom_30 = y_vel[np.logical_not(np.isnan(co_ords[:,0]))]
z_vel = np.array(fh.variables['mesh2d_ucz'][timestep,:,depthLayer])
z_vel_zoom_30 = z_vel[np.logical_not(np.isnan(co_ords[:,0]))]
z_vel_zoom_30[4106] = z_vel_zoom_30[4105]


temp = np.array(fh.variables['mesh2d_tem1'][timestep,:,depthLayer])
temp_zoom_30 = temp[np.logical_not(np.isnan(co_ords[:,0]))]
temp_zoom_30[4106] = temp_zoom_30[4105]
temp_zoom_30[3988] = temp_zoom_30[3987]
temp_zoom_30[5806] = temp_zoom_30[5805]

ts = str(time[timestep])
ts = ts[:10]

## Interpolates model output onto user defined grid
grid_x, grid_y = np.mgrid[np.min(co_ords_zoom[:,0]):np.max(co_ords_zoom[:,0]):2000j, np.min(co_ords_zoom[:,1]):np.max(co_ords_zoom[:,1]):2000j]

temp_grid_30 = griddata(co_ords_zoom, temp_zoom_30, (grid_x, grid_y), method='linear')
temp_grid_30 = np.flip(temp_grid_30, 1)
temp_grid_30 = np.transpose(temp_grid_30)

# grid_y_vel, grid_x_vel = np.mgrid[np.min(y):np.max(y):100j, np.min(x):np.max(x):100j]

# ## Normalises current speed if required
# normalized_x_vel = x_vel / np.sqrt(x_vel**2 + y_vel**2)
# normalized_y_vel = y_vel / np.sqrt(x_vel**2 + y_vel**2)

########################################

# normalized_current_speed = np.sqrt(normalized_x_vel**2 + normalized_y_vel**2)

## Specify model depth layer
#depthLayer = 44 ## 30m depth
#depthLayer = 41 ## 60m depth
depthLayer = 39 ## 90m depth

timestep = 18

co_ords_zoom = co_ords[np.logical_not(np.isnan(co_ords[:,0]))]

##Extracts variables from map file
x_vel = np.array(fh.variables['mesh2d_ucx'][timestep,:,depthLayer])
x_vel_zoom_90 = x_vel[np.logical_not(np.isnan(co_ords[:,0]))]
y_vel = np.array(fh.variables['mesh2d_ucy'][timestep,:,depthLayer])
y_vel_zoom_90 = y_vel[np.logical_not(np.isnan(co_ords[:,0]))]
z_vel = np.array(fh.variables['mesh2d_ucz'][timestep,:,depthLayer])
z_vel_zoom_90 = z_vel[np.logical_not(np.isnan(co_ords[:,0]))]
z_vel_zoom_90[4106] = z_vel_zoom_90[4105]


temp = np.array(fh.variables['mesh2d_tem1'][timestep,:,depthLayer])
temp_zoom_90 = temp[np.logical_not(np.isnan(co_ords[:,0]))]
temp_zoom_90[4106] = temp_zoom_90[4105]
temp_zoom_90[4107] = temp_zoom_90[4106]

ts = str(time[timestep])
ts = ts[:10]

## Interpolates model output onto user defined grid
grid_x, grid_y = np.mgrid[np.min(co_ords_zoom[:,0]):np.max(co_ords_zoom[:,0]):2000j, np.min(co_ords_zoom[:,1]):np.max(co_ords_zoom[:,1]):2000j]

temp_grid_90 = griddata(co_ords_zoom, temp_zoom_90, (grid_x, grid_y), method='linear')
temp_grid_90 = np.flip(temp_grid_90, 1)
temp_grid_90 = np.transpose(temp_grid_90)

# grid_y_vel, grid_x_vel = np.mgrid[np.min(y):np.max(y):100j, np.min(x):np.max(x):100j]

# ## Normalises current speed if required
# normalized_x_vel = x_vel / np.sqrt(x_vel**2 + y_vel**2)
# normalized_y_vel = y_vel / np.sqrt(x_vel**2 + y_vel**2)

# normalized_current_speed = np.sqrt(normalized_x_vel**2 + normalized_y_vel**2)

#########################################

img= cv2.imread(r'Reef_Map_R0.png')
RGB_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

src= cv2.imread(r'Sodwana01_edited_R1.bmp')
tmp = cv2.cvtColor(src, cv2.COLOR_BGR2GRAY)
_,alpha = cv2.threshold(tmp,0,255,cv2.THRESH_BINARY)
b, g, r = cv2.split(src)
rgba = [b,g,r, alpha]
dst = cv2.merge(rgba,4)

########################################

#%%

fig, ax = plt.subplots(1, 3, figsize=(15,6.5))

#im1 = ax[0].scatter(co_ords_zoom[:,0], co_ords_zoom[:,1], vmin = 20.5, vmax = 22.5, c = temp_zoom_30, cmap='jet')
im2 = ax[0].imshow(RGB_img,extent=[32.639,  32.815, -27.57, -27.37])
im3 = ax[0].imshow(dst,extent=[32.647,  32.767, -27.559, -27.405], zorder=2)
im1 = ax[0].imshow(temp_grid_30, vmin = 20.5, vmax = 22.5, cmap='jet', extent=[np.min(co_ords_zoom[:,0]),  np.max(co_ords_zoom[:,0]), np.min(co_ords_zoom[:,1]), np.max(co_ords_zoom[:,1])])
ax[0].set_xlabel("Longitude [deg]", fontsize = 12)
ax[0].set_ylabel("Latitude [deg]", fontsize = 12)
ax[0].set_xlim([lon_min,  lon_max])
ax[0].set_ylim([lat_min,  lat_max])

cbar0 = fig.colorbar(im1, ax=ax[0], orientation='horizontal', pad = 0.095)
cbar0.ax.tick_params(labelsize=7)
cbar0.ax.set_xlabel('Temperature [$^\circ$C]', fontsize = 12)


im2 = ax[1].imshow(RGB_img,extent=[32.639,  32.815, -27.57, -27.37])
im3 = ax[1].imshow(dst,extent=[32.647,  32.767, -27.559, -27.405], zorder=2)
im1 = ax[1].imshow(temp_grid_60, vmin = 16.5, vmax = 18.5, cmap='jet', extent=[np.min(co_ords_zoom[:,0]),  np.max(co_ords_zoom[:,0]), np.min(co_ords_zoom[:,1]), np.max(co_ords_zoom[:,1])])
#im1 = ax[1].scatter(co_ords_zoom[:,0], co_ords_zoom[:,1], vmin = 16.5, vmax = 18.5, c = temp_zoom_60, cmap='jet')
ax[1].set_xlabel("Longitude [deg]", fontsize = 12)
ax[1].set_xlim([lon_min,  lon_max])
ax[1].set_ylim([lat_min,  lat_max])

cbar1 = fig.colorbar(im1, ax=ax[1], orientation='horizontal', pad = 0.095)
cbar1.ax.tick_params(labelsize=7)
cbar1.ax.set_xlabel('Temperature [$^\circ$C]', fontsize = 12)

#im1 = ax[2].scatter(co_ords_zoom[:,0], co_ords_zoom[:,1], vmin = 15.5, vmax = 17.5, c = temp_zoom_90, cmap='jet')
im2 = ax[2].imshow(RGB_img,extent=[32.639,  32.815, -27.57, -27.37])
im3 = ax[2].imshow(dst,extent=[32.647,  32.767, -27.559, -27.405], zorder=2)
im1 = ax[2].imshow(temp_grid_90, vmin = 15.5, vmax = 17.5, cmap='jet', extent=[np.min(co_ords_zoom[:,0]),  np.max(co_ords_zoom[:,0]), np.min(co_ords_zoom[:,1]), np.max(co_ords_zoom[:,1])])
ax[2].set_xlabel("Longitude [deg]", fontsize = 12)
ax[2].set_xlim([lon_min,  lon_max])
ax[2].set_ylim([lat_min,  lat_max])

cbar2 = fig.colorbar(im1, ax=ax[2], orientation='horizontal', pad = 0.095)
cbar2.ax.tick_params(labelsize=7)
cbar2.ax.set_xlabel('Temperature [$^\circ$C]', fontsize = 12)

#quiv = ax.quiver(x, y, normalized_x_vel, normalized_y_vel, scale = 50, zorder=1)

plt.savefig('Canyon_Temperature.pdf')

















