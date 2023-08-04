## This program takes co-ordinates in LAT LON from an input file and generates a CSV that can be read in by Delft3D as a temperature boundary input file.

import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import xarray as xr
#from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy.interpolate import griddata
import cv2
import matplotlib.colors
from matplotlib.colors import LinearSegmentedColormap

## load data file names    
all_files = glob.glob('../2004_02/01b_cdb.dsproj_data/FlowFM/input/output/FlowFM_map.nc')
all_files.sort()
datafile = all_files[0]

###########################################

fh = xr.open_dataset(datafile)
time = np.array(fh.variables['time'][:])
x = np.array(fh.variables['mesh2d_face_x'][:])
y = np.array(fh.variables['mesh2d_face_y'][:])

co_ords = np.array([x,y])
co_ords = np.transpose(co_ords)
co_ords = pd.DataFrame(co_ords)
co_ords.loc[co_ords.isnull().any(axis=1), :] = np.nan
co_ords = np.array(co_ords)
co_ords_zoom = co_ords[np.logical_not(np.isnan(co_ords[:,0]))]

x_mask = np.array(x)
x_mask[x_mask > 33] = np.nan
x_mask[x_mask < 32] = np.nan
x_mask[x_mask > 0] = 1
y_mask = np.array(y)
y_mask[y_mask > -27] = np.nan
y_mask[y_mask < -28] = np.nan
y_mask[y_mask < 0] = 1        
mask = x_mask * y_mask
x_reduced = np.array(x) * mask
x_reduced = x_reduced[~np.isnan(x_reduced)]
y_reduced = np.array(y) * mask
y_reduced = y_reduced[~np.isnan(y_reduced)]
co_ords_reduced = np.array([x_reduced,y_reduced])
co_ords_reduced = np.transpose(co_ords_reduced)
grid_x_red, grid_y_red = np.mgrid[np.min(x_reduced):np.max(x_reduced):1000j, np.min(y_reduced):np.max(y_reduced):1000j]

depthLayer = 44

############################################

time_step_01 = 16
date_01 = time[time_step_01]
ts_01 = str(time[time_step_01])
ts_01 = ts_01[:10]

tempBottom_01 = np.zeros(len(fh.variables['mesh2d_tem1'][time_step_01,:,49]))
xVelBottom_01 = np.zeros(len(fh.variables['mesh2d_ucx'][time_step_01,:,49]))
yVelBottom_01 = np.zeros(len(fh.variables['mesh2d_ucy'][time_step_01,:,49]))

tempSlice_01 = np.array(fh.variables['mesh2d_tem1'][time_step_01,:,depthLayer])
xVelSlice_01 = np.array(fh.variables['mesh2d_ucx'][time_step_01,:,depthLayer])
yVelSlice_01 = np.array(fh.variables['mesh2d_ucy'][time_step_01,:,depthLayer])


depth_level = 1
while depth_level < 5:
    
    ## Finds layer below one I'm interested in - this is why there is a (-1) after
    ## minusing the depth layer
    realValueArray = np.array(fh.variables['mesh2d_ucx'][time_step_01,:,49-depth_level-1])
    landIdentifier = np.array(realValueArray)
    landIdentifier[landIdentifier < 0] = 0
    landIdentifier[landIdentifier >= 0] = 0
    landIdentifier = np.nan_to_num(landIdentifier, nan=1)
    
    tempSliceAtCurrentLayer = np.array(fh.variables['mesh2d_tem1'][time_step_01,:,49-depth_level])
    tempSliceAtCurrentLayer = np.nan_to_num(tempSliceAtCurrentLayer, nan=0) 
    tempSliceAtCurrentLayer = tempSliceAtCurrentLayer * landIdentifier
    tempBottom_01 = tempBottom_01 + tempSliceAtCurrentLayer
    
    xVelSliceAtCurrentLayer = np.array(fh.variables['mesh2d_ucx'][time_step_01,:,49-depth_level])
    xVelSliceAtCurrentLayer = np.nan_to_num(xVelSliceAtCurrentLayer, nan=0) 
    xVelSliceAtCurrentLayer = xVelSliceAtCurrentLayer * landIdentifier
    xVelBottom_01 = xVelBottom_01 + xVelSliceAtCurrentLayer
    
    yVelSliceAtCurrentLayer = np.array(fh.variables['mesh2d_ucy'][time_step_01,:,49-depth_level])
    yVelSliceAtCurrentLayer = np.nan_to_num(yVelSliceAtCurrentLayer, nan=0) 
    yVelSliceAtCurrentLayer = yVelSliceAtCurrentLayer * landIdentifier
    yVelBottom_01 = yVelBottom_01 + yVelSliceAtCurrentLayer
		 
    depth_level+=1

tempBottomReduced_01 = np.array(tempBottom_01) * mask
k = 0
while k < len(tempBottomReduced_01):
    if tempBottomReduced_01[k] == 0:
        tempBottomReduced_01[k] = tempSlice_01[k]
    k+=1 
tempBottomReduced_01 = tempBottomReduced_01[~np.isnan(tempBottomReduced_01)]
tempBottomGrid_01 = griddata(co_ords_reduced, tempBottomReduced_01, (grid_x_red, grid_y_red), method='linear')
tempBottomGrid_01 = np.flip(tempBottomGrid_01, 1)
tempBottomGrid_01 = np.transpose(tempBottomGrid_01)

xVelBottomReduced_01 = np.array(xVelBottom_01) * mask
k = 0
while k < len(xVelBottomReduced_01):
    if xVelBottomReduced_01[k] == 0:
        xVelBottomReduced_01[k] = xVelSlice_01[k]
    k+=1   
xVelBottomReduced_01 = xVelBottomReduced_01[~np.isnan(xVelBottomReduced_01)]
xVelBottomGrid_01 = griddata(co_ords_reduced, xVelBottomReduced_01, (grid_x_red, grid_y_red), method='linear')

yVelBottomReduced_01 = np.array(yVelBottom_01) * mask
k = 0
while k < len(yVelBottomReduced_01):
    if yVelBottomReduced_01[k] == 0:
        yVelBottomReduced_01[k] = yVelSlice_01[k]
    k+=1  
yVelBottomReduced_01 = yVelBottomReduced_01[~np.isnan(yVelBottomReduced_01)]
yVelBottomGrid_01 = griddata(co_ords_reduced, yVelBottomReduced_01, (grid_x_red, grid_y_red), method='linear')

############################################

time_step_02 = 18
date_02 = time[time_step_02]
ts_02 = str(time[time_step_02])
ts_02 = ts_02[:10]

tempBottom_02 = np.zeros(len(fh.variables['mesh2d_tem1'][time_step_02,:,49]))
xVelBottom_02 = np.zeros(len(fh.variables['mesh2d_ucx'][time_step_02,:,49]))
yVelBottom_02 = np.zeros(len(fh.variables['mesh2d_ucy'][time_step_02,:,49]))

tempSlice_02 = np.array(fh.variables['mesh2d_tem1'][time_step_02,:,depthLayer])
xVelSlice_02 = np.array(fh.variables['mesh2d_ucx'][time_step_02,:,depthLayer])
yVelSlice_02 = np.array(fh.variables['mesh2d_ucy'][time_step_02,:,depthLayer])

depth_level = 1
while depth_level < 5:
    
    ## Finds layer below one I'm interested in - this is why there is a (-1) after
    ## minusing the depth layer
    realValueArray = np.array(fh.variables['mesh2d_ucx'][time_step_02,:,49-depth_level-1])
    landIdentifier = np.array(realValueArray)
    landIdentifier[landIdentifier < 0] = 0
    landIdentifier[landIdentifier >= 0] = 0
    landIdentifier = np.nan_to_num(landIdentifier, nan=1)
    
    tempSliceAtCurrentLayer = np.array(fh.variables['mesh2d_tem1'][time_step_02,:,49-depth_level])
    tempSliceAtCurrentLayer = np.nan_to_num(tempSliceAtCurrentLayer, nan=0) 
    tempSliceAtCurrentLayer = tempSliceAtCurrentLayer * landIdentifier
    tempBottom_02 = tempBottom_02 + tempSliceAtCurrentLayer
    
    xVelSliceAtCurrentLayer = np.array(fh.variables['mesh2d_ucx'][time_step_02,:,49-depth_level])
    xVelSliceAtCurrentLayer = np.nan_to_num(xVelSliceAtCurrentLayer, nan=0) 
    xVelSliceAtCurrentLayer = xVelSliceAtCurrentLayer * landIdentifier
    xVelBottom_02 = xVelBottom_02 + xVelSliceAtCurrentLayer
    
    yVelSliceAtCurrentLayer = np.array(fh.variables['mesh2d_ucy'][time_step_02,:,49-depth_level])
    yVelSliceAtCurrentLayer = np.nan_to_num(yVelSliceAtCurrentLayer, nan=0) 
    yVelSliceAtCurrentLayer = yVelSliceAtCurrentLayer * landIdentifier
    yVelBottom_02 = yVelBottom_02 + yVelSliceAtCurrentLayer
	    
    depth_level+=1

tempBottomReduced_02 = np.array(tempBottom_02) * mask
k = 0
while k < len(tempBottomReduced_02):
    if tempBottomReduced_02[k] == 0:
        tempBottomReduced_02[k] = tempSlice_02[k]
    k+=1  
tempBottomReduced_02 = tempBottomReduced_02[~np.isnan(tempBottomReduced_02)]
tempBottomGrid_02 = griddata(co_ords_reduced, tempBottomReduced_02, (grid_x_red, grid_y_red), method='linear')
tempBottomGrid_02 = np.flip(tempBottomGrid_02, 1)
tempBottomGrid_02 = np.transpose(tempBottomGrid_02)

xVelBottomReduced_02 = np.array(xVelBottom_02) * mask
k = 0
while k < len(xVelBottomReduced_02):
    if xVelBottomReduced_02[k] == 0:
        xVelBottomReduced_02[k] = xVelSlice_02[k]
    k+=1   
xVelBottomReduced_02 = xVelBottomReduced_02[~np.isnan(xVelBottomReduced_02)]
xVelBottomGrid_02 = griddata(co_ords_reduced, xVelBottomReduced_02, (grid_x_red, grid_y_red), method='linear')

yVelBottomReduced_02 = np.array(yVelBottom_02) * mask
k = 0
while k < len(yVelBottomReduced_02):
    if yVelBottomReduced_02[k] == 0:
        yVelBottomReduced_02[k] = yVelSlice_02[k]
    k+=1  
yVelBottomReduced_02 = yVelBottomReduced_02[~np.isnan(yVelBottomReduced_02)]
yVelBottomGrid_02 = griddata(co_ords_reduced, yVelBottomReduced_02, (grid_x_red, grid_y_red), method='linear')

############################################

time_step_03 = 20
date_03 = time[time_step_03]
ts_03 = str(time[time_step_03])
ts_03 = ts_03[:10]

tempBottom_03 = np.zeros(len(fh.variables['mesh2d_tem1'][time_step_03,:,49]))
xVelBottom_03 = np.zeros(len(fh.variables['mesh2d_ucx'][time_step_03,:,49]))
yVelBottom_03 = np.zeros(len(fh.variables['mesh2d_ucy'][time_step_03,:,49]))

tempSlice_03 = np.array(fh.variables['mesh2d_tem1'][time_step_03,:,depthLayer])
xVelSlice_03 = np.array(fh.variables['mesh2d_ucx'][time_step_03,:,depthLayer])
yVelSlice_03 = np.array(fh.variables['mesh2d_ucy'][time_step_03,:,depthLayer])

depth_level = 1
while depth_level < 5:
    
    ## Finds layer below one I'm interested in - this is why there is a (-1) after
    ## minusing the depth layer
    realValueArray = np.array(fh.variables['mesh2d_ucx'][time_step_03,:,49-depth_level-1])
    landIdentifier = np.array(realValueArray)
    landIdentifier[landIdentifier < 0] = 0
    landIdentifier[landIdentifier >= 0] = 0
    landIdentifier = np.nan_to_num(landIdentifier, nan=1)
    
    tempSliceAtCurrentLayer = np.array(fh.variables['mesh2d_tem1'][time_step_03,:,49-depth_level])
    tempSliceAtCurrentLayer = np.nan_to_num(tempSliceAtCurrentLayer, nan=0) 
    tempSliceAtCurrentLayer = tempSliceAtCurrentLayer * landIdentifier
    tempBottom_03 = tempBottom_03 + tempSliceAtCurrentLayer
    
    xVelSliceAtCurrentLayer = np.array(fh.variables['mesh2d_ucx'][time_step_03,:,49-depth_level])
    xVelSliceAtCurrentLayer = np.nan_to_num(xVelSliceAtCurrentLayer, nan=0) 
    xVelSliceAtCurrentLayer = xVelSliceAtCurrentLayer * landIdentifier
    xVelBottom_03 = xVelBottom_03 + xVelSliceAtCurrentLayer
    
    yVelSliceAtCurrentLayer = np.array(fh.variables['mesh2d_ucy'][time_step_03,:,49-depth_level])
    yVelSliceAtCurrentLayer = np.nan_to_num(yVelSliceAtCurrentLayer, nan=0) 
    yVelSliceAtCurrentLayer = yVelSliceAtCurrentLayer * landIdentifier
    yVelBottom_03 = yVelBottom_03 + yVelSliceAtCurrentLayer
	
    depth_level+=1

tempBottomReduced_03 = np.array(tempBottom_03) * mask
k = 0
while k < len(tempBottomReduced_03):
    if tempBottomReduced_03[k] == 0:
        tempBottomReduced_03[k] = tempSlice_03[k]
    k+=1
tempBottomReduced_03 = tempBottomReduced_03[~np.isnan(tempBottomReduced_03)]
tempBottomGrid_03 = griddata(co_ords_reduced, tempBottomReduced_03, (grid_x_red, grid_y_red), method='linear')
tempBottomGrid_03 = np.flip(tempBottomGrid_03, 1)
tempBottomGrid_03 = np.transpose(tempBottomGrid_03)

xVelBottomReduced_03 = np.array(xVelBottom_03) * mask
k = 0
while k < len(xVelBottomReduced_03):
    if xVelBottomReduced_03[k] == 0:
        xVelBottomReduced_03[k] = xVelSlice_03[k]
    k+=1   
xVelBottomReduced_03 = xVelBottomReduced_03[~np.isnan(xVelBottomReduced_03)]
xVelBottomGrid_03 = griddata(co_ords_reduced, xVelBottomReduced_03, (grid_x_red, grid_y_red), method='linear')

yVelBottomReduced_03 = np.array(yVelBottom_03) * mask
k = 0
while k < len(yVelBottomReduced_03):
    if yVelBottomReduced_03[k] == 0:
        yVelBottomReduced_03[k] = yVelSlice_03[k]
    k+=1  
yVelBottomReduced_03 = yVelBottomReduced_03[~np.isnan(yVelBottomReduced_03)]
yVelBottomGrid_03 = griddata(co_ords_reduced, yVelBottomReduced_03, (grid_x_red, grid_y_red), method='linear')

####################################

time_step_04 = 22
date_04 = time[time_step_04]
ts_04 = str(time[time_step_04])
ts_04 = ts_04[:10]

tempBottom_04 = np.zeros(len(fh.variables['mesh2d_tem1'][time_step_04,:,49]))
xVelBottom_04 = np.zeros(len(fh.variables['mesh2d_ucx'][time_step_04,:,49]))
yVelBottom_04 = np.zeros(len(fh.variables['mesh2d_ucy'][time_step_04,:,49]))

tempSlice_04 = np.array(fh.variables['mesh2d_tem1'][time_step_04,:,depthLayer])
xVelSlice_04 = np.array(fh.variables['mesh2d_ucx'][time_step_04,:,depthLayer])
yVelSlice_04 = np.array(fh.variables['mesh2d_ucy'][time_step_04,:,depthLayer])

depth_level = 1
while depth_level < 5:
    
    ## Finds layer below one I'm interested in - this is why there is a (-1) after
    ## minusing the depth layer
    realValueArray = np.array(fh.variables['mesh2d_ucx'][time_step_04,:,49-depth_level-1])
    landIdentifier = np.array(realValueArray)
    landIdentifier[landIdentifier < 0] = 0
    landIdentifier[landIdentifier >= 0] = 0
    landIdentifier = np.nan_to_num(landIdentifier, nan=1)
    
    tempSliceAtCurrentLayer = np.array(fh.variables['mesh2d_tem1'][time_step_04,:,49-depth_level])
    tempSliceAtCurrentLayer = np.nan_to_num(tempSliceAtCurrentLayer, nan=0) 
    tempSliceAtCurrentLayer = tempSliceAtCurrentLayer * landIdentifier
    tempBottom_04 = tempBottom_04 + tempSliceAtCurrentLayer
    
    xVelSliceAtCurrentLayer = np.array(fh.variables['mesh2d_ucx'][time_step_04,:,49-depth_level])
    xVelSliceAtCurrentLayer = np.nan_to_num(xVelSliceAtCurrentLayer, nan=0) 
    xVelSliceAtCurrentLayer = xVelSliceAtCurrentLayer * landIdentifier
    xVelBottom_04 = xVelBottom_04 + xVelSliceAtCurrentLayer
    
    yVelSliceAtCurrentLayer = np.array(fh.variables['mesh2d_ucy'][time_step_04,:,49-depth_level])
    yVelSliceAtCurrentLayer = np.nan_to_num(yVelSliceAtCurrentLayer, nan=0) 
    yVelSliceAtCurrentLayer = yVelSliceAtCurrentLayer * landIdentifier
    yVelBottom_04 = yVelBottom_04 + yVelSliceAtCurrentLayer
	
    depth_level+=1
    
tempBottomReduced_04 = np.array(tempBottom_04) * mask
k = 0
while k < len(tempBottomReduced_04):
    if tempBottomReduced_04[k] == 0:
        tempBottomReduced_04[k] = tempSlice_04[k]
    k+=1    
tempBottomReduced_04 = tempBottomReduced_04[~np.isnan(tempBottomReduced_04)]
tempBottomGrid_04 = griddata(co_ords_reduced, tempBottomReduced_04, (grid_x_red, grid_y_red), method='linear')
tempBottomGrid_04 = np.flip(tempBottomGrid_04, 1)
tempBottomGrid_04 = np.transpose(tempBottomGrid_04)

xVelBottomReduced_04 = np.array(xVelBottom_04) * mask
k = 0
while k < len(xVelBottomReduced_04):
    if xVelBottomReduced_04[k] == 0:
        xVelBottomReduced_04[k] = xVelSlice_04[k]
    k+=1   
xVelBottomReduced_04 = xVelBottomReduced_04[~np.isnan(xVelBottomReduced_04)]
xVelBottomGrid_04 = griddata(co_ords_reduced, xVelBottomReduced_04, (grid_x_red, grid_y_red), method='linear')

yVelBottomReduced_04 = np.array(yVelBottom_04) * mask
k = 0
while k < len(yVelBottomReduced_04):
    if yVelBottomReduced_04[k] == 0:
        yVelBottomReduced_04[k] = yVelSlice_04[k]
    k+=1  
yVelBottomReduced_04 = yVelBottomReduced_04[~np.isnan(yVelBottomReduced_04)]
yVelBottomGrid_04 = griddata(co_ords_reduced, yVelBottomReduced_04, (grid_x_red, grid_y_red), method='linear')


#########################################

img= cv2.imread(r'Reefs_Edited_02.png')
RGB_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

src= cv2.imread(r'Sodwana01_edited_R1.bmp')
tmp = cv2.cvtColor(src, cv2.COLOR_BGR2GRAY)
_,alpha = cv2.threshold(tmp,0,255,cv2.THRESH_BINARY)
b, g, r = cv2.split(src)
rgba = [b,g,r, alpha]
dst = cv2.merge(rgba,4)

########################################

#%%

fig, ax = plt.subplots(1,4, figsize=(12,5))

im1 = ax[0].imshow(tempBottomGrid_01, vmin = 22, vmax = 26, cmap='jet', extent=[np.min(x_reduced), np.max(x_reduced), np.min(y_reduced), np.max(y_reduced)]);
#im2 = ax.imshow(RGB_img,extent=[32.65260892313,  32.76582166177, -27.54139396108, -27.40630195089])
im2 = ax[0].imshow(dst,extent=[32.647,  32.767, -27.559, -27.405], zorder=2)
ax[0].set_xlabel("Longitude [deg]", fontsize = 12)
ax[0].set_ylabel("Latitude [deg]", fontsize = 12)
ax[0].set_xlim([32.65,  32.77])
ax[0].set_ylim([-27.55, -27.41])
ax[0].set_title(ts_01)
quiv = ax[0].quiver(grid_x_red[::7,::5], grid_y_red[::7,::5], xVelBottomGrid_01[::7,::5], yVelBottomGrid_01[::7,::5], color = 'k', scale = 2)

im1 = ax[1].imshow(tempBottomGrid_02, vmin = 22, vmax = 26, cmap='jet', extent=[np.min(x_reduced), np.max(x_reduced), np.min(y_reduced), np.max(y_reduced)]);
im2 = ax[1].imshow(dst,extent=[32.647,  32.767, -27.559, -27.405], zorder=2)
ax[1].set_xlabel("Longitude [deg]", fontsize = 12)
ax[1].set_xlim([32.65,  32.77])
ax[1].set_ylim([-27.55, -27.41])
ax[1].set_yticks([])
ax[1].set_title(ts_02)
quiv = ax[1].quiver(grid_x_red[::7,::5], grid_y_red[::7,::5], xVelBottomGrid_02[::7,::5], yVelBottomGrid_02[::7,::5], color = 'k', scale = 2)

im1 = ax[2].imshow(tempBottomGrid_03, vmin = 22, vmax = 26, cmap='jet', extent=[np.min(x_reduced), np.max(x_reduced), np.min(y_reduced), np.max(y_reduced)]);
im2 = ax[2].imshow(dst,extent=[32.647,  32.767, -27.559, -27.405], zorder=2)
ax[2].set_xlabel("Longitude [deg]", fontsize = 12)
ax[2].set_xlim([32.65,  32.77])
ax[2].set_ylim([-27.55, -27.41])
ax[2].set_yticks([])
ax[2].set_title(ts_03)
quiv = ax[2].quiver(grid_x_red[::7,::5], grid_y_red[::7,::5], xVelBottomGrid_03[::7,::5], yVelBottomGrid_03[::7,::5], color = 'k', scale = 2)

im1 = ax[3].imshow(tempBottomGrid_04, vmin = 22, vmax = 26, cmap='jet', extent=[np.min(x_reduced), np.max(x_reduced), np.min(y_reduced), np.max(y_reduced)]);
im2 = ax[3].imshow(dst,extent=[32.647,  32.767, -27.559, -27.405], zorder=2)
ax[3].set_xlabel("Longitude [deg]", fontsize = 12)
ax[3].set_xlim([32.65,  32.77])
ax[3].set_ylim([-27.55, -27.41])
ax[3].set_yticks([])
ax[3].set_title(ts_04)
quiv = ax[3].quiver(grid_x_red[::7,::5], grid_y_red[::7,::5], xVelBottomGrid_04[::7,::5], yVelBottomGrid_04[::7,::5], color = 'k', scale = 2)

#plt.tight_layout()

cbar0 = fig.colorbar(im1, ax=ax, orientation='horizontal', pad = 0.2)
cbar0.ax.tick_params(labelsize=9)
cbar0.ax.set_xlabel('Temperature [deg C]', fontsize = 11)
