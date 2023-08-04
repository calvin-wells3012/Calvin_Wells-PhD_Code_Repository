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
import csv
from numpy import savetxt


data = pd.read_csv('*.csv')
lat = np.array(data['Lat'])
lon = np.array(data['Lon'])
seasons = data['Season']

spring_lat = np.zeros(len(seasons))
summer_lat = np.zeros(len(seasons))
autumn_lat = np.zeros(len(seasons))
winter_lat = np.zeros(len(seasons))

spring_lon = np.zeros(len(seasons))
summer_lon = np.zeros(len(seasons))
autumn_lon = np.zeros(len(seasons))
winter_lon = np.zeros(len(seasons))

i = 0
while i < len(seasons):
    if seasons[i] == 'Summer':
        summer_lat[i] = lat[i]
    else:
        summer_lat[i] = np.nan

    if seasons[i] == 'Spring':
        spring_lat[i] = lat[i]
    else:
        spring_lat[i] = np.nan
        
        
    if seasons[i] == 'Winter':
        winter_lat[i] = lat[i]
    else:
        winter_lat[i] = np.nan
        
    if seasons[i] == 'Autumn':
        autumn_lat[i] = lat[i]
    else:
        autumn_lat[i] = np.nan
    i+=1
    
j = 0
while j < len(seasons):
    if seasons[j] == 'Summer':
        summer_lon[j] = lon[j]
    else:
        summer_lon[j] = np.nan

    if seasons[j] == 'Spring':
        spring_lon[j] = lon[j]
    else:
        spring_lon[j] = np.nan
        
        
    if seasons[j] == 'Winter':
        winter_lon[j] = lon[j]
    else:
        winter_lon[j] = np.nan
        
    if seasons[j] == 'Autumn':
        autumn_lon[j] = lon[j]
    else:
        autumn_lon[j] = np.nan
    j+=1
    
background = np.zeros([65,69])

##Plotting    
fig01 = plt.figure(figsize=(11,7))
m = Basemap(projection='cyl', lat_0 = 0, lon_0 = 0, resolution = 'l', llcrnrlat=-34, urcrnrlat = -18, llcrnrlon = 30, urcrnrlon = 48)
border_color = 'black'    
m.drawcoastlines()    
plt.imshow(background, vmin = -3, vmax = 3, cmap='bwr', extent=[30, 48, -34, -18], interpolation='bilinear');
m.fillcontinents(color='grey')
plt.plot(32.65, -27.53, 'ro', markersize=5, color = 'deeppink')
plt.xlabel('Longitude [deg]')
plt.ylabel('Latitude [deg]')
plt.xticks(np.arange(30, 48, step=1))
plt.yticks(np.arange(-34, -18, step=1))
plt.plot(autumn_lon, autumn_lat, color = 'gold', label = 'Autumn')
plt.plot(winter_lon, winter_lat, color = 'blue', label = 'Winter')
plt.plot(spring_lon, spring_lat, color = 'limegreen', label = 'Spring')
plt.plot(summer_lon, summer_lat, color = 'red', label = 'Summer')
legend_x = 1
legend_y = 0.1
plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
#plt.plot(lon, lat)
plt.tight_layout()
plt.savefig('Cyclonic_Eddy_Tracks_Seasonal.PDF')  





