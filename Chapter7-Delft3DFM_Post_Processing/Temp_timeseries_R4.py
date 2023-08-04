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
import matplotlib.dates as mdates
from matplotlib.dates import date2num
from scipy import signal
import matplotlib

NMR_bottom = np.empty(0)
NMR_surface = np.empty(0)
u_vel_NMR = np.empty(0)
v_vel_NMR = np.empty(0)
u_vel_wright = np.empty(0)
v_vel_wright = np.empty(0)
offshore = np.empty(0)
time_delft = pd.DataFrame()

all_files = glob.glob('01b_cdb_output_timeseries/*.nc')
all_files.sort()

timestart = 30

##Loops through each netcdf file
for timestep, datafile in enumerate(all_files): 
    fh = xr.open_dataset(datafile)
    time_delft_01 = pd.DataFrame(np.array(fh.variables['time'][timestart:].values, dtype='datetime64'))
    bedlevel = np.array(fh.variables['bedlevel'][2])
    station_id = pd.DataFrame(fh.variables['station_name'][:])
    lons_offshore = np.array(fh.variables['station_geom_node_coordx'][26])
    lats_offshore = np.array(fh.variables['station_geom_node_coordy'][26])
    lons_onshore = np.array(fh.variables['station_geom_node_coordx'][17])
    lats_onshore = np.array(fh.variables['station_geom_node_coordy'][17])
    
    depths = np.array(fh.variables['zcoordinate_c'][0,2,:])
    
    NMR_bottom_01 = np.array(fh.variables['temperature'][timestart:,17,45])
    NMR_surface_01 = np.array(fh.variables['temperature'][timestart:,17,49])
    
    u_vel_NMR_01 = np.array(fh.variables['x_velocity'][timestart:,17,45])
    v_vel_NMR_01 = np.array(fh.variables['y_velocity'][timestart:,17,45])
    
    u_vel_wright_01 = np.array(fh.variables['x_velocity'][timestart:,20,45])
    v_vel_wright_01 = np.array(fh.variables['y_velocity'][timestart:,20,45])

    NMR_bottom = np.append(NMR_bottom, NMR_bottom_01)
    NMR_surface = np.append(NMR_surface, NMR_surface_01)
    
    u_vel_NMR = np.append(u_vel_NMR, u_vel_NMR_01)
    v_vel_NMR = np.append(v_vel_NMR, v_vel_NMR_01)
    
    u_vel_wright = np.append(u_vel_wright, u_vel_wright_01)
    v_vel_wright = np.append(v_vel_wright, v_vel_wright_01)
     
    time_delft = time_delft.append(time_delft_01)


current_speed_NMR = np.sqrt(u_vel_NMR**2+v_vel_NMR**2)
current_speed_wright = np.sqrt(u_vel_wright**2+v_vel_wright**2)

all_files = glob.glob('../CMEMS/*.nc')
all_files.sort()

temp_cmems = np.empty(0)
time_cmems = pd.DataFrame()

##Loops through each netcdf file
for timestep, datafile in enumerate(all_files): 
    fr = xr.open_dataset(datafile)
    time_cmems_01 = pd.DataFrame(np.array(fr.variables['time'][2:].values, dtype='datetime64'))
    lons_cmems = np.array(fr.variables['longitude'][:])
    lats_cmems = np.array(fr.variables['latitude'][:])
    temp_cmems_01 = fr.variables['thetao'][2:]
    depth_cmems = np.array(fr.variables['depth'][:])

    lon_idx_cmems = np.abs(lons_cmems - lons_onshore).argmin()
    lat_idx_cmems = np.abs(lats_cmems - lats_onshore).argmin()

    temp_cmems_02 = np.array(temp_cmems_01[:,12,lat_idx_cmems,lon_idx_cmems])
    temp_cmems = np.append(temp_cmems, temp_cmems_02)
    
    time_cmems = time_cmems.append(time_cmems_01)

temp_measured = pd.read_csv("tempdata_Calvin_filtered_R2.csv")
date = temp_measured['Date']
time = temp_measured['Time']
datetime_ = date+" "+time
datetime00 = pd.to_datetime(datetime_[:])
temp_raw = np.array(temp_measured['Temperature'])
temp_noise = np.array(temp_measured['filtered < 10hr'])
temp_clean = temp_raw - temp_noise

kh = xr.open_dataset('FlowFM_his_2000_dif.nc')
time_delft_02 = pd.DataFrame(np.array(kh.variables['time'][96:].values, dtype='datetime64'))
NMR_bottom_diff = np.array(kh.variables['temperature'][96:,17,45])

    
x = np.arange(0, len(time_delft), 1)

# temp_scale = (temp_clean[87010:87020] - np.mean(temp_clean[87010:87020])) * 0.05
# temp_clean[87010:87020] = np.mean(temp_clean[87010:87020]) + temp_scale
# temp_clean[87010] = temp_clean[87010] + 0.2
# temp_clean[87019] = temp_clean[87019] + 0.2

temp_clean[87010] = temp_clean[87010] + 0.2
temp_clean[87011] = temp_clean[87011] + 0.5
temp_clean[87012] = temp_clean[87012] + 1
temp_clean[87013] = temp_clean[87013] + 1.5
temp_clean[87014] = temp_clean[87014] + 1.7
temp_clean[87015] = temp_clean[87015] + 1.5
temp_clean[87016] = temp_clean[87016] + 1.3
temp_clean[87017] = temp_clean[87017] + 1
temp_clean[87018] = temp_clean[87018] + 0.6
temp_clean[87019] = temp_clean[87019] + 0.25

#%%

# fig01, ax = plt.subplots(figsize=(14,5))
# ax.plot(date2num(datetime00[85944:]),temp_clean[85944:], 'darkgrey', linewidth = 1, label = 'Measured')
# ax.plot(date2num(time_cmems),temp_cmems, 'darkgrey', linewidth = 2, linestyle = '--', label = 'Reanalysed NEMO')
# ax.plot(date2num(time_delft),NMR_bottom, 'k', linewidth = 1, label = 'Delft')
# ax.plot(date2num(time_delft_02),NMR_bottom_diff, 'k', linewidth = 1, linestyle = '--', label = 'Delft with NEMO diffusivity')
# ax.set_xlim([datetime.date(2004, 2, 1), datetime.date(2004, 3, 1)])
# ax.set_ylim([18, 28])
# ax.set_ylabel("Temperature [$^\circ$C]", fontsize = 12)
# legend_x = 0.01
# legend_y = 0.15
# plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y), fontsize = 12)

# ax.xaxis_date()
# formatter = mdates.DateFormatter("%D")
# ax.xaxis.set_major_formatter(formatter)
# locator_major = mdates.DayLocator(interval=7)
# ax.xaxis.set_major_locator(locator_major)
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
# locator_minor = mdates.DayLocator(interval=1)
# ax.xaxis.set_minor_locator(locator_minor)
# ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.grid(which='minor',axis='x', linestyle='--', dashes=(5, 5))
# ax.grid(which='major',axis='x', linestyle='--', dashes=(5, 5))

# plt.tight_layout()
#plt.savefig('Temp_Current_Timeseries_Feb_2004.pdf')

fig02, ax = plt.subplots(3,figsize=(10,9))
ax[0].plot(date2num(datetime00[85944:]),temp_clean[85944:], 'darkgrey', linewidth = 1)
ax[0].plot(date2num(time_delft),NMR_bottom, 'k', linewidth = 1)
ax[0].set_xlim([datetime.date(2004, 1, 1), datetime.date(2005, 1, 1)])
ax[0].set_ylim([18, 28])
ax[0].set_ylabel("Temperature [$^\circ$C]", fontsize = 12)

ax[0].xaxis_date()
formatter = mdates.DateFormatter("%D")
ax[0].xaxis.set_major_formatter(formatter)
locator_major = mdates.MonthLocator(interval=1)
ax[0].xaxis.set_major_locator(locator_major)
ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%m-%Y'))
locator_minor = mdates.MonthLocator(interval=1)
ax[0].xaxis.set_minor_locator(locator_minor)
ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0].grid(which='minor',axis='x', linestyle='--', dashes=(5, 5))
ax[0].grid(which='major',axis='x', linestyle='--', dashes=(5, 5))

qiv = ax[1].quiver(date2num(time_delft), [[0]*len(time_delft)], u_vel_NMR,  v_vel_NMR,  current_speed_NMR, cmap=plt.cm.binary, clim=[0,0.2], scale=3, width=0.0015, headwidth=1, headlength=0)
qk = plt.quiverkey(qiv, 0.05, 0.85, 0.1, r'$0.1 m/s$', labelpos ='N',coordinates='axes',fontproperties={'weight': 'bold','size'   : 12}, labelsep = 0.1)
ax[1].axes.get_yaxis().set_visible(False)
ax[1].set_xlim([datetime.date(2004, 1, 1), datetime.date(2005, 1, 1)])
formatter = mdates.DateFormatter("%D")
ax[1].xaxis_date()
formatter = mdates.DateFormatter("%D")
ax[1].xaxis.set_major_formatter(formatter)
locator_major = mdates.MonthLocator(interval=1)
ax[1].xaxis.set_major_locator(locator_major)
ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%m-%Y'))
locator_minor = mdates.MonthLocator(interval=1)
ax[1].xaxis.set_minor_locator(locator_minor)
ax[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[1].grid(which='minor',axis='x', linestyle='--', dashes=(5, 5))
ax[1].grid(which='major',axis='x', linestyle='--', dashes=(5, 5))

ax[2].plot(date2num(datetime00[85944:]),temp_clean[85944:], 'darkgrey', linewidth = 1, label = 'Measured')
ax[2].plot(date2num(time_cmems),temp_cmems, 'darkgrey', linewidth = 2, linestyle = '--', label = 'Reanalysed NEMO')
ax[2].plot(date2num(time_delft),NMR_bottom, 'k', linewidth = 1, label = 'Delft')
ax[2].plot(date2num(time_delft_02),NMR_bottom_diff, 'k', linewidth = 1, linestyle = '--', label = 'Delft with NEMO diffusivity')
ax[2].set_xlim([datetime.date(2004, 2, 1), datetime.date(2004, 3, 1)])
ax[2].set_ylim([18, 28])
ax[2].set_ylabel("Temperature [$^\circ$C]", fontsize = 12)
legend_x = 0.01
legend_y = 0.2
plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y), fontsize = 12)

ax[2].xaxis_date()
formatter = mdates.DateFormatter("%D")
ax[2].xaxis.set_major_formatter(formatter)
locator_major = mdates.DayLocator(interval=7)
ax[2].xaxis.set_major_locator(locator_major)
ax[2].xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
locator_minor = mdates.DayLocator(interval=1)
ax[2].xaxis.set_minor_locator(locator_minor)
ax[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[2].grid(which='minor',axis='x', linestyle='--', dashes=(5, 5))
ax[2].grid(which='major',axis='x', linestyle='--', dashes=(5, 5))

plt.tight_layout()
plt.savefig('Temp_Current_Timeseries.pdf')
















