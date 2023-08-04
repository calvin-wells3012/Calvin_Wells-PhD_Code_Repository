## This program plots the depth averaged (first 30m from surface) current speed spacially and creates a video

import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import fftpack
import math
from scipy import signal

start = 8100
end = 397800

## load data file names   
all_files = glob.glob('../../Raw/01/*.CSV')
all_files.sort()

##Set up variables

x_tilt = np.empty(0)
z_tilt = np.empty(0)
x_mag = np.empty(0)
y_mag = np.empty(0)
z_mag = np.empty(0)
day = np.empty(0)
month = np.empty(0)
year = np.empty(0)
hour = np.empty(0)
minute = np.empty(0)
sec = np.empty(0)

x_tilt_mean = np.empty(0)
z_tilt_mean = np.empty(0)
x_mag_mean = np.empty(0)
y_mag_mean = np.empty(0)
z_mag_mean = np.empty(0)

temp = np.empty(0)

##Loops through files and extracts measured data

for timestep, datafile in enumerate(all_files): 
    data = pd.read_csv(datafile, header = None)
    
    ##Tilt data
    x_tilt_01 = np.array(data[1])
    z_tilt_01 = np.array(data[2])
    x_tilt = np.append(x_tilt, x_tilt_01)
    z_tilt = np.append(z_tilt, z_tilt_01)
    
    ##Temp data
    temp_01 = np.array(data[6])
    temp = np.append(temp, temp_01)
    
    ##Magnetometer data
    x_mag_01 = np.array(data[3])
    y_mag_01 = np.array(data[4])
    z_mag_01 = np.array(data[5])
    x_mag_01_mean = np.mean(x_mag_01)
    y_mag_01_mean = np.mean(y_mag_01)
    z_mag_01_mean = np.mean(z_mag_01)
    x_mag = np.append(x_mag, x_mag_01)
    y_mag = np.append(y_mag, y_mag_01)
    z_mag = np.append(z_mag, z_mag_01)
    x_mag_mean = np.append(x_mag_mean, x_mag_01_mean)
    y_mag_mean = np.append(y_mag_mean, y_mag_01_mean)
    z_mag_mean = np.append(z_mag_mean, z_mag_01_mean)
    
    ##Date and time data
    date_time_01 = list(data[0])
    day_01 = np.zeros(len(date_time_01))
    month_01 = np.zeros(len(date_time_01))
    year_01 = np.zeros(len(date_time_01))
    hour_01 = np.zeros(len(date_time_01))
    min_01 = np.zeros(len(date_time_01))
    sec_01 = np.zeros(len(date_time_01))
    
    ##Error with was date reads so had to do additional processing
    j = 0
    while j < len(date_time_01):
        date_time_split =  date_time_01[j].split('/')
        day_02 = int(date_time_split[0])
        month_02 = int(date_time_split[1])
        date_time_split =  date_time_split[2].split(' ')
        year_02 = int(date_time_split[0])
        time_split = date_time_split[1].split(':')
        hour_02 = int(time_split[0])
        min_02 = int(time_split[1])
        sec_02 = int(time_split[2])
        
        day_01[j] = day_02
        month_01[j] = month_02
        year_01[j] = year_02  
        hour_01[j] = hour_02
        min_01[j] = min_02
        sec_01[j] = sec_02 
        
        j+=1
        
    day = np.append(day, day_01)
    month = np.append(month, month_01)
    year = np.append(year, year_01)
    hour = np.append(hour, hour_01)
    minute = np.append(minute, min_01)
    sec = np.append(sec, sec_01)


date_time_00 = [None] * len(day)
#
n = 0
while n < len(day):
    date_time_00[n] = datetime.datetime(int(year[n]), int(month[n]), int(day[n]), int(hour[n]), int(minute[n]), int(sec[n]))
#    if month[n] < month[n-1]:
#        month[n] = month[n-1]
#        date_time_00[n] = date_time_00[n-300] + datetime.timedelta(hours=1)
    n+=1

x_tilt = x_tilt[start:end]
z_tilt = z_tilt[start:end]
x_mag = x_mag[start:end]
y_mag = y_mag[start:end]
z_mag = z_mag[start:end]
temp = temp[start:end]
    
##Converting tilt data to velocity and current speed and direction
x_vel = np.zeros(len(x_tilt))
z_vel = np.zeros(len(x_tilt))

x_tilt_rad = np.radians(x_tilt)
z_tilt_rad = np.radians(z_tilt)
x_dist = np.sin(x_tilt_rad)
z_dist = np.sin(z_tilt_rad)
dist = np.sqrt(x_dist**2 + z_dist**2)
tilt = np.degrees(np.arcsin(dist))
  
##Magnetometer calibration and outlier filtering

k = 0
while k < len(x_mag):
    if x_mag[k] > -10:
        x_mag[k] = x_mag[k-1]     
    if x_mag[k] < -70:
        x_mag[k] = x_mag[k-1]
        
    if y_mag[k] > -55:
        y_mag[k] = y_mag[k-1]     
    if y_mag[k] < -100:
        y_mag[k] = y_mag[k-1]
        
    if z_mag[k] > 45:
        z_mag[k] = z_mag[k-1]     
    if z_mag[k] < -30:
        z_mag[k] = z_mag[k-1]
    k+=1

def sphereFit(spX,spY,spZ):
    #   Assemble the A matrix
    spX = np.array(spX)
    spY = np.array(spY)
    spZ = np.array(spZ)
    A = np.zeros((len(spX),4))
    A[:,0] = spX*2
    A[:,1] = spY*2
    A[:,2] = spZ*2
    A[:,3] = 1

    #   Assemble the f matrix
    f = np.zeros((len(spX),1))
    f[:,0] = (spX*spX) + (spY*spY) + (spZ*spZ)
    C, residules, rank, singval = np.linalg.lstsq(A,f)

    #   solve for the radius
    t = (C[0]*C[0])+(C[1]*C[1])+(C[2]*C[2])+C[3]
    radius = np.sqrt(t)

    return radius, C[0], C[1], C[2]

r, x0, y0, z0 = sphereFit(x_mag,y_mag,z_mag)
x_mag_cal = x_mag - x0
y_mag_cal = y_mag - y0
z_mag_cal = z_mag - z0

x_tilt_rad = np.radians(x_tilt)
z_tilt_rad = np.radians(z_tilt)

roll = z_tilt_rad * -1
pitch = x_tilt_rad * -1

xh =  x_mag_cal*np.cos(pitch) + z_mag_cal*np.sin(roll)*np.sin(pitch) - y_mag_cal*np.cos(roll)*np.sin(pitch)
zh =  z_mag_cal*np.cos(roll) + y_mag_cal*np.sin(roll)

angle_comp = np.degrees(np.arctan2(zh, xh))
angle_comp = angle_comp + 90 

current_dir_tilt = np.zeros(len(x_tilt))
j = 0
while j < len(x_tilt):
    current_dir_tilt[j] = np.degrees(math.atan2(x_dist[j], -1 * z_dist[j]))
    if current_dir_tilt[j] < 0:
        current_dir_tilt[j] = current_dir_tilt[j]+360
    j+=1

current_dir_tilt_cor = current_dir_tilt + angle_comp

m = 0
while m < len(current_dir_tilt_cor):
    if current_dir_tilt_cor[m] < 0:
        current_dir_tilt_cor[m] = current_dir_tilt_cor[m]+360
    if current_dir_tilt_cor[m] > 360:
        current_dir_tilt_cor[m] = current_dir_tilt_cor[m]-360
    m+=1

z_dist_cor = dist*np.cos(np.radians(current_dir_tilt_cor)) * -1
x_dist_cor = dist*np.sin(np.radians(current_dir_tilt_cor))
dist_cor = np.sqrt(x_dist_cor**2 + z_dist_cor**2)
tilt_cor = np.degrees(np.arcsin(dist_cor))

x_tilt_cor = np.degrees(np.arcsin(x_dist_cor))
z_tilt_cor = np.degrees(np.arcsin(z_dist_cor))

## Filter out any signals shorter than 4 seconds (noise and eddy shedding)
b, a = signal.butter(4, 0.25, 'low', fs=1)
x_tilt_filtered = signal.filtfilt(b, a, x_tilt_cor)
z_tilt_filtered = signal.filtfilt(b, a, z_tilt_cor)

## Butterwrth filter to remove wave signal (signals shorter than 60s)
d, c = signal.butter(4, 0.0167, 'low', fs=1)
x_tilt_smooth = signal.filtfilt(d, c, x_tilt_filtered)
z_tilt_smooth = signal.filtfilt(d, c, z_tilt_filtered)

x_dist_filtered = np.sin(np.radians(x_tilt_filtered))
z_dist_filtered = np.sin(np.radians(z_tilt_filtered))
dist_filtered = np.sqrt(x_dist_filtered**2 + z_dist_filtered**2)
tilt_filtered = np.degrees(np.arcsin(dist_filtered))

x_dist_smooth = np.sin(np.radians(x_tilt_smooth))
z_dist_smooth = np.sin(np.radians(z_tilt_smooth))
dist_smooth = np.sqrt(x_dist_smooth**2 + z_dist_smooth**2)
tilt_smooth = np.degrees(np.arcsin(dist_smooth))

current_speed = np.zeros(len(tilt_filtered))
i = 0
while i < len(tilt_cor):
    current_speed[i] = 0.0000038592*abs(tilt_filtered[i])**3 - 0.0003105897*abs(tilt_filtered[i])**2 + 0.0149545288*abs(tilt_filtered[i]) + 0.0010922498 
    i+=1

current_dir = np.zeros(len(tilt_filtered))
j = 0
while j < len(x_tilt):
    current_dir[j] = np.degrees(math.atan2(x_tilt_filtered[j], -1 * z_tilt_filtered[j]))
    if current_dir[j] < 0:
        current_dir[j] = current_dir[j]+360
    j+=1

current_speed_smooth = np.zeros(len(tilt_filtered))
i = 0
while i < len(tilt_cor):
    current_speed_smooth[i] = 0.0000038592*abs(tilt_smooth[i])**3 - 0.0003105897*abs(tilt_smooth[i])**2 + 0.0149545288*abs(tilt_smooth[i]) + 0.0010922498 
    i+=1

current_dir_smooth = np.zeros(len(tilt_filtered))
j = 0
while j < len(x_tilt):
    current_dir_smooth[j] = np.degrees(math.atan2(x_tilt_smooth[j], -1 * z_tilt_smooth[j]))
    if current_dir_smooth[j] < 0:
        current_dir_smooth[j] = current_dir_smooth[j]+360
    j+=1

x_vel = current_speed*np.sin(np.radians(current_dir))    
z_vel = current_speed*np.cos(np.radians(current_dir)) * -1

x_vel_smooth = current_speed_smooth*np.sin(np.radians(current_dir_smooth))    
z_vel_smooth = current_speed_smooth*np.cos(np.radians(current_dir_smooth)) * -1

date_time_00 = date_time_00[start:end]
x_vel_mean = np.mean(x_vel[:].reshape(-1, 300), axis=1)
z_vel_mean = np.mean(z_vel[:].reshape(-1, 300), axis=1)
temp_mean = np.mean(temp[:].reshape(-1, 300), axis=1)
current_speed_mean = np.sqrt(x_vel_mean**2+z_vel_mean**2)

x_vel_smooth_mean = np.mean(x_vel_smooth[:].reshape(-1, 300), axis=1)
z_vel_smooth_mean = np.mean(z_vel_smooth[:].reshape(-1, 300), axis=1)
current_speed_smooth_mean = np.sqrt(x_vel_smooth_mean**2+z_vel_smooth_mean**2)
#date_time_hr = date_time_00[::300]

date_time_hr = [None] * len(temp_mean)
start_time = date_time_00[0]
k = 0
while k < len(date_time_hr):
    date_time_hr[k] = start_time + datetime.timedelta(hours=k)
    k+=1

current_dir_mean = np.zeros(len(z_vel_mean))
k = 0
while k < len(z_vel_mean):
    current_dir_mean[k] = np.degrees(math.atan2(x_vel_mean[k], -1 * z_vel_mean[k]))
    if current_dir_mean[k] < 0:
        current_dir_mean[k] = current_dir_mean[k]+360
    k+=1
    
current_dir_smooth_mean = np.zeros(len(z_vel_mean))
k = 0
while k < len(z_vel_smooth_mean):
    current_dir_smooth_mean[k] = np.degrees(math.atan2(x_vel_smooth_mean[k], -1 * z_vel_smooth_mean[k]))
    if current_dir_smooth_mean[k] < 0:
        current_dir_smooth_mean[k] = current_dir_smooth_mean[k]+360
    k+=1

combined_timeseries_mean = pd.DataFrame({"Date Time" : date_time_hr, "x_vel [m/s]" : x_vel_smooth_mean, "z_vel [m/s]" : z_vel_smooth_mean, "Current Speed [m/s]" : current_speed_smooth_mean, "Current Direction [deg]" : current_dir_smooth_mean, "Temperature [deg]" : temp_mean})
combined_timeseries_mean.to_csv('TCM_Timeseries_Averaged.csv')

    
    
    
    
    