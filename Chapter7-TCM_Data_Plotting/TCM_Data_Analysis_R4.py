import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime as dtime
import matplotlib
import matplotlib.dates as mdates
from matplotlib.dates import date2num
import datetime

data_01_bikini = pd.read_csv("../../2021-10-25/TCM/Processed/10_Bikini/TCM_Timeseries_Averaged.csv")
data_02_bikini = pd.read_csv("../../2022-02-06/TCM/Processed/07/TCM_Timeseries_Averaged.csv")
data_03_bikini = pd.read_csv("../../2022-06-04/TCM/Processed/10_Bikini/TCM_Timeseries_Averaged.csv")
data_bikini = data_01_bikini.append(data_02_bikini)
data_bikini = data_bikini.append(data_03_bikini)

data_01_antons = pd.read_csv("../../2021-10-25/TCM/Processed/01_Antons/TCM_Timeseries_Averaged.csv")
data_02_antons = pd.read_csv("../../2022-02-06/TCM/Processed/11_Antons/TCM_Timeseries_Averaged.csv")
data_03_antons = pd.read_csv("../../2022-06-04/TCM/Processed/12_Antons/TCM_Timeseries_Averaged.csv")
data_antons = data_01_antons.append(data_02_antons)
data_antons = data_antons.append(data_03_antons)

time_01 = data_01_antons['Date Time']
time_02 = data_02_antons['Date Time']
time_03 = data_03_antons['Date Time']
time = time_01.append(time_02)
time = time.append(time_03)

temp_01 = np.array(data_01_antons['Temperature [deg]'])
temp_01[-1] = np.nan
temp_02 = np.array(data_02_antons['Temperature [deg]'])
temp_02[-1] = np.nan
temp_03 = np.array(data_03_antons['Temperature [deg]'])
temp_03[:12] = np.nan
temp = np.append(temp_01, temp_02)
temp = np.append(temp, temp_03)

# time_01 = data_01_bikini['Date Time']
# time_01 = time_01[9:-14]
# temp_01_bikini = np.array(data_01_bikini['Temperature [deg]'])
# temp_01_bikini = temp_01_bikini[9:-14]
# u_vel_01_bikini = np.array(data_01_bikini['x_vel [m/s]'])
# u_vel_01_bikini = u_vel_01_bikini[9:-14]
# v_vel_01_bikini = np.array(data_01_bikini['z_vel [m/s]'])
# v_vel_01_bikini = v_vel_01_bikini[9:-14]
# date_01_dly = time_01.iloc[::24,]
# temp_01_bikini_dly = np.mean(temp_01_bikini.reshape(-1, 24), axis=1)
# u_vel_01_bikini_dly = np.mean(u_vel_01_bikini.reshape(-1, 24), axis=1)
# v_vel_01_bikini_dly = np.mean(v_vel_01_bikini.reshape(-1, 24), axis=1)

# time_02 = data_02_bikini['Date Time']
# time_02 = time_02[15:-10]
# temp_02_bikini = np.array(data_02_bikini['Temperature [deg]'])
# temp_02_bikini = temp_02_bikini[15:-10]
# u_vel_02_bikini = np.array(data_02_bikini['x_vel [m/s]'])
# u_vel_02_bikini = u_vel_02_bikini[15:-10]
# v_vel_02_bikini = np.array(data_02_bikini['z_vel [m/s]'])
# v_vel_02_bikini = v_vel_02_bikini[15:-10]
# date_02_dly = time_02.iloc[::24,]
# temp_02_bikini_dly = np.mean(temp_02_bikini.reshape(-1, 24), axis=1)
# u_vel_02_bikini_dly = np.mean(u_vel_02_bikini.reshape(-1, 24), axis=1)
# v_vel_02_bikini_dly = np.mean(v_vel_02_bikini.reshape(-1, 24), axis=1)

# time_03 = data_03_bikini['Date Time']
# time_03 = time_03[5:-6]
# temp_03_bikini = np.array(data_03_bikini['Temperature [deg]'])
# temp_03_bikini = temp_03_bikini[5:-6]
# u_vel_03_bikini = np.array(data_03_bikini['x_vel [m/s]'])
# u_vel_03_bikini = u_vel_03_bikini[5:-6]
# v_vel_03_bikini = np.array(data_03_bikini['z_vel [m/s]'])
# v_vel_03_bikini = v_vel_03_bikini[5:-6]
# date_03_dly = time_03.iloc[::24,]
# temp_03_bikini_dly = np.mean(temp_03_bikini.reshape(-1, 24), axis=1)
# u_vel_03_bikini_dly = np.mean(u_vel_03_bikini.reshape(-1, 24), axis=1)
# v_vel_03_bikini_dly = np.mean(v_vel_03_bikini.reshape(-1, 24), axis=1)

# date_dly = date_01_dly.append(date_02_dly)
# date_dly = date_dly.append(date_03_dly)

# temp_bikini_dly = np.append(temp_01_bikini_dly, temp_02_bikini_dly)
# temp_bikini_dly = np.append(temp_bikini_dly, temp_03_bikini_dly)

# u_vel_bikini_dly = np.append(u_vel_01_bikini_dly, u_vel_02_bikini_dly)
# u_vel_bikini_dly = np.append(u_vel_bikini_dly, u_vel_03_bikini_dly)

# v_vel_bikini_dly = np.append(v_vel_01_bikini_dly, v_vel_02_bikini_dly)
# v_vel_bikini_dly = np.append(v_vel_bikini_dly, v_vel_03_bikini_dly)

#%%

fig, ax = plt.subplots(4,figsize = (10,12))
ax[0].plot(pd.to_datetime(np.array(data_antons['Date Time'])),temp, c = 'k', linewidth = 0.7)
ax[0].set_xlim([datetime.date(2021, 8, 15), datetime.date(2022, 4, 23)])
ax[0].set_ylim([18, 28])
ax[0].set_ylabel('Temperature [$^\circ$C]', fontsize = 12)
ax[0].xaxis_date()
formatter = mdates.DateFormatter("%D")
ax[0].xaxis.set_major_formatter(formatter)
locator_major = mdates.MonthLocator(interval=1)
ax[0].xaxis.set_major_locator(locator_major)
ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
locator_minor = mdates.MonthLocator(interval=1)
ax[0].xaxis.set_minor_locator(locator_minor)
ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0].grid(which='minor',axis='x', linestyle='--', dashes=(5, 5))
ax[0].grid(which='major',axis='x', linestyle='--', dashes=(5, 5))

qiv = ax[1].quiver(date2num(pd.to_datetime(data_bikini['Date Time'])), [[0]*len(data_bikini)], data_bikini['x_vel [m/s]'],  -1*data_bikini['z_vel [m/s]'],  data_bikini['Current Speed [m/s]'], cmap=plt.cm.binary, clim=[0,0.2], scale=3, width=0.0015, headwidth=1, headlength=0)
qk = plt.quiverkey(qiv, 0.05, 0.85, 0.1, r'$0.1 m/s$', labelpos ='N',coordinates='axes',fontproperties={'weight': 'bold','size'   : 12}, labelsep = 0.1)
ax[1].set_xlim([datetime.date(2021, 8, 15), datetime.date(2022, 4, 23)])
#ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[1].axes.get_yaxis().set_visible(False)
plt.tight_layout()
ax[1].xaxis_date()
formatter = mdates.DateFormatter("%D")
ax[1].xaxis.set_major_formatter(formatter)
locator_major = mdates.MonthLocator(interval=1)
ax[1].xaxis.set_major_locator(locator_major)
ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
locator_minor = mdates.MonthLocator(interval=1)
ax[1].xaxis.set_minor_locator(locator_minor)
ax[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[1].grid(which='minor',axis='x', linestyle='--', dashes=(5, 5))
ax[1].grid(which='major',axis='x', linestyle='--', dashes=(5, 5))

ax[2].plot(pd.to_datetime(np.array(data_antons['Date Time'])),np.array(data_antons['Temperature [deg]']), c = 'k', linewidth = 0.8)
ax[2].set_xlim([datetime.date(2021, 12, 9), datetime.date(2022, 1, 9)])
ax[2].set_ylim([19, 27])
ax[2].set_ylabel('Temperature [$^\circ$C]', fontsize = 12)
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

qiv = ax[3].quiver(date2num(pd.to_datetime(data_bikini['Date Time'])), [[0]*len(data_bikini)], data_bikini['x_vel [m/s]'],  -1*data_bikini['z_vel [m/s]'],  data_bikini['Current Speed [m/s]'], cmap=plt.cm.binary, clim=[0,0.2], scale=2, width=0.0015, headwidth=1, headlength=0)
qk = plt.quiverkey(qiv, 0.05, 0.85, 0.1, r'$0.1 m/s$', labelpos ='N',coordinates='axes',fontproperties={'weight': 'bold','size'   : 12}, labelsep = 0.1)
ax[3].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[3].axes.get_yaxis().set_visible(False)
plt.tight_layout()
ax[3].set_xlim([datetime.date(2021, 12, 9), datetime.date(2022, 1, 9)])
ax[3].xaxis_date()
formatter = mdates.DateFormatter("%D")
ax[3].xaxis.set_major_formatter(formatter)
locator_major = mdates.DayLocator(interval=7)
ax[3].xaxis.set_major_locator(locator_major)
ax[3].xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
locator_minor = mdates.DayLocator(interval=1)
ax[3].xaxis.set_minor_locator(locator_minor)
ax[3].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[3].grid(which='minor',axis='x', linestyle='--', dashes=(5, 5))
ax[3].grid(which='major',axis='x', linestyle='--', dashes=(5, 5))

plt.savefig('TCM_Measurements.pdf')

# fig, ax = plt.subplots(2,figsize = (9,6))
# ax[0].plot(pd.to_datetime(np.array(data_bikini['Date Time'])),np.array(data_bikini['Temperature [deg]']), c = 'k', linewidth = 0.8)
# ax[0].set_xlim([datetime.date(2022, 3, 18), datetime.date(2022, 4, 15)])
# ax[0].set_ylim([18, 28])
# ax[0].set_ylabel('Temperature [deg C]', fontsize = 12)
# #ax[0].axes.get_xaxis().set_visible(False)

# qiv = ax[1].quiver(date2num(pd.to_datetime(data_bikini['Date Time'])), [[0]*len(data_bikini)], data_bikini['x_vel [m/s]'],  -1*data_bikini['z_vel [m/s]'],  data_bikini['Current Speed [m/s]'], cmap=plt.cm.binary, clim=[0,0.2], scale=2, width=0.0015, headwidth=1, headlength=0)
# qk = plt.quiverkey(qiv, 0.05, 0.8, 0.1, r'$0.1 m/s$', labelpos ='N',coordinates='axes',fontproperties={'weight': 'bold','size'   : 12}, labelsep = 0.1)
# ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
# ax[1].axes.get_yaxis().set_visible(False)
# plt.tight_layout()
# ax[1].set_xlim([datetime.date(2022, 3, 18), datetime.date(2022, 4, 15)])














