import numpy as np
import pandas as pd
import datetime
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

data_filtered = pd.read_csv("*.csv")
anomaly_df = pd.read_csv("*.csv")

##Organising temperature time series data
date = data_filtered['Date']
time = data_filtered['Time']

datetime_ = date+" "+time
datetime00 = pd.to_datetime(datetime_[:])

date_int = np.empty(len(datetime00))
i=0
while i < len(datetime00):
    time_int = int(datetime00[i].strftime('%Y%m%d%H'))
    date_int[i] = time_int
    i+=1

temp_raw = np.array(data_filtered['Temperature'])
temp_noise = np.array(data_filtered['filtered < 10hr'])
temp_clean = temp_raw - temp_noise

date_anomalies = anomaly_df['Date']
time_anomalies = anomaly_df['Time']
datetime_anomalies = date_anomalies+" "+time_anomalies
datetime_anomalies = pd.to_datetime(datetime_anomalies[:])

date_anomalies_int = np.empty(len(datetime_anomalies))
j=0
while j < len(datetime_anomalies):
    time_anomalies_int = int(datetime_anomalies[j].strftime('%Y%m%d%H'))
    date_anomalies_int[j] = time_anomalies_int
    j+=1
 
index_list = np.empty(len(datetime_anomalies))
n=0
while n< len(datetime_anomalies):
    anomaly_step = datetime_anomalies[n]
    anomaly_step = int(anomaly_step.strftime('%Y%m%d%H'))
    index = np.abs(date_int - anomaly_step).argmin()
    index_list[n] = index
    n+=1
    
temp_anomaly = np.empty(len(datetime_anomalies))
k=0
while k < len(datetime_anomalies):
    x = int(index_list[k])
    temp_x = temp_clean[x]
    temp_anomaly[k] = temp_x
    k+=1    

anomaly_local = pd.read_csv("*.csv")

date_anomalies_local = anomaly_local['Date']
time_anomalies_local = anomaly_local['Time']
datetime_anomalies_local = date_anomalies_local+" "+time_anomalies_local
datetime_anomalies_local = pd.to_datetime(datetime_anomalies_local[:])

date_anomalies_int_local = np.empty(len(datetime_anomalies_local))
j=0
while j < len(datetime_anomalies_local):
    time_anomalies_int_local = int(datetime_anomalies_local[j].strftime('%Y%m%d%H'))
    date_anomalies_int_local[j] = time_anomalies_int_local
    j+=1
 
index_list_local = np.empty(len(datetime_anomalies_local))
n=0
while n< len(datetime_anomalies_local):
    anomaly_step_local = datetime_anomalies_local[n]
    anomaly_step_local = int(anomaly_step_local.strftime('%Y%m%d%H'))
    index_local = np.abs(date_int - anomaly_step_local).argmin()
    index_list_local[n] = index_local
    n+=1
    
temp_anomaly_local = np.empty(len(datetime_anomalies_local))
k=0
while k < len(datetime_anomalies_local):
    x = int(index_list_local[k])
    temp_x_local = temp_clean[x]
    temp_anomaly_local[k] = temp_x_local
    k+=1   

###Plotting figures

temp_clean[167495:167550] = 26

fig01, ax = plt.subplots(figsize=(12,5))
ax.plot(datetime00,temp_clean, 'k', linewidth=0.8)
ax.plot(datetime_anomalies,temp_anomaly, 'r', linestyle = ' ', marker = 'o', label = 'Anomalies')
ax.plot(datetime_anomalies_local,temp_anomaly_local, 'b', linestyle = ' ', marker = 'o', label = 'Anomalies')
formatter = mdates.DateFormatter("%Y")
ax.xaxis.set_major_formatter(formatter)
locator = mdates.YearLocator()
ax.xaxis.set_major_locator(locator)
plt.setp( ax.xaxis.get_majorticklabels(), rotation=90 )
plt.ylabel('Temperature [deg C]', fontsize = 14)
plt.xlabel('Year', fontsize = 14)
plt.ylim((17,30))
ax.set_xlim([datetime.date(1994, 1, 1), datetime.date(2015, 1, 1)])
plt.tight_layout()
plt.savefig('NMR_temp_anomaly_Remote_Local.pdf')


