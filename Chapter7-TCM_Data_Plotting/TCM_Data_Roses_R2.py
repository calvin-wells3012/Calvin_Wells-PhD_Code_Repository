import numpy as np
import matplotlib.pyplot as plt
from windrose import WindroseAxes
import pandas as pd
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)

data_TMR = pd.read_csv('TMR_Timeseries.csv')
direction_TMR = np.array(data_TMR['Current Direction [deg]'])
speed_TMR = np.array(data_TMR['Current Speed [m/s]'])
direction1 = np.nan_to_num(direction_TMR, nan=0)
speed1 = np.nan_to_num(speed_TMR, nan=0)

data_bikini = pd.read_csv('Bikini_Timeseries.csv')
direction_bikini = np.array(data_bikini['Current Direction [deg]'][:4200])
speed_bikini = np.array(data_bikini['Current Speed [m/s]'][:4200])
direction2 = np.nan_to_num(direction_bikini, nan=0)
speed2 = np.nan_to_num(speed_bikini, nan=0)

data_rehan = pd.read_csv('Rehan_Timeseries.csv')
direction_rehan = np.array(data_rehan['Current Direction [deg]'])
speed_rehan = np.array(data_rehan['Current Speed [m/s]'])
direction3 = np.nan_to_num(direction_rehan, nan=0)
speed3 = np.nan_to_num(speed_rehan, nan=0)

# Create a WindroseAxes object
# ax = WindroseAxes.from_ax()
# cmap = plt.cm.jet
# ax.bar(direction1, speed1, normed=True, opening=0.8, edgecolor='white', cmap=cmap)
# ax.legend(title='Wind Speed (m/s)')


# Create a figure with three subplots
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(12, 6.5), subplot_kw=dict(projection='windrose'))

bins_range = np.arange(0,0.2,0.02)

# Create a WindroseAxes object for each subplot
ax1 = WindroseAxes.from_ax(ax=axs[0])
ax2 = WindroseAxes.from_ax(ax=axs[1])
ax3 = WindroseAxes.from_ax(ax=axs[2])

# Plot the wind data for each subplot
cmap = plt.cm.jet
ax1.bar(direction1, speed1, normed=True, opening=0.8, edgecolor='white',bins=bins_range, cmap=cmap)
ax2.bar(direction2, speed2, normed=True, opening=0.8, edgecolor='white',bins=bins_range, cmap=cmap)
ax3.bar(direction3, speed3, normed=True, opening=0.8, edgecolor='white',bins=bins_range, cmap=cmap)

ax1.set_yticks(np.arange(0, 45, step=5))
ax1.set_yticklabels(np.arange(0, 45, step=5))
ax2.set_yticks(np.arange(0, 45, step=5))
ax2.set_yticklabels(np.arange(0, 45, step=5))
ax3.set_yticks(np.arange(0, 45, step=5))
ax3.set_yticklabels(np.arange(0, 45, step=5))

# Set titles for each subplot
ax1.set_title('Two-Mile Reef high rugosity')
ax2.set_title('Two-Mile Reef low rugosity')
ax3.set_title('Two-Mile Reef deep')

im = plt.imread('Legend.PNG') # insert local path of the image.
newax = fig.add_axes([0.2,0.01,0.55,0.3], anchor='S', zorder=1)
plt.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
newax.imshow(im)


# Adjust the layout and spacing
fig.tight_layout()

# Show the plot
plt.show()

plt.savefig('Current_Roses_R0.pdf')