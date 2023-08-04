import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import xarray as xr
from scipy import interpolate

from parcels import FieldSet, ParticleSet, JITParticle, ScipyParticle, AdvectionRK4_3D, Variable, Field,GeographicPolar,Geographic, plotTrajectoriesFile, DiffusionUniformKh, ErrorCode
from datetime import timedelta as delta

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from copy import copy
import cmocean
from datetime import timedelta

file_path = "../../Raw/Anomalies/*.nc"
model = xr.open_dataset(file_path)

# --------- Define meshgrid coordinates to plot velocity field with matplotlib pcolormesh ---------
latmin = 0
latmax = 133
lonmin = 0
lonmax = 121

# Velocity nodes
lon_vals, lat_vals = np.meshgrid(model['longitude'], model['latitude'])
lons_plot = lon_vals[latmin:latmax,lonmin:lonmax]
lats_plot = lat_vals[latmin:latmax,lonmin:lonmax]

dlon = 1/12
dlat = 1/12

# Centers of the gridcells formed by 4 nodes = velocity nodes + 0.5 dx
x = model['longitude'][:-1]+np.diff(model['longitude'])/2
y = model['latitude'][:-1]+np.diff(model['latitude'])/2
lon_centers, lat_centers = np.meshgrid(x, y)

color_land = copy(plt.get_cmap('Reds'))(0)
color_ocean = copy(plt.get_cmap('Reds'))(128)

def make_landmask(fielddata):
    datafile = Dataset(fielddata)

    landmask = datafile.variables['uo'][0, 0]
    landmask = np.ma.masked_invalid(landmask)
    landmask = landmask.mask.astype('int')

    return landmask

landmask = make_landmask(file_path)

# Interpolate the landmask to the cell centers - only cells with 4 neighbouring land points will be land
fl = interpolate.interp2d(model['longitude'],model['latitude'],landmask)

l_centers = fl(lon_centers[0,:],lat_centers[:,0])  

lmask = np.ma.masked_values(l_centers,1) # land when interpolated value == 1

def get_coastal_nodes(landmask):
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap -= 4*landmask
    coastal = np.ma.masked_array(landmask, mask_lap > 0)
    coastal = coastal.mask.astype('int')

    return coastal

def get_shore_nodes(landmask):
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap -= 4*landmask
    shore = np.ma.masked_array(landmask, mask_lap < 0)
    shore = shore.mask.astype('int')

    return shore

def get_coastal_nodes_diagonal(landmask):
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap += np.roll(landmask, (-1,1), axis=(0,1)) + np.roll(landmask, (1, 1), axis=(0,1))
    mask_lap += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (1, -1), axis=(0,1))
    mask_lap -= 8*landmask
    coastal = np.ma.masked_array(landmask, mask_lap > 0)
    coastal = coastal.mask.astype('int')
    
    return coastal
    
def get_shore_nodes_diagonal(landmask):
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap += np.roll(landmask, (-1,1), axis=(0,1)) + np.roll(landmask, (1, 1), axis=(0,1))
    mask_lap += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (1, -1), axis=(0,1))
    mask_lap -= 8*landmask
    shore = np.ma.masked_array(landmask, mask_lap < 0)
    shore = shore.mask.astype('int')

    return shore

coastal = get_coastal_nodes_diagonal(landmask)
shore = get_shore_nodes_diagonal(landmask)


def create_displacement_field(landmask, double_cell=False):
    shore = get_shore_nodes(landmask)
    shore_d = get_shore_nodes_diagonal(landmask) # bordering ocean directly and diagonally
    shore_c = shore_d - shore                    # corner nodes that only border ocean diagonally
    
    Ly = np.roll(landmask, -1, axis=0) - np.roll(landmask, 1, axis=0) # Simple derivative
    Lx = np.roll(landmask, -1, axis=1) - np.roll(landmask, 1, axis=1)
    
    Ly_c = np.roll(landmask, -1, axis=0) - np.roll(landmask, 1, axis=0)
    Ly_c += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (-1,1), axis=(0,1)) # Include y-component of diagonal neighbours
    Ly_c += - np.roll(landmask, (1,-1), axis=(0,1)) - np.roll(landmask, (1,1), axis=(0,1))
    
    Lx_c = np.roll(landmask, -1, axis=1) - np.roll(landmask, 1, axis=1)
    Lx_c += np.roll(landmask, (-1,-1), axis=(1,0)) + np.roll(landmask, (-1,1), axis=(1,0)) # Include x-component of diagonal neighbours
    Lx_c += - np.roll(landmask, (1,-1), axis=(1,0)) - np.roll(landmask, (1,1), axis=(1,0))
    
    v_x = -Lx*(shore)
    v_y = -Ly*(shore)
    
    v_x_c = -Lx_c*(shore_c)
    v_y_c = -Ly_c*(shore_c)
    
    v_x = v_x + v_x_c
    v_y = v_y + v_y_c

    magnitude = np.sqrt(v_y**2 + v_x**2)
    # the coastal nodes between land create a problem. Magnitude there is zero
    # I force it to be 1 to avoid problems when normalizing.
    ny, nx = np.where(magnitude == 0)
    magnitude[ny, nx] = 1

    v_x = -1 * v_x/magnitude ##  Multiply by -1 for backtracking
    v_y = -1 * v_y/magnitude

    return v_x, v_y

v_x, v_y = create_displacement_field(landmask)

def distance_to_shore(landmask, dx=1):
    ci = get_coastal_nodes(landmask) # direct neighbours
    dist = ci*dx                     # 1 dx away
    
    ci_d = get_coastal_nodes_diagonal(landmask) # diagonal neighbours
    dist_d = (ci_d - ci)*np.sqrt(2*dx**2)       # sqrt(2) dx away
        
    return dist+dist_d

d_2_s = distance_to_shore(landmask)

SMOCfile = "../../Raw/Anomalies/01.nc"
fh = xr.open_dataset("../../Raw/Anomalies/01.nc")
time =  fh.variables['time'][:]
filenames = {'U': "../../Raw/Anomalies/01.nc",
             'V': "../../Raw/Anomalies/01.nc",
             'W': "../../Particle_Tracking_R0/Vertical_Velocity_Anomalies/CMEMS_vertical_velocity_01.nc",
             'T': "../../Raw/Anomalies/01.nc"}

variables = {'U': 'uo',
              'V': 'vo',
              'W': 'wo',
              'T': 'thetao'}

dimensions = {'lat': 'latitude',
              'lon': 'longitude',
              'depth': 'depth',
              'time': 'time'}

fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)

class DisplacementParticle(JITParticle):
    dU = Variable('dU')
    dV = Variable('dV')
    d2s = Variable('d2s', initial=1e3)
    t = Variable('t', initial=fieldset.T)
    
def set_displacement(particle, fieldset, time):
    particle.d2s = fieldset.distance2shore[time, particle.depth,
                               particle.lat, particle.lon]
    if  particle.d2s < 0.5:
        dispUab = fieldset.dispU[time, particle.depth, particle.lat,
                               particle.lon]
        dispVab = fieldset.dispV[time, particle.depth, particle.lat,
                               particle.lon]
        particle.dU = dispUab
        particle.dV = dispVab
    else:
        particle.dU = 0.
        particle.dV = 0.
    
def displace(particle, fieldset, time):    
    if  particle.d2s < 0.5:
        particle.lon += particle.dU*particle.dt
        particle.lat += particle.dV*particle.dt

u_displacement = v_x
v_displacement = v_y
fieldset.add_field(Field('dispU', data=u_displacement[latmin:latmax,lonmin:lonmax],
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))
fieldset.add_field(Field('dispV', data=v_displacement[latmin:latmax,lonmin:lonmax],
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))
fieldset.dispU.units = GeographicPolar()
fieldset.dispV.units = Geographic()

fieldset.add_field(Field('landmask', landmask[latmin:latmax,lonmin:lonmax],
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))
fieldset.add_field(Field('distance2shore', d_2_s[latmin:latmax,lonmin:lonmax],
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))

fieldset.add_constant_field('Kh_zonal', 933, mesh='spherical')
fieldset.add_constant_field('Kh_meridional', 933, mesh='spherical')

pset = ParticleSet.from_list(fieldset=fieldset, pclass=DisplacementParticle, time=fieldset.U.grid.time[len(time)-2], depth = np.full(500,20),
                             lon=np.full(500,32.75), lat=np.full(500,-27.41667))

def SampleT(particle, fieldset, time):  # Custom function that samples fieldset.P at particle location
    particle.t = fieldset.T[time, particle.depth, particle.lat, particle.lon]


def DeleteParticle(particle, fieldset, time):
    particle.delete()

#kernels = pset.Kernel(AdvectionRK4_3D) + pset.Kernel(SampleT) + pset.Kernel(DiffusionUniformKh) 
kernels = pset.Kernel(displace) + pset.Kernel(AdvectionRK4_3D) + pset.Kernel(DiffusionUniformKh) + pset.Kernel(set_displacement) + pset.Kernel(SampleT)

output_file = pset.ParticleFile(name="01.nc", outputdt=delta(hours=24))

pset.execute(kernels,
              runtime=timedelta(days=40),
              dt=-timedelta(hours=0.1),
              output_file=output_file,
              recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
output_file.close()

















