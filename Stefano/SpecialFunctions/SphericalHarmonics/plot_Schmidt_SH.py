from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import lpmv # associate legendre functions, in vector form
import math
from SH_library import SchmidtSH
from SH_library import DthSchmidtSH

# Calculate Schmidt semi-normalized SH
# Based on Connerney, 2015 (Treatise on Geophysics) and
# Encyclopedia of Geomagnetism, page 382
# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.


plt.figure(1, figsize=(10,8))

map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.25)

# draw the edge of the map projection region (the projection limb)

# draw lat/lon grid lines every 30 degrees.
map1.drawmeridians(np.arange(0,360,30))
map1.drawparallels(np.arange(-90,90,30))

# make up some data on a regular lat/lon grid.
nlats = 73; nlons = 145; delta = 2.*np.pi/(nlons-1)

#original definitions of lats and lons

lats = np.linspace(np.pi/2,-np.pi/2,num=nlats)
lons = np.linspace(-np.pi,np.pi,num=nlons)
lats, lons = np.meshgrid(lats,lons)
lats = np.transpose(lats)
lons = np.transpose(lons)

# define order and degree
m=3
l=4
# adjust latitude
theta = -lats+np.pi/2

# Spherical harmonic
SH1 = SchmidtSH(l,m,theta,lons,'c')
# derivative
dSH1dth = DthSchmidtSH(l,m,theta,lons,'c')
# compute native map projection coordinates of lat/lon grid.
x, y = map1(lons*180./np.pi, lats*180./np.pi)

# contour data over the map.
cs = map1.contourf(x,y,SH1, cmap='seismic')
plt.colorbar(cs,fraction=0.03, pad=0.04)
plt.show()

# double check on the zero meridian
plt.figure(2)
plt.plot(theta[:,72],SH1[:,72], label = 'SH at phi=0') # plot at the zero meridian

expected = np.sqrt(7*5/float(2))/float(2) * np.multiply(np.sin(theta[:,72])**3,np.cos(theta[:,72])) 

plt.plot(theta[:,72],expected, label = 'expected') # plot at the zero meridian
plt.plot(theta[:,72],dSH1dth[:,72], label = 'dSHdth at phi=0') # plot at the zero meridian
plt.legend()
plt.show()