#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 14:21:00 2019

@author: earsmaf

"""

import numpy as np
import math
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


# Create equispoaced grid on a sphere
# based on a document by Markus Deserno (How to generate equidistributed points on the surface of a sphere)

Nth = 4

r_a = 6371e3 # radius of the earth

dth = np.pi/(Nth-1) # distance in lat
npoint = 0
x=[]
y=[]
z=[]
lats=[]
lons=[]
for ith in range(Nth):
    rs = r_a * np.sin(ith*dth) # radius of the circle    
    zs = r_a * np.cos(ith*dth) # and its height above equator
    if (ith==0) or (ith==Nth-1):
        nphi = 0
        dphi = 0
        x.append(0) 
        y.append(0)
        z.append(zs)
        lats.append(ith*dth)
        lons.append(0)
    else:
        nphi = math.ceil(2*np.pi*rs/dth/r_a)
        dphi=2*np.pi/nphi; # distance in lon    
        for iphi in range(int(nphi)):
            x.append(rs*np.cos(iphi*dphi))
            y.append(rs*np.sin(iphi*dphi))
            z.append(zs)
            lats.append(ith*dth)
            lons.append(iphi*dphi)
            
x = np.asarray(x)
y = np.asarray(y)
z = np.asarray(z)
lats = np.asarray(lats)
lons = np.asarray(lons)


theta = -lats+np.pi/2

#plot
theta = np.rad2deg(theta)
lons = np.rad2deg(lons)

#lons, lats = np.meshgrid(lons, lats)

fig = plt.figure(figsize=(11, 6))
#    plt.cla()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
#ax.gridlines(linewidth=1, alpha=0.5, linestyle='-')

ax.coastlines()
ax.set_global()

# geomagnetic poles location
plt.scatter(lons,theta,
           s=20, c='b',edgecolor='k',alpha=1, zorder=3,marker='o',
           transform=ccrs.PlateCarree())

plt.show()