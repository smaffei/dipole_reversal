import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


folder = '/nfs/see-fs-01_users/earsmaf/SWIGS/OptimalFlow/dipole_CHAOS6_2019.0007_max_F_Leeds/'
input_file = folder+'FLOW_VECTORS_RANDOM.DAT'
output_file = folder+'flow.pdf'

X,Y,U,V = np.loadtxt(input_file,unpack=True)

# display max flow on input grid
USQ = []
for i in range(len(U)):
    USQ.append( ( U[i]**2 + V[i]**2 )**0.5 )
print ('Largest value of |u| on grid is ', max(USQ))
long = X
lat = Y

plt.figure( figsize=(10,8))

#Mollweide projection
map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# Orthographic projection
map = Basemap(projection='ortho',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.50)


# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0,360,30),linewidth=1)
map.drawparallels(np.arange(-90,90,30),linewidth=1)
# make up some data on a regular lat/lon grid.

uproj,vproj,xx,yy = \
map.rotate_vector(np.array(V),-np.array(U),np.array(long),np.array(lat),returnxy=True)


for i in range(0,len(xx)):
    Q = map.quiver(xx[i],yy[i],uproj[i],vproj[i],angles='uv',pivot='tail',scale=400,width=0.002,headwidth=6)
plt.rcParams['font.size'] = 18
plt.quiverkey(Q,0.5,-0.1, 20,'20 km/yr',coordinates='axes')

#plt.show()
plt.savefig(output_file, bbox_inches='tight',pad_inches=0.7,dpi=400)
