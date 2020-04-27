import numpy as np
import matplotlib.pyplot as plt

secinyr=31536000.0
r2d    =180.0/np.pi

# theta is colatitude in radians
# time t in advection units
t,theta = np.loadtxt('gauss_coeffs_surface_dipmom', skiprows=1, usecols=(1,6), unpack='true')
  
tilt= theta * r2d
dt  = t[1:] -  t[0:-1]

#dtheta = np.gradient(theta,dt,edge_order=2) # Would not work with my python :(
dtilt = np.gradient(tilt[1:],dt)

plt.xlabel("time")
plt.plot(t,tilt)
plt.show()

plt.plot(t[1:],dtilt)
plt.show()

print(np.max(dtilt))
