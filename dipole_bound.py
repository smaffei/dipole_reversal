import sys
sys.path.insert(0,'/nfs/see-fs-01_teaching/ee12sg/anaconda3/envs/dipole_reversal/lib/python3.7/site-packages')
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import pyshtools


def dipole_bound(lmax_u, lmax_B, gauss_coeffs, target_rms, toroidal_only=False, save_files=False):

    """
    Function to call FORTRAN program dipole_bound_modified.F90

    Parameters
    ----------
        'lmax_u' : int.
            Maximum spherical harmonic degree, l, of the optimised flow.
        'lmax_B' : int.
            Maximum spherical harmonic degree, l, of the input magnetic field.
        'gauss_coeffs' : numpy array.
            Array of gauss coefficients for magnetic field. Must be a 2D array
            with columns relating to the sin and cosine terms in the format:
                g10 (sin), h10 (cosine)
                g11 (sin), h11 (cosine)
                g20 (sin), h20 (cosine)
                g21 (sin), h21 (cosine)
                etc...
        'target_rms' : float.
            RMS velocity in km/yr of the optimised flow.

        'toroidal_only' : bool.
            if True the optimised flow will have zero poloidal component (default: False)

        'save_files' : bool.
            if True the output files from dipole_bound.F90 will not be deleted (default: False)


    Returns
    ---------

        dipole_rate (float)
            Optimised rate of change of g10 gauss coefficient in nT/yr

        pol_ke (float)
            Kinetic energy in poloidal component of flow

        tor_ke (float)
            Kinetic energy in toroidal component of flow

    """

    if toroidal_only:
        restriction = 2
    else:
        restriction = 0

    if save_files:
        save_flag = 1
    else:
        save_flag = 0


    shape = gauss_coeffs.shape
    assert shape[1] == 2, 'gauss_coeffs array not correct shape\nRequires 2 columns instead of {}'.format(shape[2])


    harmonics = np.ones(shape)

    l, m = 1,0
    for i in range(shape[0]):
        harmonics[i,0] = l
        harmonics[i,1] = m

        if m < l:
            m+=1
        else:
            l+=1
            m = 0

        if l == lmax_B+1:
            break

    harmonics = harmonics[:i+1,:]
    gauss_coeffs = gauss_coeffs[:harmonics.shape[0],:]

    try:

        filename_gauss = write_gauss_coeffs(harmonics, gauss_coeffs)
        filename_input = write_input_file(lmax_u, lmax_B, target_rms, restriction, save_flag, filename_gauss)

    except:
        breakpoint()


    output = subprocess.run(['./modified_code/dipole_bound < {}'.format(filename_input)], shell=True, capture_output=True)

    stdout = output.stdout.decode('utf-8')
    stderr = output.stderr.decode('utf-8')

    if len(stderr) > 0:
        raise ValueError(stderr)


    os.remove(filename_gauss)
    os.remove(filename_input)

    files = ['BR_DOT_ES.DAT', 'BR_ES.DAT', 'BR_DOT_ES_CENTRED.DAT', 'BR_ES_CENTRED.DAT', 'FLOW_L_SPEC.DAT',
             'FLOW_VECTORS_CENTRED.DAT', 'FLOW_VECTORS_RANDOM.DAT', 'HARMONICS.DAT', 'OPTIMAL_FLOW.DAT']


    if not save_files:
        for f in files:
            if os.path.isfile(f):
                os.remove(f)
    else:
        if not os.path.isdir('code_output'):
            os.mkdir('./code_output')
        for f in files:
            if os.path.isfile(f):
                os.rename(f,'./code_output/'+f)

    dipole_rate, pol_ke, tor_ke = [float(x) for x in stdout.split()]

    return dipole_rate, pol_ke, tor_ke

def write_input_file(lmax_u, lmax_B_obs, target_rms, restriction, save_flag, filename_gauss):
    inputs = [lmax_u,lmax_B_obs,target_rms,restriction,save_flag]

    filename = 'temporary_input_file.txt'

    with open(filename, 'a') as f:
        for i in inputs:
            f.write('{:.0f}\n'.format(i))
        f.write(filename_gauss)

    return filename

def write_gauss_coeffs(harmonics, gauss_coeffs):
    filename='temp_file_gauss.txt'
    np.savetxt(filename, np.hstack([harmonics,gauss_coeffs]), fmt=('%d','%d','%.5e','%.5e'), delimiter='\t', header='Gauss Coefficients')
    return filename

def read_BR(directory='code_output/'):

    files = ['BR_DOT_ES.DAT', 'BR_ES.DAT']

    Br_dot_data = np.genfromtxt(directory+files[0])

    lats = np.unique(Br_dot_data[:,1])[::-1]
    lons = np.unique(Br_dot_data[:,0])[::-1]-180

    LATS, LONS = np.meshgrid(lats,lons)

    a = lons.size
    b = lats.size

    shape = (a,b)

    Br_dot = Br_dot_data[:,2].reshape(shape)


    Br_data = np.genfromtxt(files[1])
    Br = Br_data[:,2].reshape(shape)


    #copy data at longitude = 180 to -180
    LONS = np.vstack((LONS[-1,:]+360,LONS))
    LATS = np.vstack((LATS[-1,:],LATS))
    Br = np.vstack((Br[-1,:],Br))
    Br_dot = np.vstack((Br_dot[-1,:],Br_dot))


    return LONS,LATS,Br,Br_dot

def read_vectors(directory='code_output/'):

    file = 'FLOW_VECTORS_CENTRED.DAT'

    data = np.genfromtxt(directory+file)

    lats = np.unique(data[:,1])
    lons = np.unique(data[:,0])+180

    LATS, LONS = np.meshgrid(lats,lons)

    a = lats.size
    b = lons.size

    shape = (a,b)
    v = data[:,2].reshape(shape).transpose()
    u = data[:,3].reshape(shape).transpose()

    return LONS,LATS,u,v

def setup_map_figure(*args,figsize=(10,4),**kwargs):

    #Set up figure
    plt.figure(*args,figsize=figsize,**kwargs)

    #Add features
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.set_global()

    ax.add_feature(cartopy.feature.LAND,facecolor=[0.8,0.8,0.8], edgecolor='face')
    #ax.add_feature(cartopy.feature.OCEAN,facecolor=[0.95,0.95,0.95])

    return ax

def plot_vectors(ax,lats,lons,u,v, linewidth=2, *args, **kwargs):

    speed = np.sqrt(u**2 + v**2)
    lw = linewidth*speed/np.max(speed)

    #ax.quiver(lons, lats, vector_u, vector_v, color=[0.8,0,0], transform=ccrs.PlateCarree(), width=width, scale=scale)
    im = ax.streamplot(lons, lats, u, v, *args, color=speed, cmap='jet',
                       transform=ccrs.PlateCarree(), linewidth=lw, **kwargs)

    cbar = plt.colorbar(im.lines)
    cbar.set_label('Speed [km/yr]')
    #ax.coastlines()

    return im, cbar

def plot_contourf(ax,lats,lons,values,*args,**kwargs):

    im = ax.contourf(lons,lats,values,*args,transform=ccrs.PlateCarree(),**kwargs)
    cbar = plt.colorbar(im)
    ax.coastlines()

    return im, cbar

def power_spectra(lmax,gauss_coeffs):

    R = np.zeros(lmax)
    c = 0
    for l in range(1,lmax+1):
        m = 0
        for m in range(l+1):
            R[l-1] += (l+1)*np.sum(gauss_coeffs[c,:]**2)
            c += 1
    return R

def extrapolate_model(lmax,lmax_new,gauss_coeffs):

    # Current spectra
    R = power_spectra(lmax,gauss_coeffs)
        
    fit = np.polyfit(range(1,lmax+1),np.log10(R),1);

    #  Calculate the predicted power for higher degrees based on this for i = 5:10
    R_ex = np.zeros(10)
    R_ex[:4] = R[:]

    for i in range(lmax+1,lmax_new+1):
        R_ex[i-1] = 10**np.polyval(fit,i)

    f = np.sum(R)/np.sum(R_ex)

    R_ex = f*R_ex


    coeffs_ex = np.sqrt(f)*gauss_coeffs[:]
    for l in range(lmax+1,lmax_new+1):
        extrap = []
        for m in range(l+1):
            extrap.append(np.random.rand(2))
            if m == 0:
                extrap[-1][1] = 0

        extrap = np.array(extrap)
        f = R_ex[l-1]/((l+1)*np.sum(extrap**2))

        coeffs_ex = np.vstack((coeffs_ex, np.sqrt(f)*extrap))


    return coeffs_ex

def evaluate_B(gauss_coeffs, lat, lon):

    if not type(lat) == np.ndarray:
        lat = np.array([lat])
        lon = np.array([lon])

    theta = (np.pi/180)*(90-lat)
    phi = (np.pi/180)*lon


    flag = True
    l, m = 0, 0
    c = 1
    while flag:
        m += 1
        if m > l:
            l += 1
            m = 0

        if c == gauss_coeffs.shape[0]:
            flag = False
        else:
            c += 1

    lmax = l
    a = 6371e3
    r = 6371e3

    #Evaluate Legendre polynomials and their derivatives at these locations
    leg = np.zeros((theta.size,(lmax+1)*(lmax+2)//2))
    leg_deriv = leg.copy()

    for i in range(theta.size):
        try:
            leg[i,:], leg_deriv[i,:] = pyshtools.legendre.PlmSchmidt_d1(lmax, np.cos(theta[i]))
            leg_deriv[i,:] = leg_deriv[i,:]*-np.sin(theta[i])
        except:
            breakpoint()

    leg = leg[:,1:]
    leg_deriv = leg_deriv[:,1:]

    Br = np.zeros((theta.size,phi.size))
    Bth = Br.copy()
    Bph = Br.copy()
    c=0
    
    for l in range(1,lmax+1):
        for m in range(l+1):

            f1 = gauss_coeffs[c,0]*np.cos(m*phi) + gauss_coeffs[c,1]*np.sin(m*phi)
            t1,t2 = np.meshgrid(f1,leg[:,c])
            d = t1*t2
            Br = Br + (l+1)*(a/r)**(l+2)*d

            t1,t2 = np.meshgrid(f1,leg_deriv[:,c])
            d = t1*t2
            Bth = Bth - (a/r)**(l+2)*d

            f1 = m*(-gauss_coeffs[c,0]*np.sin(m*phi) + gauss_coeffs[c,1]*np.cos(m*phi))
            t1,t2 = np.meshgrid(f1,leg[:,c]/np.sin(theta))
            d = t1*t2
            Bph = Bph - (a/r)**(l+2) *d

            c = c+1


    return Br, Bth, Bph
