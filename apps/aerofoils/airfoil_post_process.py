""" djl: 07/2022: Post-processing of airfoil coefficients. Adapted from an original version by Dr. Markus Zauner."""

from opensbli.postprocess.plot_functions import *
from opensbli.core.block import SimulationBlock
import numpy as np
import h5py
import os

## Useful functions, turn into a class later
def instantaneous_q(fname, block_number):
    # Read the data only on the airfoil surface
    data_file, block_name, dsets, shape = PP.read_block(fname, block_number)
    if conservative:
        rho = PP.read_full_dset(data_file, block_name, "rho_B%d" % block_number, partial_slice = surface)
        rhou = PP.read_full_dset(data_file, block_name, "rhou0_B%d" % block_number, partial_slice = surface)
        rhov = PP.read_full_dset(data_file, block_name, "rhou1_B%d" % block_number, partial_slice = surface)
        rhow = PP.read_full_dset(data_file, block_name, "rhou2_B%d" % block_number, partial_slice = surface)
        E = PP.read_full_dset(data_file, block_name, "rhoE_B%d" % block_number, partial_slice = surface)
        u = rhou/rho
        v = rhov/rho
        w = rhow/rho
        p = (gamma - 1)*(E - 0.5*(u**2+v**2+w**2)*rho)
    else:
        rho = PP.read_full_dset(data_file, block_name, "rho_B%d" % block_number, partial_slice = surface)
        u = PP.read_full_dset(data_file, block_name, "u0_B%d" % block_number, partial_slice = surface)
        v = PP.read_full_dset(data_file, block_name, "u1_B%d" % block_number, partial_slice = surface)
        w = PP.read_full_dset(data_file, block_name, "u2_B%d" % block_number, partial_slice = surface)
        E = PP.read_full_dset(data_file, block_name, "Et_B%d" % block_number, partial_slice = surface)
        p = (gamma - 1)*(E - 0.5*(u**2+v**2+w**2))*rho

    # The rest
    a = np.sqrt(gamma*p/rho)
    M = np.sqrt(u**2 + v**2 + w**2)/a
    T = gamma*(Minf**2)*p/rho
    mu = compute_viscosity(T)
    return rho, u, v, w, E, p, T, M, mu

def time_averaged_q(fname, block_number):
    # Read the data only on the airfoil surface
    data_file, block_name, dsets, shape = PP.read_block(fname, block_number)
    if conservative:
        rho = PP.read_full_dset(data_file, block_name, "rhomean_B%d" % block_number, partial_slice = surface)
        rhou = PP.read_full_dset(data_file, block_name, "rhou0mean_B%d" % block_number, partial_slice = surface)
        rhov = PP.read_full_dset(data_file, block_name, "rhou1mean_B%d" % block_number, partial_slice = surface)
        rhow = PP.read_full_dset(data_file, block_name, "rhou2mean_B%d" % block_number, partial_slice = surface)
        E = PP.read_full_dset(data_file, block_name, "rhoEmean_B%d" % block_number, partial_slice = surface) # WARNING: check E is correct here
        u = rhou/rho
        v = rhov/rho
        w = rhow/rho
        p = (gamma - 1)*(E - 0.5*(u**2+v**2+w**2)*rho)
    else:
        rho = PP.read_full_dset(data_file, block_name, "rhomean_B%d" % block_number, partial_slice = surface)
        rhou = PP.read_full_dset(data_file, block_name, "rhou0mean_B%d" % block_number, partial_slice = surface)
        rhov = PP.read_full_dset(data_file, block_name, "rhou1mean_B%d" % block_number, partial_slice = surface)
        rhow = PP.read_full_dset(data_file, block_name, "rhou2mean_B%d" % block_number, partial_slice = surface)
        E = PP.read_full_dset(data_file, block_name, "E_mean_B%d" % block_number, partial_slice = surface)
        u = rhou/rho
        v = rhov/rho
        w = rhow/rho
        p = (gamma - 1)*(E - 0.5*(u**2+v**2+w**2))*rho

    # The rest
    a = np.sqrt(gamma*p/rho)
    M = np.sqrt(u**2 + v**2 + w**2)/a
    T = gamma*(Minf**2)*p/rho
    mu = compute_viscosity(T)
    return rho, u, v, w, E, p, T, M, mu

def compute_wall_normal_derivative(variable):
    ny = Ny
    # Ly = np.max(y)
    Ly = Ly
    delta = Ly/(ny-1.0)
    # D11 = D11[0:6, :]
    var = variable[0:6, :]
    coeffs = np.array([-1.83333333333334, 3.00000000000002, -1.50000000000003, 0.333333333333356, -8.34617916606957e-15, 1.06910884386911e-15])
    coeffs = coeffs.reshape([6, 1])
    df_dy = sum(var*coeffs)/delta
    return df_dy

def calculate_angles(D00, D01, D10, D11, detJ):
    """ Calculated on a single y plane."""
    D00, D01, D10, D11, detJ = D00[0,:], D01[0,:], D10[0,:], D11[0,:], detJ[0,:]
    Nx = np.size(D00)
    # det = (D00*D11-D01*D10)
    aidet = 1.0/detJ
    # dydx = D10/D00
    dydx = D00 / (-1*D01)
    # Storage arrays
    ss, dss, th_, S_, S2_ = np.zeros(Nx), np.zeros(Nx), np.zeros(Nx), np.zeros(Nx), np.zeros(Nx)
    # Calculate the angles
    for i in range(0,Nx):
        th_[i]=np.arctan(dydx[i])
        if(D11[i]*detJ[i]>0): # dxdxi
            S_[i] = 1
            # print("S positive", "x= %.3f" % x[0,i])
        else:
            # print("S negative", "x= %.3f" % x[0,i])
            S_[i] = -1
        if(D10[i]*detJ[i]>0): # dydxi
            S2_[i] = 1
        else:
            S2_[i] = -1
        if(i>0):
            ss[i] = ss[i-1]+((x[0,i]-x[0,i-1])**2+(y[0,i]-y[0,i-1])**2)**0.5
            dss[i] = np.fabs(ss[i]-ss[i-1])
    # exit()
    return ss, dss, th_, S_, S2_, aidet, dydx

def compute_viscosity(T):
    # mu = (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT))
    mu = T**0.7
    return mu


def compute_wall_normal_derivative(variable):
    delta = Ly/(Ny-1.0)
    var = variable[0:6, :]
    coeffs = np.array([-1.83333333333334, 3.00000000000002, -1.50000000000003, 0.333333333333356, -8.34617916606957e-15, 1.06910884386911e-15])
    coeffs = coeffs.reshape([6, 1])
    df_dy = sum(var*coeffs)/delta
    return df_dy

def surface_derivatives(fname, block_number, D00, D01, D10, D11, detJ):
    """ Uses the metric relations to find the derivatives on the airfoil surface."""
    ## WARNING: Can change this to not use loops later
    D00, D01, D10, D11, detJ = D00[0,:], D01[0,:], D10[0,:], D11[0,:], detJ[0,:]
    Nx = np.size(D00)
    # corrf = 1.0/(Ny-1)
    corrf = 1.0
    # Read the file again
    data_file, block_name, dsets, shape = PP.read_block(fname, block_number)
    rho = PP.read_full_dset(data_file, block_name, "rhomean_B%d" % block_number, partial_slice=np.s_[:,0:6+2*nhalo,startx:])
    rhou = PP.read_full_dset(data_file, block_name, "rhou0mean_B%d" % block_number, partial_slice=np.s_[:,0:6+2*nhalo,startx:])
    rhov = PP.read_full_dset(data_file, block_name, "rhou1mean_B%d" % block_number, partial_slice=np.s_[:,0:6+2*nhalo,startx:])
    print(rho.shape)
    u, v = rhou/rho, rhov/rho
    
    dudeta, dvdeta = np.zeros((Nz, Nx)), np.zeros((Nz, Nx))

    for k in range(0, Nz):
        dudeta[k,:] = compute_wall_normal_derivative(u[k,:])
        dvdeta[k,:] = compute_wall_normal_derivative(v[k,:])

    dudeta = dudeta / corrf
    dvdeta = dvdeta / corrf


    dudy = np.zeros((Nz, Nx)) # Nz, Nx ordering
    dvdx = np.zeros((Nz, Nx))

    for k in range(0,Nz):
        dudy[k,:] = (dudeta[k,:] / (detJ[:]*D00[:]))
        dvdx[k,:] = (dvdeta[k,:] / (-detJ[:]*D01[:]))
    return dudy, dvdx

def aerodynamic_coefficients(dudy, dvdx, p, mu):
    Nz = p.shape[0]
    Nx = p.shape[-1]
    # print(p.shape
    # exit()
    pinf = 1.0 / (gamma*Minf*Minf)
    tau_wall, Cp = np.zeros((Nz,Nx)), np.zeros((Nz,Nx))
    Cl, Cdp, Cdf = np.zeros(Nz), np.zeros(Nz), np.zeros(Nz)
    i = 0
    for k in range(0,Nz):
        tau_wall [k,i]=S_[i]*mu[k,i]*(dudy[k,i]*np.abs(np.cos(th_[i]))-dvdx[k,i]*np.abs(np.sin(th_[i])))
        Cp[k,i] =(p[k,i]-pinf)/(0.5*gamma*Minf*Minf)

    for i in range(1,Nx):
        for k in range(0,Nz):
            dlts=ss[i]-ss[i-1]
            fa = -S_[i]*(p[k,i]-pinf)*np.abs(np.cos(th_[i]))
            fb = -S_[i-1]*(p[k, i-1]-pinf)*np.abs(np.cos(th_[i-1]))
            # Lift coefficient
            Cl[k]=Cl[k]+0.5*(fa+fb)*dlts  / (0.5*1*1**2) # changeme
            # Drag 

            fa = S2_[i]*p[k,i]*abs(np.sin(th_[i]))
            fb = S2_[i]*p[k,i-1]*abs(np.sin(th_[i-1]))
            Cdp[k] = Cdp[k]+0.5*(fa+fb)*dlts  / (0.5*1*1**2) # changeme

            fa = S_[i]*mu[k,i]*(dudy[k,i]*np.abs(np.cos(th_[i]))-dvdx[k,i]*np.abs(np.sin(th_[i])))
            fb = S_[i]*mu[k,i-1]*( dudy[k,i-1]*np.abs(np.cos(th_[i-1]))-dvdx[k,i-1]*np.abs(np.sin(th_[i-1])) )
            Cdf[k] = Cdf[k]+0.5*(fa+fb)*dlts
            # Wall shear stress
            tau_wall [k,i] = S_[i]*mu[k,i]*(dudy[k,i]*np.abs(np.cos(th_[i]))-dvdx[k,i]*np.abs(np.sin(th_[i])))
            # Pressure coefficient
            Cp[k,i] =(p[k,i]-pinf)/(0.5*gamma*Minf*Minf)
    return Cl, Cdp, Cdf, Cp, tau_wall


# Step 1 read a file
inPath = './'
outPath = './file_output/'

startx = 0
try:
    os.mkdir(outPath)
except FileExistsError:
    pass

# Utility function
PP = OpenSBLIPreProcess()

block_number = 0
# Read the grid
grid_file, block_name, dsets, shape = PP.read_block('./data.h5', block_number)
halos = np.abs(grid_file[block_name][dsets[0]].attrs['d_p'])
nhalo = halos[0]
# Get the coordinates on the surface
yloc = 0
surface = np.s_[:,yloc + nhalo,startx:]
x = PP.read_full_dset(grid_file, block_name, 'x0_B%d' % block_number, partial_slice=surface)

y = PP.read_full_dset(grid_file, block_name, 'x1_B%d' % block_number, partial_slice=surface)
z = PP.read_full_dset(grid_file, block_name, 'x2_B%d' % block_number, partial_slice=surface)
# Number of grid points per direction (without halos) on this surface slice
Nx, Ny, Nz = shape[-1], shape[1], shape[0]
# for i in range(Ny):
#     print(y[0,i])
# exit()
# # Domain size
Lx, Ly, Lz = 0, 10.0, 0
# Grid spacing
dx = 0
dy = 0
dz = z[1,0] - z[0,0]
# Load metrics on the surface
D00 = PP.read_full_dset(grid_file, block_name, 'D00_B%d' % block_number, partial_slice=surface)
D01 = PP.read_full_dset(grid_file, block_name, 'D01_B%d' % block_number, partial_slice=surface)
D10 = PP.read_full_dset(grid_file, block_name, 'D10_B%d' % block_number, partial_slice=surface)
D11 = PP.read_full_dset(grid_file, block_name, 'D11_B%d' % block_number, partial_slice=surface)
detJ = PP.read_full_dset(grid_file, block_name, 'detJ_B%d' % block_number, partial_slice=surface)
# Calculate angles over the surface
# Calculate angles around the surface
ss, dss, th_, S_, S2_, aidet, dydx = calculate_angles(D00, D01, D10, D11, detJ)


input_files = ['./restart.h5']
input_files = ['./stats_output.h5']

time_averaged = True
conservative = False
gamma = 1.4
Minf = 0.72
Re = 5.0e5

for idx, fname in enumerate(input_files):
    # Read flow data
    if time_averaged:
        rho, u, v, w, E, p, T, M, mu = time_averaged_q(fname, block_number)
    else:
        rho, u, v, w, E, p, T, M, mu = instantaneous_q(fname, block_number)
    # Get velocity derivatives at the surface
    dudy, dvdx = surface_derivatives(fname, block_number, D00, D01, D10, D11, detJ)


    # Calculate aerodynamic coefficients
    Cl, Cdp, Cdf, Cp, tau_wall = aerodynamic_coefficients(dudy, dvdx, p, mu)
    print("Cl: {:}, Cdp: {:}, Cdf: {:}".format(np.mean(Cl), np.mean(Cdp), np.mean(Cdf)))

    # Test plot
    # for i in Cp:
    #     print(i)
    # exit()
    # Average over the span
    Cp = np.mean(Cp, axis=0)
    # Calculate skin friction
    Cf = tau_wall/(0.5*Re)
    Cf = np.mean(Cf, axis=0)

    print(Cf)
    print(np.min(Cf), np.max(Cf))
    plt.plot(x[startx,:], Cf)
    # plt.gca().invert_yaxis()
    plt.show()

    plt.clf()
    plt.plot(x[startx,:], Cp)
    plt.gca().invert_yaxis()
    plt.xlabel('x')
    plt.ylabel('Cp')
    plt.show()