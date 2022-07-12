from opensbli.postprocess.plot_functions import *
from opensbli.core.block import SimulationBlock
import numpy as np
import h5py
import os
""" Extracts a 2D slice from a 3D datafile and writes it back out as an HDF5 file for Tecplot. DJL: 06/22."""

# Step 1 read a file
inPath = './'
outPath = './file_output/'

try:
    os.mkdir(outPath)
except FileExistsError:
    pass


class AirfoilPostProcess(object):
    """ Routines for post-processing monitors during the simulation of airfoils."""
    def __init__(self, restarts, stats, grid, multi_block=False):
        # directions = ['x', 'y', 'z']
        # self.ndim = block.ndim
        # self.block = block
        self.conservative = False
        self.gamma = 1.4
        self.Minf = 0.72
        self.Re = 5.0e5
        self.yloc = 1 # 1 point off the airfoil surface
        self.restarts = restarts
        self.stats = stats
        self.grid = grid
        if multi_block:
            self.nblocks = multi_block.nblocks
        else:
            self.nblocks = 1
        # Process everything
        self.main()
        return


    def read_metrics(self, block_no):
        a = OpenSBLIPreProcess()
        # Grid read
        a.read_grid(block_no)
        self.x, self.y, self.z = a.x, a.y, a.z
        self.shape = a.shape
        # Assuming uniform spacing in spanwise z direction
        self.dz = a.z[1,10,10] - a.z[0,10,10]
        # Metrics read
        a.read_file(self.grid, block_no)
        D00, D01, D10, D11 = a.read_full_dset('D00_B%d' % block_no)[:,self.yloc,:], a.read_full_dset('D01_B%d' % block_no)[:,self.yloc,:], a.read_full_dset('D10_B%d' % block_no)[:,self.yloc,:], a.read_full_dset('D11_B%d' % block_no)[:,self.yloc,:]
        self.metrics = {'D00': D00, 'D01': D01, 'D10' : D10, 'D11': D11}
        self.PP = a
        return

    def calculate_angles(self, block_no):
        """ Calculated on a single y plane."""
        D00, D01, D10, D11 = self.metrics['D00'][0,:], self.metrics['D01'][0,:], self.metrics['D10'][0,:], self.metrics['D11'][0,:]
        det = self.dz*(D00*D11-D01*D10)
        self.aidet = 1.0/det
        dydx = D10/D00
        # Calculate the angles around the surface
        nxt = self.shape[-1]
        x, y = self.x[:,self.yloc,:], self.y[:,self.yloc,:]
        # Storage arrays
        ss, dss, th_, S_, S2_ = np.zeros(nxt), np.zeros(nxt), np.zeros(nxt), np.zeros(nxt), np.zeros(nxt)
        # Calculate the angles
        for i in range(0,nxt):
            th_[i]=np.arctan(dydx[i])
            if(D00[i]>0):
                S_[i] = 1
            else:
                S_[i] = -1
            if(D10[i]>0):
                S2_[i] = 1
            else:
                S2_[i] = -1
            if(i>0):
                ss[i] = ss[i-1]+((x[0,i]-x[0,i-1])**2+(y[0,i]-y[0,i-1])**2)**0.5
                dss[i] = np.fabs(ss[i]-ss[i-1])
        self.ss, self.th_ = ss, th_
        return

    def extract_flow_variables(self, block_no):
        # Read the data only on the airfoil surface
        self.PP.read_file(self.restarts, block_no)
        if self.conservative:
            rho = self.PP.read_full_dset("rho_B%d" % block_no)[:,self.yloc,:]
            rhou = self.PP.read_full_dset("rhou0_B%d" % block_no)[:,self.yloc,:]
            rhov = self.PP.read_full_dset("rhou1_B%d" % block_no)[:,self.yloc,:]
            rhow = self.PP.read_full_dset("rhou2_B%d" % block_no)[:,self.yloc,:]
            E = self.PP.read_full_dset("rhoE_B%d" % block_no)[:,self.yloc,:]
            u = rhou/rho
            v = rhov/rho
            w = rhow/rho
            p = (self.gamma - 1)*(E - 0.5*(u**2+v**2+w**2)*rho)
        else:
            rho = self.PP.read_full_dset("rho_B%d" % block_no)[:,self.yloc,:]
            u = self.PP.read_full_dset("u0_B%d" % block_no)[:,self.yloc,:]
            v = self.PP.read_full_dset("u1_B%d" % block_no)[:,self.yloc,:]
            w = self.PP.read_full_dset("u2_B%d" % block_no)[:,self.yloc,:]
            E = self.PP.read_full_dset("Et_B%d" % block_no)[:,self.yloc,:]
            p = (self.gamma - 1)*(E - 0.5*(u**2+v**2+w**2))*rho

        # The rest
        a = np.sqrt(self.gamma*p/rho)
        M = np.sqrt(u**2 + v**2 + w**2)/a
        T = self.gamma*(self.Minf**2)*p/rho
        mu = self.compute_viscosity(T)
        return rho, u, v, w, E, p, T, M, mu

    def compute_wall_normal_derivative(self, variable):
        ny = np.size(self.y[:, 0])
        Ly = self.Ly
        delta = Ly/(ny-1.0)
        # D11 = self.D11[0:6, :]
        var = variable[0:6, :]
        coeffs = np.array([-1.83333333333334, 3.00000000000002, -1.50000000000003, 0.333333333333356, -8.34617916606957e-15, 1.06910884386911e-15])
        coeffs = coeffs.reshape([6, 1])
        df_dy = sum(var*coeffs)/delta
        return df_dy

    def compute_viscosity(self, T):
        # mu = (T**(1.5)*(1.0+self.SuthT/self.RefT)/(T+self.SuthT/self.RefT))
        mu = T**0.7
        return mu

    def compute_skin_friction(self, u, mu):
        # Wall viscosity all x points
        mu_wall = mu[0, :]
        dudy = self.compute_wall_derivative(u)
        tau_wall = dudy*mu_wall
        Cf = tau_wall/(0.5*self.Re)
        return Cf

    def derivative_x(self, var):

        return

    def aerodynamic_coefficients(self, dudy, dvdx, p, mu):
        Nx, Nz = self.shape[-1], self.shape[0]
        th_ = self.th_

        tau_wall = np.zeros((Nx,Nz))
        i=0
        for k in range(0,Nz):
            tau_wall [i,k]=S_[i]*mu[i,k]*(dudy[i,k]*np.abs(np.cos(th_[i]))-dvdx[i,k]*np.abs(np.sin(th_[i])))

        for i in range(1,Nx):
            for k in range(0,Nz):
                dlts=ss[i]-ss[i-1]
                fa = -S_[i]*(P[i,k]-pinf)*np.abs(np.cos(th_[i]))
                fb = -S_[i-1]*(P[i-1,k]-pinf)*np.abs(np.cos(th_[i-1]))
                Cl[ind,k]=Cl[ind,k]+0.5*(fa+fb)*dlts

                fa = S2_[i]*P[i,k]*abs(np.sin(th_[i]))
                fb = S2_[i]*P[i-1,k]*abs(np.sin(th_[i-1]))
                Cdp[ind,k] = Cdp[ind,k]+0.5*(fa+fb)*dlts

                fa = S_[i]*mu[i,k]*(dudy[i,k]*np.abs(np.cos(th_[i]))-dvdx[i,k]*np.abs(np.sin(th_[i])))
                fb = S_[i]*mu[i-1,k]*( dudy[i-1,k]*np.abs(np.cos(th_[i-1]))-dvdx[i-1,k]*np.abs(np.sin(th_[i-1])) )
                Cdf[ind,k] = Cdf[ind,k]+0.5*(fa+fb)*dlts
                tau_wall [i,k]=S_[i]*mu[i,k]*(dudy[i,k]*np.abs(np.cos(th_[i]))-dvdx[i,k]*np.abs(np.sin(th_[i])))
        return

    def surface_derivatives(self, u, v):
        """ Uses the metric relations to find the derivatives on the airfoil surface."""
        ## WARNING: Can change this to not use loops later
        D00, D01, D10, D11 = self.metrics['D00'][0,:], self.metrics['D01'][0,:], self.metrics['D10'][0,:], self.metrics['D11'][0,:]
        # qin=np.zeros((nxt,nyt))
        # qin[:,0]=q[:,0,5]/corrf
        # dudxi=d1xi_2(qin,nxt,nyt,bc4)
        corrf = 1.0

        Nz, Nx = self.shape[0], self.shape[-1]
        print(Nz, Nx)
        dudeta = u / corrf
        dvdeta = v / corrf
        dudy = np.zeros((Nz, Nx)) # Nz, Nx ordering
        dvdx = np.zeros((Nz, Nx))

        for k in range(0,Nz):
            dudy[k,:] = self.aidet[:]*(dudeta[k,:]*D00[:]*self.dz)
            dvdx[k,:] = self.aidet[:]*(-dvdeta[k,:]*D10[:]*self.dz)
        return dudy, dvdx

    def main(self):
        for block_no in range(self.nblocks):
            self.data = {}
            # Read the grid and metrics
            mets = self.read_metrics(block_no)
            # Calculate angles around the surface
            angles = self.calculate_angles(block_no)
            # Read flow data
            rho, u, v, w, E, p, T, M, mu = self.extract_flow_variables(block_no)
            dudy, dvdx = self.surface_derivatives(u, v)
            # Calculate aerodynamic coefficients
            self.aerodynamic_coefficients(dudy, dvdx, p, mu)



        return








# already_processed = sorted(glob.glob(outPath + '/opensbli_output_*.h5'))
# already_processed = [re.findall("\d+", s)[0].lstrip('0') for s in already_processed]
# # Extract data
# a = OpenSBLIPreProcess()
# fnames, iters = a.find_files(inPath)

blocks = 1

restarts = './restart.h5'
stats = './stat_output.h5'
grid = './data.h5'

a = AirfoilPostProcess(restarts, stats, grid)






# for i, f in enumerate(fnames):
#   if iters[i] not in already_processed:
#       # Write to a new 2D HDF5 file
#       h5f = h5py.File(outPath + 'opensbli_output_%s.h5' % iters[i], 'w')
#       for block_no in range(blocks):
#           a.read_grid(block_no)
#           zloc = int(a.shape[0]/2.0)
#           zloc = 0
#           add_grid = True
#           data = a.read_file(f, block_no)
#           # Extract a single slice
#           data_2D = {}
#           # Extract all the variables at this slice
#           for dset_name in a.dsets:
#               data_2D[dset_name] = a.read_full_dset(dset_name)[zloc,:,:]
#           if add_grid:
#               data_2D['x0_B%d' % block_no] = a.x[zloc,:,:]
#               data_2D['x1_B%d' % block_no] = a.y[zloc,:,:]

#           b = SimulationBlock(2, block_number=block_no)
#           g1 = h5f.create_group(b.blockname)
#           nhalos = [5, 5]
#           halo = [[-i for i in nhalos], nhalos]
#           # apply_group_attributes(g1, b)
#           for dset_name in data_2D.keys():
#               dset = g1.create_dataset('%s' % (dset_name), data=data_2D[dset_name])
#               # set_hdf5_metadata(dset, halos=halo, npoints=[OPS_shape[0], OPS_shape[1]], block=b)
#       h5f.close()
