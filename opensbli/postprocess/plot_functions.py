import numpy as np
import h5py
import matplotlib.pyplot as plt
import glob, os, time, re

class OpenSBLIPreProcess(object):
    """ Commonly used plotting routines in OpenSBLI."""
    def __init__(self):
        # self.file_name = file_name
        self.open_files = []
        self.gamma = 1.4
        return

    def read_file(self, file_name, blocknumber, remove_halos=True):
        print("Reading from file: %s" % file_name)
        f = h5py.File(file_name, 'r')
        block_name = list(f.keys())[blocknumber]
        dsets = list(f[block_name].keys())
        print("Found %d datasets: %s, with dimensions: %s" % (len(dsets), dsets, f[block_name][dsets[0]].shape[::-1]))
        self.f, self.block_name, self.dsets = f, block_name, dsets
        self.nhalos = np.abs(f[block_name][dsets[0]].attrs['d_m'])
        self.shape = list(f[block_name][dsets[0]].shape)
        self.shape = tuple([x-10 for x in self.shape])
        return

    def read_grid(self, blocknumber, file_name='./data.h5', remove_halos=True):
        print("Reading from file: %s" % file_name)
        f = h5py.File(file_name, 'r')
        block_name = list(f.keys())[blocknumber]
        dsets = list(f[block_name].keys())
        print("Found %d datasets: %s, with dimensions: %s" % (len(dsets), dsets, f[block_name][dsets[0]].shape[::-1]))
        self.nhalos = np.abs(f[block_name][dsets[0]].attrs['d_m'])
        self.shape = list(f[block_name][dsets[0]].shape)
        self.shape = tuple([x-10 for x in self.shape])
        self.x = self.remove_halos(f[block_name]['x0'+'_B%d' % blocknumber])
        self.y = self.remove_halos(f[block_name]['x1'+'_B%d' % blocknumber])
        self.z = self.remove_halos(f[block_name]['x2'+'_B%d' % blocknumber])
        self.blocknumber = blocknumber
        return

    def find_files(self, directory):
        """ Finds a list of OpenSBLI HDF5 files from a specified directory."""
        file_list = sorted(glob.glob(directory + '/opensbli_output_*.h5'))
        iteration_numbers = [re.findall("\d+", s)[0].lstrip('0') for s in file_list]
        print("Found {:} OpenSBLI output files:".format(len(file_list)))
        for f in file_list:
            print(f)
        return file_list, iteration_numbers

    def NaN_check(self, dset):
        data = self.f[self.block_name][dset]
        NaN = np.isnan(np.sum(data))
        if NaN:
            print("NaN detected in dataset: %s" % dset)
        else:
            print("No NaN detected.")
        return

    def MinMax(self):
        for dset in self.dsets:
            data = self.read_full_dset(dset)
        return

    def remove_halos(self, dataset):
        size = dataset.shape
        read_start = [abs(d) for d in self.nhalos]
        read_end = [s-abs(d) for d, s in zip(self.nhalos, size)]
        if len(read_end) == 1:
            read_data = dataset[read_start[0]:read_end[0]]
        elif len(read_end) == 2:
            read_data = dataset[read_start[0]:read_end[0], read_start[1]:read_end[1]]
        else:
            read_data = dataset[read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
        return read_data

    def read_full_dset(self, dset, remove_halos=True, min_max=True):
        if remove_halos:
            data = self.remove_halos(self.f[self.block_name][dset])
        else:
            data = self.f[self.block_name][dset]
        print("Reading dataset: {:}, Min: {:.3f}, Max: {:.3f}".format(dset, np.min(data), np.max(data)))
        return data


    def x_slice(self, dset, location, remove_halos=True):
        """ Only used for 3D arrays, slicing in 3 directions without loading the entire array."""
        data = self.f[self.block_name][dset][:,:,location]
        if remove_halos:
            data = self.remove_halos(data)
        return data

    def y_slice(self, dset, location, remove_halos=True):
        """ Only used for 3D arrays, slicing in 3 directions without loading the entire array."""
        data = self.f[self.block_name][dset][:,location,:]
        if remove_halos:
            data = self.remove_halos(data)
        return data

    def z_slice(self, dset, location, remove_halos=True):
        """ Only used for 3D arrays, slicing in 3 directions without loading the entire array."""
        data = self.f[self.block_name][dset][location,:,:]
        if remove_halos:
            data = self.remove_halos(data)
        return data

    def add_flow_attributes(self):
        """ Add the numerical values of flow paramters to the HDF5 file."""
        return


    def kinetic_energy(self, conservative=True):
        gamma = 1.4
        if conservative:
            rho, rhou, rhov, rhow, rhoE = self.read_full_dset('rho_B%d' % self.blocknumber), self.read_full_dset('rhou0_B%d' % self.blocknumber), self.read_full_dset('rhou1_B%d' % self.blocknumber), self.read_full_dset('rhou2_B%d' % self.blocknumber), self.read_full_dset('rhoE_B%d' % self.blocknumber)
            u, v, w = rhou/rho, rhov/rho, rhow/rho
            KE = np.sum(0.5*rho*(u**2 + v**2 + w**2)) / (self.shape[0]*self.shape[1]*self.shape[2])
        return KE

    def pressure(self, conservative=True):
        rho, rhou, rhov, rhow, rhoE = self.read_full_dset('rho_B%d' % self.blocknumber), self.read_full_dset('rhou0_B%d' % self.blocknumber), self.read_full_dset('rhou1_B%d' % self.blocknumber), self.read_full_dset('rhou2_B%d' % self.blocknumber), self.read_full_dset('rhoE_B%d' % self.blocknumber)
        u, v, w = rhou/rho, rhov/rho, rhow/rho
        p = (self.gamma-1)*(rhoE - 0.5*rho*(u**2 + v**2 + w**2))
        return p

    def GlobalMach(self, output='Max', conservative=True):
        rho, rhou, rhov, rhow, rhoE = self.read_full_dset('rho_B%d' % self.blocknumber), self.read_full_dset('rhou0_B%d' % self.blocknumber), self.read_full_dset('rhou1_B%d' % self.blocknumber), self.read_full_dset('rhou2_B%d' % self.blocknumber), self.read_full_dset('rhoE_B%d' % self.blocknumber)
        u, v, w = rhou/rho, rhov/rho, rhow/rho
        p = (self.gamma-1)*(rhoE - 0.5*rho*(u**2 + v**2 + w**2))
        c = np.sqrt(self.gamma*p/rho)
        Mach = np.sqrt(u**2 + v**2 + w**2) / c
        if output == 'Max':
            argmax = np.argmax(Mach)
            Mach = np.max(Mach)
            print("Maximum Mach number is: ", Mach)
            print("Maximum is located at (x,y,z):", np.ravel(self.x)[argmax], np.ravel(self.y)[argmax], np.ravel(self.z)[argmax])
        return Mach


class OpenSBLIPlot(object):
    """ Commonly used plotting routines in OpenSBLI."""
    def __init__(self):
        # self.file_name = file_name
        self.open_files = []
        self.nlevels = 30
        return

    def create_figure(self):
        fig, ax = plt.subplots()
        self.fig, self.ax = fig, ax
        return fig, ax

    def get_figure(self):
        return self.fig, self.ax

    def set_labels(self, ax, direction):
        if direction == 'xy':
            ax.set_xlabel('x')
            ax.set_ylabel('y')
        elif direction == 'xz':
            ax.set_xlabel('x')
            ax.set_ylabel('z')
        elif direction == 'zy':
            ax.set_xlabel('z')
            ax.set_ylabel('y')
        else:
            raise NotImplementedError("Direction should be xy, xz, or zy.")
        return ax


    def simple_imshow_plot(self, data, var_label):
        import cmocean
        cmap = cmocean.cm.balance
        self.ax.imshow(data, cmap=cmap)
        self.fig
        return

    def simple_contour_plot(self, dset, coordinates, count_plot=0):
        import cmocean
        cmap = cmocean.cm.balance
        x, y = coordinates[0], coordinates[1]

        if count_plot == 0:
            dset_im = self.ax.contourf(x, y, dset.T, self.nlevels, cmap=cmap)
        return dset_im

    def line_plot(self, xvar, yvar):
        self.ax.plot(xvar, yvar)
        return

    def save_figure(self, fname, dpi=300):
        plt.savefig('%s.png' % fname, dpi=dpi, bbox_inches='tight')
        return



# class OpenSBLIProcessNS(object):
#     """ Post processing calculations from an OpenSBLI output file, using NumPy."""
#     def __init__(self):
#         # self.file_name = file_name
#         self.open_files = []
#         self.nlevels = 30
#         return

#     def 

#     def kinetic_energy(self, conservative=True):
#         if conservative:


#         return KE
