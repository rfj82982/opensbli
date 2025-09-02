""" Vortex plotting script in 3D."""
import numpy as np
import numexpr as ne
import h5py
import matplotlib.pyplot as plt
import os, glob, re
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import multiprocessing
import time
import matplotlib.cm as mpl_cm
import matplotlib.ticker as ticker
import warnings
import shutil


figure_dpi = 300
def get_levels(data,nvar):
    """ Sets contour levels continuous around zero."""
    print("Variable: {}, Min: {:.3g}, Max: {:.3g}".format(nvar, np.min(data), np.max(data)))
    abs_max = np.max(np.sqrt(data*data))
    if nvar == 'passive-scalar' or nvar == 'Mach':
        vmin = 0
        vmax = abs_max
        cpal = "jet"
    elif nvar == 'T':
        vmin = 0.95
        vmax = 1.05
        cpal = 'jet'
    elif nvar == 'omega_z':
        vmin = -10
        vmax = 10
        cpal = "jet_r"
    else:
        cpal = "RdBu_r"
        vmin = -abs_max
        vmax =  abs_max
    vstep =   1
    levels = np.arange(vmin, vmax+vstep, vstep)
    cmap_cont = mpl_cm.get_cmap(cpal)
    return vmin, vmax, cmap_cont

def set_label(f, f_idx, nvar):
    """ Format the text to be added to the plot."""
    if nvar == "rho":
        label = "Density"
    elif nvar == "vortz":
        label = "Vorticity"
    elif nvar == "p":
        label = 'Pressure'
    elif nvar == 'u':
        label = 'u Velocity'
    elif nvar == 'v':
        label = 'v Velocity'
    elif nvar == 'w':
        label = 'w Velocity'
    elif nvar == 'Mach':
        label = 'Mach'
    #label += ' : Mode: {} : Freq: {:.3f}'.format(f_idx, f)
    label += ' : $St=$ {:.3f}'.format(f)
    return label

def format_side(ax, CS, label):
    """ Labels and formatting options for side view plots."""
    add_color_bar = True
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    ax.set_ylim([-10, 10])
    plt.gca().set_aspect('equal')
    if add_color_bar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(CS, cax=cax)
        cbar.ax.set_ylabel(r'$%s$' % label)
    plt.tight_layout()
    return

def create_fig():
    """ Control figure size."""
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6.4, 4.8))
    return fig, ax

def read_dataset(file, dataset):
    group = file["opensbliblock00"]
    d_m = group["%s" % (dataset)].attrs['d_m']
    size = group["%s" % (dataset)].shape
    start=[abs(d) for d in d_m]
    end=[s-abs(d-1) for d, s in zip(d_m, size)]
    read_data=group["%s" % (dataset)][start[0]:end[0],start[1]:end[1],start[2]:end[2]]
    return read_data

def plot_instant(packed):
    index, data, nvar = packed[0], packed[1], packed[2]
    if nvar == 'p':
        data = data - 1.0/(gama*Minf**2.0)
    # Set label for text
    # label = set_label(f, f_idx, nvar)
    vmin, vmax, cmap_cont = get_levels(data, nvar)
    # print("Plotting mode: {} with frequency: {:.3f}, phase: {} and ID: {}, data min/max: {:.3f}, {:.3f}".format(mode_idx, f, phase, f_idx, np.min(data_phased), np.max(data_phased)))

    fig, ax = create_fig()
    CS = ax.pcolormesh(x0, x1, data, cmap=cmap_cont, vmin=vmin,vmax=vmax, shading='auto')
    format_side(ax, CS, nvar)
    plt.savefig(output_dir + 'vortex_core_{}_{}'.format(nvar, index), dpi=figure_dpi, bbox_inches='tight')
    plt.clf()
    return

def find_files(directory):
    """ Finds a list of OpenSBLI HDF5 files from a specified directory."""
    file_list = sorted(glob.glob(directory + '/opensbli_output_*.h5'))
    iteration_numbers = [re.findall("\d+", s)[0].lstrip('0') for s in file_list]
    print("Found {:} OpenSBLI output files:".format(len(file_list)))
    print(file_list[0], ",......,", file_list[-1])
    return file_list, iteration_numbers

print('Reading data')

input_dir = './'
output_dir = './images/'
GIF_dir = './GIFs/'
try:
    os.mkdir(output_dir)
except FileExistsError:
    pass
try:
    os.mkdir(GIF_dir)
except FileExistsError:
    pass

fnames, iterations = find_files(input_dir)
simulation_times = []

for index, fname in enumerate(fnames):
    ff=h5py.File(fname, 'r')
    if index == 0:
        # Get the constants
        dt, Minf, Re, gama = ff['dt'][0], ff['Minf'][0], ff['Re'][0], ff['gama'][0]
        Nx, Ny, Nz = ff['block0np0'][0], ff['block0np1'][0], ff['block0np2'][0]
        y_0, y_1, Lx = ff['y_0'][0], ff['y_1'][0], ff['Lx'][0]
        stretch = ff['stretch'][0]
        # Shear-layer parameters.

        # Mixing layer configuration
        print("Processing case with parameters: dt: {:.4f}, Re: {:.4f}, Minf: {:.4f}".format(dt, Re, Minf))
        print("Grid info: (Nx, Ny, Nz) = ({:d}, {:d}, {:d}), Stretching factor: {:.4f}".format(Nx, Ny, Nz, stretch))

        # Get the grid spacing within the shear-layer
        zloc = int(Nz/2.0)
        dx = read_dataset(ff,'x0_B0')[zloc,int(Ny/2),1] - read_dataset(ff,'x0_B0')[zloc,int(Ny/2),0]
        dy = np.abs(read_dataset(ff,'x1_B0')[zloc,int(Ny/2),1] - read_dataset(ff,'x1_B0')[zloc,int(Ny/2)+1,1])
        dz = read_dataset(ff,'x2_B0')[1,int(Ny/2),0] - read_dataset(ff,'x2_B0')[0,int(Ny/2),0]
        print("Grid resolution at centreline: delta_x: {:3f}, delta_y: {:3f}, delta_z: {:3f}".format(dx, dy, dz))
        print("Shear-layer parameters: y_0: {:.4f}, y_1: {:.3f}, L_x: {:.1f}, y_0/L_x: {:.4f}".format(y_0, y_1, Lx, y_0/Lx))
        edge_index = np.abs(read_dataset(ff,'x1_B0')[zloc,:,1] - 2.5).argmin()
        print(edge_index)
        dy_edge = np.abs(read_dataset(ff,'x1_B0')[zloc,edge_index+1,1] - read_dataset(ff,'x1_B0')[zloc,edge_index,1])
        print(dy_edge)
        #exit()
        # Load a plane of the grid for plotting
        x0 = read_dataset(ff,'x0_B0')[zloc,:,:]
        x1 = read_dataset(ff,'x1_B0')[zloc,:,:]
        x2 = read_dataset(ff,'x2_B0')[zloc,:,:]
        # exit()

    # Check if file already processed
    file_check = output_dir + 'vortex_core_omega_z_%d.png' % index
    if not os.path.isfile(file_check):
        # Get simulation time
        simulation_times.append(ff['simulation_time'][0])
        print("Processing file: {} at time {:.3f}".format(fname, simulation_times[-1]))


        rho = read_dataset(ff, 'rho_B0')[zloc,:,:]
        rhou = read_dataset(ff, 'rhou0_B0')[zloc,:,:]
        rhov = read_dataset(ff, 'rhou1_B0')[zloc,:,:]
        rhow = read_dataset(ff, 'rhou2_B0')[zloc,:,:]
        rhoE = read_dataset(ff, 'rhoE_B0')[zloc,:,:]
        # Post process quantities
        # wx = read_dataset(ff, 'wx_B0')[zloc,:,:]
        # wy = read_dataset(ff, 'wy_B0')[zloc,:,:]
        wz = read_dataset(ff, 'wz_B0')[zloc,:,:]
        # divV = read_dataset(ff, 'divV_B0')[zloc,:,:]
        dudy = read_dataset(ff, 'dudy_B0')[zloc,:,:]

        print('Calculating variables')
        u = ne.evaluate('rhou/rho')
        v = ne.evaluate('rhov/rho')
        w = ne.evaluate('rhow/rho')
        e = ne.evaluate('rhoE/rho-0.5*(u**2+v**2+w**2)')
        p = ne.evaluate('(gama-1.0)*rho*e')
        # T = ne.evaluate('e*Minf**2*gama*(gama-1.0)')
        #f = ne.evaluate('rhof/rho')
        a = np.sqrt(gama*p/rho)
        Mach = np.sqrt(u**2 + v**2 + w**2) / a
        # mu = T**0.7

        variables = ['Mach', 'p', 'u', 'v', 'w', 'omega_z', 'dudy']
        nvariables = len(variables)
        # raw_data = [rho, Mach, T, p, u, v, w, wx, wy, wz, divV, dudy, mu*(wx**2+wy**2+wz**2)]
        raw_data = [Mach, p, u, v, w, wz, dudy]
        #raw_data = [u, v, wz]
        # raw_data = [dset[zloc,:,:] for dset in raw_data]
        index = [index for _ in range(len(variables))]
        packed_data = zip(index, raw_data, variables)
        # Plot images in parallel
        with multiprocessing.Pool(processes=nvariables) as pool:
            pool.map(plot_instant, packed_data)
        # # Wait for all tasks to finish processing
        barrier = multiprocessing.Barrier(1) # wait one second
    else:
        print("Skipping file: {}".format(index))


make_GIFs = True
# # Animate the images into a gif
if make_GIFs:
    os.chdir(GIF_dir)
        
    print("Making animated GIFs")
    for var_name in variables:
        print("Making GIF for: {}".format(var_name))
        # for i in range(0,len(fnames)):
        os.system('convert -delay 20 -dispose previous -loop 0 $(ls ../images/vortex_core_'+var_name+'_*.png | sort -V) vortex_core_'+var_name+'_'+'.gif')
        # shutil.move("vortex_core_"+var_name+"_"+".gif", )
