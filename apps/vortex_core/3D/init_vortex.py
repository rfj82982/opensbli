import numpy as np
from numpy import sinh
from opensbli import *
import h5py

# Random number generation for the initial condition
def apply_group_attributes(group, block):
    group.attrs.create("dims", [block.ndim], dtype="int32")
    group.attrs.create("ops_type", u"ops_block",dtype="S8")
    group.attrs.create("index", [block.blocknumber], dtype="int32")
    return

def set_hdf5_metadata(dset, halos, npoints, block):
    """ Function to set hdf5 metadata required by OPS to a dataset. """
    d_m = halos[0]
    d_p = halos[1]

    dset.attrs.create("d_p", d_p, dtype="int32")
    dset.attrs.create("d_m", d_m, dtype="int32")
    dset.attrs.create("dim", [1], dtype="int32")
    dset.attrs.create("ops_type", u"ops_dat",dtype="S10")
    dset.attrs.create("block_index", [block.blocknumber], dtype="int32")
    dset.attrs.create("base", [0 for i in range(block.ndim)], dtype="int32")
    dset.attrs.create("type", u"double",dtype="S15")
    dset.attrs.create("block", u"%s" % block.blockname,dtype="S25")
    dset.attrs.create("size", npoints, dtype="int32")
    return

def output_hdf5(array, array_name, halos, npoints, block):
    """ Creates an HDF5 file for reading in data to a simulation, 
    sets the metadata required by the OPS library. """
    if not isinstance(array, list):
        array = [array]
    if not isinstance(array_name, list):
        array_name = [array_name]
    assert len(array) == len(array_name)
    with h5py.File('data.h5', 'w') as hf:
        # Create a group
        if (isinstance(block, MultiBlock)):
            all_blocks = block.blocks
        else:
            all_blocks = [block]
        for b in all_blocks:
            g1 = hf.create_group(b.blockname)
            # Loop over all the dataset inputs and write to the hdf5 file
            for ar, name in zip(array, array_name):
                g1.attrs.create("dims", [b.ndim], dtype="int32")
                g1.attrs.create("ops_type", u"ops_block",dtype="S9")
                g1.attrs.create("index", [b.blocknumber], dtype="int32")
                block_dset_name = b.location_dataset(name).base
                dset = g1.create_dataset('%s' % (block_dset_name), data=ar)
                set_hdf5_metadata(dset, halos, npoints, b)
    return

## User inputs
# Problem parameters
U_0 = 1.0 # Freestream velocity
y_1 = 0.1 # Amplitude of initial wave
# Length in width and axial direction
Lx = 10.0 # Width of vortex sheet
Lz = 1.0 # Axial length
stretch = 4.0 # Stretch factor in y
Ly = 40.0 # Length in y
######### Set the number of points within the initial vorticity thickness: uniform Delta x ~ Delta y ~ Delta z < y_0*grid_factor ###########
grid_factor = 1/20 # grid spacing is this factor smaller than y_0
### Other conditions to set
Minf = 0.4
# Select which way to set the thickness
vary_y0 = False

# Step 1: Select y_0 for fixed Lx, or set y_0 based on the ratio of Reynolds numbers
if vary_y0:
    print("Manually setting the y_0 vorticity thickness for a fixed width Lx.")
    y_0 = 1
    Re_0 = y_0*U_0 # divided by nu
    Re = 2*Lx*U_0 # divided by nu
    Ratio = Re/Re_0
else: # Set based on ratio of Reynolds number definitions
    print("Selecting the y_0 value based on the ratio of the Reynolds number definitions.")
    # y_0 * U_0 / nu
    Re_0 = 10**4
    # 2*U_0*Lx
    Re = 10**6
    Ratio = Re/Re_0
    y_0 = 2*Lx / Ratio

# Ratio = 2*U_0*Lx / ( y_0*U_0 ) = 2*Lx/y_0 -> y_0 = 2*Lx / Ratio
print("Re ratio is: {}".format(Ratio))
# Reynolds number to set in the code
nu_factor = np.round(1.0 / (y_0*U_0 / Re_0), 0)
print("Initial vorticity thickness y_0 is: {}, for 1/nu (Re in the code) of: {}".format(y_0, nu_factor))
# dx = Lx / (Nx - 1)
# Width of the domain
Nx = round(Lx / (y_0*grid_factor) + 1)
dx = Lx / (Nx - 1)
print("Nx and dx are: {}, {}".format(Nx, dx))
# Length of the vortex
Nz = round(Lz / (y_0*grid_factor) + 1)
dz = Lz / (Nz - 1)
print("Nz and dz are: {}, {}".format(Nz, dz))
# Calculate Ny we need to maintain the same resolution to edge of y0/2

# Find y value at jth grid point
y_func = lambda Ny, j : 0.5*Ly*sinh(stretch*(j-(Ny-1)/2)/((Ny-1)/2))/sinh(stretch)
target_dy = dx
print("Target dy is: {}".format(target_dy))
# Find which Ny value to use for fixed Ly, stretch
Ntest = [i for i in range(3, 10000, 2)]
for Ny in Ntest:
    centreline = int(Ny/2)
    j = int(centreline + int(0.5*(1/grid_factor)))
    dy = y_func(Ny, j+1) - y_func(Ny, j)
    # print("Ny = {}, Centreline = {}, j = {}, dy = {}.".format(Ny, centreline, j, dy))
    if np.isclose(dy, dx, rtol=dx):
        Ny += 2
        print("Selected dy: {:.5f} which requires Ny: {}".format(dy, Ny))
        break
print("Final grid distribution is: (Nx, Ny, Nz) = {}, for total N = {:.3e} grid points.".format((Nx, Ny, Nz), Nx*Ny*Nz))

# Substitute these values to the correct place in the C code
subs_dict = {
'Lx' : Lx,
'Ly' : Ly,
'Lz' : Lz,
'stretch' : stretch,
'block0np0' : Nx,
'block0np1' : Ny,
'block0np2' : Nz,
'y_0'  : y_0,
'y_1' : y_1,
'Re' : nu_factor,
'U_0' : U_0,
'Minf' : Minf,
}

print("Vortex input parameter summary: {}".format(subs_dict))
# Generate random numbers input file
ndim = 3
nhalo = 5
name = 'random_nums'
fname = "data.h5"
h5f = h5py.File(fname, 'w')
npoints = [Nx, Ny, Nz]
print("Generating random numbers of size: [Nx, Ny, Nz] = {}".format(npoints))
npoints = npoints[::-1]
size_including_halo = [npoints[i] + nhalo*2 for i in range(ndim)]
# Generate the initial white noise seeding
random_numbers = np.random.rand(*size_including_halo) - 0.5
print("Random numbers: Min: {:.3f}, Max: {:.3f}, Mean: {:.3f}".format(np.min(random_numbers), np.max(random_numbers), np.mean(random_numbers)))
halo = [[-nhalo for _ in range(ndim)]] + [[nhalo for _ in range(ndim)]]
# Make an OpenSBLI block
b = SimulationBlock(ndim, block_number=0)
g1 = h5f.create_group(b.blockname)
apply_group_attributes(g1, b)
block_dset_name = b.location_dataset("%s" % name).base
dset = g1.create_dataset('%s' % (block_dset_name), data=random_numbers)
npoints = npoints[::-1]
set_hdf5_metadata(dset, halos=halo, npoints=[npoints[0], npoints[1], npoints[2]], block=b)
h5f.close()

# Update the simulation parameters in the C code
# Read the number of grid points being used in the C code
file_path = "./opensbli.cpp"
with open(file_path) as f:
    lines = f.readlines()

for k, v in subs_dict.items():
    for no, line in enumerate(lines[0:100]):
        if '%s = ' % k in line:
            split_line = line.split()
            old_value = split_line[-1]
            # Update the value
            split_line[-1] = str(v) + ';' + '\n'
            lines[no] = ' '.join(split_line)
            # print("Quantity: {} - Before: {} After: {}".format(k, old_value, split_line[-1]))
            break # Update only first instance of the constant
# Update the simulation code file
with open(file_path, 'w') as output_file:
    output_file.writelines(lines)

