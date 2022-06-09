import h5py
import numpy as np
from opensbli import *
import os


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


def create_z_coordinates(x, y, full_z):
    inner_shape = x.shape[::-1]
    full_shape = full_z.shape[::-1]
    # Uniform spacing in the span
    dz = Lz / (float(nz))
    print("Grid spacing in z is: {:.5f}".format(dz))
    z_coordinates = [k*dz for k in range(full_z.shape[0] - 2*nhalo)]
    zm = [z_coordinates[0] - k*dz for k in reversed(range(1, nhalo+1))]
    zp = [z_coordinates[-1] + k*dz for k in range(1, nhalo+1)]
    z_coordinates = np.around(np.array(zm + z_coordinates + zp), decimals=10)
    print("Z coordinates including halo points are:", z_coordinates)
    for k in range(full_z.shape[0]):
        z = np.full((full_shape[0], full_shape[1]), z_coordinates[k])
        full_z[k, :, :] = np.transpose(z)
    return full_z


def fill_halo_coordinates(x, y):
    # Create an array with zeros padded around the data
    shape = [nz] + list(x.shape) 
    print(shape)
    new_shape = tuple([shape[i]+ 2*nhalo for i in range(ndim)])
    print("Reversed shape for C-style indexing", new_shape)
    # Full arrays with halo points on the outside
    full_x = np.zeros(new_shape)
    full_y = np.zeros(new_shape)
    full_z = np.zeros(new_shape) 

    # Fill out the interior data
    # Create slices of the interior points to reuse (Nz, Ny, Nz) (they have been transposed into C style)
    x_slice = np.s_[nhalo:new_shape[2] -nhalo]
    y_slice = np.s_[nhalo:new_shape[1] -nhalo]
    # z_slice = np.s_[nhalo:new_shape[0] -nhalo]

    # Filling out the full data
    for k in range(new_shape[0]):        
        full_x[k, y_slice, x_slice] = x
        full_y[k, y_slice, x_slice] = y

    full_z = create_z_coordinates(x, y, full_z)

    for k in range(new_shape[0]):
        # xm, x halos
        full_x[k, y_slice, 4] = x[:,0]
        full_x[k, y_slice, 3] = x[:,1]
        full_x[k, y_slice, 2] = x[:,2]
        full_x[k, y_slice, 1] = x[:,3]
        full_x[k, y_slice, 0] = x[:,4]
        # xm, y halos
        full_y[k, y_slice, 4] = y[:,0]
        full_y[k, y_slice, 3] = y[:,1]
        full_y[k, y_slice, 2] = y[:,2]
        full_y[k, y_slice, 1] = y[:,3]
        full_y[k, y_slice, 0] = y[:,4]
        # xp, x halos
        full_x[k, y_slice, -5] = initial_x[:, 1]
        full_x[k, y_slice, -4] = initial_x[:, 2]
        full_x[k, y_slice, -3] = initial_x[:, 3]
        full_x[k, y_slice, -2] = initial_x[:, 4]
        full_x[k, y_slice, -1] = initial_x[:, 5]
        # xp, y halos
        full_y[k, y_slice, -5] = initial_y[:, 1]
        full_y[k, y_slice, -4] = initial_y[:, 2]
        full_y[k, y_slice, -3] = initial_y[:, 3]
        full_y[k, y_slice, -2] = initial_y[:, 4]
        full_y[k, y_slice, -1] = initial_y[:, 5]

        #ym x halos
        dx = full_x[k, 6, :] - full_x[k, 5, :]
        full_x[k, 4, :] = full_x[k, 5, :] - dx
        full_x[k, 3, :] = full_x[k, 4, :] - dx
        full_x[k, 2, :] = full_x[k, 3, :] - dx
        full_x[k, 1, :] = full_x[k, 2, :] - dx
        full_x[k, 0, :] = full_x[k, 1, :] - dx

        #ym, y halos
        dy = full_y[k, 6, :] - full_y[k, 5, :]
        full_y[k, 4, :] = full_y[k, 5, :] - dy
        full_y[k, 3, :] = full_y[k, 4, :] - dy
        full_y[k, 2, :] = full_y[k, 3, :] - dy
        full_y[k, 1, :] = full_y[k, 2, :] - dy
        full_y[k, 0, :] = full_y[k, 1, :] - dy




    return full_x, full_y, full_z

ndim = 3
# Get the header
fname = '../naca0012_rounded-te_coarse.dat'
f = open(fname)
size = f.readline().split('\n', 1)[0]
nx = int(size.split(' ')[0])
ny = int(size.split(' ')[1])
data = np.genfromtxt(fname, skip_header=1)
initial_x = data[:,0].reshape(nx, ny)
initial_y = data[:,1].reshape(nx, ny)

# Transpose the data to C style
initial_x = np.transpose(initial_x)
initial_y = np.transpose(initial_y)

# Number of halo points
nhalo = 5
nz = 50
Lz = 0.05

# full_x, full_y = np.zeros((nz+2*nhalo, ny+2*nhalo, nx+2*nhalo)), np.zeros((nz+2*nhalo, ny+2*nhalo, nx+2*nhalo))

full_x, full_y, full_z = fill_halo_coordinates(initial_x, initial_y)

print("Full 2D slice size with halos is (Nx, Ny):", full_x.shape)



# # Fill coordinates in the halos over the interace
# # x negative halos, increasing x, copy from the other side over the periodic interface
# full_x[y_slice, 4] = initial_x[:, -2]
# full_x[y_slice, 3] = initial_x[:, -3]
# full_x[y_slice, 2] = initial_x[:, -4]
# full_x[y_slice, 1] = initial_x[:, -5]
# full_x[y_slice, 0] = initial_x[:, -6]

# # y negative halos, copy from the other side over the periodic interface
# full_y[y_slice, 4] = initial_y[:, -2]
# full_y[y_slice, 3] = initial_y[:, -3]
# full_y[y_slice, 2] = initial_y[:, -4]
# full_y[y_slice, 1] = initial_y[:, -5]
# full_y[y_slice, 0] = initial_y[:, -6]



# # x positive halos, increasing x, copy from the other side over the periodic interface
# full_x[y_slice, -5] = initial_x[:, 1]
# full_x[y_slice, -4] = initial_x[:, 2]
# full_x[y_slice, -3] = initial_x[:, 3]
# full_x[y_slice, -2] = initial_x[:, 4]
# full_x[y_slice, -1] = initial_x[:, 5]

# # y positive halos, increasing x, copy from the other side over the periodic interface
# full_y[y_slice, -5] = initial_y[:, 1]
# full_y[y_slice, -4] = initial_y[:, 2]
# full_y[y_slice, -3] = initial_y[:, 3]
# full_y[y_slice, -2] = initial_y[:, 4]
# full_y[y_slice, -1] = initial_y[:, 5]

# Add the z direction


### FINISHED, write the file
# Create the HDF5 file for reading into OpenSBLI
nhalos = [5, 5, 5]
# Output grid file name
fname = "data.h5"
h5f = h5py.File(fname, 'w')


# (Nx, Ny, Nz)
OPS_shape = list(initial_x.shape)[::-1] + [nz]

# Make an OpenSBLI block
b = SimulationBlock(ndim, block_number=0)
g1 = h5f.create_group(b.blockname)
halo = [[-nhalo, -nhalo, -nhalo], [nhalo, nhalo, nhalo]]
apply_group_attributes(g1, b)
print("OpenSBLI block shape without halo points: %s" % OPS_shape)
output = [full_x, full_y, full_z]
for i, out in enumerate(output):
    block_dset_name = b.location_dataset("x%d" % i).base
    dset = g1.create_dataset('%s' % (block_dset_name), data=out)
    set_hdf5_metadata(dset, halos=halo, npoints=OPS_shape, block=b)


h5f.close()
