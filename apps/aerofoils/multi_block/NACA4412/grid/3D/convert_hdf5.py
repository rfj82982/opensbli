# Converts the SBLI 2d block grid into HDF5 format. For reading into OpenSBLI
# Structured mesh, sharp trailing edge version.

import numpy as np
import h5py
from opensbli import *

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


# Specify the input grid files
input_files = ["../Bl1.dat", "../Bl2.dat","../Bl3.dat"]
# Number of halo points to add on each side of each direction (default 5)
nhalo = [5, 5, 5] # was 5, 5, 2
# Output grid file name
fname = "data.h5"
h5f = h5py.File(fname, 'w')

# Number of points in the periodic span.
nz = 5
# Grid spacing
dz = 0.02

sharp_TE = True

total_grid_points = 0

# Loop over all of the grid points
for block_number, block in enumerate(input_files):
    print("\n\n\nReading from %s." % block)
    f = open(block)
    nx,ny = map(int, f.readlines()[0].split())
    print("Nx, Ny from the input file for block %d" % block_number, nx, ny)
    f.close()
    # Read the data
    x,y,z = np.loadtxt(block, skiprows =1, unpack=True)
    x = x.reshape(nx, ny)
    y = y.reshape(nx, ny)
    
    # Sharp trailing edge -> take away one point at the start and end of the grid in x, from block 2
    if block_number == 1:
        if sharp_TE:
            print("Taking off 2 columns in x direction for the sharp trailing edge.")
            x = x[1:-1,:]
            y = y[1:-1,:]
    shape = list(x.shape) +[nz]
    total = shape[0]*shape[1]*shape[2]
    print("Block %d has %e grid points." % (block_number, int(total)))
    total_grid_points +=  total
    print("Original 3D shape: %s" % shape)
    new_shape = tuple(reversed([shape[i]+ 2*nhalo[i] for i in range(3)]))
    print("Reversed shape for C-style indexing", new_shape)
    #exit()
    newx = np.zeros(new_shape)
    newy = np.zeros(new_shape)
    newz = np.zeros(new_shape)
    for k in range(nz + 2*nhalo[2]):
        zloc = dz * float(k - nhalo[2])
        # print(zloc)
        z = np.full(x.shape, zloc)
        #print z.shape
        newx[k,nhalo[1]:new_shape[1] -nhalo[1], nhalo[0]:new_shape[2] -nhalo[0]] = np.transpose(x)
        newy[k,nhalo[1]:new_shape[1] -nhalo[1], nhalo[0]:new_shape[2] -nhalo[0]] = np.transpose(y)
        newz[k,nhalo[1]:new_shape[1] -nhalo[1], nhalo[0]:new_shape[2] -nhalo[0]] = np.transpose(z)

    # Make an OpenSBLI block
    b = SimulationBlock(3, block_number=block_number)
    g1 = h5f.create_group(b.blockname)
    halo = [[-i for i in nhalo], nhalo]
    apply_group_attributes(g1, b)
    block_dset_name = b.location_dataset("x0").base
    print("OpenSBLI block shape without halo points: %s" % shape)

    # Create x coordinates
    dset = g1.create_dataset('%s' % (block_dset_name), data=newx)
    set_hdf5_metadata(dset, halos=halo, npoints=[shape[0], shape[1], nz], block=b)
    # Create y coordinates
    block_dset_name = b.location_dataset("x1").base
    dset = g1.create_dataset('%s' % (block_dset_name), data=newy)
    set_hdf5_metadata(dset, halos=halo, npoints=[shape[0], shape[1], nz], block=b)
    # Create z coordinates
    block_dset_name = b.location_dataset("x2").base
    dset = g1.create_dataset('%s' % (block_dset_name), data=newz)
    set_hdf5_metadata(dset, halos=halo, npoints=[shape[0], shape[1], nz], block=b)
    
    print("Length in x for block %d:" % block_number, abs(np.amin(x) - np.amax(x)))
    print("Length in y for block %d:" % block_number, abs(np.amin(y) - np.amax(y)))
    print("Length in z for block %d:" % block_number, zloc)
print("Total grid points: %g" % total_grid_points)

h5f.close()
