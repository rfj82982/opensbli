import numpy as np
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



ndim = 3
nhalo = 5
name = 'random_nums'
fname = "data.h5"
h5f = h5py.File(fname, 'w')

# Read the number of grid points being used in the C code
file_path = "./opensbli.cpp"
with open(file_path) as f:
    lines = f.readlines()
npoints = []
for no, line in enumerate(lines):
    for i in range(ndim):
        check_string = "block0np%d = " % i
        if check_string in line:
            npoints += [int(line.split('=')[-1].split(';')[0])]

npoints = npoints[::-1]
print("Generating random numbers of size: [Nx, Ny, Nz] = {}".format(npoints))
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


