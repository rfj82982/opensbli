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

already_processed = sorted(glob.glob(outPath + '/opensbli_output_*.h5'))
already_processed = [re.findall("\d+", s)[0].lstrip('0') for s in already_processed]
# Extract data
a = OpenSBLIPreProcess()
a.read_grid()
fnames, iters = a.find_files(inPath)

# Slice location
zloc = int(a.shape[0]/2.0)
add_grid = True

for i, f in enumerate(fnames):
	if iters[i] not in already_processed:
		data = a.read_file(f)
		# Extract a single slice
		data_2D = {}
		# Extract all the variables at this slice
		for dset_name in a.dsets:
			data_2D[dset_name] = a.read_full_dset(dset_name)[zloc,:,:]
		if add_grid:
			data_2D['x0_B0'] = a.x[zloc,:,:]
			data_2D['x1_B0'] = a.y[zloc,:,:]
		# Write to a new 2D HDF5 file
		h5f = h5py.File(outPath + 'opensbli_output_%s.h5' % iters[i], 'w')
		b = SimulationBlock(2, block_number=0)
		g1 = h5f.create_group(b.blockname)
		nhalos = [5, 5]
		halo = [[-i for i in nhalos], nhalos]
		# apply_group_attributes(g1, b)
		for dset_name in data_2D.keys():
			dset = g1.create_dataset('%s' % (dset_name), data=data_2D[dset_name])
			# set_hdf5_metadata(dset, halos=halo, npoints=[OPS_shape[0], OPS_shape[1]], block=b)
		h5f.close()
