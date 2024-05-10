import numpy as np
from opensbli import *
# Random number generation for the initial condition
from opensbli.utilities.helperfunctions import output_hdf5
ndim = 3
block = SimulationBlock(ndim=ndim)
name = 'random_nums'

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

print("Generating random numbers of size: [Nx, Ny, Nz] = {}".format(npoints))
halos = [(-5, 5) for _ in range(ndim)]
size_including_halo = [npoints[i] + sum(np.absolute(halos[i])) for i in range(ndim)]
# Generate the initial white noise seeding
random_numbers = np.random.rand(*size_including_halo) - 0.5
print("Random numbers: Min: {:.3f}, Max: {:.3f}, Mean: {:.3f}".format(np.min(random_numbers), np.max(random_numbers), np.mean(random_numbers)))
output_hdf5(random_numbers, name, halos, npoints, block)