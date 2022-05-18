""" File to change the array declarations in OpenSBLI from initialising to zero, to restarting from a restart file.
Updated to work with multi-block 05/2022 (djl)."""
import shutil

# Number of simulation blocks in the problem
nblocks = 3
restart_grid = False
# Solution vector to restart
Q = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
if restart_grid:
    coords = ['x0', 'x1', 'x2']
# Backup the original array declarations
filename = 'defdec_data_set.h'
original = 'init_defdec.h'
restart_file_name = 'restart.h5'
shutil.copyfile(filename, original)


for i in range(nblocks):
    arrays = ['%s_B%d' % (x, i) for x in Q]
    print(arrays)
    for variable in arrays:
        f = open(filename, 'r')
        file_data = f.read()
        f.close()
        text_to_search = 'ops_decl_dat(opensbliblock0{:}, 1, size, base, halo_m, halo_p, value, "double", "{:}");'.format(i, variable)
        replacement_text = 'ops_decl_dat_hdf5(opensbliblock0{:}, 1, "double", "{:}", "{:}");'.format(i, variable, restart_file_name)
        newdata = file_data.replace(text_to_search, replacement_text)
        f = open(filename, 'w')
        f.write(newdata)
        f.close()
    arrays = []