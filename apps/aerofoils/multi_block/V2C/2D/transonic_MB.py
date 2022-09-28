#!/usr/bin/env python
from opensbli import *
from sympy import sin, exp, pi, tan, cos
import copy
from opensbli.multiblock.algorithm import TraditionalAlgorithmRKMB
from opensbli.postprocess.airfoil import *
from sympy.functions.elementary.piecewise import Piecewise, ExprCondPair
import os
# Disable the gmpy library for this case to avoid deepcopy issues
os.environ['MPMATH_NOGMPY'] = '1'

import itertools
def create_exchange_calls_codes(multiblock_descriptor, dsets):
    kernels = []
    for block in multiblock_descriptor.blocks:
        arrays = [block.location_dataset(a) for a in flatten(dsets)]
        kernels += block.apply_interface_bc(arrays, multiblock_descriptor, full_halo_swap=True)
    return kernels

ndim = 2
nblocks = 3
# Set non-conservative LHS to reduce array storage
conservative = True
multi_block = MultiBlock(ndim, nblocks, conservative=conservative)
SimulationDataType.set_datatype(Double)

# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "RefT", "SuthT"]
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)
#Create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs
Einstein_expansion = EinsteinEquation()
Einstein_expansion.optional_subs_dict = optional_subs_dict
metric_vel = "Eq(U_i, D_i_j*u_j)"
eqns = Einstein_expansion.expand(metric_vel, ndim, coordinate_symbol, [], constants)
for eq in eqns:
    Einstein_expansion.optional_subs_dict[eq.lhs] = eq.rhs


# NS = NS_Split('Kennedy_Gruber', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='constant')
NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic', energy_formulation='enthalpy', debug=False)

mass, momentum, energy = NS.mass, NS.momentum, NS.energy
# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)

# Formulas for the variables used in the equations
constituent = ConstituentRelations()
if conservative:
    velocity = "Eq(u_i, rhou_i/rho)"
    eqns = Einstein_expansion.expand(velocity, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)
    pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
    enthalpy = "Eq(H, (rhoE + p) / rho)"
else:
    pressure = "Eq(p, (gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j)))"
    enthalpy = "Eq(H, Et + p / rho)"

temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

eqns = Einstein_expansion.expand(pressure, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
eqns = Einstein_expansion.expand(temperature, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
eqns = Einstein_expansion.expand(viscosity, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
eqns = Einstein_expansion.expand(enthalpy, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)

# Transform the equations into curvilinear form
simulation_eq.apply_metrics(metriceq)
# Specify the numerical schemes
schemes = {}
rk = RungeKuttaLS(3, formulation='SSP')
schemes[rk.name] = rk
# cent = Central(4)
cent = StoreSome(4, 'u0 u1 T')
schemes[cent.name] = cent
multi_block.set_discretisation_schemes(schemes)

# Initial conditions
## Need to change for swept cases
conserve_vector = flatten(simulation_eq.time_advance_arrays)
d, u0, u1, p = symbols("d, u0:2, p", **{'cls':GridVariable})
gama, Minf = symbols("gama, Minf", **{'cls':ConstantObject})
initial_equations = []
initial_equations += [Eq(d, 1.0)]
initial_equations += [Eq(u0, 1.0)]
initial_equations += [Eq(u1, 0.0)]
initial_equations += [Eq(p, 1.0/(gama*Minf**2.0))]

# Set the q vector values for initial condition and farfield boundaries
q_vector = flatten(simulation_eq.time_advance_arrays)
if conservative:
    initial_equations += [Eq(q_vector[0], d)]
    initial_equations += [Eq(q_vector[1], d*u0)]
    initial_equations += [Eq(q_vector[2], d*u1)]
    initial_equations += [Eq(q_vector[3], p/(gama-1.0) + 0.5* d *(u0**2+u1**2))]
else:
    initial_equations += [Eq(q_vector[0], d)]
    initial_equations += [Eq(q_vector[1], u0)]
    initial_equations += [Eq(q_vector[2], u1)]
    initial_equations += [Eq(q_vector[3], p/(d*(gama-1.0)) + 0.5*(u0**2+u1**2))]

initial = GridBasedInitialisation()
initial.add_equations(copy.deepcopy(initial_equations))
multi_block.set_equations([initial])

# block 0 boundary conditions
mb_bcs = {0:None, 1:None, 2:None}
# Boundary conditions for block 0
# The boundary conditions are [InterfaceBC, outflow] in x0 direction and [SharedInterfaceBC, Inflow]  in x1 direction 
# Matching boundaries are located at are [1,0,0] and [2, 1, 0]
block0_bc = []
direction = 0
side = 0
block0_bc.append(InterfaceBC(direction=0, side=0,  halos=[-2,2], name="block0_to_block1", match=(1, 0, 0, True)))
block0_bc.append(ExtrapolationBC(direction=0, side=1, order=0))
block0_bc.append(SharedInterfaceBC(direction=1, side=0,  halos=[-4,4], name="block0_to_block2", match=(2, 1, 0, True)))
block0_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
mb_bcs[0] = block0_bc

# Boundary conditions for block 1
#The boundary conditions are [InterfaceBC, InterfaceBC] in x0 direction and [wall, Inflow]  in x1 direction 
# Matching boundaries are located at are [0,0,0] and [2, 0, 0]
block1_bc = []
block1_bc.append(InterfaceBC(direction=0, side=0,  halos=[-2,2], name="block1_to_block0", match=(0, 0, 0, True)))
block1_bc.append(InterfaceBC(direction=0, side=1,  halos=[-2,2], name="block1_to_block2", match=(2, 0, 0, False)))
# Wall temperature is required for halo points
Twall = ConstantObject('Twall')
Twall.value = 1.0
if conservative:
    wall_energy = [Eq(q_vector[-1], Twall*q_vector[0]/((gama-1.0)*gama*Minf*Minf))]
else:
    wall_energy = [Eq(DataObject('Et'), Twall/((gama-1.0)*gama*Minf*Minf))]
block1_bc.append(IsothermalWallBC(direction=1, side=0, corners=False, equations=wall_energy, multi_block=True))
block1_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
mb_bcs[1] = block1_bc

# Boundary conditions for block 2
# The boundary conditions are [InterfaceBC, outflow] in x0 direction and  SharedInterfaceBC, Inflow]  in x1 direction 
# Matching boundaries are located at are [1,0,1] and [0, 1, 0]
block2_bc = []
block2_bc.append(InterfaceBC(direction=0, side=0,  halos=[-2,2], name="block2_to_block1", match=(1, 0, 1, False)))
block2_bc.append(ExtrapolationBC(direction=0, side=1, order=0))
block2_bc.append(SharedInterfaceBC(direction=1, side=0,  halos=[-4,4], name="block2_to_block0", match=(0, 1, 0, True)))
block2_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
mb_bcs[2] = block2_bc
# Set the multi block boundary conditions
multi_block.set_block_boundaries(mb_bcs)

# Set the equations on the blocks
multi_block.set_equations([simulation_eq, constituent, metriceq])

# Add statsistics gathering
stats = True
if stats:
    # Create the statistics equations, this shows another way of writing the equations
    from airfoil_stats import favre_averaged_stats
    q_vector = flatten(simulation_eq.time_advance_arrays)
    stat_equation_classes, stats_arrays = favre_averaged_stats(ndim, q_vector, conservative=conservative)
else:
    stat_equation_classes, stats_arrays = [], []

multi_block.set_equations(stat_equation_classes)

# Add filters to each block
filters = {0:[], 1:[], 2:[]}
for no, block in enumerate(multi_block.blocks):
    if no == 0 or no == 1 or no == 2: # Main aerofoil block, C-mesh. Don't filter near the aerofoil
        filters[no] += [WENOFilter(block, order=7, metrics=metriceq, dissipation_sensor='Ducros', airfoil=True, flux_type='LLF').equation_classes]

# Add DRP filters for freestream
for no, block in enumerate(multi_block.blocks):
    filters[no] += [ExplicitFilter(block, [0,1], width=9, filter_type='DRP', optimized=False, sigma=0.2, wall_control=True, multi_block=multi_block).equation_classes]

# Add a binomial filter on the outlet boundary to kill reflections
for no, block in enumerate(multi_block.blocks):
    i, j = block.grid_indexes[0], block.grid_indexes[1]
    if no == 0:
        grid_condition = i >= 794
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
        grid_condition = j >= 475
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
    elif no == 1:
        grid_condition = j >= 475
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
    elif no == 2:
        grid_condition = i >= 794
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
        grid_condition = j >= 475
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
multi_block.set_filters(filters)

# HDF5 input/output
kwargs = {'iotype': "Write"}
q_hdf5 = iohdf5(save_every=1000, **kwargs)
q_hdf5.add_arrays(simulation_eq.time_advance_arrays)
q_hdf5.add_arrays([DataObject('kappa'), DataObject('Mach_sensor'), DataObject('q0'), DataObject('q1'), DataObject('q2'), DataObject('q3')])
# Read in the grid file
kwargs = {'iotype': "Read"}
x,y = symbols("x0, x1", **{'cls':DataObject})
grid_hdf5 = iohdf5(**kwargs)
grid_hdf5.add_arrays([x, y])
# Stats HDF5 and write metrics to the grid file
metrics_hdf5 = iohdf5(arrays=metriceq.grid_der_wks, **{'position': "init", 'iotype': 'Write', 'name': "metrics.h5"})
# HDF5 output of statistics arrays
kwargs = {'iotype': "Write", 'name': "stats_output.h5"}
stats_hdf5 = iohdf5(arrays=stats_arrays, **kwargs)
# Set the I/O on the block
multi_block.setio([q_hdf5, grid_hdf5, metrics_hdf5, stats_hdf5])
# Perform the discretization
multi_block.discretise()


# Add the wake treatment kernels
wake_ker = generate_wake_kernel(q_vector, multi_block, wall_energy[0])
# Sponge zones for outer boundaries
# Outlet
outlet_sponge_block0 = generate_outlet_sponge(q_vector, multi_block.get_block(0), Lx=4.5, npoints=62)
outlet_sponge_block2 = generate_outlet_sponge(q_vector, multi_block.get_block(2), Lx=4.5, npoints=62)
# Farfield
farfield_sponge_block0 = generate_farfield_sponge(q_vector, multi_block.get_block(0), Ly=7.5, npoints=62)
farfield_sponge_block1 = generate_farfield_sponge(q_vector, multi_block.get_block(1), Ly=7.5, npoints=62)
farfield_sponge_block2 = generate_farfield_sponge(q_vector, multi_block.get_block(2), Ly=7.5, npoints=62)

# Add wake exchanges and kernels to block2 boundary conditions
b = multi_block.get_block(2)
for no, eq in enumerate(b.list_of_equation_classes):
    # Add the sponge kernels after updating the residuals for all blocks, before time advancement
    if isinstance(eq, SimulationEquations):
        eq.Kernels += [outlet_sponge_block0, outlet_sponge_block2, farfield_sponge_block0, farfield_sponge_block1, farfield_sponge_block2]
        eq.boundary_kernels += wake_ker

# Make some full swaps for interfaces before filtering
filter_swaps = create_exchange_calls_codes(multi_block, simulation_eq.time_advance_arrays)
for block in multi_block.blocks:
    for no, eq in enumerate(block.list_of_equation_classes):
        if isinstance(eq, UserDefinedEquations):
            if eq.full_swap:
                eq.Kernels += filter_swaps

# Create the OPS C code
alg = TraditionalAlgorithmRKMB(multi_block)
OPSC(alg, OPS_diagnostics=1)
# NaN check and iteration counter
print_iteration_ops(NaN_check='rho', every=100, nblocks=nblocks)
# Substitute simulation parameter values
constants = ['gama', 'Minf', 'Pr', 'Re', 'dt', 'niter', 'sigma_filt', 'SuthT', 'RefT', 'stat_frequency']
values = ['1.4', '0.70', '0.72', '5.0e5', '3.0e-5', '1000000', '0.01', '110.4', '268.67', '10']
# Block 0
constants += ['block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'inv_rfact0_block0', 'inv_rfact1_block0']
values += ['1099', '980', '11.5/(block0np0 - 1.0)', '22.5/(block0np1 - 1.0)', '1.0/Delta0block0', '1.0/Delta1block0']
# Block 1
constants += ['block1np0', 'block1np1', 'Delta0block1', 'Delta1block1', 'inv_rfact0_block1', 'inv_rfact1_block1']
values += ['1495', '980', '22.5/(block1np0 - 1.0)', '22.5/(block1np1 - 1.0)', '1.0/Delta0block1', '1.0/Delta1block1']
# Block 2
constants += ['block2np0', 'block2np1', 'Delta0block2', 'Delta1block2', 'inv_rfact0_block2', 'inv_rfact1_block2']
values += ['1099', '980', '11.5/(block2np0 - 1.0)', '22.5/(block2np1 - 1.0)', '1.0/Delta0block2', '1.0/Delta1block2']
substitute_simulation_parameters(constants, values)
