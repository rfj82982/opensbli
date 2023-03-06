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

ndim = 3
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
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True), (False, False)], 2, latex_debug=False)
#Create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs
Einstein_expansion = EinsteinEquation()
Einstein_expansion.optional_subs_dict = optional_subs_dict
metric_vel = "Eq(U_i, D_i_j*u_j)"
eqns = Einstein_expansion.expand(metric_vel, ndim, coordinate_symbol, [], constants)
for eq in eqns:
    Einstein_expansion.optional_subs_dict[eq.lhs] = eq.rhs

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
rk = RungeKuttaLS(4)
schemes[rk.name] = rk
# cent = Central(4)
cent = StoreSome(4, 'u0 u1 u2 T')
schemes[cent.name] = cent
multi_block.set_discretisation_schemes(schemes)

# Initial conditions
## Need to change for swept cases
d, u0, u1, u2, p = symbols("d, u0:3, p", **{'cls':GridVariable})
gama, Minf = symbols("gama, Minf", **{'cls':ConstantObject})
initial_equations = []
initial_equations += [Eq(d, 1.0)]
initial_equations += [Eq(u0, 1.0)]
initial_equations += [Eq(u1, 0.0)]
initial_equations += [Eq(u2, 0.0)]
initial_equations += [Eq(p, 1.0/(gama*Minf**2.0))]

# Set the q vector values for initial condition and farfield boundaries
q_vector = flatten(simulation_eq.time_advance_arrays)
if conservative:
    initial_equations += [Eq(q_vector[0], d)]
    initial_equations += [Eq(q_vector[1], d*u0)]
    initial_equations += [Eq(q_vector[2], d*u1)]
    initial_equations += [Eq(q_vector[3], d*u2)]
    initial_equations += [Eq(q_vector[4], p/(gama-1.0) + 0.5* d *(u0**2+u1**2+ u2**2))]
else:
    initial_equations += [Eq(q_vector[0], d)]
    initial_equations += [Eq(q_vector[1], u0)]
    initial_equations += [Eq(q_vector[2], u1)]
    initial_equations += [Eq(q_vector[3], u2)]
    initial_equations += [Eq(q_vector[4], p/(d*(gama-1.0)) + 0.5*(u0**2+u1**2+ u2**2))]

temp_x2 = [Eq(GridVariable('temp'), DataObject('x2'))]
initial = GridBasedInitialisation()
initial.add_equations(copy.deepcopy(initial_equations) + temp_x2)
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
block0_bc.append(SharedInterfaceBC(direction=1, side=0, halos=[-4,4], name="block0_to_block2", match=(2, 1, 0, True)))
block0_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
block0_bc.append(PeriodicBC(direction=2, side=0, halos=[-2,2], corners=False))
block0_bc.append(PeriodicBC(direction=2, side=1, halos=[-2,2], corners=False))
mb_bcs[0] = block0_bc

# Boundary conditions for block 1
#The boundary conditions are [InterfaceBC, InterfaceBC] in x0 direction and [wall, Inflow]  in x1 direction 
# Matching boundaries are located at are [0,0,0] and [2, 0, 0]
block1_bc = []
block1_bc.append(InterfaceBC(direction=0, side=0, halos=[-2,2], name="block1_to_block0", match=(0, 0, 0, True)))
block1_bc.append(InterfaceBC(direction=0, side=1, halos=[-2,2], name="block1_to_block2", match=(2, 0, 0, False)))

## Wall condition
# Isothermal wall in x1 direction
gama, Minf, Twall = symbols('gama Minf Twall', **{'cls': ConstantObject})
# Boundary-layer tripping
tripped = True
direction, side = 1, 0
if tripped:
    Amp, sigma, xts, xtp = symbols('tripA tripSigma xts xtp', **{'cls':ConstantObject})
    # Time dependence
    # current_iter = multi_block.get_block(nblocks-1).get_temporal_schemes[0].iteration_number # Current iteration number
    current_iter = Globalvariable("iter", integer=True)
    dt, omega0, omega1, omega2 = symbols('dt omega_0 omega_1 omega_2', **{'cls': ConstantObject})
    t = dt*current_iter 
    # Spatial Modes
    k0, k1, k2 = symbols('k_0 k_1 k_2', **{'cls': ConstantObject})
    phi0, phi1, phi2 = symbols('phi_0 phi_1 phi_2', **{'cls': ConstantObject})
    # Coordinate arrays
    x0, z0 = DataObject('x0'), DataObject('x2')
    conditional_expressions = []
    # Suction side trip
    SS_trip = Amp*(exp(-(x0 - xts)**2  / (2*sigma**2))*(sin(k0**2 * z0)*sin(omega0*t + phi0) + sin(k1**2 * z0)*sin(omega1*t + phi1) + sin(k2**2 * z0)*sin(omega2*t + phi2)))
    # Pressure side trip
    PS_trip = Amp*(exp(-(x0 - xtp)**2  / (2*sigma**2))*(sin(k0**2 * z0)*sin(omega0*t + phi0) + sin(k1**2 * z0)*sin(omega1*t + phi1) + sin(k2**2 * z0)*sin(omega2*t + phi2)))
    # Index in x direction, assuming anti-clockwise grid configuration here
    idx = multi_block.get_block(1).grid_indexes[0] # x index on block 1 (airfoil block)
    expr_condition_pairs = Piecewise((SS_trip, idx > ConstantObject('block1np0')/2), (PS_trip, idx < ConstantObject('block1np0')/2),  (0, True))
    v = OpenSBLIEq(GridVariable('v'), expr_condition_pairs)
    rhov_wall = OpenSBLIEq(DataObject('rhou1'), DataObject('rho')*GridVariable('v'))
    # Adding the non-zero V component to the calculation of rhoE at the wall
    rhoE_wall = OpenSBLIEq(DataObject('rhoE'), DataObject('rho')*Twall/(gama*(gama-1.0)*Minf**2.0) + 0.5*DataObject('rho')*GridVariable('v')**2)
    # Equations to set rhov and rhoE on the wall
    wall_eqns = [rhov_wall, rhoE_wall]
    wall_energy = [rhoE_wall]
    block1_bc.append(ForcingStripBC(direction, 0, v, wall_eqns, corners=False, multi_block=True))
else:
    # Energy on the wall is set
    if conservative:
        wall_energy = [Eq(q_vector[-1], q_vector[0]*Twall / (gama * Minf**2.0 * (gama - S.One)))]
    else:
        wall_energy = [Eq(q_vector[-1], Twall / (gama * Minf**2.0 * (gama - S.One)))]
    direction = 1
    lower_wall_eq = wall_energy[:]

    block1_bc.append(IsothermalWallBC(direction=1, side=0, corners=False, equations=lower_wall_eq, multi_block=True))


## Farfield
block1_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
block1_bc.append(PeriodicBC(direction=2, side=0, halos=[-2,2], corners=False))
block1_bc.append(PeriodicBC(direction=2, side=1, halos=[-2,2], corners=False))
mb_bcs[1] = block1_bc

# Boundary conditions for block 2
# The boundary conditions are [InterfaceBC, outflow] in x0 direction and  SharedInterfaceBC, Inflow]  in x1 direction 
# Matching boundaries are located at are [1,0,1] and [0, 1, 0]
block2_bc = []
block2_bc.append(InterfaceBC(direction=0, side=0,  halos=[-2,2], name="block2_to_block1", match=(1, 0, 1, False)))
block2_bc.append(ExtrapolationBC(direction=0, side=1, order=0))
block2_bc.append(SharedInterfaceBC(direction=1, side=0,  halos=[-4,4], name="block2_to_block0", match=(0, 1, 0, True)))
block2_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
block2_bc.append(PeriodicBC(direction=2, side=0, halos=[-2,2], corners=False))
block2_bc.append(PeriodicBC(direction=2, side=1, halos=[-2,2], corners=False))
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
shock_filters = []
for no, block in enumerate(multi_block.blocks):
    if no == 0 or no == 1 or no == 2: # Main aerofoil block, C-mesh. Don't filter near the aerofoil
        WF = WENOFilter(block, order=5, metrics=metriceq, dissipation_sensor='Ducros', airfoil=True, flux_type='LLF')
        shock_filters.append(WF)
        filters[no] += [WF.equation_classes]

# Add DRP filters for freestream
for no, block in enumerate(multi_block.blocks):
    filters[no] += [ExplicitFilter(block, [0,1,2], width=9, filter_type='DRP', optimized=False, sigma=0.3333333, wall_control=True, multi_block=multi_block).equation_classes]

# Add a binomial filter on the outlet boundary to kill reflections
for no, block in enumerate(multi_block.blocks):
    i, j, k = block.grid_indexes[0], block.grid_indexes[1], block.grid_indexes[2]
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
x,y,z = symbols("x0, x1, x2", **{'cls':DataObject})
kwargs = {'iotype': "Write"}
q_hdf5 = iohdf5(save_every=1000, **kwargs)
q_hdf5.add_arrays(simulation_eq.time_advance_arrays)
q_hdf5.add_arrays([DataObject('kappa')])
# Read in the grid file
kwargs = {'iotype': "Read"}
grid_hdf5 = iohdf5(**kwargs)
grid_hdf5.add_arrays([x, y, z])
# Stats HDF5 and write metrics to the grid file
metrics_hdf5 = iohdf5(arrays=metriceq.grid_der_wks, **{'position': "init", 'iotype': 'Write', 'name': "metrics.h5"})
# HDF5 output of statistics arrays
kwargs = {'iotype': "Write", 'name': "stats_output.h5"}
stats_hdf5 = iohdf5(arrays=stats_arrays, **kwargs)
# Set the I/O on the block
multi_block.setio([q_hdf5, grid_hdf5, metrics_hdf5, stats_hdf5])
# Perform the discretization
multi_block.discretise()

# Add a periodic boundary condition call for WENO filters
for i, block in enumerate(multi_block.blocks):
    shock_filters[i].update_periodic_boundary(block, halos=[-4,4])

# Add the wake treatment kernels
wake_ker = generate_wake_kernel(q_vector, multi_block, wall_energy[0])
# Sponge zones for outer boundaries
# Outlet
outlet_sponge_block0 = generate_outlet_sponge(q_vector, multi_block.get_block(0), Lx=4.5, npoints=12)
outlet_sponge_block2 = generate_outlet_sponge(q_vector, multi_block.get_block(2), Lx=4.5, npoints=12)
# Farfield
farfield_sponge_block0 = generate_farfield_sponge(q_vector, multi_block.get_block(0), Ly=7.5, npoints=12)
farfield_sponge_block1 = generate_farfield_sponge(q_vector, multi_block.get_block(1), Ly=7.5, npoints=12)
farfield_sponge_block2 = generate_farfield_sponge(q_vector, multi_block.get_block(2), Ly=7.5, npoints=12)

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

# Add some full [-5,5] halo swaps over the periodic directions only when the filter is called
def create_periodic_BCs(multi_block, dsets):
    kernels = []
    for block in multi_block.blocks:
        arrays = [block.location_dataset(a) for a in flatten(dsets)]
        for direction in [2]:
            for side in [0,1]:
                BC = PeriodicBC(direction, side, halos=[-4,4], corners=False)
                kernels += [BC.apply(arrays, block)]
    return kernels

# Periodic boundary condition for DRP filters
if conservative:
    dsets = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
else:
    dsets = ['rho', 'u0', 'u1', 'u2', 'Et']
DRP_periodic = create_periodic_BCs(multi_block, dsets)
for block in multi_block.blocks:
    for no, eq in enumerate(block.list_of_equation_classes):
        if isinstance(eq, UserDefinedEquations):
            if eq.order == 0:
                eq.Kernels += DRP_periodic

# Create the OPS C code
alg = TraditionalAlgorithmRKMB(multi_block)
OPSC(alg, OPS_diagnostics=1)
# NaN check and iteration counter
print_iteration_ops(NaN_check='rho', every=100, nblocks=nblocks)
# Substitute simulation parameter values
constants = ['gama', 'Minf', 'Pr', 'Re', 'dt', 'niter', 'sigma_filt', 'SuthT', 'RefT', 'stat_frequency', 'shock_factor', 'Twall']
values = ['1.4', '0.70', '0.72', '5.0e5', '4.0e-5', '1000000', '0.01', '110.4', '268.67', '100', '1.0', '1.0']
# Block 0
constants += ['block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0', 'inv_rfact0_block0', 'inv_rfact1_block0', 'inv_rfact2_block0']
values += ['1099', '980', '50', '11.5/(block0np0 - 1.0)', '22.5/(block0np1 - 1.0)', '0.05/block0np2', '1.0/Delta0block0', '1.0/Delta1block0', '1.0/Delta2block0']
# Block 1
constants += ['block1np0', 'block1np1', 'block1np2', 'Delta0block1', 'Delta1block1', 'Delta2block1', 'inv_rfact0_block1', 'inv_rfact1_block1', 'inv_rfact2_block1']
values += ['1495', '980', '50', '2.0461756979465546/(block1np0 - 1.0)', '22.5/(block1np1 - 1.0)', '0.05/block1np2', '1.0/Delta0block1', '1.0/Delta1block1', '1.0/Delta2block1']
# Block 2
constants += ['block2np0', 'block2np1', 'block2np2', 'Delta0block2', 'Delta1block2', 'Delta2block2', 'inv_rfact0_block2', 'inv_rfact1_block2', 'inv_rfact2_block2']
values += ['1099', '980', '50', '11.5/(block2np0 - 1.0)', '22.5/(block2np1 - 1.0)', '0.05/block2np2', '1.0/Delta0block2', '1.0/Delta1block2', '1.0/Delta2block2']

# Add forcing modes
constants += ['tripA', 'tripSigma', 'xts', 'xtp', 'omega_0', 'omega_1', 'omega_2', 'k_0', 'k_1', 'k_2', 'phi_0', 'phi_1', 'phi_2']
values += ['0.05', '0.00833', '0.1', '0.1', '26', '88', '200', '120*M_PI', '160*M_PI', '160*M_PI', '0.0', 'M_PI', '-M_PI/2']

substitute_simulation_parameters(constants, values)
