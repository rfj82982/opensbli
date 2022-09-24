#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters, debug_equation
from sympy import pi, sin, cos, Abs, sqrt

ndim = 3
stats = True
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True), (False, False)], 2)
#Create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs

# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]
# symbol for the coordinate system in the equations
conservative = False
# NS = NS_Split('Feiereisen', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic')
NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic', energy_formulation='enthalpy', debug=False)

mass, momentum, energy = NS.mass, NS.momentum, NS.energy
# pprint(-1*debug_equation(ndim, energy, 0))

# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)

einstein_eq = EinsteinEquation()
einstein_eq.optional_subs_dict = optional_subs_dict

metric_vel = "Eq(U_i, D_i_j*u_j)"
eqns = einstein_eq.expand(metric_vel, ndim, coordinate_symbol, [], constants)
for eq in eqns:
    einstein_eq.optional_subs_dict[eq.lhs] = eq.rhs

# Constituent relations
if conservative:
    pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
    velocity = "Eq(u_i, rhou_i/rho)"
    enthalpy = "Eq(H, (rhoE + p) / rho)"
else:
    pressure = "Eq(p, rho*(gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j)))"
    enthalpy = "Eq(H, Et + p / rho)"

temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object
# Expand momentum and add the expanded equations to the constituent relations
if conservative:
    velocity = "Eq(u_i, rhou_i/rho)"
    eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)

# Expand pressure and add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand temperature and add the expanded equations to the constituent relations
eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand enthalpy and add the expanded equations to the constituent relations
eqns = einstein_eq.expand(enthalpy, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# # Expand viscosity and add the expanded equations to the constituent relations
eqns = einstein_eq.expand(viscosity, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)


# Create a simulation block
block = SimulationBlock(ndim, block_number=0, conservative=conservative)
simulation_eq.apply_metrics(metriceq)

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Initial conditions as strings
u0 = "Eq(GridVariable(u0), 1.0)"
u1 = "Eq(GridVariable(u1), 0.0)"
u2 = "Eq(GridVariable(u2), 0.0)"
p = "Eq(GridVariable(p), 1/(gama*Minf*Minf))"
r = "Eq(GridVariable(r), gama*Minf*Minf*p)"

rho = "Eq(DataObject(rho), r)"
rhou0 = "Eq(DataObject(u0), u0)"
rhou1 = "Eq(DataObject(u1), u1)"
rhou2 = "Eq(DataObject(u2), u2)"
rhoE = "Eq(DataObject(Et), p/(r*(gama-1)) + 0.5*(u0**2+ u1**2 + u2**2))"
temp_eq = "Eq(GridVariable('temp'), DataObject(x2))"
eqns = [u0, u1, u2, p, r, rho, rhou0, rhou1, rhou2, rhoE, temp_eq]

# parse the initial conditions
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
# Low storage optimisation for the central scheme
fns = 'u0 u1 u2 T'
cent = StoreSome(4, fns)
# cent = Central(4)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(4)
schemes[rk.name] = rk

# Create boundaries, one for each side per dimension
q_vector = flatten(simulation_eq.time_advance_arrays)
boundaries = []
direction = 0
# Apply a periodic boundary over the shared mesh line
boundaries += [PeriodicBC(direction, 0, halos=[-2,2])]
boundaries += [PeriodicBC(direction, 1, halos=[-2,2])]
# Isothermal wall in x1 direction
gama, Minf, Twall = symbols('gama Minf Twall', **{'cls': ConstantObject})
# Energy on the wall is set
wall_energy = [Eq(q_vector[-1], Twall / (gama * Minf**2.0 * (gama - S.One)))]
direction = 1
lower_wall_eq = wall_energy[:]
boundaries += [IsothermalWallBC(direction, 0, lower_wall_eq)]
# Far field boundary
direction, side = 1,1
boundaries += [DirichletBC(direction, side, initial_equations)]
# Periodic span
direction = 2
boundaries += [PeriodicBC(direction, 0, halos=[-2,2])]
boundaries += [PeriodicBC(direction, 1, halos=[-2,2])]

# set the boundaries for the block
block.set_block_boundaries(boundaries)

# Set the equations to be solved on the block
if stats:
    # Create the statistics equations, this shows another way of writing the equations
    from airfoil_stats import favre_averaged_stats
    q_vector = flatten(simulation_eq.time_advance_arrays)
    stat_equation_classes, stats_arrays = favre_averaged_stats(ndim, q_vector, conservative=conservative)
else:
    stat_equation_classes, stats_arrays = [], []

block.set_equations([constituent, simulation_eq, initial, metriceq] + stat_equation_classes)

# Set the IO class to write out arrays
h5 = iohdf5(save_every=1000, **{'iotype': "Write"})
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('kappa'), DataObject('Mach_sensor')]) # shock sensor array
# Read grid file
h5_read = iohdf5(**{'iotype': "Read"})
h5_read.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('x2')])
# HDF5 output of statistics arrays
kwargs = {'iotype': "Write", 'name': "stats_output.h5"}
stats_hdf5 = iohdf5(arrays=stats_arrays, **kwargs)
# Write grid metrics to a file
metrics_hdf5 = iohdf5(arrays=metriceq.grid_der_wks, **{'position': "init", 'iotype': 'Write', 'name': "metrics.h5"})
block.setio([h5, h5_read, stats_hdf5, metrics_hdf5])

# Various filters and shock capturing
j = block.grid_indexes[1]
grid_condition = j >= 642
BF = BinomialFilter(block, order=4, directions=3, grid_condition=grid_condition, sigma=0.2)
block.set_equations(BF.equation_classes)

DRP = ExplicitFilter(block, [0,1,2], width=11, filter_type='DRP', optimized=True, Mach_sensor = True, sigma=0.3, wall_control=True, multi_block=None)
block.set_equations(DRP.equation_classes)

# WENO filter for shock-capturing
WF = WENOFilter(block, order=5, metrics=metriceq, dissipation_sensor='Ducros', flux_type='LLF', airfoil=True)
block.set_equations(WF.equation_classes)

# set the discretisation schemes
block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
block.discretise()

WF.update_periodic_boundary(block, halos=[-3,4])

# Add some full [-5,5] halo swaps over the periodic directions only when the filter is called
def create_exchange_calls_codes(block, dsets):
    kernels = []
    arrays = [block.location_dataset(a) for a in flatten(dsets)]
    for direction in [0,2]:
        for side in [0,1]:
            BC = PeriodicBC(direction, side, halos=[-5,5])
            kernels += [BC.apply(arrays, block)]
    return kernels

# Make some full swaps for interfaces before filtering
filter_swaps = create_exchange_calls_codes(block, ['rho', 'u0', 'u1', 'u2', 'Et'])
for no, eq in enumerate(block.list_of_equation_classes):
    if isinstance(eq, UserDefinedEquations):
        if eq.full_swap:
            eq.Kernels += filter_swaps

alg = TraditionalAlgorithmRK(block)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg, OPS_diagnostics=1)
# Simulation parameters
constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0', 'Twall', 'stat_frequency', 'RefT', 'SuthT', 'inv_rfact0', 'inv_rfact1', 'inv_rfact2']
values = ['5.0e5', '1.4', '0.70', '0.71', '3.0e-5', '500000000', '3001', '647', '50', '49.0/(block0np0-1)', '49.0/(block0np1-1)', '0.05/(block0np2-1)', '1.0', '10', '273.15', '110.4', '(block0np0-1)/49.0', '(block0np1-1)/49.0', '(block0np2)/0.05' ]
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho')
