#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import pi, sin, cos, Abs, sqrt

# Problem dimension
ndim = 2
# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "RefT", "SuthT"]
# # symbol for the coordinate system in the equations
coordinate_symbol = "x"
# symbol for the coordinate system in the equations
conservative = True
NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic', energy_formulation='enthalpy')

mass, momentum, energy = NS.mass, NS.momentum, NS.energy
# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)

# Constituent relations used in the system
if conservative:
    pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
    velocity = "Eq(u_i, rhou_i/rho)"
    enthalpy = "Eq(H, (rhoE + p) / rho)"
else:
    pressure = "Eq(p, rho*(gama-1)*(E - (1/2)*(KD(_i,_j)*u_i*u_j)))"
    enthalpy = "Eq(H, E + p / rho)"

temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)
optional_subs_dict = metriceq.metric_subs

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object
einstein_eq = EinsteinEquation()
einstein_eq.optional_subs_dict = optional_subs_dict

# Expand momentum add the expanded equations to the constituent relations
if conservative:
    eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)
# Expand pressure add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand enthalpy add the expanded equations to the constituent relations
eqns = einstein_eq.expand(enthalpy, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand temperature add the expanded equations to the constituent relations
eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand viscosity add the expanded equations to the constituent relations
eqns = einstein_eq.expand(viscosity, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Create a simulation block
block = SimulationBlock(ndim, block_number=0, conservative=conservative)


# Transform to curvilinear
simulation_eq.apply_metrics(metriceq)

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Initial conditions as strings
u0 = "Eq(GridVariable(u0),1.0)"
u1 = "Eq(GridVariable(u1), 0.0,)"
p = "Eq(GridVariable(p), 1/(gama*Minf*Minf))"
r = "Eq(GridVariable(r), gama*Minf*Minf*p)"

rho = "Eq(DataObject(rho), r)"
rhou0 = "Eq(DataObject(rhou0), r*u0)"
rhou1 = "Eq(DataObject(rhou1), r*u1)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r *(u0**2+ u1**2))"
eqns = [u0, u1, p, r, rho, rhou0, rhou1, rhoE]

# parse the initial conditions
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
# Low storage optimisation for the central scheme
fns = 'u0 u1 T'
cent = StoreSome(4, fns)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(4)
schemes[rk.name] = rk

# Create boundaries, one for each side per dimension
q_vector = flatten(simulation_eq.time_advance_arrays)
boundaries = []
direction = 0
# Apply a periodic boundary over the shared mesh line
boundaries += [PeriodicBC(direction, 0)]
boundaries += [PeriodicBC(direction, 1)]
# Isothermal wall in x1 direction
gama, Minf, Twall = symbols('gama Minf Twall', **{'cls': ConstantObject})
# Energy on the wall is set
if conservative:
    wall_energy = [Eq(q_vector[-1], Twall*q_vector[0] / (gama * Minf**2.0 * (gama - S.One)))]
else:
    wall_energy = [Eq(q_vector[-1], Twall / (gama * Minf**2.0 * (gama - S.One)))]

direction = 1
lower_wall_eq = wall_energy[:]
boundaries += [IsothermalWallBC(direction, 0, lower_wall_eq, scheme=ReducedAccess())]
# Far field boundary
direction, side = 1,1
boundaries += [DirichletBC(direction, side, initial_equations, scheme=ReducedAccess())]
# set the boundaries for the block
block.set_block_boundaries(boundaries)

# Set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=5000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('kappa')])
kwargs = {'iotype': "Read"}
h5_read = iohdf5(**kwargs)
h5_read.add_arrays([DataObject('x0'), DataObject('x1')])
block.setio([h5, h5_read])

# Add SFD filtering
# SFD = SFD(block, chifilt=0.1, omegafilt=1.0/0.75)

j = block.grid_indexes[1]
grid_condition = j >= 595
BF = BinomialFilter(block, order=6, directions=[0,1], grid_condition=grid_condition, sigma=0.1)


# Set the equations to be solved on the block
block.set_equations([constituent, simulation_eq, initial, metriceq])
ShockFilter = WENOFilter(block, order=5, metrics=metriceq, dissipation_sensor='Ducros', Mach_correction=False, flux_type='LLF')

block.set_equations(ShockFilter.equation_classes + BF.equation_classes)

DRP = ExplicitFilter(block, [0,1], width=11, filter_type='DRP', optimized=True, sigma=0.2, wall_control=True, multi_block=None)
block.set_equations(DRP.equation_classes)

# set the discretisation schemes
block.set_discretisation_schemes(schemes)
# Discretise the equations on the block
block.discretise()
# Full 5 swaps for the filter over the interface
# Add some full [-5,5] halo swaps over the periodic directions only when the filter is called
def create_exchange_calls_codes(block, dsets):
    kernels = []
    arrays = [block.location_dataset(a) for a in flatten(dsets)]
    for direction in [0]:
        for side in [0,1]:
            BC = PeriodicBC(direction, side, full_depth=True)
            kernels += [BC.apply(arrays, block)]
    return kernels

# Make some full swaps for interfaces before filtering
if conservative:
    filter_swaps = create_exchange_calls_codes(block, ['rho', 'rhou0', 'rhou1', 'rhoE'])
else:
    filter_swaps = create_exchange_calls_codes(block, ['rho', 'u0', 'u1', 'E'])

for no, eq in enumerate(block.list_of_equation_classes):
    if isinstance(eq, UserDefinedEquations):
        if eq.full_swap:
            eq.Kernels += filter_swaps

alg = TraditionalAlgorithmRK(block)
# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)
# Write the code for the algorithm
OPSC(alg)
# Simulation parameters
constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'Twall', 'SuthT', 'RefT']
values = ['300.0', '1.4', '1.5', '0.71', '0.0002', '5000000', '471', '602', '1.0/(block0np0-1)', '1.0/(block0np1-1)', '1.0', '110.4', '273.15']
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho', every=100)
