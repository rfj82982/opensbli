#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import pi, sin, cos, Abs, sqrt

ndim = 2
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)
#Create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs

#Define the compresible Navier-Stokes equations in Einstein notation, by default the scheme is Central no need to
#Specify the schemes
mass = "Eq(Der(rho,t), - Skew(rho*u_j,x_j))"
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j))"
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j))"
# Substitutions used in the equations
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
substitutions = [stress_tensor, heat_flux]
# Constants that are used
constants = ["Re", "Pr", "gama", "Minf"]
# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, T**0.7)"

einstein_eq = EinsteinEquation()
einstein_eq.optional_subs_dict = optional_subs_dict

metric_vel = "Eq(U_i, D_i_j*u_j)"
eqns = einstein_eq.expand(metric_vel, ndim, coordinate_symbol, substitutions, constants)
for eq in eqns:
    einstein_eq.optional_subs_dict[eq.lhs] = eq.rhs

# Change the symbol to xi as we will be using metrics
simulation_eq = SimulationEquations()
# Perform the expansion
eqns = einstein_eq.expand(mass, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)
eqns = einstein_eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)
eqns = einstein_eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object
# Expand momentum and add the expanded equations to the constituent relations
eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
# Expand pressure and add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
# Expand temperature and add the expanded equations to the constituent relations
eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
# Expand viscosity and add the expanded equations to the constituent relations
eqns = einstein_eq.expand(viscosity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)


# Create a simulation block
block = SimulationBlock(ndim, block_number=0)
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
boundaries += [PeriodicBC(direction, 0, full_swap=True)]
boundaries += [PeriodicBC(direction, 1, full_swap=True)]
# Isothermal wall in x1 direction
gama, Minf, Twall = symbols('gama Minf Twall', **{'cls': ConstantObject})
# Energy on the wall is set
wall_energy = [Eq(q_vector[-1], Twall*q_vector[0] / (gama * Minf**2.0 * (gama - S.One)))]
direction = 1
lower_wall_eq = wall_energy[:]
boundaries += [IsothermalWallBC(direction, 0, lower_wall_eq)]
# Far field boundary
direction, side = 1,1
boundaries += [DirichletBC(direction, side, initial_equations)]
# set the boundaries for the block
block.set_block_boundaries(boundaries)

# Set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=1000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('kappa')])
kwargs = {'iotype': "Read"}
h5_read = iohdf5(**kwargs)
h5_read.add_arrays([DataObject('x0'), DataObject('x1')])
block.setio([h5, h5_read])

# Boundary filtering
j = block.grid_indexes[1]
grid_condition = j >= 185
BF = BinomialFilter(block, order=6, grid_condition=grid_condition, sigma=0.1)

# Set the equations to be solved on the block
block.set_equations([constituent, simulation_eq, initial, metriceq])
block.set_equations(BF.equation_classes)
DRP = ExplicitFilter(block, [0,1], width=11, filter_type='DRP', optimized=True, sigma=0.2, wall_control=True, multi_block=None)
block.set_equations(DRP.equation_classes)
# WENO filter for shock-capturing
WF = WENOFilter(block, order=7, metrics=metriceq, dissipation_sensor='Ducros', Mach_correction=False, flux_type='LLF')
block.set_equations(WF.equation_classes)
# set the discretisation schemes
block.set_discretisation_schemes(schemes)
# Discretise the equations on the block
block.discretise()

alg = TraditionalAlgorithmRK(block)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg)
# Simulation parameters
constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'Twall']
values = ['2.1e5', '1.4', '0.7', '0.71', '0.00002', '500000000', '2301', '192', '37.6887/(block0np0-1)', '36.9844/(block0np1-1)', '1.0']
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho')
