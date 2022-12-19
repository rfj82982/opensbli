#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

# Number of dimensions of the system to be solved
ndim = 3
# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf"]
coordinate_symbol = "x"
# symbol for the coordinate system in the equations
conservative = True
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
# Constituent relations
if conservative:
    pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
    velocity = "Eq(u_i, rhou_i/rho)"
    enthalpy = "Eq(H, (rhoE + p) / rho)"
else:
    pressure = "Eq(p, rho*(gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j)))"
    enthalpy = "Eq(H, Et + p / rho)"

temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.4042)/(T+0.40417)))" ## Modified sutherland law

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

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Initial conditions as strings
x0 = "Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])"
x1 = "Eq(GridVariable(x1), block.deltas[1]*block.grid_indexes[1])"
x2 = "Eq(GridVariable(x2), block.deltas[2]*block.grid_indexes[2])"

u0 = "Eq(GridVariable(u0),sin(x0)*cos(x1)*cos(x2))"
u1 = "Eq(GridVariable(u1),-cos(x0)*sin(x1)*cos(x2))"
u2 = "Eq(GridVariable(u2), 0.0)"
p = "Eq(GridVariable(p), 1.0/(gama*Minf*Minf)+ (1.0/16.0) * (cos(2.0*x0)+cos(2.0*x1))*(2.0 + cos(2.0*x2)))"
r = "Eq(GridVariable(r), gama*Minf*Minf*p)"

if conservative:
    rho = "Eq(DataObject(rho), r)"
    rhou0 = "Eq(DataObject(rhou0), r*u0)"
    rhou1 = "Eq(DataObject(rhou1), r*u1)"
    rhou2 = "Eq(DataObject(rhou2), r*u2)"
    rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r *(u0**2+ u1**2 + u2**2))"
else:
    rho = "Eq(DataObject(rho), r)"
    rhou0 = "Eq(DataObject(u0), u0)"
    rhou1 = "Eq(DataObject(u1), u1)"
    rhou2 = "Eq(DataObject(u2), u2)"
    rhoE = "Eq(DataObject(Et), p/(r*(gama-1)) + 0.5*(u0**2+ u1**2 + u2**2))"    

eqns = [x0, x1, x2, u0, u1, u2, p, r, rho, rhou0, rhou1, rhou2, rhoE]

# parse the initial conditions
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
fns = 'u0 u1 u2'
cent = StoreSome(4, fns)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(3)
schemes[rk.name] = rk

boundaries = []
# Create boundaries, one for each side per dimension, so in total 6 BC's for 3D'
for direction in range(ndim):
    boundaries += [PeriodicBC(direction, 0, halos=[-2,2])]
    boundaries += [PeriodicBC(direction, 1, halos=[-2,2])]

# set the boundaries for the block
block.set_block_boundaries(boundaries)
# set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=10000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
block.setio(copy.deepcopy(h5))
# set the equations to be solved on the block

# Dispersion relation preserving filters
# DRP = ExplicitFilter(block, [0,1,2], width=11, filter_type='DRP', optimized=True, sigma=0.2, wall_control=False, multi_block=None)
# block.set_equations(DRP.equation_classes)

# WENO filter for shock-capturing
WF = WENOFilter(block, order=7, dissipation_sensor='Ducros', flux_type='LLF', airfoil=False)
block.set_equations(WF.equation_classes)

block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])
# set the discretisation schemes
block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
block.discretise()

WF.update_periodic_boundary(block, halos=[-5,5])

# create an algorithm from the discretised computations
alg = TraditionalAlgorithmRK(block)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Add the simulation monitor to the algorithm
alg = TraditionalAlgorithmRK(block)

# Write the code for the algorithm
OPSC(alg, OPS_diagnostics=2, OPS_V2=True)

# NaN check and iteration counter
constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0', 'shock_factor', 'inv_rfact0_block0', 'inv_rfact1_block0', 'inv_rfact2_block0']
values = ['1600.0', '1.4', '1.25', '0.71', '0.0003385', '59085', '256', '256', '256', '2*M_PI/block0np0', '2*M_PI/block0np1', '2*M_PI/block0np2', '1', '1.0/Delta0block0', '1.0/Delta1block0', '1.0/Delta2block0']
substitute_simulation_parameters(constants, values)
