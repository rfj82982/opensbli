#!/usr/bin/env python
# --------------------------------------------------------------------------------------------------------------------------------------------
# flatplate, simple 2D flat plate M2 simulation without shock capturing and central scheme with feieresen
#
# author. gnsa1e21, 2022
# university of southampton
# --------------------------------------------------------------------------------------------------------------------------------------------

from opensbli import *
import copy
from opensbli.utilities.flat_init import Initialise_Flatplate 
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import sin, cos, sinh, tanh, exp, pi, log


# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# define equations																														
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

ndim = 2

# Define the compresible Navier-Stokes equations in Einstein notation
# Feiereisen quadratic skew-symmetric formulation, no change in continuity
mass = "Eq(Der(rho, t), - Der(rhou_j, x_j))"

# Feiereisen quadratic skew-symmetric momentum
QSSFm = "(1/2) * (Conservative(rhou_i*u_j, x_j) + rhou_j* Der(u_i,x_j) + u_i * Der(rhou_j,x_j))"
momentum = "Eq(Der(rhou_i, t), - Der(p, x_i) + Der(tau_i_j, x_j) - KD(_i,_j)*c_j )"

QSSFe = "(1/2) * (Conservative(rhoE*u_j, x_j) + rhou_j*Conservative(rhoE/rho, x_j) + (rhoE/rho) * Der(rhou_j, x_j))"
energy = "Eq(Der(rhoE, t), - %s - Conservative(p*u_j, x_j) - Dot(c_j, u_j) + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j) )" % (QSSFe)

stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"

# Substitutions
substitutions = [stress_tensor, heat_flux]
constants = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

# Instatiate equation classes
eq = EinsteinEquation()

# Create SimulationEquations and Constituent relations, add the expanded equations
simulation_eq = SimulationEquations()
constituent = ConstituentRelations()

# Expand momentum add the expanded equations to the simulation equations
expanded_Feiereisen = eq.expand(QSSFm, ndim, coordinate_symbol, substitutions, constants)
eqns = eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
# Substract the inviscid part ot the RHS
for no, value in enumerate(eqns):
    eqns[no] = Eq(eqns[no].lhs,  eqns[no].rhs - expanded_Feiereisen[no])
simulation_eq.add_equations(eqns)

base_eqns = [mass, energy]
constituent_eqns = [velocity, pressure, speed_of_sound, temperature, viscosity]
# Expand the base equations
for i, base in enumerate(base_eqns):
    base_eqns[i] = eq.expand(base, ndim, coordinate_symbol, substitutions, constants)
# Expand the constituent relations
for i, CR in enumerate(constituent_eqns):
    constituent_eqns[i] = eq.expand(CR, ndim, coordinate_symbol, substitutions, constants)

block = SimulationBlock(ndim, block_number=0)

# Create metrics before the scheme selection
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False)], 2)


for eqn in base_eqns:
    simulation_eq.add_equations(eqn)

for eqn in constituent_eqns:
    constituent.add_equations(eqn)

# Grid is stretched normal to the wall
simulation_eq.apply_metrics(metriceq)


# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# assign central scheme																														
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
# Low storage optimisation for the central scheme
fns = 'u0 u1 T'
cent = StoreSome(4, fns)
# cent = Central(4)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(3)
schemes[rk.name] = rk
# Set the discretisation schemes to be used (a python dictionary)
block.set_discretisation_schemes(schemes)

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# boundary conditions																														
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------


local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

x_loc = parse_expr("Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)

rho = parse_expr("Eq(DataObject(rho), d)", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), d*u0)", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), d*u1)", local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), p/(gama-1) + 0.5* d *(u0**2+u1**2))", local_dict=local_dict)

boundaries = [[0, 0] for t in range(ndim)]
# Left pressure extrapolation at x= 0, inlet conditions
direction = 0
side = 0
boundaries[direction][side] = InletTransferBC(0, 0, scheme=ReducedAccess())
# Right extrapolation at outlet
direction = 0
side = 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0, scheme=ReducedAccess())
# boundaries[direction][side] = ZeroGradientOutletBC(0, 1)
# Bottom no-slip isothermal wall
direction = 1
side = 0
wall_const = ["Minf", "Twall"]
for con in wall_const:
    local_dict[con] = ConstantObject(con)
# Isothermal wall condition
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
wall_eqns = [rhoE_wall]
boundaries[direction][side] = IsothermalWallBC(1, 0, wall_eqns, scheme=ReducedAccess())
# Top dirichlet shock generator condition
direction = 1
side = 1

boundaries[direction][side] = ZeroGradientOutletBC(1, 1)

block.set_block_boundaries(boundaries)

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# initial condition																														
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

# Perform initial condition
# Reynolds number, Mach number and free-stream temperature for the initial profile, additional Twall added
Re, xMach, Tinf, Twall = 950.0, 2.0, 288.0, 1.67619431
## Ensure the grid size passed to the initialisation routine matches the grid sizes used in the simulation parameters
polynomial_directions = [(False, DataObject('x0')), (True, DataObject('x1'))]
n_poly_coefficients = 50
grid_const = ["Lx1", "by"]
for con in grid_const:
    local_dict[con] = ConstantObject(con)
gridx0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
gridx1 = parse_expr("Eq(DataObject(x1), Lx1*sinh(by*block.deltas[1]*block.grid_indexes[1]/Lx1)/sinh(by))", local_dict=local_dict)
coordinate_evaluation = [gridx0, gridx1]

# Isothermal wall with a specified wall temperature (works for a larger range),
# initial = Initialise_Flatplate(polynomial_directions, n_poly_coefficients, Re, xMach, Tinf, Twall=Twall, adiabaticwall_condition=True, coordinate_evaluations=coordinate_evaluation)

# Adiabtic wall, the code produce default conditions for an adiabatic wall.
initial = Initialise_Flatplate(polynomial_directions, n_poly_coefficients, Re, xMach, Tinf, coordinate_evaluations=coordinate_evaluation)

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# read/write options																														
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=1000000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('D11')])
block.setio(h5)

# Add SFD filtering
# SFD = SFD(block, chifilt=0.1, omegafilt=1.0/0.75)
j = block.grid_indexes[1]
grid_condition = j >= 169
F = BinomialFilter(block, order=10, grid_condition=grid_condition)

# Set equations on the block and discretise
block.set_equations([constituent, simulation_eq, initial, metriceq])
block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
# Substitute simulation parameter values

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# define consstants
#                   Note: str(xMach), str(Twall), str(Re) and str(Tinf) staying consistent with the definitions above																														#
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

constants = ['gama', 'Minf', 'Pr', 'Re', 'Twall', 'dt', 'niter', 'block0np0', 'block0np1',
                 'Delta0block0', 'Delta1block0', 'SuthT', 'RefT', 'eps', 'Lx1', 'by', 'epsilon']
values = ['1.4', str(xMach), '0.72', str(Re), str(Twall), '0.04', '25000000', '500', '250',
              '400.0/(block0np0-1)', '115.0/(block0np1-1)', '110.4', str(Tinf), '1e-15', '115.0', '5.0', '1.0e-30']
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho', every=1000)
