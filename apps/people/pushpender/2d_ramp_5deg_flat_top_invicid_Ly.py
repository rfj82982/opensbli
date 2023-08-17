#!/usr/bin/env python
# --------------------------------------------------------------------------------------------------------------------------------------------
# flat13 
#
# author. 
# university of southampton
# --------------------------------------------------------------------------------------------------------------------------------------------

from opensbli import *
import copy
from opensbli.utilities.katzer_init import Initialise_Katzer
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import sin, cos, tan, sinh, cosh, tanh, exp, pi, log

#############################################################################################################################################
#																																			#
# Inputs			 																														#
#																																			#
#############################################################################################################################################
 
input_dict = {
    "gama"                 : "1.4",
    "Minf"                 : "1.5", 
    "Pr"                   : "0.72",
    "Re"                   : "1.0e50",
    "Twall"                : "1.3809973268575328",
    "dt"                   : "0.005", 
    "niter"                : "250000", 
    "block0np0"            : "500", 
    "block0np1"            : "250",    
    "Delta0block0"         : "400.0/(block0np0-1)",
    "Delta1block0"         : "115.0/(block0np1-1)",
    "SuthT"                : "110.4",
    "RefT"                 : "202.17",
    "eps"                  : "1e-15",
    "TENO_CT"              : "1e-5",
    "Lx"                   : "400.0",
    "Ly"                   : "115.0",
    "tramp"                : "tan(5.0*M_PI/180.0)",     # tan of ramp angle
    "aramp"                : "1.0",                     # parameter for ramp grid
    "xramp"                : "200.0",                   # x-location of ramp start
    "Rcorn"                : "0.01",                    # radius of corner (0.01 for sharp, 1.0 rounded but y goes slightly negative)
    "Rcurv"                : "0.1",                     # rate of increase of curvature factor
    "b"                    : "5.0",
    "A"                    : "0.05",
    "omega"                : "0.5",
    "xF"                   : "400.0",
    "yF"                   : "155.0",
    "teno_a1"              : "10.5",
    "teno_a2"              : "4.5",
    "epsilon"              : "1.0e-30",
}

constants = input_dict.keys()
values = input_dict.values()

#############################################################################################################################################
#																																			#
# Governing equations 																														#
#																																			#
#############################################################################################################################################

ndim = 2
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
# added metric class to genenerate curvilinear transformation
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, True), (True, True)], 2)

# Define the compresible Navier-Stokes equations in Einstein notation.
sc1 = "**{\'scheme\':\'Weno\'}"

a = "Conservative(detJ * rho*U_j,xi_j,%s)" % sc1
mass = "Eq(Der(rho,t), - %s/detJ)" % (a)
a = "Conservative(detJ * (rhou_i*U_j + p*D_j_i), xi_j , %s)" % sc1
momentum = "Eq(Der(rhou_i,t) , -  %s/detJ)" % (a)
a = "Conservative(detJ * (p+rhoE)*U_j,xi_j, %s)" % sc1
energy = "Eq(Der(rhoE,t), - %s/detJ)" % (a)

# auxilliary equations for NS
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"

# Substitutions
substitutions = [stress_tensor, heat_flux]
constants_var = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]
# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

# Instatiate equation classes
eq = EinsteinEquation()

# viscous momentum and energy components
visc_momentum = "Eq(Der(rhou_i, t), Der(tau_i_j,x_j))"
visc_momentum = eq.expand(visc_momentum, ndim, coordinate_symbol, substitutions, constants)
visc_momentum = [metriceq.apply_transformation(v) for v in visc_momentum]

visc_energy = "Eq(Der(rhoE, t), -Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j))"
visc_energy = eq.expand(visc_energy, ndim, coordinate_symbol, substitutions, constants)
visc_energy = metriceq.apply_transformation(visc_energy)

# create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs
eq.optional_subs_dict = optional_subs_dict

metric_vel = "Eq(U_i, D_i_j*u_j)"
# eqns = eq.expand(metric_vel, ndim, coordinate_symbol, substitutions, constants)
# for eq in eqns:
#     optional_subs_dict[eq.lhs] = eq.rhs
# Create SimulationEquations and Constituent relations, add the expanded equations
simulation_eq = SimulationEquations()
constituent = ConstituentRelations()

# change coordinate symbol to curvilinear
coordinate_symbol = "xi"
base_eqns = [mass, momentum, energy]
constituent_eqns = [velocity, pressure, speed_of_sound, temperature, viscosity, metric_vel]
# Expand the base equations
for i, base in enumerate(base_eqns):
    base_eqns[i] = eq.expand(base, ndim, coordinate_symbol, substitutions, constants_var)
    if base==momentum:
        for no, b in enumerate(base_eqns[i]):
            base_eqns[i][no] = Eq(base_eqns[i][no].lhs, base_eqns[i][no].rhs + visc_momentum[no].rhs)
    else:
        if base==energy:
            base_eqns[i] = Eq(base_eqns[i].lhs, base_eqns[i].rhs + visc_energy.rhs)

# Expand the constituent relations
for i, CR in enumerate(constituent_eqns):
    constituent_eqns[i] = eq.expand(CR, ndim, coordinate_symbol, substitutions, constants_var)

for eqn in base_eqns: 
    simulation_eq.add_equations(eqn)

for eqn in constituent_eqns:
    constituent.add_equations(eqn)

# reset coordinate symbol
coordinate_symbol = "x"

#############################################################################################################################################
#																																			#
# Schemes																																	#
#																																			#
#############################################################################################################################################

schemes = {}
weno_order = 5
# averaging procedure to be used for the eigen system evaluation
Avg = SimpleAverage([0, 1])
# LF scheme
LF = LFWeno(weno_order, formulation='Z', averaging=Avg)
# add to schemes
schemes[LF.name] = LF

#fns = 'u0 u1 T'
# cent = StoreSome(4, fns)

cent = Central(4)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(3)
schemes[rk.name] = rk

################################################################################

#############################################################################################################################################
#																																			#
# Grid and intial conditions																												#
#																																			#
#############################################################################################################################################

# Perform initial condition
## Ensure the grid size passed to the initialisation routine matches the grid sizes used in the simulation parameters
grid_const = ["Lx", "Ly", "tramp","Rcorn","Rcurv","xramp","b"]

block = SimulationBlock(ndim, block_number=0)
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
for con in grid_const:
    local_dict[con] = ConstantObject(con)

gridx0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)

gridx1=parse_expr("Eq(DataObject(x1),Ly*sinh(b*block.deltas[1]*block.grid_indexes[1]/115.0)/sinh(b)-(sinh(b*block.deltas[1]*block.grid_indexes[1]/115.0)/sinh(b)-1.0)*0.5*tramp*(DataObject(x0)-xramp)*(1.0+tanh((block.deltas[0]*block.grid_indexes[0]-xramp)/(Rcorn+Rcurv*Lx*block.deltas[1]*block.grid_indexes[1]/115.0))))", local_dict=local_dict)

coordinate_evaluation = [gridx0, gridx1]

d_in = parse_expr("Eq(GridVariable(d), 1.0)", local_dict=local_dict)
u0_in = parse_expr("Eq(GridVariable(u0), 1.0)", local_dict=local_dict)
u1_in = parse_expr("Eq(GridVariable(u1), 0.0)", local_dict=local_dict)
p_in = parse_expr("Eq(GridVariable(p), 1.0/(gama*Minf**2.0))", local_dict=local_dict)  # p = 1/(gama*Minf^2)


rho = parse_expr("Eq(DataObject(rho), d)", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), d*u0)", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), d*u1)", local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), p/(gama-1) + 0.5* d *(u0**2+u1**2))", local_dict=local_dict)

initial = GridBasedInitialisation()
initial.add_equations([gridx0, gridx1, d_in, u0_in, u1_in, p_in, rho, rhou0, rhou1, rhoE])


#############################################################################################################################################
#																																			#
# Boundary conditions																														#
#																																			#
#############################################################################################################################################

# Set the discretisation schemes to be used (a python dictionary)
block.set_discretisation_schemes(schemes)

boundaries = [[0, 0] for t in range(ndim)]

# Left Dirchlet at x= 0, inlet conditions
# Inflow conditions
direction = 0
side = 0
inlet_eq = [d_in, u0_in, u1_in, p_in, rho, rhou0, rhou1, rhoE]
boundaries[direction][side] = DirichletBC(direction, side, inlet_eq)


# Right extrapolation at outlet
direction = 0
side = 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0, scheme=ReducedAccess())

# Bottom inviscid wall
direction = 1
side = 0
boundaries[direction][side] = SymmetryBC(direction, side)

# top: extrapolation on top
direction = 1
side = 1
boundaries[direction][side] = ZeroGradientOutletBC(1, 1)

# # top: extrapolation on top
# direction = 1
# side = 1
# boundaries[direction][side] = ExtrapolationBC(direction, side, order=0)

block.set_block_boundaries(boundaries)



#############################################################################################################################################
#																																			#
# Data i/o 																																	#
#																																			#
#############################################################################################################################################

kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=10, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('D11')])
h5.add_arrays([DataObject('p')]) # save pressure
block.setio(h5)

# Set equations on the block and discretise
block.set_equations([constituent, simulation_eq, initial, metriceq])
block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho')
