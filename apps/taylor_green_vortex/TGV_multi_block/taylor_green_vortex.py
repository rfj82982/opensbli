#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

def TGV_initial_condition(block_number):
    # Create a simulation block
    block = SimulationBlock(ndim, block_number=block_number)
    # Local dictionary for parsing the expressions
    local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
    # Set the local coordinates
    if block_number == 0:
        x0 = "Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])"
        x1 = "Eq(DataObject(x1), block.deltas[1]*block.grid_indexes[1])"
    elif block_number == 1 :
        x0 = "Eq(DataObject(x0), M_PI + block.deltas[0]*block.grid_indexes[0])"
        x1 = "Eq(DataObject(x1), block.deltas[1]*block.grid_indexes[1])"
    elif block_number == 2:
        x0 = "Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])"
        x1 = "Eq(DataObject(x1), M_PI + block.deltas[1]*block.grid_indexes[1])"
    elif block_number == 3:
        x0 = "Eq(DataObject(x0), M_PI + block.deltas[0]*block.grid_indexes[0])"
        x1 = "Eq(DataObject(x1), M_PI + block.deltas[1]*block.grid_indexes[1])"
    else:
        raise NotImplementedError("The multi-block case has been defined for 4 blocks.")
    # z coordinate
    x2 = "Eq(DataObject(x2), block.deltas[2]*block.grid_indexes[2])"       
    coords = [parse_expr(eq, local_dict=local_dict) for eq in [x0, x1, x2]]

    # Initial conditions as strings
    u0 = "Eq(GridVariable(u0),sin(DataObject(x0))*cos(DataObject(x1))*cos(DataObject(x2)))"
    u1 = "Eq(GridVariable(u1),-cos(DataObject(x0))*sin(DataObject(x1))*cos(DataObject(x2)))"
    u2 = "Eq(GridVariable(u2), 0.0)"
    p = "Eq(GridVariable(p), 1.0/(gama*Minf*Minf)+ (1.0/16.0) * (cos(2.0*DataObject(x0))+cos(2.0*DataObject(x1)))*(2.0 + cos(2.0*DataObject(x2))))"
    r = "Eq(GridVariable(r), gama*Minf*Minf*p)"
    # Conservative form
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
    # Parse the initial conditions
    vortex_condition = [parse_expr(eq, local_dict=local_dict) for eq in [u0, u1, u2, p, r, rho, rhou0, rhou1, rhou2, rhoE]]
    
    init_class = GridBasedInitialisation()
    init_class.add_equations(coords + vortex_condition)
    return [init_class]

def TGV_boundaries(block_number, match_conditions):
    xm = InterfaceBC(direction=0, side=0,  match=match_conditions[block_number][0])
    xp = InterfaceBC(direction=0, side=1,  match=match_conditions[block_number][1])
    ym = InterfaceBC(direction=1, side=0,  match=match_conditions[block_number][2])
    yp = InterfaceBC(direction=1, side=1,  match=match_conditions[block_number][3])
    zm = PeriodicBC(direction=2, side=0)
    zp = PeriodicBC(direction=2, side=1)
    return [xm, xp, ym, yp, zm, zp]

# Number of dimensions of the system to be solved
ndim = 3
nblocks = 4
multi_block = MultiBlock(ndim, nblocks)
SimulationDataType.set_datatype(Double)

# Number of dimensions of the system to be solved
ndim = 3
# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "mu"]
# # symbol for the coordinate system in the equations
coordinate_symbol = "x"
# symbol for the coordinate system in the equations
conservative = False
# NS = NS_Split('Kennedy_Gruber', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='constant')
NS = NS_Split('Feiereisen', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='constant')

mass, momentum, energy = NS.mass, NS.momentum, NS.energy
# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)

# Constituent relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
if conservative:
    pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
    velocity = "Eq(u_i, rhou_i/rho)"
else:
    pressure = "Eq(p, rho*(gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j)))"

temperature = "Eq(T, p*gama*Minf*Minf/(rho))"

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object
einstein_eq = EinsteinEquation()

# Expand momentum add the expanded equations to the constituent relations
if conservative:
    eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)

# Expand pressure add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)

# Expand temperature add the expanded equations to the constituent relations
eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)

# Set the initial conditions on each of the blocks
mb_initial_conditions = {0:None, 1:None, 2:None, 3:None}
for i in range(nblocks):
    init_eq = TGV_initial_condition(i)
    mb_initial_conditions[i] = init_eq
multi_block.set_initial_conditions(mb_initial_conditions)

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
fns = 'u0 u1 u2 T'
cent = Central(4)#, fns)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKutta(3)
schemes[rk.name] = rk

# Create boundaries, one for each side per dimension, so in total 6 BC's for 3D'
mb_bcs = {0:None, 1:None, 2:None, 3:None}
# Matching conditions
match_conditions = {0: None, 1:None, 2:None, 3:None}
match_conditions[0] = [(1, 0, 1, True), (1, 0, 0, True), (2, 1, 1, True), (2, 1, 0, True)]
match_conditions[1] = [(0, 0, 1, True), (0, 0, 0, True), (3, 1, 1, True), (3, 1, 0, True)]
match_conditions[2] = [(3, 0, 1, True), (3, 0, 0, True), (0, 1, 1, True), (0, 1, 0, True)]
match_conditions[3] = [(2, 0, 1, True), (2, 0, 0, True), (1, 1, 1, True), (1, 1, 0, True)]

for i in range(nblocks):
    mb_bcs[i] = TGV_boundaries(i, match_conditions)
# set the boundaries for the block
multi_block.set_block_boundaries(mb_bcs)

# Input/output arguments
x,y,z = symbols("x0, x1, x2", **{'cls':DataObject})
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=10000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays + [x, y, z])
multi_block.setio([h5])


# Add DRP filtering on each block
filter_list = []
for no, block in enumerate(multi_block.blocks):
    filter_list += [ExplicitFilter(block, [0,1,2], width=11, filter_type='DRP', optimized=True, sigma=0.2, wall_control=False, multi_block=multi_block).equation_classes]
multi_block.set_filters(filter_list)

# set the equations to be solved on the block
multi_block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq)])
# set the discretisation schemes
multi_block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
multi_block.discretise()

# create an algorithm from the discretised computations
alg = TraditionalAlgorithmRKMB(multi_block)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg, OPS_diagnostics=5, OPS_V2=True)

# NaN check and iteration counter
print_iteration_ops(NaN_check='rho_B0')

constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter']
values = ['1600.0', '1.4', '0.1', '0.71', '0.003385', '100']
constants += ['block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0']
constants += ['block1np0', 'block1np1', 'block1np2', 'Delta0block1', 'Delta1block1', 'Delta2block1']
constants += ['block2np0', 'block2np1', 'block2np2', 'Delta0block2', 'Delta1block2', 'Delta2block2']
constants += ['block3np0', 'block3np1', 'block3np2', 'Delta0block3', 'Delta1block3', 'Delta2block3']
values += ['128', '128', '128', 'M_PI/block0np0', 'M_PI/block0np1', '2*M_PI/block0np2']
values += ['128', '128', '128', 'M_PI/block1np0', 'M_PI/block1np1', '2*M_PI/block1np2']
values += ['128', '128', '128', 'M_PI/block2np0', 'M_PI/block2np1', '2*M_PI/block2np2']
values += ['128', '128', '128', 'M_PI/block3np0', 'M_PI/block3np1', '2*M_PI/block3np2']
substitute_simulation_parameters(constants, values)
