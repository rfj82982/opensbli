#!/usr/bin/env python
# --------------------------------------------------------------------------------------------------------------------------------------------
# flatplate.py 
#               2D flat plate M6 adiabatic wall (Re_deltastar0 = 8200) simulation with central differencing
#
# author. 
# university of southampton
# --------------------------------------------------------------------------------------------------------------------------------------------
from opensbli import *
import copy
from opensbli.utilities.flat_init import Initialise_Flatplate 
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import sin, cos, sinh, tanh, exp, pi, log
import numpy as np

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# define equations		
# 
#  should i include time step and iteration number here? food for thought
#           time_step, niterations = 0.01, 1200000
# 																																		
# --------------------------------------------------------------------------------------------------------------------------------------------
# number of simulation dimensions
ndim = 2
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = 'x'

# split scheme, kgp or feiereisen
split_scheme = 'normal'

# define compressible navier stokes equations in einstein notation
mass     = 'Eq(Der(rho,t)    , - Conservative(rhou_j,x_j))'
momentum = 'Eq(Der(rhou_i,t) , - Conservative(rhou_i*u_j,x_j)  - KD(_i,_j)*Der(p,x_j)        + Der(tau_i_j,x_j))'
energy   = 'Eq(Der(rhoE,t)   , - Conservative(p*u_j,x_j)       - Conservative(rhoE*u_j, x_j) - Der(q_j,x_j)   + Der(u_i*tau_i_j ,x_j) )'

constants = ['Re', 'Pr', 'gama', 'Minf', 'SuthT', 'RefT']

# kennedy gruber scheme from splitforms
# NS = NS_Split('KGP', ndim, constants, coordinate_symbol='x', conservative=True, viscosity='dynamic', energy_formulation='enthalpy', debug=False)
# mass, momentum, energy = NS.mass, NS.momentum, NS.energy

# substitutions 
stress_tensor = 'Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j) + Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))'
heat_flux     = 'Eq(q_j    , (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))'

substitutions = [stress_tensor, heat_flux]

# constituent equations
velocity       = 'Eq(u_i , rhou_i/rho)'
pressure       = 'Eq(p   , (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))'
speed_of_sound = 'Eq(a   , (gama*p/rho)**0.5)'
temperature    = 'Eq(T   , p*gama*Minf*Minf/(rho))'

# viscocity sutherland's law
viscosity      = 'Eq(mu  , (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))'

constituent_eqns = [velocity, pressure, speed_of_sound, temperature, viscosity]

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# add equations to the kernels																													
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

# instantiate insteinEquation class for expanding the einstein indices in the equations
eq = EinsteinEquation()

base_eqns = [mass, momentum, energy]

# expand the base equations
for i, base in enumerate(base_eqns):
    base_eqns[i] = eq.expand(base, ndim, coordinate_symbol, substitutions, constants)
    
# expand the constituent relations
for i, CR in enumerate(constituent_eqns):
    constituent_eqns[i] = eq.expand(CR, ndim, coordinate_symbol, substitutions, constants)

block = SimulationBlock(ndim, block_number=0)

# create metrics before the scheme selection
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False)], 2)

# Create SimulationEquations and Constituent relations, add the expanded equations
simulation_eq = SimulationEquations()
constituent = ConstituentRelations()

for eqn in base_eqns:
    simulation_eq.add_equations(eqn)

for eqn in constituent_eqns:
    constituent.add_equations(eqn)

# Grid is stretched normal to the wall
simulation_eq.apply_metrics(metriceq)


# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# create a schemes dictionary to be used for discretisation																													
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
# Low storage optimisation for the central scheme
fns                 = 'u0 u1 T'
cent                = StoreSome(4,fns)
# cent = Central(4)
schemes[cent.name]  = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk                  = RungeKuttaLS(4)
schemes[rk.name]    = rk

# Set the discretisation schemes to be used (a python dictionary)
block.set_discretisation_schemes(schemes)

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# grid generation and initial conditions																															
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

local_dict = {'block': block, 'GridVariable': GridVariable, 'DataObject': DataObject}

# perform initial condition
# reynolds number, mach number and free-stream temperature for the initial profile
Re, xMach, Tinf, Twall = 8200.0, 6.0, 288.0, 7.021

# ensure the grid size passed to the initialisation routine matches the grid sizes used in the simulation parameters
polynomial_directions = [(False, DataObject('x0')), (True, DataObject('x1'))]
n_poly_coefficients = 50
grid_const = ['Lx1', 'by']
for con in grid_const:
    local_dict[con] = ConstantObject(con)

gridx0 = parse_expr('Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])', local_dict=local_dict)
gridx1 = parse_expr('Eq(DataObject(x1), Lx1*sinh(by*block.deltas[1]*block.grid_indexes[1]/Lx1)/sinh(by))', local_dict=local_dict)
coordinate_evaluation = [gridx0, gridx1]

# intiialise flat plate compressible similarity equations, for isothermal wall, add Twall = Twall into the function
initial = Initialise_Flatplate(polynomial_directions, n_poly_coefficients, Re, xMach, Tinf, coordinate_evaluations=coordinate_evaluation)

# --------------------------------------------------------------------------------------------------------------------------------------------																																			
# 
# boundary conditions		
#               (0,0) - inflow
#               (0,1) - outflow
#               (1,0) - wall 
#               (1,1) - farfield																																																												
#
#  --------------------------------------------------------------------------------------------------------------------------------------------

boundaries = []

boundaries += [InletTransferBC(direction=0, side=0, scheme=Carpenter())]
boundaries += [ExtrapolationBC(direction=0, side=1, order=0, scheme=Carpenter())]

wall_const = ['Minf', 'Twall']
for con in wall_const:
    local_dict[con] = ConstantObject(con)

# isothermal wall condition
rhoE_wall = parse_expr('Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))', local_dict=local_dict)
wall_eqns = [rhoE_wall]

boundaries += [IsothermalWallBC(1, 0, wall_eqns, scheme=Carpenter())]
boundaries += [ZeroGradientOutletBC(direction=1, side=1)]

# # forcing strip wall implementation ------------------------------------- 2d strip wall
# uref, dt = symbols('uref dt')
# current_iter = Globalvariable('iter', integer=True)
# x0 = DataObject('x0')
# t = dt*current_iter
# trip_eqn = exp(-((x0 - 20))**2.0/ ( 0.425*2/(0.1)))  * 1e-3 * (   sin(2*np.pi*0.12 * t )) 

# wall_normal_velocity = OpenSBLIEq(DataObject('rhou1'), (DataObject('rho'))*trip_eqn)
# boundaries += [ForcingStripBC(1, 0,  wall_normal_velocity, wall_eqns, scheme=ReducedAccess())]

block.set_block_boundaries(boundaries)

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# read/write definitions and output arrays																												
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

# monitor points
arrays = [	'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p',\
			'p', 'p', 'p' ]
arrays = [block.location_dataset('%s' % dset) for dset in arrays]
indices = [ (110, 0), (115, 0), (120, 0), \
            (125, 0), (130, 0), (135, 0), \
            (140, 0), (145, 0), (150, 0), \
		    (155, 0), (160, 0), (165, 0), \
			(170, 0), (175, 0), (180, 0), \
			(185, 0), (190, 0), (195, 0), \
            \
            (200, 0), (205, 0), (210, 0), \
            (215, 0), (220, 0), (225, 0), \
		    (230, 0), (235, 0), (240, 0), \
			(245, 0), (250, 0), (255, 0), \
			(260, 0), (265, 0), (270, 0), \
			(275, 0), (280, 0), (285, 0), \
			(290, 0), (295, 0), (300, 0), \
            \
			(305, 0), (310, 0), (315, 0), \
            (320, 0), (325, 0), (330, 0), \
		    (335, 0), (340, 0), (345, 0), \
			(350, 0), (355, 0), (360, 0), \
			(365, 0), (370, 0), (375, 0), \
			(380, 0), (385, 0), (390, 0), \
            (395, 0), (400, 0), (405, 0), \
            \
            (410, 0), (415, 0), (420, 0), \
			(425, 0), (430, 0), (435, 0), \
			(440, 0), (445, 0), (450, 0), \
			(455, 0), (460, 0), (465, 0), \
			(470, 0), (475, 0), (480, 0), \
			(485, 0), (490, 0), (495, 0), \
			\
			(500, 0), (505, 0), (510, 0), \
			(515, 0), (520, 0), (525, 0), \
			(530, 0), (535, 0), (540, 0), \
			(545, 0), (550, 0), (555, 0), \
			(560, 0), (565, 0), (570, 0), \
			(575, 0), (580, 0), (585, 0), \
			(590, 0), (595, 0), (600, 0), \
			\
			(605, 0), (610, 0), (615, 0), \
			(620, 0), (625, 0), (630, 0), \
			(635, 0), (640, 0), (645, 0), \
			(650, 0), (655, 0), (660, 0), \
			(665, 0), (670, 0), (675, 0), \
			(680, 0), (685, 0), (690, 0), \
            (695, 0), (700, 0), (705, 0), \
            \
    		(710, 0), (715, 0), (720, 0), \
            (725, 0), (730, 0), (735, 0), \
            (740, 0), (745, 0), (750, 0), \
		    (755, 0), (760, 0), (765, 0), \
			(770, 0), (775, 0), (780, 0), \
			(785, 0), (790, 0), (795, 0), \
            \
            (800, 0), (805, 0), (810, 0), \
            (815, 0), (820, 0), (825, 0), \
		    (830, 0), (835, 0), (840, 0), \
			(845, 0), (850, 0), (855, 0), \
			(860, 0), (865, 0), (870, 0), \
			(875, 0), (880, 0), (885, 0), \
			(890, 0), (895, 0), (900, 0), \
            \
			(905, 0), (910, 0), (915, 0), \
            (920, 0), (925, 0), (930, 0), \
		    (935, 0), (940, 0), (945, 0), \
			(950, 0), (955, 0), (960, 0), \
			(965, 0), (970, 0), (975, 0), \
			(980, 0), (985, 0), (990, 0), \
            \
			(995, 0), (1000, 0), (1005, 0), \
            (1010, 0), (1015, 0), (1020, 0), \
			(1025, 0), (1030, 0), (1035, 0), \
			(1040, 0), (1045, 0), (1050, 0), \
			(1055, 0), (1060, 0), (1065, 0), \
			(1070, 0), (1075, 0), (1080, 0), \
			(1085, 0), (1090, 0), (1095, 0), \
			\
			(1100, 0), (1105, 0), (1110, 0), \
			(1115, 0), (1120, 0), (1125, 0), \
			(1130, 0), (1135, 0), (1140, 0), \
			(1145, 0), (1150, 0), (1155, 0), \
			(1160, 0), (1165, 0), (1170, 0), \
			(1175, 0), (1180, 0), (1185, 0), \
			(1190, 0), (1195, 0), (1200, 0), \
			\
			(1205, 0), (1210, 0), (1215, 0), \
			(1220, 0), (1225, 0), (1230, 0), \
			(1235, 0), (1240, 0), (1245, 0), \
			(1250, 0), (1255, 0), (1260, 0), \
			(1265, 0), (1270, 0), (1275, 0), \
			(1280, 0), (1285, 0), (1290, 0), \
			(1295, 0), (1300, 0), (1305, 0), \
			\
			(1310, 0), (1315, 0), (1320, 0), \
			(1325, 0), (1330, 0), (1335, 0), \
			(1340, 0), (1345, 0), (1350, 0), \
			(1355, 0), (1360, 0), (1365, 0), \
			(1370, 0), (1375, 0), (1380, 0), \
			(1385, 0), (1390, 0), (1395, 0), \
			\
			(1400, 0), (1405, 0), (1410, 0), \
			(1415, 0), (1420, 0), (1425, 0), \
			(1430, 0), (1435, 0), (1440, 0), \
			(1445, 0), (1450, 0), (1455, 0), \
			(1460, 0), (1465, 0), (1470, 0), \
			(1475, 0), (1480, 0), (1485, 0), \
			(1490, 0), (1495, 0), (1500, 0), \
			\
			(1505, 0), (1510, 0), (1515, 0), \
			(1520, 0), (1525, 0), (1530, 0), \
			(1535, 0), (1540, 0), (1545, 0), \
			(1550, 0), (1555, 0), (1560, 0), \
			(1565, 0), (1570, 0), (1575, 0), \
			(1580, 0), (1585, 0), (1590, 0), \
			(1595, 0), (1600, 0), (1605, 0), \
			\
			(1610, 0), (1615, 0), (1620, 0), \
			(1625, 0), (1630, 0), (1635, 0), \
			(1640, 0), (1645, 0), (1650, 0), \
			(1655, 0), (1660, 0), (1665, 0), \
			(1670, 0), (1675, 0), (1680, 0), \
			(1685, 0), (1690, 0), (1695, 0), \
			\
			(1700, 0), (1705, 0), (1710, 0), \
			(1715, 0), (1720, 0), (1725, 0), \
			(1730, 0), (1735, 0), (1740, 0), \
			(1745, 0), (1750, 0), (1755, 0), \
			(1760, 0), (1765, 0), (1770, 0), \
			(1775, 0), (1780, 0), (1785, 0), \
			(1790, 0), (1795, 0), (1800, 0), \
			\
			(1805, 0), (1810, 0), (1815, 0), \
			(1820, 0), (1825, 0), (1830, 0), \
			(1835, 0), (1840, 0), (1845, 0), \
			(1850, 0), (1855, 0), (1860, 0), \
			(1865, 0), (1870, 0), (1875, 0), \
			(1880, 0), (1885, 0), (1890, 0), \
			(1895, 0), (1900, 0), (1905, 0), \
			\
			(1910, 0), (1915, 0), (1920, 0), \
			(1925, 0), (1930, 0), (1935, 0), \
			(1940, 0), (1945, 0), (1950, 0), \
			(1955, 0), (1960, 0), (1965, 0), \
			(1970, 0), (1975, 0), (1980, 0), \
			(1985, 0), (1990, 0), (1995, 0)]

SM = SimulationMonitor(arrays, indices, block, print_frequency=100, fp_precision=12, output_file='baseflow_m6.log')

# write hdf5 files
kwargs = {'iotype': 'Write'}
h5 = iohdf5(save_every=500000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('D11')])
block.setio(h5)

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# add to opensbli 																											
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

# set equations on the block and discretise
block.set_equations([constituent, simulation_eq, initial, metriceq])
# set the equations to be solved on the block
block.discretise()

alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)
SimulationDataType.set_datatype(Double)
OPSC(alg)

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# define constants																											
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

# physical constants
physical_constants   = ['gama',     'Minf',    'Pr',    'Re',    'Twall', 'SuthT',     'RefT']
physical_values      = [ '1.4', str(xMach),  '0.72', str(Re), str(Twall), '110.4',  str(Tinf)]
substitute_simulation_parameters(physical_constants, physical_values)

# numerical constants
numerical_constants  = [  'dt',   'niter', 'block0np0', 'block0np1',         'Delta0block0',          'Delta1block0',    'Lx1',   'by']
numerical_values     = ['0.01', '1200000',     '500*4',     '250*2',  '600.0/(block0np0-1)',   '100.0/(block0np1-1)',  '100.0',  '5.0']
substitute_simulation_parameters(numerical_constants, numerical_values)

print_iteration_ops(every=100, NaN_check='rho_B0')