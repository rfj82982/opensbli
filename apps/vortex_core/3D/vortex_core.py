#!/usr/bin/env python

# Import all the functions from opensbli
from opensbli import *
from sympy import sin, cos, sinh, tanh, exp, pi
#import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

simulation_parameters = {
'Re'        :   '500000.0',   
'gama'      :   '1.4',   
'Minf'      :   '0.5',   
'Pr'        :   '0.72',   
'Sc'        :   '1.0',   
'dt'        :   '0.0025',   
'niter'     :   '200000',   
'block0np0'     :   '1500',   
'block0np1'     :   '521',   
'block0np2'     :   '125',
'Lx'        : '10.0',
'Ly'        :   '80.0',
'Lz'        : '1.0',
'Delta0block0'      :   'Lx/(block0np0)',   
'Delta1block0'      :   'Ly/(block0np1-1)',   
'Delta2block0'      :   'Lz/(block0np2)',
'stretch'       :   '4.0',
'y_0'        : '0.1',
'y_1'        : '0.1',
'U_0'        : '1.0',
'RN_amplitude' : '1.0',
}

viscosity_relation = 'constant'
# Instantiate EinsteinEquation class for expanding the Einstein indices in the equations
einstein_eq = EinsteinEquation()
# symbol for the coordinate system in the equations
coordinate_symbol = "x"
# Constants that are used
constants = ["Re", "Pr", "gama", "Minf"]
# Number of dimensions of the system to be solved
ndim = 3
stats = False
# Define the compresible Navier-Stokes equations in Einstein notation, by default the scheme is Central no need to
# NS = NS_Split('Feiereisen', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=True, viscosity='dynamic', energy_formulation='enthalpy', debug=False)
NS = NS_Split('KEEP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=True, viscosity=viscosity_relation, energy_formulation='enthalpy', debug=False)
mass, momentum, energy = NS.mass, NS.momentum, NS.energy
# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)

# Constituent relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
enthalpy = "Eq(H, (rhoE + p) / rho)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
internal_energy = "Eq(e, p / (rho*(gama-1)))"
# viscosity = "Eq(mu, T**0.7)"

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object

CRs = [velocity, pressure, enthalpy, temperature, internal_energy]
for eqn in CRs:
    eqns = einstein_eq.expand(eqn, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# initial conditions
dx, dy, dz = block.deltas
i, j, k = block.grid_indexes
gama, Minf = symbols('gama Minf', **{'cls': ConstantObject})
# Mixing layer parameters
y_0, y_1, U_0, RN_amplitude = symbols('y_0 y_1 U_0 RN_amplitude', **{'cls': ConstantObject})
q_vector = flatten(simulation_eq.time_advance_arrays)
x, y, z = symbols('x0:%d' % ndim, **{'cls': DataObject})
nx, ny, nz, stretch = symbols('block0np0 block0np1 block0np2 stretch', **{'cls': ConstantObject})
Lx, Ly, Lz = symbols('Lx Ly Lz', **{'cls': ConstantObject})
grid_equations= []
# Stretched in y between slip conditions
stretch_eqn=0.5*Ly*sinh(stretch*(j-(ny-1)/2)/((ny-1)/2))/sinh(stretch)
grid_equations += [Eq(x, -Lx/2.0 + i*dx), Eq(y,stretch_eqn), Eq(z, -Lz/2.0 + k*dz)]

initial_equations = []
rho, u, v, w, p, T = symbols('rho, u, v, w, p, T', **{'cls': GridVariable})

initial_equations += [Eq(u, U_0*tanh((y - y_1*cos(2*pi*(x+Lx/4)/Lx))/y_0))]
# initial_equations += [Eq(vpert,0.1*cos((2.0*pi*(x+Lx/4))/Lx)*exp(-y**2.0/10.0))]
conditions = Piecewise((RN_amplitude*DataObject('random_nums'), Abs(y) < y_0/2), (0, True))
initial_equations += [Eq(v, conditions)]
initial_equations += [Eq(w, 0.0)]
# initial_equations += [Eq(T,1.0+Minf**2*(gama-1.0)/2.0*(1.0-u**2))]
initial_equations += [Eq(T,1.0)]
initial_equations += [Eq(p,1.0/(gama*Minf**2.0))]
initial_equations += [Eq(rho,gama*Minf**2*p/T)]

for eqn in initial_equations:
    pprint(eqn)
# exit()
initial_equations += [Eq(q_vector[0],rho)]
initial_equations += [Eq(q_vector[1],rho*u)]
initial_equations += [Eq(q_vector[2],rho*v)]
initial_equations += [Eq(q_vector[3],rho*w)]
initial_equations += [Eq(q_vector[4],p/(gama-1.0)+0.5*rho*(u**2+v**2+w**2))]

# initial conditions
initial = GridBasedInitialisation()
initial.add_equations(grid_equations + initial_equations)

# metrics
metriceq =  MetricsEquation()
#metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False), (False, False)], 2)
simulation_eq.apply_metrics(metriceq)

# Create a schemes dictionary to be used for discretisation
schemes = {}
fns = 'u0 u1 u2'
cent = StoreSome(4, fns)
# cent = Central(4)
schemes[cent.name] = cent
rk = RungeKuttaLS(4)
schemes[rk.name] = rk

# Create boundaries, one for each side per dimension, so in total 6 BC's for 3D'
boundaries = []
direction=0
boundaries += [PeriodicBC(direction, side=0, halos=[-2,2])]
boundaries += [PeriodicBC(direction, side=1, halos=[-2,2])]
direction=1
boundaries += [SymmetryBC(direction, 0)]
boundaries += [SymmetryBC(direction, 1)]
direction=2
boundaries += [PeriodicBC(direction, side=0, halos=[-2,2])]
boundaries += [PeriodicBC(direction, side=1, halos=[-2,2])]
block.set_block_boundaries(boundaries)

# Post processing
# Velocity in 3D
vel = symbols("u0:%d"%ndim,  **{'cls':DataObject})
# Vorticity-z
wx, wy, wz = symbols("wx wy wz",  **{'cls':DataObject})
dudy_lhs = symbols("dudy", **{'cls' : DataObject})
# coordinates
coord = symbols("x0:%d"%ndim,  **{'cls':CoordinateObject})
# Matrix of derivatives
der_matrix = Matrix(ndim,ndim,[CentralDerivative(u,x) for u in vel for x in coord])
post = UserDefinedEquations()
post.kernel_merge = True

post.algorithm_place = InTheSimulation(frequency=100)
post.computation_name = 'Vortex core post-processing'
post.order = 10000000 # appear at the end of the kernels at the end of the time-loop

dudy = OpenSBLIEq(dudy_lhs, der_matrix[0,1])
dudy = metriceq.apply_transformation(dudy)
dudy_max = ReductionMax('dudy_max')
vorticity_thickness = OpenSBLIEq(dudy_max, dudy_lhs)
metriceq.apply_transformation(vorticity_thickness)

# # # X vorticity
vortx = der_matrix[2,1] - der_matrix[1,2]
# vortx = metriceq.apply_transformation(vortx)
vortx = Eq(wx, vortx)
post.add_equations(vortx)
# # Y vorticity
vorty = der_matrix[0,2] - der_matrix[2,0]
# vorty = metriceq.apply_transformation(vorty)
vorty = Eq(wy, vorty)
post.add_equations(vorty)
# # Z vorticity
vortz = der_matrix[1,0] - der_matrix[0,1]
vortz = metriceq.apply_transformation(vortz)
vortz = Eq(wz, vortz)
# # # Dilatation
divV = symbols("divV", **{'cls':DataObject})
dil = Eq(divV, der_matrix[0,0] + der_matrix[1,1] + der_matrix[2,2])
dil = metriceq.apply_transformation(dil)

# # Evaluate quantities required for dissipation measures
if viscosity_relation == 'constant':
    mu = 1
else:
    mu = DataObject('mu')
rho_m, KE, eps_D, eps_S = ReductionSum('rhom'), ReductionSum('KE'), ReductionSum('dilatation_dissipation'), ReductionSum('enstrophy_dissipation')
rho_eqn = OpenSBLIEq(rho_m, rho_m + DataObject('rho'))
ke_eqn = OpenSBLIEq(KE, KE + 0.5*DataObject('rho')*sum([u**2 for u in vel]))
# dilatation_eqn = OpenSBLIEq(eps_D, eps_D + Rational(4,3)*mu*divV**2)
# enstrophy_eqn = OpenSBLIEq(eps_S, eps_S + mu*(wx**2 + wy**2 + wz**2))
# post.add_equations([dudy, vorticity_thickness, dil, vortz, rho_eqn, ke_eqn, dilatation_eqn, enstrophy_eqn])
post.add_equations([dudy, vorticity_thickness, vortz, rho_eqn, ke_eqn])

# Dispersion relation preserving filters
DRP = ExplicitFilter(block, [0,1,2], width=11, filter_type='DRP', frequency=25000000, optimized=True, sigma=0.1, airfoil=False, multi_block=None)

# set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=500, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays + [x, y, z] + [DataObject('dudy'), DataObject('wz')])
block.setio(copy.deepcopy(h5))

# Read in random numbers
kwargs = {'iotype': "Read"}
h5_read = iohdf5(**kwargs)
h5_read.add_arrays([DataObject('random_nums')])
block.setio([h5_read])

# set the equations to be solved on the block
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial, metriceq, post] + DRP.equation_classes)
# set the discretisation schemes
block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
block.discretise()

#Add some full [-5,5] halo swaps over the periodic directions only when the filter is called
def create_exchange_calls_codes(block, dsets):
    kernels = []
    arrays = [block.location_dataset(a) for a in flatten(dsets)]
    for direction in [0,2]:
        for side in [0,1]:
            BC = PeriodicBC(direction, side, halos=[-5, 5], full_depth=True)
            kernels += [BC.apply(arrays, block)]
    return kernels

# Make some full swaps for interfaces before filtering
filter_swaps = create_exchange_calls_codes(block, ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE'])
for no, eq in enumerate(block.list_of_equation_classes):
    if isinstance(eq, UserDefinedEquations):
        if eq.full_swap:
            eq.Kernels += filter_swaps

arrays = ['KE', 'rhom']
probe_locations = [(None), (None)]
SM = SimulationMonitor(arrays, probe_locations, block, output_file='vortex_history.log', print_frequency=100)

# create an algorithm from the discretised computations
alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg)
# Add the simulation constants to the OPS C code
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check='rho')
