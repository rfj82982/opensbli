#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

simulation_parameters = {
'Re'        :   '1600.0',
'gama'      :   '1.4',
'Minf'      :   '1.25',
'Pr'        :   '0.71',
'dt'        :   '0.0005',
'niter'     :   '40000',
'block0np0'     :   '256',
'block0np1'     :   '256',
'block0np2'     :   '256',
'Delta0block0'      :   '2*M_PI/block0np0',
'Delta1block0'      :   '2*M_PI/block0np1',
'Delta2block0'      :   '2*M_PI/block0np2',
'shock_factor'      :   '1',
'inv_rfact0_block0'     :   '1.0/Delta0block0',
'inv_rfact1_block0'     :   '1.0/Delta1block0',
'inv_rfact2_block0'     :   '1.0/Delta2block0',
'TENO_CT'       :   '1e-6',
'eps'       :   '1.0e-30',
'teno_a1'       :   '10.5',
'teno_a2'       :   '4.5',
}

# Number of dimensions of the system to be solved
ndim = 3
# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf"]
coordinate_symbol = "x"
# symbol for the coordinate system in the equations
conservative = True
teno = True
einstein_eq = EinsteinEquation()
# Central scheme plus WENO filtering, otherwise pure TENO
if not teno:
    NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic', energy_formulation='enthalpy', debug=False)
    mass, momentum, energy = NS.mass, NS.momentum, NS.energy
    # Expand the simulation equations, for this create a simulation equations class
    simulation_eq = SimulationEquations()
    simulation_eq.add_equations(mass)
    simulation_eq.add_equations(momentum)
    simulation_eq.add_equations(energy)

else:
    sc1 = "**{\'scheme\':\'Teno\'}"
    # Define the compresible Navier-Stokes equations in Einstein notation.
    mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))" % sc1
    momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s) + Der(tau_i_j,x_j) )" % sc1
    energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )" % sc1
    stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
    heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
    # Substitutions
    substitutions = [stress_tensor, heat_flux]
    constants = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]
    base_eqns = [mass, momentum, energy]
    # Expand the base equations
    for i, base in enumerate(base_eqns):
        base_eqns[i] = einstein_eq.expand(base, ndim, coordinate_symbol, substitutions, constants)
    simulation_eq = SimulationEquations()
    for eqn in base_eqns:
        simulation_eq.add_equations(eqn)

# Constituent relations
pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
velocity = "Eq(u_i, rhou_i/rho)"
enthalpy = "Eq(H, (rhoE + p) / rho)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.4042)/(T+0.40417)))" ## Modified sutherland law
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object
for input_eqn in [velocity, pressure, temperature, enthalpy, viscosity, speed_of_sound]:
    eqns = einstein_eq.expand(input_eqn, ndim, coordinate_symbol, [], constants)
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

store_grid = True
if store_grid:
    eqns += ["Eq(DataObject(x0), GridVariable(x0))"]
    eqns += ["Eq(DataObject(x1), GridVariable(x1))"]
    eqns += ["Eq(DataObject(x2), GridVariable(x2))"]

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
rk = RungeKuttaLS(3, formulation='SSP')
schemes[rk.name] = rk

if teno:
    teno_order = 6
    Avg = RoeAverage([0, 1])
    # Adaptive changing of TENO CT coefficients paired with a shock-sensor
    adaptive = False    
    if adaptive:
        SS = ShockSensor()
        shock_sensor, sensor_array = SS.ducros_equations(block, coordinate_symbol)
        # Add shock Ducros sensor to constituent relations
        constituent.add_equations(shock_sensor)
        LF = LFTeno(teno_order, averaging=Avg, flux_type='LLF', formulation='adaptive', sensor=sensor_array, store_sensor=True)
    else:
        # LF = LFTeno(teno_order, averaging=Avg, flux_type='LLF')
        LF = HLLCTeno(teno_order, averaging=Avg, flux_type='HLLC-LM')
    schemes[LF.name] = LF

boundaries = []
# Create boundaries, one for each side per dimension, so in total 6 BC's for 3D'
if not teno:
    halos = [-2,2]
else:
    halos = [-4,4]
# Set periodic boundary with desired halo depth
for direction in range(ndim):
    boundaries += [PeriodicBC(direction, 0, halos=halos)]
    boundaries += [PeriodicBC(direction, 1, halos=halos)]

q_vector = simulation_eq.time_advance_arrays
# set the boundaries for the block
block.set_block_boundaries(boundaries)
# set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=10000, **kwargs)
h5.add_arrays(q_vector + [DataObject('x0'), DataObject('x1'), DataObject('x2')])

if not teno:
    h5.add_arrays([DataObject('kappa')] + [DataObject('q0')]) 
# Write out a grid file
if store_grid:
    h5_grid = iohdf5(**kwargs)
    h5_grid.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('x2')])
block.setio(copy.deepcopy(h5))

# Dispersion relation preserving filters
# DRP = ExplicitFilter(block, [0,1,2], width=11, filter_type='DRP', optimized=True, sigma=0.2, wall_control=False, multi_block=None)
# block.set_equations(DRP.equation_classes)

## Post-processing for TGV case, kinetic energy and enstrophy reductions
# Velocity in 3D
vel = symbols("u0:%d"%ndim,  **{'cls':DataObject})
# Vorticity-z
wx, wy, wz = symbols("wx wy wz",  **{'cls':GridVariable})
# coordinates
coord = symbols("x0:%d"%ndim,  **{'cls':CoordinateObject})
# Matrix of derivatives
der_matrix = Matrix(ndim,ndim,[CentralDerivative(u,x) for u in vel for x in coord])
post = UserDefinedEquations()
post.kernel_merge = True
post.algorithm_place = InTheSimulation(frequency=100)
post.computation_name = 'Taylor-Green vortex post-processing'
post.order = 10000000 # appear at the end of the kernels at the end of the time-loop
# # X vorticity
vortx = der_matrix[2,1] - der_matrix[1,2]
# vortx = metriceq.apply_transformation(vortx)
vortx = Eq(wx, vortx)
post.add_equations(vortx)
# Y vorticity
vorty = der_matrix[0,2] - der_matrix[2,0]
# vorty = metriceq.apply_transformation(vorty)
vorty = Eq(wy, vorty)
post.add_equations(vorty)
# Z vorticity
vortz = der_matrix[1,0] - der_matrix[0,1]
# vortz = metriceq.apply_transformation(vortz)
vortz = Eq(wz, vortz)
post.add_equations(vortz)
# # Dilatation
divV = symbols("divV", **{'cls':GridVariable})
dil = Eq(divV, der_matrix[0,0] + der_matrix[1,1] + der_matrix[2,2])
post.add_equations(dil)

# Evaluate quantities required for dissipation measures
# rhom, KE, eps_D, eps_S = ReductionVariable('rhom', **{'reduction_type' : 'sum'}), ReductionVariable('KE', **{'reduction_type' : 'sum'}), ReductionVariable('dilatation_dissipation', **{'reduction_type' : 'sum'}), ReductionVariable('enstrophy_dissipation', **{'reduction_type' : 'sum'})
rhom, KE, eps_D, eps_S = ReductionSum('rhom'), ReductionSum('KE'), ReductionSum('dilatation_dissipation'), ReductionSum('enstrophy_dissipation')
rho_eqn = OpenSBLIEq(rhom, rhom + DataObject('rho'))
ke_eqn = OpenSBLIEq(KE, KE + 0.5*DataObject('rho')*sum([u**2 for u in vel]))
dilatation_eqn = OpenSBLIEq(eps_D, eps_D + Rational(4,3)*DataObject('mu')*divV**2)
enstrophy_eqn = OpenSBLIEq(eps_S, eps_S + DataObject('mu')*(wx**2 + wy**2 + wz**2))
post.add_equations([rho_eqn, ke_eqn, dilatation_eqn, enstrophy_eqn])

block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial, post])
# set the discretisation schemes
block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
block.discretise()
# Apply a periodic BC for WENO filter
if not teno:
    WF.update_periodic_boundary(block, halos=[-5,5])

# Simulation monitor
arrays = ['KE', 'dilatation_dissipation', 'enstrophy_dissipation', 'rhom']
probe_locations = [(None), (None), (None), (None)]
SM = SimulationMonitor(arrays, probe_locations, block, output_file='TGV.log', print_frequency=100)
# Add the simulation monitor to the algorithm
alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg, OPS_diagnostics=2, OPS_V2=True)
# Add the simulation constants to the OPS C code
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check='rho')
