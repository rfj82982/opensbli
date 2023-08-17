#!/usr/bin/env python
# --------------------------------------------------------------------------------------------------------------------------------------------
# mixflat_N2 - pure N2 frozen flow, should be a dimensional extension of flatplate
#
# author. gnsa1e21, 2023
# university of southampton
# --------------------------------------------------------------------------------------------------------------------------------------------

# Import all the functions from opensbli
from opensbli import *
from sympy import sin, cos, sinh, tanh, exp, pi, log
#import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from opensbli.utilities.flatmix_init import Initialise_Flatmix
# from opensbli.utilities.mixflat_init import Initialise_MixFlat

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# define equations																														
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------
# Number of dimensions of the system to be solved
ndim = 2
stats = False
# Define the compresible Navier-Stokes equations in Einstein notation
mass_O = "Eq(Der(rhoO,t), - Skew(rhoO*u_j,x_j) + Der(mu/(Re*Sc)*Der(yO,x_j),x_j)  )" #+wdotO
mass_O2 = "Eq(Der(rhoO2,t), - Skew(rhoO2*u_j,x_j) + Der(mu/(Re*Sc)*Der(yO2,x_j),x_j)  )" #+wdotO2
mass_N = "Eq(Der(rhoN,t), - Skew(rhoN*u_j,x_j) + Der(mu/(Re*Sc)*Der(yN,x_j),x_j)  )" #+wdotN
mass_N2 = "Eq(Der(rhoN2,t), - Skew(rhoN2*u_j,x_j) + Der(mu/(Re*Sc)*Der(yN2,x_j),x_j)  )" #+wdotN2
mass_NO = "Eq(Der(rhoNO,t), - Skew(rhoNO*u_j,x_j) + Der(mu/(Re*Sc)*Der(yNO,x_j),x_j)  )" #+wdotNO
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j))"
evib = "Eq(Der(rhoev,t), - Skew(rhoev*u_j,x_j) + (rhoO2*eveqO2+rhoN2*eveqN2+rhoNO*eveqNO - rho*ev)/tau +Der(mu/(Re*Sc)*(evO2*Der(yO2,x_j)+evN2*Der(yN2,x_j)+evNO*Der(yNO,x_j)),x_j) - Der(qv_j,x_j) )" # + wdotO2*evO2+wdotN2*evN2+wdotNO*evNO
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) +Der(mu/(Re*Sc)*Rhat*T*(5.0/(2.0*MO)*Der(yO,x_j)+5.0/(2.0*MN)*Der(yN,x_j)+7.0/(2.0*MO2)*Der(yO2,x_j)+7.0/(2.0*MN2)*Der(YN2,x_j)+7.0/(2.0*MNO)*Der(YNO,x_j)),x_j) + Der(mu/(Re*Sc)*4.1868e6*(dhO/MO*Der(YO,x_j)+dhN/MN*Der(YN,x_j)+dhNO/MNO*Der(YNO,x_j)),x_j) + Der(mu/(Re*Sc)*(evO2*Der(yO2,x_j)+evN2*Der(yN2,x_j)+evNO*Der(yNO,x_j)),x_j) - Der(q_j,x_j) - Der(qv_j,x_j) + Der(u_i*tau_i_j ,x_j))"
scalar = "Eq(Der(rhof,t), - Skew(rhof*u_j,x_j) + Der(mu/(Re*Sc)*Der(f,x_j),x_j))" # non-reacting scalar is useful as a (diffusing) marker of original fluid regions

# Substitutions used in the equations
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))" 
# note that we need to keep Re (set=1) as this is used in scheme.py to split terms into convective and viscous (Q for DL - does this apply to conduction and species diffusion too?)
heat_flux = "Eq(q_j, -(kappa/Re)*Der(T,x_j))"
heat_flux_vib = "Eq(qv_j, -(kappav/Re)*Der(Tv,x_j))" 
evibration = "Eq(ev, rhoev/rho)"
#evO2 = "Eq(evO2, thetavO2*Rhat/(MO2*(exp(thetavO2/Tv)-1.0)))"
#evN2 = "Eq(evN2, thetavN2*Rhat/(MN2*(exp(thetavN2/Tv)-1.0)))"
#evNO = "Eq(evNO, thetavNO*Rhat/(MNO*(exp(thetavNO/Tv)-1.0)))"
density = "Eq(rho, (rhoO+rhoO2+rhoN+rhoN2+rhoNO))"
molesum = "Eq(ysum, (rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO))"
molesumM = "Eq(ysumM, (rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))"
hformation = "Eq(dhf, 4.1868e6*(dhO*rhoO/MO+dhN*rhoN/MN+dhNO*rhoNO/MNO))"
timeconst = "Eq(tau, 0.1 * (rhoO2/MO2+rhoN2/MN2+rhoNO/MNO)*101325.0/(p*(rhoO2/(MO2*ptauO2)+rhoN2/(MN2*ptauN2)+rhoNO/(MNO*ptauNO))))"
thetavset = "Eq(thetavnum, (thetavO2*rhoO2/MO2+thetavN2*rhoN2/MN2+thetavNO*rhoNO/MNO))"


# make substitutions
substitutions = [stress_tensor, heat_flux, heat_flux_vib, evibration, density, molesum, molesumM, hformation, timeconst, thetavset]


# Constants that are used
constants = ["Re", "Sc", "uref", "pref", "rhoref", "pexp", "Rhat", "MO", "MO2", "MN", "MN2", "MNO","rhoEref"]
constants=constants+["dhO", "dhN", "dhNO", "thetavO2", "thetavN2", "thetavNO"]


# symbol for the coordinate system in the equations
coordinate_symbol = "x"

# Variable relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
mixturefraction = "Eq(f, rhof/rho)"
# energyam = "Eq(E, rhoE/rho)"				# Added to involve rhoE in the equaation
pressure = "Eq(p, Rhat*T*(rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO))"
temperature = "Eq(T, (rhoE -rhoev - dhf - rho*(1./2.)*(KD(_i,_j)*u_i*u_j))/(Rhat*(3.0/2.0*(rhoO/MO+rhoN/MN)+5.0/2.0*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))) )"
#temperature = "Eq(T, (rhoE - 4.1868e6*(dhO*rhoO/MO+dhN*rhoN/MN+dhNO*rhoNO/MNO)- rho*(1./2.)*(KD(_i,_j)*u_i*u_j))/(Rhat*(3/2*(rhoO/MO+rhoN/MN)+5/2*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))) )"

# vibrational terms
tempv = "Eq(Tv, thetavnum/(ysumM*log(1.0+thetavnum*Rhat/(rhoev))))" # method to find Tv based on a mole-weighted thetav (compare with N-R or Cv-based method later)
#tempc = "Eq(Tw, rhothetav/(rho*log(1.0+rhothetav*Rhat*ysum/(rho**2*ev))))"
#Tchem = "Eq(Tc, Tv**0.3)"
#Tchem = "Eq(Tc, (rhothetav/(rho*log(1.0+rhothetav*Rhat*ysum/(rho**2*ev))))**0.3*((rhoE -rhoev - dhf - rho*(1./2.)*(KD(_i,_j)*u_i*u_j))/(Rhat*(3.0/2.0*(rhoO/MO+rhoN/MN)+5.0/2.0*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))))**0.7)" # OpenSBLI seems to prevent this being shortened - it is only (Tv**0.3)*(T**0.7)
evequilO2 = "Eq(eveqO2, thetavO2*Rhat/(MO2*(exp(thetavO2/T)-1.0)))"
evequilN2 = "Eq(eveqN2, thetavN2*Rhat/(MN2*(exp(thetavN2/T)-1.0)))"
evequilNO = "Eq(eveqNO, thetavNO*Rhat/(MNO*(exp(thetavNO/T)-1.0)))"
evvO2 = "Eq(evO2, thetavO2*Rhat/(MO2*(exp(thetavO2/Tv)-1.0)))"
evvN2 = "Eq(evN2, thetavN2*Rhat/(MN2*(exp(thetavN2/Tv)-1.0)))"
evvNO = "Eq(evNO, thetavNO*Rhat/(MNO*(exp(thetavNO/Tv)-1.0)))"
timefactorO2 = "Eq(ptauO2, (rhoO/MO*exp(129.0*(T**(-1.0/3.0)-0.0271)-18.42)+rhoO2/MO2*exp(129.0*(T**(-1.0/3.0)-0.0300)-18.42)+rhoN/MN*exp(129.0*(T**(-1.0/3.0)-0.0265)-18.42)+rhoN2/MN2*exp(129.0*(T**(-1.0/3.0)-0.0295)-18.42)+rhoNO/MNO*exp(129.0*(T**(-1.0/3.0)-0.0298)-18.42))/ysum)" 
timefactorN2 = "Eq(ptauN2, (rhoO/MO*exp(220.0*(T**(-1.0/3.0)-0.0268)-18.42)+rhoO2/MO2*exp(220.0*(T**(-1.0/3.0)-0.0295)-18.42)+rhoN/MN*exp(220.0*(T**(-1.0/3.0)-0.0262)-18.42)+rhoN2/MN2*exp(220.0*(T**(-1.0/3.0)-0.0290)-18.42)+rhoNO/MNO*exp(220.0*(T**(-1.0/3.0)-0.0293)-18.42))/ysum)" 
timefactorNO = "Eq(ptauNO, (rhoO/MO*exp(168.0*(T**(-1.0/3.0)-0.0270)-18.42)+rhoO2/MO2*exp(168.0*(T**(-1.0/3.0)-0.0298)-18.42)+rhoN/MN*exp(168.0*(T**(-1.0/3.0)-0.0264)-18.42)+rhoN2/MN2*exp(168.0*(T**(-1.0/3.0)-0.0293)-18.42)+rhoNO/MNO*exp(168.0*(T**(-1.0/3.0)-0.0295)-18.42))/ysum)" 

# fluid properties
#viscosity = "Eq(mu, 4.644e-7*T**0.65)" # Power laws from Hirshel (good for 300<T<2000)
#conductivity = "Eq(kappa, 3.4957e-4*T**0.75)" 
viscosity = "Eq(mu, ((yO2+yN2+yNO)*0.1*exp(-11.2202)*T**(0.021823*log(T)+0.34357)+(yO+yN)*0.1*exp(-11.7344)*T**(0.022652*log(T)+0.342509))*(1.0-exp(-0.010568*T)) )" # NDS simplication of Blottner/Sutherland
conductivity = "Eq(kappa, (1410.0*(yO2+yN2+yNO)*0.1*exp(-11.2202)*T**(0.021823*log(T)+0.34357)+2210.0*(yO+yN)*0.1*exp(-11.7344)*T**(0.022652*log(T)+0.342509))*(1.0-exp(-0.010568*T)) )" # NDS simplification of Blottner/Sutherland/Eucken (see Viscosity_model3.m and Viscosity_model_v2_optimise.m)
conductivity_vib = "Eq(kappav, (286.7*(yO2+yN2+yNO)*0.1*exp(-11.2202)*T**(0.021823*log(T)+0.34357)+519.6*(yO+yN)*0.1*exp(-11.7344)*T**(0.022652*log(T)+0.342509))*(1.0-exp(-0.010568*T)) )" # NDS simplification of Blottner/Sutherland/Eucken (see Viscosity_model3.m and Viscosity_model_v2_optimise.m) vibration is just taken as one factor of Rhat/Mhat(needs a closer look maybe)

# chemistry
molefractionO = "Eq(yO, rhoO/(MO*ysum))"
molefractionO2 = "Eq(yO2, rhoO2/(MO2*ysum))"
molefractionN = "Eq(yN, rhoN/(MN*ysum))"
molefractionN2 = "Eq(yN2, rhoN2/(MN2*ysum))"
molefractionNO = "Eq(yNO, rhoNO/(MNO*ysum))"


# Instantiate EinsteinEquation class for expanding the Einstein indices in the equations
eq = EinsteinEquation()

# Expand the simulation equations
simulation_eq = SimulationEquations()
base_eqns = [mass_O, mass_O2, mass_N, mass_N2, mass_NO, momentum, evib, energy, scalar]
for i, base in enumerate(base_eqns):
	base_eqns[i]=eq.expand(base, ndim, coordinate_symbol, substitutions, constants)

for eqn in base_eqns:
	simulation_eq.add_equations(eqn)

# Expand the constituent relations
constituent = ConstituentRelations()  
constituent_eqns = [velocity, pressure, temperature, mixturefraction, viscosity, conductivity, conductivity_vib, molefractionO, molefractionO2, molefractionN, molefractionN2, molefractionNO, evequilO2, evequilN2, evequilNO, evvO2, evvN2, evvNO, timefactorO2, timefactorN2, timefactorNO, tempv] # Added energy am
for i, CR in enumerate(constituent_eqns):
	constituent_eqns[i] = eq.expand(CR, ndim, coordinate_symbol, substitutions, constants)

for eqn in constituent_eqns:
	constituent.add_equations(eqn)

# Write the expanded equations to a Latex file with a given name and titile
latex = LatexWriter()
latex.open('equations.tex', "Einstein Expansion of the simulation equations")
latex.write_string("Simulation equations\n")
for index, eq in enumerate(flatten(simulation_eq.equations)):
    latex.write_expression(eq)

latex.write_string("Constituent relations\n")
for index, eq in enumerate(flatten(constituent.equations)):
    latex.write_expression(eq)

latex.close()

block = SimulationBlock(ndim, block_number=0)

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}


# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# grid generation and initial conditions																															
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

dx, dy = block.deltas
x, y = symbols('x0:%d' % ndim, **{'cls': DataObject})
i, j = block.grid_indexes
nx, ny = symbols('block0np0 block0np1', **{'cls': ConstantObject}) 
Lx1, by = symbols('Lx1 by', **{'cls': ConstantObject}) 
Lx=nx*dx
q_vector=flatten(simulation_eq.time_advance_arrays)
grid_equations= []
x0=i*dx
x1=Lx1*sinh(by*j*dy/Lx1)/sinh(by)

# Added in magic to make it work - Rhys Nov 2022
initial_equations = []
uref, pref, rhoref,rhoEref, Rhat, MO, MO2, MN, MN2, MNO, dhO, dhN, dhNO, thetavO2, thetavN2, thetavNO, Tref = symbols('uref pref rhoref rhoEref Rhat MO MO2 MN MN2 MNO dhO dhN dhNO thetavO2 thetavN2 thetavNO Tref', **{'cls': ConstantObject})
rhoO, rhoO2, rhoN, rhoN2, rhoNO, u, v, p, T ,f, ev, evequilO2, evequilN2, evequilNO,evO2, evN2, evNO, Tv = symbols('rhoO, rhoO2, rhoN, rhoN2, rhoNO, u, v, p, T, f, ev, evequilO2, evequilN2, evequilNO,evO2, evN2, evNO, Tv', **{'cls': GridVariable})

# Mach 2.0 non-catalytic wall case
Re, xMach, Tinf, Twall, Sc = 950.0, 2.0, 288.0, 288.0*1.71138101, 1.0
cN2, cN, cO2, cO, cNO = 0.95, 0.05, 0.0, 0.0, 0.0
adiabatic_condition, catalytic_condition = True, False
MN2, MN, MO2, MO, MNO = 28.0, 14.0, 32.0, 16.0, 30.0 # molar mass
sigOtoNe = (2.0*cO2/MO2+cO/MO+cNO/MNO)/(2.0*cN2/MN2+cN/MN+cNO/MNO)

# pref, rhoref, uref, blthickness = 100.0, 0.0011136546984423536, 713.0769048663936, 0.021092341667534868
pref, rhoref, uref, blthickness = 1000.0, ConstantObject('rhoref'), ConstantObject('uref'), ConstantObject('blthicknesss')

delta0block0m, delta1block0m = blthickness*400, blthickness*100
delta0block0, delta1block0 = str(delta0block0m) + '/(block0np0-1)', str(delta1block0m) + '/(block0np1-1)'


# # Mach 6.0 non-catalytic wall case -- N2 ----------------------------------------------------------------------------------------------
# Re, xMach, Tinf, Twall, Sc = 950.0, 6.0, 288.0, 288.0*6.868697232381297, 1.0
# cN2, cN, cO2, cO, cNO = 1.0, 0.00, 0.0, 0.0, 0.0
# adiabatic_condition, catalytic_condition = True, True
# sigOtoNe = (2.0*cO2/MO2+cO/MO+cNO/MNO)/(2.0*cN2/MN2+cN/MN+cNO/MNO)
# MN2, MN, MO2, MO, MNO = 28.0, 14.0, 32.0, 16.0, 30.0 # molar mass
# pref, rhoref, uref, blthickness = 100.0, 0.0011693374333644712, 2075.9590553653807, 0.007037143221907268

# delta0block0m, delta1block0m = blthickness*400, blthickness*100
# delta0block0, delta1block0 = str(delta0block0m) + '/(block0np0-1)', str(delta1block0m) + '/(block0np1-1)'


## Ensure the grid size passed to the initialisation routine matches the grid sizes used in the simulation parameters
polynomial_directions = [(False, DataObject('x0')), (True, DataObject('x1'))]
n_poly_coefficients = 100
grid_const = ["Lx1", "by"]
for con in grid_const:
    local_dict[con] = ConstantObject(con)
gridx0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
gridx1 = parse_expr("Eq(DataObject(x1), Lx1*sinh(by*block.deltas[1]*block.grid_indexes[1]/Lx1)/sinh(by))", local_dict=local_dict)
coordinate_evaluation = [gridx0, gridx1]

# initial = Initialise_Flatmix(polynomial_directions, n_poly_coefficients,  Re, xMach, Tinf, rhoref, uref, ydomain, blthickness, coordinate_evaluation)
initial = Initialise_Flatmix(polynomial_directions, n_poly_coefficients,  Re, xMach, Tinf, Twall, Sc, adiabatic_condition, catalytic_condition, cN2, cN, cO2, cO, cNO, pref, rhoref, uref, blthickness, coordinate_evaluation)

# metrics
metriceq =  MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False)], 2)
simulation_eq.apply_metrics(metriceq)


# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# create a schemes dictionary to be used for discretisation																													
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

schemes = {}
# low storage
fns = 'u0 u1 T'
cent = StoreSome(4,fns)
# cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# boundary conditions																													
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

boundaries = []

boundaries += [InletTransferBC(0, 0, scheme=Carpenter())]
boundaries += [ExtrapolationBC(0, 1, order=0, scheme=Carpenter())]
boundaries += [nonCatalyticIsothermalWallBC(1,0, scheme=Carpenter())]
boundaries += [ZeroGradientOutletBC(1, 1)]

block.set_block_boundaries(boundaries)

# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# read/write definitions and output arrays																												
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=1000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays + [x, y])
h5.add_arrays([DataObject('D11')])
block.setio(h5)

# set monitor points
arrays = ['u1','p']
arrays = [block.location_dataset('%s' % dset) for dset in arrays]
indices = [(0, '(block0np1-1)/2'), ('block0np0/4', '(block0np1-1)/2')]
SM = SimulationMonitor(arrays, indices, block, print_frequency=20,fp_precision=12, output_file='output.log')

# set the equations to be solved on the block
block.set_equations([constituent, simulation_eq, initial, metriceq])
# set the discretisation schemes
block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
block.discretise()

# create an algorithm from the discretised computations
alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
# OPSC(alg)
OPSC(alg) # ,OPS_V2=True


# --------------------------------------------------------------------------------------------------------------------------------------------
#																																			
# define constants, in line with intialisation above																												
#																																			
# --------------------------------------------------------------------------------------------------------------------------------------------

# mixlayer conditions
# physical_constants = ['Re', 'Sc', 'uref', 'pref', 'rhoref', 'pexp', 'Twall', 'Twn']
# physical_values = ['1.0', str(Sc), str(uref), str(pref), str(rhoref), '0.0', str(Twall), '0.0']

physical_constants = ['Re', 'Sc', 'pexp',  'Twn']
physical_values = ['1.0', str(Sc), '0.0', '0.0']


substitute_simulation_parameters(physical_constants, physical_values)

gas_data = ['Rhat', 'MO', 'MO2', 'MN', 'MN2', 'MNO', 'dhO', 'dhN', 'dhNO', 'thetavO2', 'thetavN2', 'thetavNO']
physical_values = ['8314.3', '16.0', '32.0', '14.0', '28.0', '30.0', '59.544', '112.951', '21.6009', '2270.0', '3390.0', '2740.0']
substitute_simulation_parameters(gas_data, physical_values)

# concentration_data = ['cOe', 'cNe', 'cO2e', 'cN2e', 'cNOe', 'sigOtoNe']
concentration_data = ['cN2e', 'cNe', 'cO2e', 'cOe', 'cNOe', 'sigOtoNe']
physical_values    = [  str(cN2), str(cN), str(cO2), str(cO), str(cNO), str(sigOtoNe)]
substitute_simulation_parameters(concentration_data, physical_values)

numerical_constants= ['dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'Lx1', 'by']
numerical_values=['1e-9', '4000000', '500', '250', delta0block0, delta1block0, str(delta1block0m), '5.0']

substitute_simulation_parameters(numerical_constants, numerical_values)

reaction_constants = ['Cf1','nf1','thetaf1','B11','B12','B13','B14','B15']
reaction_values = ['2.900e23','-2.0','5.975e4','2.855','0.988','-6.181','-0.023','-0.001']
substitute_simulation_parameters(reaction_constants, reaction_values)
reaction_constants = ['Cf2','nf2','thetaf2','B21','B22','B23','B24','B25']
reaction_values = ['9.680e22','-2.0','5.975e4','2.855','0.988','-6.181','-0.023','-0.001']
substitute_simulation_parameters(reaction_constants, reaction_values)
reaction_constants = ['Cf3','nf3','thetaf3','B31','B32','B33','B34','B35']
reaction_values = ['1.600e22','-1.6','1.132e5','1.858','-1.325','-9.856','-0.174','0.008']
substitute_simulation_parameters(reaction_constants, reaction_values)
reaction_constants = ['Cf4','nf4','thetaf4','B41','B42','B43','B44','B45']
reaction_values = ['4.980e22','-1.6','1.132e5','1.858','-1.325','-9.856','-0.174','0.008']
substitute_simulation_parameters(reaction_constants, reaction_values)
reaction_constants = ['Cf5','nf5','thetaf5','B51','B52','B53','B54','B55']
reaction_values = ['3.700e21','-1.6','1.132e5','1.858','-1.325','-9.856','-0.174','0.008']
substitute_simulation_parameters(reaction_constants, reaction_values)
reaction_constants = ['Cf6','nf6','thetaf6','B61','B62','B63','B64','B65']
reaction_values = ['4.980e21','-1.6','1.132e5','1.858','-1.325','-9.856','-0.174','0.008']
substitute_simulation_parameters(reaction_constants, reaction_values)
reaction_constants = ['Cf7','nf7','thetaf7','B71','B72','B73','B74','B75']
reaction_values = ['7.950e23','-2.0','7.550e4','0.792','-0.492','-6.761','-0.091','0.004']
substitute_simulation_parameters(reaction_constants, reaction_values)
reaction_constants = ['Cf8','nf8','thetaf8','B81','B82','B83','B84','B85']
reaction_values = ['8.370e12','0.0','1.945e4','-2.063','-1.480','-0.580','-0.114','0.005']
substitute_simulation_parameters(reaction_constants, reaction_values)
reaction_constants = ['Cf9','nf9','thetaf9','B91','B92','B93','B94','B95']
reaction_values = ['6.440e17','-1.0','3.837e4','1.066','-0.833','-3.095','-0.084','0.004']
substitute_simulation_parameters(reaction_constants, reaction_values)

print_iteration_ops(every=1000,NaN_check='rhoN2_B0')
