#!/usr/bin/env python

# --------------------------------------------------------------------------------------------------------------------------------------------
# mixflat_transition, 
#          - debugging mode: transition case still in development, but this works and provides with turbulent results
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
from opensbli.utilities.oblique_shock import ShockConditions

import time_averaging
# settings to turn on statistics gathering & read from restart file
stats = True
restart = True
stats_class = []
if stats:
    stats_class = time_averaging.get_stats_classes()

# Declare constant values
restart_iteration_no = symbols("restart_iteration_no", **{'cls': ConstantObject})
restart_iteration_no.datatype = Int()
CTD.add_constant([restart_iteration_no])

ndim = 3
stats = True
# Define the compresible Navier-Stokes equations in Einstein notation
mass_O = "Eq(Der(rhoO,t), - Skew(rhoO*u_j,x_j) + Der(mu/(Re*Sc)*Der(yO,x_j),x_j)  )" #+wdotO
mass_O2 = "Eq(Der(rhoO2,t), - Skew(rhoO2*u_j,x_j) + Der(mu/(Re*Sc)*Der(yO2,x_j),x_j)  )" #+wdotO2
mass_N = "Eq(Der(rhoN,t), - Skew(rhoN*u_j,x_j) + Der(mu/(Re*Sc)*Der(yN,x_j),x_j)  )" #+wdotN
mass_N2 = "Eq(Der(rhoN2,t), - Skew(rhoN2*u_j,x_j) + Der(mu/(Re*Sc)*Der(yN2,x_j),x_j)  )" #+wdotN2
mass_NO = "Eq(Der(rhoNO,t), - Skew(rhoNO*u_j,x_j) + Der(mu/(Re*Sc)*Der(yNO,x_j),x_j)  )" #+wdotNO

momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j))"

# evib =   "Eq(Der(rhoev,t), - Skew(rhoev*u_j,x_j) + (rhoO2*eveqO2+rhoN2*eveqN2+rhoNO*eveqNO - rho*ev)/tau +Der(mu/(Re*Sc)*(evO2*Der(yO2,x_j)+evN2*Der(yN2,x_j)+evNO*Der(yNO,x_j)),x_j) )" # - Der(qv_j,x_j) + wdotO2*evO2+wdotN2*evN2+wdotNO*evNO
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) +Der(mu/(Re*Sc)*Rhat*T*(5.0/(2.0*MO)*Der(yO,x_j)+5.0/(2.0*MN)*Der(yN,x_j)+7.0/(2.0*MO2)*Der(yO2,x_j)+7.0/(2.0*MN2)*Der(YN2,x_j)+7.0/(2.0*MNO)*Der(YNO,x_j)),x_j) + Der(mu/(Re*Sc)*4.1868e6*(dhO/MO*Der(YO,x_j)+dhN/MN*Der(YN,x_j)+dhNO/MNO*Der(YNO,x_j)),x_j)  - Der(q_j,x_j) - Der(qv_j,x_j)  + Der(u_i*tau_i_j ,x_j))" # + Der(mu/(Re*Sc)*(evO2*Der(yO2,x_j)+evN2*Der(yN2,x_j)+evNO*Der(yNO,x_j)),x_j) - Der(q_j,x_j) - Der(qv_j,x_j) 
scalar = "Eq(Der(rhof,t), - Skew(rhof*u_j,x_j) + Der(mu/(Re*Sc)*Der(f,x_j),x_j))" # non-reacting scalar is useful as a (diffusing) marker of original fluid regions

# Substitutions used in the equations
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))" 
# note that we need to keep Re (set=1) as this is used in scheme.py to split terms into convective and viscous (Q for DL - does this apply to conduction and species diffusion too?)
heat_flux = "Eq(q_j, -(kappa/Re)*Der(T,x_j))"
heat_flux_vib = "Eq(qv_j, -(kappav/Re)*Der(Tv,x_j))" 
# heat_flux_vib = "Eq(qv_j, -(kappav/Re)*Der(Tv,x_j))" 
# evibration = "Eq(ev, rhoev/rho)"
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
substitutions = [stress_tensor, heat_flux, heat_flux_vib, density, molesum, molesumM, hformation, timeconst, thetavset]



# Constants that are used
constants = ["Re", "Sc", "uref", "pref", "rhoref", "pexp", "Rhat", "MO", "MO2", "MN", "MN2", "MNO","rhoEref"]
constants=constants+["dhO", "dhN", "dhNO", "thetavO2", "thetavN2", "thetavNO"]


# symbol for the coordinate system in the equations
coordinate_symbol = "x"

# Variable relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
mixturefraction = "Eq(f, rhof/rho)"
# energyam = "Eq(E, rhoE/rho)"				# Added to involve rhoE in the equaation
pressure = "Eq(p, Rhat*T*(rhoO/MO+rhoO2/MO2+rhoN2/MN2+rhoN/MN+rhoNO/MNO))" # 
temperature = "Eq(T, (rhoE -rhoev - dhf - rho*(1./2.)*(KD(_i,_j)*u_i*u_j))/(Rhat*(3.0/2.0*(rhoO/MO+rhoN/MN)+5.0/2.0*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))) )" # 
#temperature = "Eq(T, (rhoE - 4.1868e6*(dhO*rhoO/MO+dhN*rhoN/MN+dhNO*rhoNO/MNO)- rho*(1./2.)*(KD(_i,_j)*u_i*u_j))/(Rhat*(3/2*(rhoO/MO+rhoN/MN)+5/2*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))) )"

# vibrational terms
tempv = "Eq(Tv, (thetavnum/(ysumM*log(1.0+thetavnum*Rhat/(1)))))" # method to find Tv based on a mole-weighted thetav (compare with N-R or Cv-based method later)
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
base_eqns = [mass_O, mass_O2, mass_N, mass_N2, mass_NO, momentum, energy, scalar]
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

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)
# Local dictionary for the block
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Metric transformation
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False), (False, False)], 2)
# # Create the Ducros equations for the shock sensor
# SS = ShockSensor()
# shock_sensor, sensor_array = SS.ducros_equations(block, coordinate_symbol, metriceq)
# store_sensor = True
# teno_order = 6
# Avg = RoeAverage([0, 1])
# LLF = LLFTeno(teno_order, formulation='adaptive', averaging=Avg, sensor=sensor_array, store_sensor=store_sensor)
schemes = {}
# schemes[LLF.name] = LLF
fns = 'u0 u1 u2 T'
cent = StoreSome(4, fns)
schemes[cent.name] = cent
rk = RungeKuttaLS(3)
schemes[rk.name] = rk
# metrics
# metriceq =  MetricsEquation()
# metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False)], 2)
simulation_eq.apply_metrics(metriceq)

# # Create a schemes dictionary to be used for discretisation
# schemes = {}
# # low storage
# fns = 'u0 u1 T'
# cent = StoreSome(4,fns)
# #cent = Central(4)
# schemes[cent.name] = cent
# rk = RungeKutta(3)
# schemes[rk.name] = rk

# block.set_discretisation_schemes(schemes)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create the body force after specifying a time stepping scheme
# A, B, xF, yF, bta, omega, dt = symbols('A B xF yF bta omega dt', **{'cls': ConstantObject})
# forcing_const = ["A", "B", "xF", "yF", "bta", "omega", "dt"]

# current_iter = block.get_temporal_schemes[0].iteration_number
# x0, x1, x2 = symbols('x0 x1 x2', **{'cls': DataObject})
# # damping = (1 - exp(-x1 / 0.5))
# body_force = Eq(DataObject('BF'), A*exp(-((x0-xF)**2+(x1-yF)**2))*cos(bta*x2)*sin(omega*dt*current_iter))

# # Add forcing as an acoustic term to the continuity equation
# for i, eq in enumerate(base_eqns):
#     if i == 0:
#         base_eqns[i] = Eq(eq.lhs, eq.rhs + body_force.lhs)

# # Create SimulationEquations and Constituent relations, add the expanded equations
# simulation_eq = SimulationEquations()
# constituent = ConstituentRelations()
# for eqn in base_eqns:
#     simulation_eq.add_equations(eqn)
# for eqn in constituent_eqns:
#     constituent.add_equations(eqn)

# Apply metric transformation to the simulation equations
# simulation_eq.apply_metrics(metriceq)


# Define a shock sensor for the TENO schemee
# Add the shock sensor to constituent relations
# constituent.add_equations(shock_sensor)
# Add the body forcing term to the constituent relations
# constituent.add_equations(body_force)

# ---------------------------------------------------------------------------------------------------------------
# Boundary conditions onwards
boundaries = [[0, 0] for t in range(ndim)]
# Left pressure extrapolation at x= 0, inlet conditions
direction, side = 0, 0
boundaries[direction][side] = InletTransferBC(direction, side, scheme=ReducedAccess())
# Right extrapolation at outlet
direction, side = 0, 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0, scheme=ReducedAccess())
# Bottom no-slip isothermal wall
local_dict['Lx1'] = ConstantObject('Lx1')
local_dict['by'] = ConstantObject('by')
local_dict['Twall'] = ConstantObject('Twall')
direction, side = 1, 0
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
# wall_eqns = [rhoE_wall]
# boundaries[direction][side] = IsothermalWallBC(direction, side, wall_eqns, scheme=ReducedAccess())


Amp, sigma, xts, xtp = symbols('tripA tripSigma xts xtp', **{'cls':ConstantObject})
dt, omega0, omega1, omega2 = symbols('dt omega_0 omega_1 omega_2', **{'cls': ConstantObject})
# Spatial Modes
k0, k1, k2 = symbols('k_0 k_1 k_2', **{'cls': ConstantObject})
phi0, phi1, phi2 = symbols('phi_0 phi_1 phi_2', **{'cls': ConstantObject})



# current_iter = symbols('current_iter', **{'cls': ConstantObject})
current_iter = Globalvariable('iter', integer=True)
t = dt*current_iter
# Coordinate arrays
x0, z0 = DataObject('x0'), DataObject('x2')
conditional_expressions = []
# trip_eqn = Amp*(exp(-(x0 - xts)**2  / (2*sigma**2))*(sin(k0 * z0)*sin(omega0*dt*current_iter + phi0) + sin(k1 * z0)*sin(omega1*dt*current_iter + phi1) + sin(k2 * z0)*sin(omega2*dt*current_iter + phi2)))
# trip_eqn = exp(-(x0 - 20)**2.0/ ( 0.425*2/(0.1*0.00025439)))  * 0.1* 585.1815* sin(0.15 * t) * sin(1.4 * z0) 10000000
# trip_eqn = exp(-((x0/0.00025439 - 20))**2.0/ ( 0.425*2/(0.1)))  * 0.15* 585.1815* sin(0.15 * 10000000  * t) * sin((0.8/0.00025439) * z0)
uref,deltastar = 585.1815, 0.00025439
trip_eqn = exp(-((x0/deltastar - 20))**2.0/ ( 0.425*2/(0.1)))  * 0.2* uref* sin(0.15 * (uref/deltastar)  * t) * sin((0.92*(1.0/deltastar)) * z0)
wall_normal_velocity = OpenSBLIEq(DataObject('rhou1'), (DataObject('rhoO') + DataObject('rhoN') + DataObject('rhoO2') + DataObject('rhoN2') + DataObject('rhoNO'))*trip_eqn)
# DataObject('rho')

local_dict['Lx1'] = ConstantObject('Lx1')
local_dict['by'] = ConstantObject('by')
local_dict['Twall'] = ConstantObject('Twall')
direction, side = 1, 0
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rhoev) + 4.1868e6*(DataObject(rhoO)*dhO/MO + DataObject(rhoN)*dhN/MN + DataObject(rhoNO)*dhNO/MNO ) + Twall*Rhat* (1.5*(DataObject(rhoO)/MO + DataObject(rhoN)/MN) + 2.5*(DataObject(rhoO2)/MO2 +  DataObject(rhoN2)/MN2 + DataObject(rhoNO)/MNO)) )", local_dict=local_dict)


# rhoO_wall = parse_expr("Eq(DataObject(rhoO), 0.0)")
# rhoN_wall = parse_expr("Eq(DataObject(rhoN), 0.0)")
# rhoO2_wall = parse_expr("Eq(DataObject(rhoO2), 0.0)")
# rhoN2_wall = parse_expr("Eq(DataObject(rhoN2), DataObject(rhoO) + DataObject(rhoN) + DataObject(rhoO2) + DataObject(rhoN2) + DataObject(rhoN)) ")
# rhoNO_wall = parse_expr("Eq(DataObject(rhoNO), 0.0)")

rhoO_wall = OpenSBLIEq(DataObject('rhoO'), 0.0)
rhoN_wall = OpenSBLIEq(DataObject('rhoN'), 0.0)
rhoO2_wall = OpenSBLIEq(DataObject('rhoO2'), 0.0)
rhoN2_wall = OpenSBLIEq(DataObject('rhoN2'),  DataObject('rhoO') + DataObject('rhoN') + DataObject('rhoO2') + DataObject('rhoN2') + DataObject('rhoNO'))
rhoNO_wall = OpenSBLIEq(DataObject('rhoNO'), 0.0)


# rhoEv_wall = parse_expr("Eq(DataObject(rhoE), Rhat)")  
Rhat, Twall, thetavO2, thetavN2, thetavNO = symbols('Rhat Twall thetavO2 thetavN2 thetavNO', **{'cls':ConstantObject})
MO, MN, MO2, MN2, MNO = symbols('MO MN MO2 MN2 MNO', **{'cls':ConstantObject})

# rhoEv_wall = parse_expr("Eq(DataObject(rhoev), Rhat*(DataObject(rhoO2)*thetavO2/(MO2*(exp(thetavO2/Twall)-1.0)) +DataObject(rhoN2)*thetavN2/(MN2*(exp(thetavN2/Twall)-1.0)) +DataObject(rhoNO)*thetavNO/(MNO*(exp(thetavNO/Twall)-1.0))))")
# rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rhoev) + 4.1868e6*(DataObject(rhoO)*dhO/MO + DataObject(rhoN)*dhN/MN + DataObject(rhoNO)*dhNO/MNO ) + Twall*Rhat* (1.5*(DataObject(rhoO)/MO + DataObject(rhoN)/MN) + 2.5*(DataObject(rhoO2)/MO2 +  DataObject(rhoN2)/MN2 + DataObject(rhoNO)/MNO)) )")
rhoEv_wall = OpenSBLIEq(DataObject('rhoev'), Rhat*(DataObject('rhoO2')*thetavO2/(MO2*(exp(thetavO2/Twall)-1.0)) +DataObject('rhoN2')*thetavN2/(MN2*(exp(thetavN2/Twall)-1.0)) +DataObject('rhoNO')*thetavNO/(MNO*(exp(thetavNO/Twall)-1.0))))

rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rhoev) + 4.1868e6*(DataObject(rhoO)*dhO/MO + DataObject(rhoN)*dhN/MN + DataObject(rhoNO)*dhNO/MNO ) + Twall*Rhat* (1.5*(DataObject(rhoO)/MO + DataObject(rhoN)/MN) + 2.5*(DataObject(rhoO2)/MO2 +  DataObject(rhoN2)/MN2 + DataObject(rhoNO)/MNO)) )", local_dict=local_dict)

# rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
wall_eqns = [rhoO_wall, rhoN_wall, rhoO2_wall, rhoN2_wall, rhoNO_wall, rhoEv_wall, rhoE_wall]
# wall_eqns = [rhoEv_wall, rhoE_wall]
# wall_energy = [rhoE_wall]
# ForcingStripBC(direction, side, v, wall_eqns)
boundaries[direction][side] = ForcingStripBC(direction, side,  wall_normal_velocity, wall_eqns)
# boundaries[direction][side] = catalyticWallBC(direction, side, scheme=ReducedAccess())


# Top dirichlet shock generator condition
direction, side = 1, 1
wave_angle = 44.6607551
xmach = 2.0
gamma = 1.4
pre_shock = (1.0, 1.0, 0.00466654053208844, (1.0/(gamma*xmach**2))/(gamma-1.0) + 0.5*(1.0*1.0**2 + 0.00466654053208844**2))
OS = ShockConditions(wave_angle, xmach, gamma)
post_shock = OS.conservative_post_shock_conditions(1.0)

rho = parse_expr("Eq(DataObject(rho), Piecewise((%.15f, DataObject(x0)>2000.0), (%.15f, True)))" % (post_shock[0], pre_shock[0]), local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), Piecewise((%.15f, DataObject(x0)>2000.0), (%.15f, True)))" % (post_shock[1], pre_shock[1]), local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), Piecewise((%.15f, DataObject(x0)>2000.0), (%.15f, True)))" % (post_shock[2], pre_shock[2]), local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), Piecewise((%.15f + 0.5*DataObject(rhou2)**2 / DataObject(rho), DataObject(x0)>2000.0), (%.15f + 0.5*DataObject(rhou2)**2 / DataObject(rho), True)))" % (post_shock[3], pre_shock[3]), local_dict=local_dict)
upper_eqns = [rho, rhou0, rhou1, rhoE]

boundaries[direction][side] = ZeroGradientOutletBC(direction, side)
# Periodic direction 2
direction = 2
for side in [0,1]:
    boundaries[direction][side] = PeriodicBC(direction, side)
block.set_block_boundaries(boundaries)

# Perform initial condition
# Re, xMach, Tinf = 950.0, 2.0, 288
# Re, xMach, Tinf, Twall = 1000.0, 2.0, 288.0, 288.0
# rhoref, uref, ydomain, blthickness = 0.11693, 585.1815, 0.00025439*110, 0.00025439
# Ensure the grid size passed to the initialisation routine matches the grid sizes used in the simulation parameters
grid_const = ["Lx1", "by"]
for con in grid_const:
    local_dict[con] = ConstantObject(con)
gridx0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
gridx1 = parse_expr("Eq(DataObject(x1), Lx1*sinh(by*block.deltas[1]*block.grid_indexes[1]/Lx1)/sinh(by))", local_dict=local_dict)
gridx2 = parse_expr("Eq(DataObject(x2), block.deltas[2]*block.grid_indexes[2])", local_dict=local_dict)
coordinate_evaluation = [gridx0, gridx1, gridx2]
polynomial_directions = [(False, DataObject('x0')), (True, DataObject('x1')), (False, DataObject('x2'))]


# Added in magic to make it work - Rhys Nov 2022
initial_equations = []
uref, pref, rhoref,rhoEref, Rhat, MO, MO2, MN, MN2, MNO, dhO, dhN, dhNO, thetavO2, thetavN2, thetavNO, Tref = symbols('uref pref rhoref rhoEref Rhat MO MO2 MN MN2 MNO dhO dhN dhNO thetavO2 thetavN2 thetavNO Tref', **{'cls': ConstantObject})
rhoO, rhoO2, rhoN, rhoN2, rhoNO, u, v, p, T ,f, ev, evequilO2, evequilN2, evequilNO,evO2, evN2, evNO, Tv = symbols('rhoO, rhoO2, rhoN, rhoN2, rhoNO, u, v, p, T, f, ev, evequilO2, evequilN2, evequilNO,evO2, evN2, evNO, Tv', **{'cls': GridVariable})

n_poly_coefficients = 50

Re, xMach, Tinf, Twall = 1000.0, 2.0, 288.0, 288.0
rhoref, uref, ydomain, blthickness = 0.11693, 585.1815, 0.00025439*110, 0.00025439


# initial = Initialise_Katzer(polynomial_directions, n_poly_coefficients, Re, xMach, Tinf, coordinate_evaluations=coordinate_evaluation)
initial = Initialise_Flatmix(polynomial_directions, n_poly_coefficients,  Re, xMach, Tinf, Twall, Sc, adiabatic_condition, catalytic_condition, cN2, cN, cO2, cO, cNO, pref, rhoref, uref, blthickness, coordinate_evaluation)

kwargs = {'iotype': "Write"}
h5 = iohdf5(arrays=simulation_eq.time_advance_arrays, save_every=5000, **kwargs)
# h5 = iohdf5(arrays=simulation_eq.time_advance_arrays, save_every=100, **kwargs)
h5.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('x2'), DataObject('D11')])
# h5.add_arrays(simulation_eq.time_advance_arrays + [x, y])
h5.add_arrays([DataObject('D11')])
# if store_sensor:
#     h5.add_arrays([DataObject('TENO')])
block.setio(h5)

if stats:
    stats_arrays = time_averaging.get_arrays()
    kwargs = {'iotype': "Write", "name": "stats_output.h5"}
    h5_stats = iohdf5(**kwargs)
    h5_stats.add_arrays(stats_arrays)
    block.setio(h5_stats)

# set monitor points
arrays = ['u1','p']
arrays = [block.location_dataset('%s' % dset) for dset in arrays]
indices = [(0, '(block0np1-1)/2'), ('block0np0/4', '(block0np1-1)/2')]
SM = SimulationMonitor(arrays, indices, block, print_frequency=20,fp_precision=12, output_file='output.log')

# Set equations on the block and discretise
block.set_equations([simulation_eq, constituent, initial, metriceq] + stats_class)
block.set_discretisation_schemes(schemes)
block.discretise()
# Create an algorithm and write the OPS C code
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
# Substitute the simulation parameters
# constants = ['block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0', 'SuthT', 'RefT', 'eps', 'Lx1', 'by', 'A', 'bta', 'omega', 'xF', 'yF', 'teno_a1', 'teno_a2', 'epsilon', 'restart_iteration_no']
# values = ['500*2', '200', '100', '375.0*2.0/(block0np0-1)', '140.0/(block0np1-1)', '27.32/(block0np2)', '110.4', '288.0', '1e-30', '140.0', '5.0', '0.0', '0.23', '0.1011', '20.0', '4.0', '9.5', '3.5', '1.0e-16', '0']
# substitute_simulation_parameters(constants, values)
# print_iteration_ops()

# numerical_constants= ['dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'Lx1', 'by']
# numerical_values=['7e-9', '4000000', '500', '250', '0.00025439*300/(block0np0-1)', '0.00025439*110/(block0np1-1)', '0.00025439*110', '5.0']

numerical_constants = ['dt', 'niter', 'block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0', 'eps', 'Lx1', 'by', 'epsilon', 'restart_iteration_no']
numerical_values=['7e-9', '300000', '500*2', '200', '100', '0.00025439*375.0*2.0/(block0np0-1)', '0.00025439*140.0/(block0np1-1)', '0.00025439*27.32/(block0np2)', '1e-30', '0.00025439*140.0', '5.0', '1.0e-16', '0']
# numerical_values=['7e-9', '500', '500*2', '200', '100', '0.00025439*375.0*2.0/(block0np0-1)', '0.00025439*140.0/(block0np1-1)', '0.00025439*27.32/(block0np2)', '1e-30', '0.00025439*140.0', '5.0', '1.0e-16', '0']

# numerical_values=['7e-9', '4000000', '500*2', '200', '100', '0.00025439*375.0*2.0/(block0np0-1)', '0.00025439*140.0/(block0np1-1)', '0.00025439*27.32/(block0np2)', '1e-30', '0.00025439*140.0', '5.0']

# mixlayer conditions
physical_constants = ['Re', 'Sc', 'uref', 'pref', 'rhoref', 'pexp', 'Twall', 'Twn']
physical_values = ['1.0', '0.71', '585.1815', '10000.0', '0.11693', '0.0', '288.0*1.67619431', '0.0']

substitute_simulation_parameters(physical_constants, physical_values)

gas_data = ['Rhat', 'MO', 'MO2', 'MN', 'MN2', 'MNO', 'dhO', 'dhN', 'dhNO', 'thetavO2', 'thetavN2', 'thetavNO']
physical_values = ['8314.3', '16.0', '32.0', '14.0', '28.0', '30.0', '59.544', '112.951', '21.6009', '2270.0', '3390.0', '2740.0']
substitute_simulation_parameters(gas_data, physical_values)

# # mach 1.2 flow
# numerical_constants= ['dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'Lx1', 'by']
# numerical_values=['7e-10', '40000', '250', '125', '0.01764/(block0np0-1)', '0.00441/(block0np1-1)', '0.00441', '5.0']

# # mach 0.8 flow
# numerical_constants= ['dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'Lx1', 'by']
# numerical_values=['7e-10', '40000', '250', '125', '0.02646/(block0np0-1)', '0.00661/(block0np1-1)', '0.00661', '5.0']

concentration_data = ['concO', 'concN', 'concO2', 'concN2', 'concNO']
physical_values    = [  '0.0',   '0.0',    '0.0',    '1.0',    '0.0']
substitute_simulation_parameters(concentration_data, physical_values)

# mach 0.6 flow
# numerical_constants= ['dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'Lx1', 'by']
# numerical_values=['7e-9', '4000000', '500', '250', '0.00025439*300/(block0np0-1)', '0.00025439*110/(block0np1-1)', '0.00025439*110', '5.0']
# numerical_values=['7e-10', '40000', '250', '125', '15.0*0.0133482/block0np0', '1.0/(block0np1-1)', '100.0/15.0*0.133482', '5.0']
# numerical_values=['7e-10', '40000', '250', '125', '0.03528/(block0np0-1)', '0.00882/(block0np1-1)', '0.00882', '5.0']


# numerical_values=['0.00000000000002', '4000', '120', '251', '15.0*0.0133482/block0np0', '1.0/(block0np1-1)', '100.0/15.0*0.133482', '3.4']
#numerical_values=['0.0000002', '6000', '120', '251', '15.0*0.0133482/block0np0', '1.0/(block0np1-1)', '100.0/15.0*0.133482', '3.4']
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

print_iteration_ops(NaN_check='rhoN2')
