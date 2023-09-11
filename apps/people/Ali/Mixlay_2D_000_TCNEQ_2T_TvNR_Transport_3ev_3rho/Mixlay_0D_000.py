#!/usr/bin/env python
# ------------------------------------------------------------------------------------
from opensbli import *
from sympy import sin, cos, sinh, tanh, exp, pi, log, Piecewise, Or
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
import chemistry
# ------------------------------------------------------------------------------------

#########################################################################################################################
#																														#
# Inputs			 																									#
#																														#
#########################################################################################################################
Chemistry_model = 'none'									# Choose from 'GGS' or 'Park01' or 'none'
# viscosity: 											mu_Yos_Gupta - mu_Yos_Brokaw_Gupta - Sandham - Blottner - Sutherland - none
# thermal conductivity: 								mu_Yos_Gupta - mu_Yos_Brokaw_Gupta - Sandham - Blottner - none
Transport_model = ['Sandham','Sandham','Y'] 	# Viscosity - Thermal conductivity - 3xTC (Y-N)
Diffusion		= 'CSc'								# Options: 'CSc' - 'Lee85'
TvNR_py			= 'N'									# Enable or disable TvNR exchange - Y (yes) or N (NO)
inf_to_zero		= 'Y'									# Enable or disable the replacement of inf to zero for 3-rho ev models

Physical = {
    'Re'			:'1.0',
	'Sc'			:'0.71',
	'uref'			:'623.005172',
	'pref'			:'37691.56',
	'rhoref'		:'0.02',
	'pexp'			:'0.0',
	'vorthick'		:'0.0047229'
}
Grid = {
	'dt'			:'0.0000001',
	'niter'			:'4000',
	'block0np0'		:'120',
	'block0np1'		:'375',
	'Delta0block0'	:'15.0 * 0.0047229 /block0np0',
	'Delta1block0'	:'1.0/(block0np1-1)',
	'Ly'			:'1000.0/15.0 * 0.0047229 ',
	'stretch'		:'3.0',
}
Gas = {
	'Rhat'			:'8314.3',
	'MO'			:'16.0',
	'MO2'			:'32.0',
	'MN'			:'14.0',
	'MN2'			:'28.0',
	'MNO'			:'30.0',
	'dhO'			:'59.544',
	'dhN'			:'112.951',
	'dhNO'			:'21.6009',
	'thetavO2'		:'2270.0',
	'thetavN2'		:'3390.0',
	'thetavNO'		:'2740.0',
	'boltz'			:'1.380649e-23',
}
set_constants = {**Physical , **Grid , **Gas}					# REF[1]: Unpacks the dictionaries in a single one

#########################################################################################################################
#																														#
# Simulation Equations																									#
#																														#
#########################################################################################################################
# Number of simulation dimensions
ndim  = 2
stats = False

# Define the compresible Navier-Stokes equations in Einstein notation
mass_O   = "Eq(Der(rhoO,t),  	- Skew(rhoO*u_j,x_j)  + Der(rho *DO* Der(XO,x_j),x_j) 		+ wdotO	  - rhoO  * (Der(DO* Der(XO,x_j),x_j)+Der(DO2* Der(XO2,x_j),x_j)+Der(DN* Der(XN,x_j),x_j)+Der(DN2* Der(XN2,x_j),x_j)+Der(DNO* Der(XNO,x_j),x_j) ))"
mass_O2  = "Eq(Der(rhoO2,t), 	- Skew(rhoO2*u_j,x_j) + Der(rho *DO2*Der(XO2,x_j),x_j)  	+ wdotO2  - rhoO2 * (Der(DO* Der(XO,x_j),x_j)+Der(DO2* Der(XO2,x_j),x_j)+Der(DN* Der(XN,x_j),x_j)+Der(DN2* Der(XN2,x_j),x_j)+Der(DNO* Der(XNO,x_j),x_j) ))"
mass_N   = "Eq(Der(rhoN,t),  	- Skew(rhoN*u_j,x_j)  + Der(rho *DN* Der(XN,x_j),x_j)  		+ wdotN   - rhoN  * (Der(DO* Der(XO,x_j),x_j)+Der(DO2* Der(XO2,x_j),x_j)+Der(DN* Der(XN,x_j),x_j)+Der(DN2* Der(XN2,x_j),x_j)+Der(DNO* Der(XNO,x_j),x_j) ))"
mass_N2  = "Eq(Der(rhoN2,t), 	- Skew(rhoN2*u_j,x_j) + Der(rho *DN2*Der(XN2,x_j),x_j)  	+ wdotN2  - rhoN2 * (Der(DO* Der(XO,x_j),x_j)+Der(DO2* Der(XO2,x_j),x_j)+Der(DN* Der(XN,x_j),x_j)+Der(DN2* Der(XN2,x_j),x_j)+Der(DNO* Der(XNO,x_j),x_j) ))"
mass_NO  = "Eq(Der(rhoNO,t), 	- Skew(rhoNO*u_j,x_j) + Der(rho *DNO*Der(XNO,x_j),x_j)  	+ wdotNO  - rhoNO * (Der(DO* Der(XO,x_j),x_j)+Der(DO2* Der(XO2,x_j),x_j)+Der(DN* Der(XN,x_j),x_j)+Der(DN2* Der(XN2,x_j),x_j)+Der(DNO* Der(XNO,x_j),x_j) ))"
momentum = "Eq(Der(rhou_i,t), 	- Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j))"
evibO2   = "Eq(Der(rhoevO2,t),  - Skew(rhoevO2*u_j,x_j) + (rhoO2*eveqO2 - rhoevO2)/  (101325.0/(p*ptauO2)) + Der(qvO2_j,x_j)  + Der(DO2*(evO2*Der(XO2,x_j)),x_j) + wdotO2*evO2  )"
evibN2   = "Eq(Der(rhoevN2,t),  - Skew(rhoevN2*u_j,x_j) + (rhoN2*eveqN2 - rhoevN2)/  (101325.0/(p*ptauN2)) + Der(qvN2_j,x_j)  + Der(DN2*(evN2*Der(XN2,x_j)),x_j) + wdotN2*evN2  )"
evibNO   = "Eq(Der(rhoevNO,t),  - Skew(rhoevNO*u_j,x_j) + (rhoNO*eveqNO - rhoevNO)/  (101325.0/(p*ptauNO)) + Der(qvNO_j,x_j)  + Der(DNO*(evNO*Der(XNO,x_j)),x_j) + wdotNO*evNO  )"
energy   = "Eq(Der(rhoE,t), 	- Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) " \
		   "+Der(Rhat*T*(5.0/(2.0*MO)*DO*Der(XO,x_j)+5.0/(2.0*MN)*DN*Der(XN,x_j)+7.0/(2.0*MO2)*DO2*Der(XO2,x_j)+7.0/(2.0*MN2)*DN2*Der(XN2,x_j)" \
		   "+7.0/(2.0*MNO)*DNO*Der(XNO,x_j)),x_j) + Der(4.1868e6*(dhO/MO*DO*Der(XO,x_j)+dhN/MN*DN*Der(XN,x_j)+dhNO/MNO*DNO*Der(XNO,x_j)),x_j) " \
		   "+ Der((evO2*DO2*Der(XO2,x_j)+evN2*DN2*Der(XN2,x_j)+evNO*DNO*Der(XNO,x_j)),x_j) + Der(q_j,x_j) + Der(qvO2_j,x_j) + Der(qvN2_j,x_j) + Der(qvNO_j,x_j) + Der(u_i*tau_i_j ,x_j))"
scalar   = "Eq(Der(rhof,t), 	- Skew(rhof*u_j,x_j) + Der(mu/(Re*Sc)*Der(f,x_j),x_j))"
dummy 	 = "Eq(Der(rhodum,t), 	- Der(rhodum*u_j,x_j) + dum*1.0e-12 +Tvref*1.0e-12)"

# ----------------------------------------------------------------------------------------------------------------------
# Constituent Equations																									#
# ----------------------------------------------------------------------------------------------------------------------

# Variable relations used in the system
velocity 		= "Eq(u_i, rhou_i/rho)"
mixturefraction = "Eq(f, rhof/rho)"
pressure 		= "Eq(p, Rhat*T*(rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO))"
temperature 	= "Eq(T, (rhoE - (rhoevO2 + rhoevN2 + rhoevNO) - dhf - rho*(1./2.)*(KD(_i,_j)*u_i*u_j)) / " \
				 		"(Rhat*(3.0/2.0*(rhoO/MO+rhoN/MN)+5.0/2.0*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))) )"

tempvO2 	 	= "Eq(TvO2,  thetavO2/(log(1.0+thetavO2*Rhat/(MO2*rhoevO2/rhoO2))))"
tempvN2 	 	= "Eq(TvN2,  thetavN2/(log(1.0+thetavN2*Rhat/(MN2*rhoevN2/rhoN2))))"
tempvNO 	 	= "Eq(TvNO,  thetavNO/(log(1.0+thetavNO*Rhat/(MNO*rhoevNO/rhoNO))))"

# Condition for TvNR
if TvNR_py == 'Y':
	tempv 		 	= 	"Eq(Tvref, 1.0 	)"
	tempvNR		 	=  	"Eq(Tv, Tvref + " \
							"((rhoevO2+rhoevN2+rhoevNO) - (  (rhoO2*thetavO2*(Rhat/MO2)/(exp(thetavO2/Tvref)-1)) +"\
															"(rhoN2*thetavN2*(Rhat/MN2)/(exp(thetavN2/Tvref)-1)) +"\
															"(rhoNO*thetavNO*(Rhat/MNO)/(exp(thetavNO/Tvref)-1))  ) ) /" \
							"(  (rhoO2*thetavO2**(2.0)*(Rhat/MO2)*exp(thetavO2/Tvref)/(Tvref**(2.0)*(exp(thetavO2/Tvref)-1)**(2.0)))  "\
							" + (rhoN2*thetavN2**(2.0)*(Rhat/MN2)*exp(thetavN2/Tvref)/(Tvref**(2.0)*(exp(thetavN2/Tvref)-1)**(2.0)))  "\
							" + (rhoNO*thetavNO**(2.0)*(Rhat/MNO)*exp(thetavNO/Tvref)/(Tvref**(2.0)*(exp(thetavNO/Tvref)-1)**(2.0))) )   )"
else:
	tempv		 	=   "Eq(Tv, thetavnum/(ysumM*log(1.0+thetavnum*Rhat/((rhoevO2 + rhoevN2 + rhoevNO)))))"
	tempvNR		 	=   "Eq(Tvref, 1.0)"


evequilO2 	 	= "Eq(eveqO2, 	thetavO2*Rhat/(MO2*(exp(thetavO2/T)-1.0)))"
evequilN2 	 	= "Eq(eveqN2, 	thetavN2*Rhat/(MN2*(exp(thetavN2/T)-1.0)))"
evequilNO 	 	= "Eq(eveqNO, 	thetavNO*Rhat/(MNO*(exp(thetavNO/T)-1.0)))"
evibrationsO2 	= "Eq(evO2, 	rhoevO2/rhoO2)"
evibrationsN2 	= "Eq(evN2, 	rhoevN2/rhoN2)"
evibrationsNO 	= "Eq(evNO, 	rhoevNO/rhoNO)"
evdum 		 	= "Eq(rhoev, 	(rhoevO2 + rhoevN2 + rhoevNO) )"
timefactorO2 	= "Eq(ptauO2, 	((rhoO/MO)/exp(129.0*(T**(-1.0/3.0)-0.0271)-18.42)" \
				   "+ (rhoO2/MO2)/exp(129.0*(T**(-1.0/3.0)-0.0300)-18.42)" \
				   "+ (rhoN/MN)/exp(129.0*(T**(-1.0/3.0)-0.0265)-18.42)" \
				   "+ (rhoN2/MN2)/exp(129.0*(T**(-1.0/3.0)-0.0295)-18.42)" \
				   "+ (rhoNO/MNO)/exp(129.0*(T**(-1.0/3.0)-0.0298)-18.42))/ysum )"
timefactorN2 	= "Eq(ptauN2, 	( (rhoO/MO)/exp(220.0*(T**(-1.0/3.0)-0.0268)-18.42)" \
				   "+ (rhoO2/MO2)/exp(220.0*(T**(-1.0/3.0)-0.0295)-18.42)" \
				   "+ (rhoN/MN)/exp(220.0*(T**(-1.0/3.0)-0.0262)-18.42)" \
				   "+ (rhoN2/MN2)/exp(220.0*(T**(-1.0/3.0)-0.0290)-18.42)" \
				   "+ (rhoNO/MNO)/exp(220.0*(T**(-1.0/3.0)-0.0293)-18.42))/ysum )"
timefactorNO 	= "Eq(ptauNO, 	( (rhoO/MO)/exp(168.0*(T**(-1.0/3.0)-0.0270)-18.42)" \
				   "+ (rhoO2/MO2)/exp(168.0*(T**(-1.0/3.0)-0.0298)-18.42)" \
				   "+ (rhoN/MN)/exp(168.0*(T**(-1.0/3.0)-0.0264)-18.42)" \
				   "+ (rhoN2/MN2)/exp(168.0*(T**(-1.0/3.0)-0.0293)-18.42)" \
				   "+ (rhoNO/MNO)/exp(168.0*(T**(-1.0/3.0)-0.0295)-18.42))/ysum )"
TtauparkO2   	= "Eq(tauparkO2, 1/( (rhoO2*Rhat/(MO2*boltz)) * (8.0*boltz*T/(3.14159265358979323846*(MO2*boltz/Rhat)))**0.5 * 3.0*1.0e-21*(50000.0/T)**2.0 )  )"
TtauparkN2   	= "Eq(tauparkN2, 1/( (rhoN2*Rhat/(MN2*boltz)) * (8.0*boltz*T/(3.14159265358979323846*(MN2*boltz/Rhat)))**0.5 * 3.0*1.0e-21*(50000.0/T)**2.0 )  )"
TtauparkNO   	= "Eq(tauparkNO, 1/( (rhoNO*Rhat/(MNO*boltz)) * (8.0*boltz*T/(3.14159265358979323846*(MNO*boltz/Rhat)))**0.5 * 3.0*1.0e-21*(50000.0/T)**2.0 )  )"

# fluid properties
viscosity 		 = "Eq(mu,  	 (XO+XO2+XN+XN2+XNO)*T*Tv*1.0e-6  )" 		# Placeholder for viscosity
conductivity 	 = "Eq(kappatr,  (XO+XO2+XN+XN2+XNO)*T*Tv*1.0e-6  )" 		# Placeholder for thermal conductivity tr
conductivity_vib = "Eq(kappavib, (XO+XO2+XN+XN2+XNO)*T*Tv*1.0e-6  )"  		# Placeholder for thermal conductivity vib

# Diffusion - Depending on the
if Diffusion == 'Lee85':
	DiffO  = "Eq(DO,  	 (XO+XO2+XN+XN2+XNO)*T*rhoref*rhoO*rhoO2*rhoN*rhoN2*rhoNO*1.0e-6  )" 		# Placeholder
	DiffO2 = "Eq(DO2,  	 (XO+XO2+XN+XN2+XNO)*T*rhoref*rhoO*rhoO2*rhoN*rhoN2*rhoNO*1.0e-6  )" 		# Placeholder
	DiffN  = "Eq(DN,  	 (XO+XO2+XN+XN2+XNO)*T*rhoref*rhoO*rhoO2*rhoN*rhoN2*rhoNO*1.0e-6  )" 		# Placeholder
	DiffN2 = "Eq(DN2,  	 (XO+XO2+XN+XN2+XNO)*T*rhoref*rhoO*rhoO2*rhoN*rhoN2*rhoNO*1.0e-6  )" 		# Placeholder
	DiffNO = "Eq(DNO,  	 (XO+XO2+XN+XN2+XNO)*T*rhoref*rhoO*rhoO2*rhoN*rhoN2*rhoNO*1.0e-6  )" 		# Placeholder
elif Diffusion == 'CSc':
	DiffO  = "Eq(DO,  	 mu/(rho *Re*Sc)  )"
	DiffO2 = "Eq(DO2,  	 mu/(rho *Re*Sc)  )"
	DiffN  = "Eq(DN,  	 mu/(rho *Re*Sc)  )"
	DiffN2 = "Eq(DN2,  	 mu/(rho *Re*Sc)  )"
	DiffNO = "Eq(DNO,  	 mu/(rho *Re*Sc)  )"


# comment out if not using 3tc
conductivity_vibO2 	= "Eq(kappavibO2, (XO+XO2+XN+XN2+XNO)*T*Tv*1.0e-6  )"  	# Placeholder for thermal conductivity vib
conductivity_vibN2 	= "Eq(kappavibN2, (XO+XO2+XN+XN2+XNO)*T*Tv*1.0e-6  )"  	# Placeholder for thermal conductivity vib
conductivity_vibNO 	= "Eq(kappavibNO, (XO+XO2+XN+XN2+XNO)*T*Tv*1.0e-6  )"  	# Placeholder for thermal conductivity vib
heat_flux_vibO2 	= "Eq(qvO2_j, 		(kappavibO2/Re)*Der(Tv,x_j))"
heat_flux_vibN2 	= "Eq(qvN2_j, 		(kappavibN2/Re)*Der(Tv,x_j))"
heat_flux_vibNO 	= "Eq(qvNO_j, 		(kappavibNO/Re)*Der(Tv,x_j))"


# chemistry
molefractionO  = "Eq(XO,  rhoO/(MO*ysum))"
molefractionO2 = "Eq(XO2, rhoO2/(MO2*ysum))"
molefractionN  = "Eq(XN,  rhoN/(MN*ysum))"
molefractionN2 = "Eq(XN2, rhoN2/(MN2*ysum))"
molefractionNO = "Eq(XNO, rhoNO/(MNO*ysum))"


# ----------------------------------------------------------------------------------------------------------------------
# Substitutions																											#
# ----------------------------------------------------------------------------------------------------------------------

# Substitutions used in the equations
stress_tensor   = "Eq(tau_i_j, 		(mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux 		= "Eq(q_j, 		    (kappatr/Re)*Der(T,x_j))"
heat_flux_vib 	= "Eq(qv_j, 		(kappavib/Re)*Der(Tv,x_j))"
density    		= "Eq(rho, 	 	 	(rhoO+rhoO2+rhoN+rhoN2+rhoNO))"
molesum    		= "Eq(ysum, 	 	(rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO))"
molesumM   		= "Eq(ysumM, 	 	(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))"
hformation 		= "Eq(dhf,       	4.1868e6*(dhO*rhoO/MO+dhN*rhoN/MN+dhNO*rhoNO/MNO))"
timeconst  		= "Eq(tau,       	(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO)*101325.0/(p*(rhoO2/(MO2*ptauO2)+rhoN2/(MN2*ptauN2)+rhoNO/(MNO*ptauNO))))"
thetavset  		= "Eq(thetavnum, 	(thetavO2*rhoO2/MO2+thetavN2*rhoN2/MN2+thetavNO*rhoNO/MNO))"
dummysub 		= "Eq(dum, 		 	TvO2+TvN2+TvNO+evO2+evN2+rhoev)"

# make substitutions
substitutions = [stress_tensor, heat_flux, heat_flux_vib,heat_flux_vibO2,heat_flux_vibN2,heat_flux_vibNO, density, molesum, molesumM, hformation, timeconst, thetavset,dummysub]
if Diffusion == 'CSc':
	substitutions += [DiffO, DiffO2, DiffN, DiffN2, DiffNO,density]
# Chemistry
substitutions += chemistry.substitutions(Chemistry_model)


# ----------------------------------------------------------------------------------------------------------------------
# Expanding the Equations																								#
# ----------------------------------------------------------------------------------------------------------------------

# Instantiate EinsteinEquation class for expanding the Einstein indices in the equations
eq = EinsteinEquation()
# symbol for the coordinate system in the equations
coordinate_symbol = "x"

# Constants that are used
constants = ["Re", "Sc", "uref", "pref", "rhoref", "pexp", "Rhat", "MO", "MO2", "MN", "MN2", "MNO",
			 "dhO", "dhN", "dhNO", "thetavO2", "thetavN2", "thetavNO","boltz","vorthick"]

# Add the constants used for the chemistry model - (their names)
constants += list(chemistry.constants(Chemistry_model))

# Expand the simulation equations
simulation_eq = SimulationEquations()
base_eqns 	  = [mass_O, mass_O2, mass_N, mass_N2, mass_NO, momentum, evibO2, evibN2,evibNO, energy, scalar,dummy]
for i, base in enumerate(base_eqns):
	base_eqns[i] = eq.expand(base, ndim, coordinate_symbol, substitutions, constants)
for eqn in base_eqns:
	simulation_eq.add_equations(eqn)

# Expand the constituent relations
constituent 	 = ConstituentRelations()
constituent_eqns = [velocity, pressure, temperature,tempvO2,tempvN2,tempvNO, mixturefraction, viscosity, conductivity,
					conductivity_vib,conductivity_vibO2,conductivity_vibN2,conductivity_vibNO,
					molefractionO, molefractionO2, molefractionN, molefractionN2, molefractionNO, evequilO2, evequilN2,
					evequilNO,evdum,TtauparkO2,TtauparkN2,TtauparkNO, timefactorO2, timefactorN2, timefactorNO,tempv,
					tempvNR,evibrationsNO,evibrationsO2,evibrationsN2]
if Diffusion == 'Lee85':
	constituent_eqns+= [DiffO, DiffO2, DiffN, DiffN2, DiffNO]
# Chemistry constiuent equations
constituent_eqns+= chemistry.constituent(Chemistry_model) 		# rate of reaction



for i, CR in enumerate(constituent_eqns):
	constituent_eqns[i] = eq.expand(CR, ndim, coordinate_symbol, substitutions, constants)
for eqn in constituent_eqns:
	constituent.add_equations(eqn)

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

#########################################################################################################################
#																														#
# Grid																													#
#																														#
#########################################################################################################################

# Set the discretisation schemes - low storage
schemes = {}
fns  				= 'u0 u1 T Tvref'
cent 				= StoreSome(4,fns)
schemes[cent.name] 	= cent
rk 					= RungeKuttaLS(4)
# rk 					= RungeKutta(3)
schemes[rk.name] 	= rk
block.set_discretisation_schemes(schemes)

# ----------------------------------------------------------------------------------------------------------------------

dx, dy  			= block.deltas
x, y 				= symbols('x0:%d' % ndim, **{'cls': DataObject})
i, j 				= block.grid_indexes
nx, ny, Ly, stretch = symbols('block0np0 block0np1 Ly stretch', **{'cls': ConstantObject})
Lx					= nx*dx
q_vector			= flatten(simulation_eq.time_advance_arrays)
grid_equations 		= []
stretch_eqn 		= 0.5*Ly*sinh(stretch*(j-(ny-1)/2)/((ny-1)/2))/sinh(stretch)
grid_equations 	   += [Eq(x, i*dx), Eq(y,stretch_eqn)]

# Metrics: Stretching or curvature of the grid
metriceq =  MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False)], 2)
simulation_eq.apply_metrics(metriceq)

# Boundary Conditions
boundaries  = []
direction   = 0
boundaries += [PeriodicBC(direction, side=0)]
boundaries += [PeriodicBC(direction, side=1)]
direction   = 1
boundaries += [SymmetryBC(direction, 0)]
boundaries += [SymmetryBC(direction, 1)]
block.set_block_boundaries(boundaries)

# # Binomial Filter
# j 		= block.grid_indexes[1]
# nfilter = 20
# grid_condition = Or(j<=nfilter,j>=ny-nfilter)
# BF 		= BinomialFilter(block, order=2, grid_condition=grid_condition ,sigma=0.20)


#########################################################################################################################
#																														#
# Initial Conditions																									#
#																														#
#########################################################################################################################
uref, pref, rhoref,rhoEref, Rhat, MO, MO2, MN, MN2, MNO, dhO, dhN, dhNO, thetavO2, thetavN2, thetavNO,boltz, vorthick \
	= symbols('uref pref rhoref rhoEref Rhat MO MO2 MN MN2 MNO dhO dhN dhNO thetavO2 thetavN2 thetavNO boltz vorthick', **{'cls': ConstantObject})
rhoO, rhoO2, rhoN, rhoN2, rhoNO,rho, u, v, p, T ,f, ev, evequilO2, evequilN2, evequilNO,evO2, evN2, evNO,evO2eq, evN2eq, evNOeq, Tv \
	= symbols('rhoO, rhoO2, rhoN, rhoN2, rhoNO,rho, u, v, p, T, f, ev, evequilO2, evequilN2, evequilNO,evO2, evN2, evNO,evO2eq, evN2eq, evNOeq, Tv', **{'cls': GridVariable})


initial_equations = [

	# Velcity/Scalar
	Eq(u,uref*tanh(2.0*y/vorthick)),
	Eq(v,0.01*uref*cos(2.0*pi*x/Lx)*exp(-(y/vorthick)**2.0/10.0)),
	Eq(f, 0.5 * (1.0 + tanh(2.0 * y / vorthick))),

	# Species/Density
	Eq(rhoO,rhoref*0.000),
	Eq(rhoO2,rhoref*0.21),
	Eq(rhoN,rhoref*0.000),
	Eq(rhoN2,rhoref*0.79),
	Eq(rhoNO,rhoref*0.000),
	Eq(rho,rhoO+rhoO2+rhoN+rhoN2+rhoNO),

	# Vibrational Energy
	# Eq(Tv,pref/(Rhat*(rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO))),
	Eq(Tv,pref/(Rhat*(rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO)) - 1000.0),
	Eq(evO2, thetavO2*Rhat/(MO2*(exp(thetavO2/Tv)-1.0))),		# Added this
	Eq(evN2, thetavN2*Rhat/(MN2*(exp(thetavN2/Tv)-1.0))),		# Added this
	Eq(evNO, thetavNO*Rhat/(MNO*(exp(thetavNO/Tv)-1.0))),		# Added this

	# q Vectors
	Eq(q_vector[0],rhoO),
	Eq(q_vector[1],rhoO2),
	Eq(q_vector[2],rhoN),
	Eq(q_vector[3],rhoN2),
	Eq(q_vector[4],rhoNO),
	Eq(q_vector[5],rho*u),
	Eq(q_vector[6],rho*v),
	Eq(q_vector[7],rhoO2*evO2),
	Eq(q_vector[8],rhoN2*evN2),
	Eq(q_vector[9],rhoNO*evNO),
	Eq(q_vector[10],pref*(3.0/2.0*(rhoO/MO+rhoN/MN)+5.0/2.0*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))/(rhoO/MO+rhoN/MN+rhoO2/MO2+rhoN2/MN2+rhoNO/MNO)+4.1868e6*(dhO*rhoO/MO+dhN*rhoN/MN+dhNO*rhoNO/MNO)+rhoO2*evO2+rhoN2*evN2+rhoNO*evNO+0.5*rho*(u**2+v**2)),
	Eq(q_vector[11],rho*f),
	Eq(q_vector[12],1.0),
]

# Parse the initial conditions
initial = GridBasedInitialisation()
initial.add_equations(grid_equations + initial_equations)


#########################################################################################################################
#																														#
# Printouts & Monitor Points / Latex																					#
#																														#
#########################################################################################################################

# set the IO class to write out arrays
kwargs  = {'iotype': "Write"}
h5 		= iohdf5(save_every=1000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays + [x, y])
h5.add_arrays([DataObject('T'),DataObject('Tv'),DataObject('p'),DataObject('evO2'),DataObject('evN2'),DataObject('evNO')])
block.setio(copy.deepcopy(h5))

# set monitor points
arrays  = ['rhoO','rhoO2','rhoN','rhoN2','rhoNO','rhoE','rhoev','T','Tv','mu']
# arrays  = ['wdotO','wdotO2','wdotN','wdotN2','wdotNO','rhoE','rhoev','T','Tv']
arrays  = [block.location_dataset('%s' % dset) for dset in arrays]
indices = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)] # ('block0np0/4', '(block0np1-1)/2')
SM 		= SimulationMonitor(arrays, indices, block, print_frequency=2,fp_precision=12, output_file='output.log')

# ----------------------------------------------------------------------------------------------------------------------

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

#########################################################################################################################
#																														#
# Miscellaneous																											#
#																														#
#########################################################################################################################
# set chemical model constants
set_constants.update(chemistry.constants(Chemistry_model))

# set the equations to be solved on the block
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial, metriceq]) # "BF.equation_classes" is for the filter

# Discretise the equations on the block
block.discretise()

# create an algorithm from the discretised computations
alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg) # ,OPS_V2=True

# Add the data from "Input" section
substitute_simulation_parameters(set_constants.keys(), set_constants.values())

print_iteration_ops(NaN_check='rhoN2_B0')


#############################################################################################################################################
#																																			#
# Refrences																																	#
#																																			#
# 	[1]: Unpacks the dictionaries in a single one																							#
#			https://stackoverflow.com/questions/13361510/typeerror-unsupported-operand-types-for-dict-items-and-dict-items					#
#############################################################################################################################################





