#!/usr/bin/env python

# Started from v1 11/1/22 - full Park kinetics
# v4 14/1/22 Added scalar, vib nonequilibrium and chemistry/vibration coupling (slows it down a lot). Still missing some diffusion and molecular depletion terms and improved, mu, kappa, D models
# v5 addition of catalytic boundary conditions

# Import all the functions from opensbli
from opensbli import *
from sympy import sin, cos, sinh, tanh, exp, pi, log
#import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

# Number of dimensions of the system to be solved
ndim = 2
stats = False
# Define the compresible Navier-Stokes equations in Einstein notation
mass_O = "Eq(Der(rhoO,t), - Skew(rhoO*u_j,x_j) + Der(mu/(Re*Sc)*Der(yO,x_j),x_j) +wdotO )"
mass_O2 = "Eq(Der(rhoO2,t), - Skew(rhoO2*u_j,x_j) + Der(mu/(Re*Sc)*Der(yO2,x_j),x_j) +wdotO2 )"
mass_N = "Eq(Der(rhoN,t), - Skew(rhoN*u_j,x_j) + Der(mu/(Re*Sc)*Der(yN,x_j),x_j) +wdotN )"
mass_N2 = "Eq(Der(rhoN2,t), - Skew(rhoN2*u_j,x_j) + Der(mu/(Re*Sc)*Der(yN2,x_j),x_j) +wdotN2 )"
mass_NO = "Eq(Der(rhoNO,t), - Skew(rhoNO*u_j,x_j) + Der(mu/(Re*Sc)*Der(yNO,x_j),x_j) +wdotNO )"
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j))"
evib = "Eq(Der(rhoev,t), - Skew(rhoev*u_j,x_j) + (rhoO2*eveqO2+rhoN2*eveqN2+rhoNO*eveqNO - rho*ev)/tau +Der(mu/(Re*Sc)*(evO2*Der(yO2,x_j)+evN2*Der(yN2,x_j)+evNO*Der(yNO,x_j)),x_j) - Der(qv_j,x_j) + (wdotO2*evO2+wdotN2*evN2+wdotNO*evNO))" 
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) +Der(mu/(Re*Sc)*Rhat*T*(5.0/(2.0*MO)*Der(yO,x_j)+5.0/(2.0*MN)*Der(yN,x_j)+7.0/(2.0*MO2)*Der(yO2,x_j)+7.0/(2.0*MN2)*Der(yN2,x_j)+7.0/(2.0*MNO)*Der(yNO,x_j)),x_j) + Der(mu/(Re*Sc)*4.1868e6*(dhO/MO*Der(yO,x_j)+dhN/MN*Der(yN,x_j)+dhNO/MNO*Der(yNO,x_j)),x_j) + Der(mu/(Re*Sc)*(evO2*Der(yO2,x_j)+evN2*Der(yN2,x_j)+evNO*Der(yNO,x_j)),x_j) - Der(q_j,x_j) - Der(qv_j,x_j) + Der(u_i*tau_i_j ,x_j))"
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
timeconst = "Eq(tau, (rhoO2/MO2+rhoN2/MN2+rhoNO/MNO)*101325.0/(p*(rhoO2/(MO2*ptauO2)+rhoN2/(MN2*ptauN2)+rhoNO/(MNO*ptauNO))))"
thetavset = "Eq(thetavnum, (thetavO2*rhoO2/MO2+thetavN2*rhoN2/MN2+thetavNO*rhoNO/MNO))"
#
Rf1="Eq(Rf1,Rf1part*(0.001*rhoO2/MO2)*(0.001*rhoN/MN+0.001*rhoO/MO))"
Rb1="Eq(Rb1,Rf1part/K1*(0.001*rhoO/MO)**2*(0.001*rhoN/MN+0.001*rhoO/MO))"
Rf2="Eq(Rf2,Rf2part*(0.001*rhoO2/MO2)*(0.001*rhoN2/MN2+0.001*rhoO2/MO2+0.001*rhoNO/MNO))"
Rb2="Eq(Rb2,Rf2part/K2*(0.001*rhoO/MO)**2*(0.001*rhoN2/MN2+0.001*rhoO2/MO2+0.001*rhoNO/MNO))"
Rf3="Eq(Rf3,Rf3part*(0.001*rhoN2/MN2)*(0.001*rhoN/MN))"
Rb3="Eq(Rb3,Rf3part/K3*(0.001*rhoN/MN)**3)"
Rf4="Eq(Rf4,Rf4part*(0.001*rhoN2/MN2)*(0.001*rhoO/MO))"
Rb4="Eq(Rb4,Rf4part/K4*(0.001*rhoN/MN)**2*(0.001*rhoO/MO))"
Rf5="Eq(Rf5,Rf5part*(0.001*rhoN2/MN2)*(0.001*rhoN2/MN2+0.001*rhoO2/MO2))"
Rb5="Eq(Rb5,Rf5part/K5*(0.001*rhoN/MN)**2*(0.001*rhoN2/MN2+0.001*rhoO2/MO2))"
Rf6="Eq(Rf6,Rf6part*(0.001*rhoN2/MN2)*(0.001*rhoNO/MNO))"
Rb6="Eq(Rb6,Rf6part/K6*(0.001*rhoN/MN)**2*(0.001*rhoNO/MNO))"
Rf7="Eq(Rf7,Rf7part*(0.001*rhoNO/MNO)*(0.001*rhoO/MO+0.001*rhoO2/MO2+0.001*rhoN/MN+0.001*rhoN2/MN2+0.001*rhoNO/MNO))"
Rb7="Eq(Rb7,Rf7part/K7*(0.001*rhoN/MN)*(0.001*rhoO/MO)*(0.001*rhoO/MO+0.001*rhoO2/MO2+0.001*rhoN/MN+0.001*rhoN2/MN2+0.001*rhoNO/MNO))"
Rf8="Eq(Rf8,Rf8part*(0.001*rhoNO/MNO)*(0.001*rhoO/MO))"
Rb8="Eq(Rb8,Rf8part/K8*(0.001*rhoO2/MO2)*(0.001*rhoN/MN))"
Rf9="Eq(Rf9,Rf9part*(0.001*rhoN2/MN2)*(0.001*rhoO/MO))"
Rb9="Eq(Rb9,Rf9part/K9*(0.001*rhoNO/MNO)*(0.001*rhoN/MN))"
#
Rf1factor="Eq(Rf1part,1000.0*Cf1*(T**nf1)*exp(-thetaf1/T))"
Rf2factor="Eq(Rf2part,1000.0*Cf2*(T**nf2)*exp(-thetaf2/T))"
Rf3factor="Eq(Rf3part,1000.0*Cf3*(T**nf3)*exp(-thetaf3/T))"
Rf4factor="Eq(Rf4part,1000.0*Cf4*(T**nf4)*exp(-thetaf4/T))"
Rf5factor="Eq(Rf5part,1000.0*Cf5*(T**nf5)*exp(-thetaf5/T))"
Rf6factor="Eq(Rf6part,1000.0*Cf6*(T**nf6)*exp(-thetaf6/T))"
Rf7factor="Eq(Rf7part,1000.0*Cf7*(T**nf7)*exp(-thetaf7/T))"
Rf8factor="Eq(Rf8part,1000.0*Cf8*(T**nf8)*exp(-thetaf8/T))"
Rf9factor="Eq(Rf9part,1000.0*Cf9*(T**nf9)*exp(-thetaf9/T))"
#Rf1factor="Eq(Rf1part,1000.0*Cf1*((Tv**pexp*T**(1.0-pexp))**nf1)*exp(-thetaf1/(Tv**pexp*T**(1.0-pexp))))" # this pexp (Park approach) slows things down considerably (2-3x)
#Rf2factor="Eq(Rf2part,1000.0*Cf2*((Tv**pexp*T**(1.0-pexp))**nf2)*exp(-thetaf2/(Tv**pexp*T**(1.0-pexp))))"
#Rf3factor="Eq(Rf3part,1000.0*Cf3*((Tv**pexp*T**(1.0-pexp))**nf3)*exp(-thetaf3/(Tv**pexp*T**(1.0-pexp))))"
#Rf4factor="Eq(Rf4part,1000.0*Cf4*((Tv**pexp*T**(1.0-pexp))**nf4)*exp(-thetaf4/(Tv**pexp*T**(1.0-pexp))))"
#Rf5factor="Eq(Rf5part,1000.0*Cf5*((Tv**pexp*T**(1.0-pexp))**nf5)*exp(-thetaf5/(Tv**pexp*T**(1.0-pexp))))"
#Rf6factor="Eq(Rf6part,1000.0*Cf6*((Tv**pexp*T**(1.0-pexp))**nf6)*exp(-thetaf6/(Tv**pexp*T**(1.0-pexp))))"
#Rf7factor="Eq(Rf7part,1000.0*Cf7*((Tv**pexp*T**(1.0-pexp))**nf7)*exp(-thetaf7/(Tv**pexp*T**(1.0-pexp))))"
#Rf8factor="Eq(Rf8part,1000.0*Cf8*((Tv**pexp*T**(1.0-pexp))**nf8)*exp(-thetaf8/(Tv**pexp*T**(1.0-pexp))))"
#Rf9factor="Eq(Rf9part,1000.0*Cf9*((Tv**pexp*T**(1.0-pexp))**nf9)*exp(-thetaf9/(Tv**pexp*T**(1.0-pexp))))"
K1factor="Eq(K1,exp(B11+B12*log(Z)+B13*Z+B14*Z**2+B15*Z**3))"
K2factor="Eq(K2,exp(B21+B22*log(Z)+B23*Z+B24*Z**2+B25*Z**3))"
K3factor="Eq(K3,exp(B31+B32*log(Z)+B33*Z+B34*Z**2+B35*Z**3))"
K4factor="Eq(K4,exp(B41+B42*log(Z)+B43*Z+B44*Z**2+B45*Z**3))"
K5factor="Eq(K5,exp(B51+B52*log(Z)+B53*Z+B54*Z**2+B55*Z**3))"
K6factor="Eq(K6,exp(B61+B62*log(Z)+B63*Z+B64*Z**2+B65*Z**3))"
K7factor="Eq(K7,exp(B71+B72*log(Z)+B73*Z+B74*Z**2+B75*Z**3))"
K8factor="Eq(K8,exp(B81+B82*log(Z)+B83*Z+B84*Z**2+B85*Z**3))"
K9factor="Eq(K9,exp(B91+B92*log(Z)+B93*Z+B94*Z**2+B95*Z**3))"
Zfactor = "Eq(Z, 10000.0/T)"
#Zfactor = "Eq(Z, 10000.0/(Tv**pexp*T**(1.0-pexp)))"

# make substitutions
substitutions = [stress_tensor, heat_flux, heat_flux_vib, evibration, density, molesum, molesumM, hformation, timeconst, thetavset]
substitutions=substitutions+[Rf1, Rb1, Rf1factor, K1factor] 
substitutions=substitutions+[Rf2, Rb2, Rf2factor, K2factor] 
substitutions=substitutions+[Rf3, Rb3, Rf3factor, K3factor] 
substitutions=substitutions+[Rf4, Rb4, Rf4factor, K4factor] 
substitutions=substitutions+[Rf5, Rb5, Rf5factor, K5factor] 
substitutions=substitutions+[Rf6, Rb6, Rf6factor, K6factor] 
substitutions=substitutions+[Rf7, Rb7, Rf7factor, K7factor] 
substitutions=substitutions+[Rf8, Rb8, Rf8factor, K8factor] 
substitutions=substitutions+[Rf9, Rb9, Rf9factor, K9factor] 
substitutions=substitutions+[Zfactor]

# Constants that are used
constants = ["Re", "Sc", "uref", "pref", "rhoref", "pexp", "Rhat", "MO", "MO2", "MN", "MN2", "MNO"]
constants=constants+["dhO", "dhN", "dhNO", "thetavO2", "thetavN2", "thetavNO"]
constants=constants+["B11", "B12", "B13", "B14", "B15", "Cf1", "nf1", "thetaf1"]
constants=constants+["B21", "B22", "B23", "B24", "B25", "Cf2", "nf2", "thetaf2"]
constants=constants+["B31", "B32", "B33", "B34", "B35", "Cf3", "nf3", "thetaf3"]
constants=constants+["B41", "B42", "B43", "B44", "B45", "Cf4", "nf4", "thetaf4"]
constants=constants+["B51", "B52", "B53", "B54", "B55", "Cf5", "nf5", "thetaf5"]
constants=constants+["B61", "B62", "B63", "B64", "B65", "Cf6", "nf6", "thetaf6"]
constants=constants+["B71", "B72", "B73", "B74", "B75", "Cf7", "nf7", "thetaf7"]
constants=constants+["B81", "B82", "B83", "B84", "B85", "Cf8", "nf8", "thetaf8"]
constants=constants+["B91", "B92", "B93", "B94", "B95", "Cf9", "nf9", "thetaf9"]

# symbol for the coordinate system in the equations
coordinate_symbol = "x"

# Variable relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
mixturefraction = "Eq(f, rhof/rho)"
pressure = "Eq(p, Rhat*T*(rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO))"
temperature = "Eq(T, (rhoE -rhoev - dhf - rho*(1./2.)*(KD(_i,_j)*u_i*u_j))/(Rhat*(3.0/2.0*(rhoO/MO+rhoN/MN)+5.0/2.0*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))) )"
#temperature = "Eq(T, (rhoE - 4.1868e6*(dhO*rhoO/MO+dhN*rhoN/MN+dhNO*rhoNO/MNO)- rho*(1./2.)*(KD(_i,_j)*u_i*u_j))/(Rhat*(3/2*(rhoO/MO+rhoN/MN)+5/2*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))) )"

# vibrational terms
tempv = "Eq(Tv, thetavnum/(ysumM*log(1.0+thetavnum*Rhat/(rho*ev))))" # method to find Tv based on a mole-weighted thetav (compare with N-R or Cv-based method later)
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
#reactionrateO = "Eq(wdotO,0.0*MO*(2.0*(Rf1-Rb1)+2.0*(Rf2-Rb2)+1.0*(Rf7-Rb7)-1.0*(Rf8-Rb8)-1.0*(Rf9-Rb9)) )"
#reactionrateO2 = "Eq(wdotO2,0.0*MO2*(-1.0*(Rf1-Rb1)-1.0*(Rf2-Rb2)+1.0*(Rf8-Rb8)) )"
#reactionrateN = "Eq(wdotN,0.0*MN*(2.0*(Rf3-Rb3)+2.0*(Rf4-Rb4)+2.0*(Rf5-Rb5)+2.0*(Rf6-Rb6)+1.0*(Rf7-Rb7)+1.0*(Rf8-Rb8)+1.0*(Rf9-Rb9)) )"
#reactionrateN2 = "Eq(wdotN2,0.0*MN2*(-1.0*(Rf3-Rb3)-1.0*(Rf4-Rb4)-1.0*(Rf5-Rb5)-1.0*(Rf6-Rb6)-1.0*(Rf9-Rb9)) )"
#reactionrateNO = "Eq(wdotNO,0.0*MNO*(-1.0*(Rf7-Rb7)-1.0*(Rf8-Rb8)+1.0*(Rf9-Rb9)) )"
reactionrateO = "Eq(wdotO,0.0*MO*(2.0*(Rf1-Rb1)+2.0*(Rf2-Rb2)+1.0*(Rf7-Rb7)-1.0*(Rf8-Rb8)-1.0*(Rf9-Rb9)) )"
reactionrateO2 = "Eq(wdotO2,0.0*MO2*(-1.0*(Rf1-Rb1)-1.0*(Rf2-Rb2)+1.0*(Rf8-Rb8)) )"
reactionrateN = "Eq(wdotN,0.0*MN*(2.0*(Rf3-Rb3)+2.0*(Rf4-Rb4)+2.0*(Rf5-Rb5)+2.0*(Rf6-Rb6)+1.0*(Rf7-Rb7)+1.0*(Rf8-Rb8)+1.0*(Rf9-Rb9)) )"
reactionrateN2 = "Eq(wdotN2,0.0*MN2*(-1.0*(Rf3-Rb3)-1.0*(Rf4-Rb4)-1.0*(Rf5-Rb5)-1.0*(Rf6-Rb6)-1.0*(Rf9-Rb9)) )"
reactionrateNO = "Eq(wdotNO,0.0*MNO*(-1.0*(Rf7-Rb7)-1.0*(Rf8-Rb8)+1.0*(Rf9-Rb9)) )"

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
constituent_eqns = [velocity, pressure, temperature, mixturefraction, viscosity, conductivity, conductivity_vib, molefractionO, molefractionO2, molefractionN, molefractionN2, molefractionNO, reactionrateO, reactionrateO2, reactionrateN, reactionrateN2, reactionrateNO, evequilO2, evequilN2, evequilNO, evvO2, evvN2, evvNO, timefactorO2, timefactorN2, timefactorNO, tempv]
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

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

#
q_vector=flatten(simulation_eq.time_advance_arrays)

# initial conditions
#dx, dy = block.deltas
#x, y = symbols('x0:%d' % ndim, **{'cls': DataObject})
#i, j = block.grid_indexes
#nx, ny, Ly, stretch = symbols('block0np0 block0np1 Ly stretch', **{'cls': ConstantObject}) 
#Lx=nx*dx
#grid_equations= []
#stretch_eqn=0.5*Ly*sinh(stretch*(j-(ny-1)/2)/((ny-1)/2))/sinh(stretch)
#grid_equations += [Eq(x, i*dx), Eq(y,stretch_eqn)]

# grid 
dx, dy = block.deltas
i, j = block.grid_indexes
nx, ny, Lr, stretch = symbols('block0np0 block0np1 Lr stretch', **{'cls': ConstantObject}) 
x0 = Eq(DataObject('x0'),(0.00008+Lr*sinh(stretch*j/(ny-1))/sinh(stretch))*cos(i*dx) , evaluate=False)
x1 = Eq(DataObject('x1'),-(0.00008+Lr*sinh(stretch*j/(ny-1))/sinh(stretch))*sin(i*dx) , evaluate=False)
#x0 = Eq(DataObject('x0'),(0.0025+Lr*sinh(stretch*j/(ny-1))/sinh(stretch))*cos(i*dx) , evaluate=False)
#x1 = Eq(DataObject('x1'),-(0.0025+Lr*sinh(stretch*j/(ny-1))/sinh(stretch))*sin(i*dx) , evaluate=False)
#pprint(x0)
#pprint(x1)

initial_equations = []
uref, pref, rhoref, Twall, Twalld, Rhat, MO, MO2, MN, MN2, MNO, dhO, dhN, dhNO, thetavO2, thetavN2, thetavNO = symbols('uref pref rhoref Twall Twalldummy Rhat MO MO2 MN MN2 MNO dhO dhN dhNO thetavO2 thetavN2 thetavNO', **{'cls': ConstantObject})
rhoO, rhoO2, rhoN, rhoN2, rhoNO, u, v, p, T ,f, ev, evequilO2, evequilN2, evequilNO = symbols('rhoO, rhoO2, rhoN, rhoN2, rhoNO, u, v, p, T, f, ev, evequilO2, evequilN2, evequilNO', **{'cls': GridVariable})

initial_equations += [Eq(u,uref)] 	# so vorticity thickness=1mm
initial_equations += [Eq(f,1.0)] # a useful smooth function (0 in lower stream 1 in upper stream)
initial_equations += [Eq(v,0.0)]
#initial_equations += [Eq(T,Tref+0.5*(1.0-u)**2.0/1005.0)]   # rough Crocco-Busemann to be going on with (need full form for cp)
#initial_equations += [Eq(T,Tref)]

initial_equations += [Eq(rhoO,rhoref*0.22883)]
initial_equations += [Eq(rhoO2,rhoref*0.00025)]
initial_equations += [Eq(rhoN,rhoref*0.26973)]
initial_equations += [Eq(rhoN2,rhoref*0.49065)]
initial_equations += [Eq(rhoNO,rhoref*0.01053)]
#initial_equations += [Eq(rhoO,rhoref*0.0)]
#initial_equations += [Eq(rhoO2,rhoref*0.2347)]
#initial_equations += [Eq(rhoN,rhoref*0.0)]
#initial_equations += [Eq(rhoN2,rhoref*0.7653)]
#initial_equations += [Eq(rhoNO,rhoref*0.0)]
initial_equations += [Eq(T,pref/(Rhat*(rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO)))]
initial_equations += [Eq(evequilO2,thetavO2*Rhat/(MO2*(exp(thetavO2/T)-1.0)))] 
initial_equations += [Eq(evequilN2,thetavN2*Rhat/(MN2*(exp(thetavN2/T)-1.0)))] 
initial_equations += [Eq(evequilNO,thetavNO*Rhat/(MNO*(exp(thetavNO/T)-1.0)))] 
#initial_equations += [Eq(ev,evequilO2+evequilN2+evequilNO)] # equilibrium ev 
#initial_equations += [Eq(ev,0.0)] 
initial_equations += [Eq(q_vector[0],rhoO)]  			
initial_equations += [Eq(q_vector[1],rhoO2)] 		
initial_equations += [Eq(q_vector[2],rhoN)] 						
initial_equations += [Eq(q_vector[3],rhoN2)]	 	
initial_equations += [Eq(q_vector[4],rhoNO)]			
initial_equations += [Eq(q_vector[5],rhoref*u)]
initial_equations += [Eq(q_vector[6],rhoref*v)]
#initial_equations += [Eq(q_vector[7],0.0)]
#initial_equations += [Eq(q_vector[8],rhoref*5.e6)]
initial_equations += [Eq(q_vector[7],rhoO2*evequilO2+rhoN2*evequilN2+rhoNO*evequilNO)]
initial_equations += [Eq(q_vector[8],pref*(3.0/2.0*(rhoO/MO+rhoN/MN)+5.0/2.0*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))/(rhoO/MO+rhoN/MN+rhoO2/MO2+rhoN2/MN2+rhoNO/MNO)+4.1868e6*(dhO*rhoO/MO+dhN*rhoN/MN+dhNO*rhoNO/MNO)+rhoO2*evequilO2+rhoN2*evequilN2+rhoNO*evequilNO+0.5*rhoref*(u**2+v**2))]
initial_equations += [Eq(q_vector[9],rhoref*f + Twall*Twalld)]

# parse the initial conditions
initial = GridBasedInitialisation()
#initial.add_equations(grid_equations + initial_equations)
initial.add_equations([x0, x1] + initial_equations)

# metrics
metriceq =  MetricsEquation()
#metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (False, False)], 2)
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)
simulation_eq.apply_metrics(metriceq)

# Create a schemes dictionary to be used for discretisation
schemes = {}
# low storage
fns = 'u0 u1 T Tv yO yO2 yN yN2 yNO f'
cent = StoreSome(4,fns)
#cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

# Create boundaries, one for each side per dimension, so in total 6 BC's for 3D'
boundaries = []
boundaries += [PeriodicBC(direction=0, side=0)]
boundaries += [PeriodicBC(direction=0, side=1)]
#boundaries += [PeriodicBC(direction=1, side=0)]
#boundaries += [PeriodicBC(direction=1, side=1)]
#direction=1
#boundaries += [SymmetryBC(direction, 0)]
#boundaries += [SymmetryBC(direction, 1)]

# Isothermal wall in x1 direction
#gama, Minf = symbols('gama Minf', **{'cls': ConstantObject})
#Twall = ConstantObject("Twall")
#wall_energy = [Eq(q_vector[8], Twall*(q_vector[0]+q_vector[1]+q_vector[2]+q_vector[3]+q_vector[4]) )]
#direction = 1
#lower_wall_eq = wall_energy[:]
#boundaries += [IsothermalWallBC(direction, 0, lower_wall_eq)]

direction=1
# Extrapolation bc (for wall, must manually modify for non-slip condition)
side=0
# boundaries += [ExtrapolationBC(direction, side, order=0)]
boundaries += [catalyticIsothermalWallBC(direction, side)]
# Far field boundary
side=1
boundaries += [DirichletBC(direction, side, initial_equations)]

block.set_block_boundaries(boundaries)


# set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=10000, **kwargs)
#block.setio(copy.deepcopy(h5))
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1')]) #,  DataObject('rho_filt'), DataObject('rhou0_filt'), DataObject('rhou1_filt'), DataObject('rhoE_filt')])
block.setio([h5])
#kwargs = {'iotype': "Read"}
#h5_read = iohdf5(**kwargs)
#h5_read.add_arrays([DataObject('x0'), DataObject('x1')])
#block.setio([h5, h5_read])

# Add SFD filtering
# SFD = SFD(block, chifilt=0.1, omegafilt=1.0/0.75)
j = block.grid_indexes[1]
grid_condition = j >= 169
F = BinomialFilter(block, order=10, grid_condition=grid_condition)

# set monitor points
arrays = ['rhoO','rhoO2','rhoN','rhoN2','rhoNO','rhoev']
arrays = [block.location_dataset('%s' % dset) for dset in arrays]
indices = [(0,0), (0,0), (0,0), (0,0), (0,0), (0,0)]
SM = SimulationMonitor(arrays, indices, block, print_frequency=200,fp_precision=12, output_file='output.log')

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
OPSC(alg)


physical_constants = ['Re', 'Sc', 'uref', 'pref', 'rhoref', 'pexp', 'Twall', 'Twalldummy']
physical_values = ['1.0', '1.0', '438.75', '500000.0', '0.16234', '0.0', '3000.0', '0.0']
#physical_values = ['1.0', '0.71', '470.3628', '10000.0', '0.0072319', '0.0']
#physical_values = ['1.0', '0.71', '0.0', '1000.0', '0.005', '0.0']
substitute_simulation_parameters(physical_constants, physical_values)

gas_data = ['Rhat', 'MO', 'MO2', 'MN', 'MN2', 'MNO', 'dhO', 'dhN', 'dhNO', 'thetavO2', 'thetavN2', 'thetavNO']
physical_values = ['8314.3', '16.0', '32.0', '14.0', '28.0', '30.0', '59.544', '112.951', '21.6009', '2270.0', '3390.0', '2740.0']
substitute_simulation_parameters(gas_data, physical_values)

numerical_constants= ['dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'Lr', 'stretch']
numerical_values=['0.0000000000125', '1000', '360', '181', '2.0*M_PI/block0np0', '1.0/(block0np1-1)', '0.00392', '6.0']
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
