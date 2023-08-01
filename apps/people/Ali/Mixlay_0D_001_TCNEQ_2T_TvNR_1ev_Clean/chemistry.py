# Coded by: A. Musawi - 2023.06.04
# Last edited: 2023.06.04



def constants(model):
	if model == 'none':
		return []
	elif model == 'GGS' :
		return constants_GGS

	elif model == 'Park01':
		return constants_P01

def substitutions(model):
	if model == 'none':
		return []
	elif model == 'GGS' :
		return [
				ORf0_GGS, ORf1_GGS, ORf2_GGS, ORf3_GGS, ORf4_GGS, ORf5_GGS, ORf6_GGS, ORf7_GGS, ORf8_GGS,  			# forward reaction
				ORb0_GGS, ORb1_GGS, ORb2_GGS, ORb3_GGS, ORb4_GGS, ORb5_GGS, ORb6_GGS, ORb7_GGS, ORb8_GGS,  				# Backward reaction
				ORb00_GGS, ORb01_GGS, ORb02_GGS, ORb03_GGS, ORb04_GGS, ORb05_GGS, ORb06_GGS, ORb07_GGS, ORb08_GGS,  	# Backward reaction coefficient
				OK0_GGS, OK1_GGS, OK2_GGS, OK3_GGS, OK4_GGS, OK5_GGS, OK6_GGS, OK7_GGS, OK8_GGS,  						# equilibrium constant
				ORf00_GGS, ORf01_GGS, ORf02_GGS, ORf03_GGS, ORf04_GGS, ORf05_GGS, ORf06_GGS, ORf07_GGS, ORf08_GGS,  	# forward reaction coefficient
				OZ_GGS,
				twotee]

	elif model == 'Park01':
		return [
				ORf0_P01,ORf1_P01,ORf2_P01,ORf3_P01,ORf4_P01,ORf5_P01,				# forward reaction
				ORb0_P01,ORb1_P01,ORb2_P01,ORb3_P01,ORb4_P01,ORb5_P01,				# Backward reaction
				ORb00_P01, ORb01_P01, ORb02_P01, ORb03_P01, ORb04_P01, ORb05_P01,		# Backward reaction coefficient
				OK0_P01, OK1_P01, OK2_P01, OK3_P01, OK4_P01, OK5_P01,  # equilibrium constant
				ORf00_P01, ORf01_P01, ORf02_P01, ORf03_P01, ORf04_P01, ORf05_P01,  # forward reaction coefficient
				OZ_P01,
				twotee]

def constituent(model):
	if model == 'none':
		return []
	elif model == 'GGS' :
		return [OwdotO_GGS,OwdotO2_GGS,OwdotN_GGS,OwdotN2_GGS,OwdotNO_GGS]

	elif model == 'Park01':
		return [OwdotO_P01,OwdotO2_P01,OwdotN_P01,OwdotN2_P01,OwdotNO_P01]



twotee			= "Eq(twoT, (Tv*T)**(0.5))"


# Chemical equations and constant for OpenSBLI cases
# Park 1989 [GGS]

OwdotO_GGS  = "Eq( wdotO  , MO  * (2.0 * (Rf0-Rb0) + 2.0*(Rf1-Rb1)+1.0*(Rf6-Rb6)-1.0*(Rf7-Rb7)-1.0*(Rf8-Rb8)) )"
OwdotO2_GGS = "Eq( wdotO2 , MO2 * (-1.0*(Rf0-Rb0)-1.0*(Rf1-Rb1)+1.0*(Rf7-Rb7)) )"
OwdotN_GGS  = "Eq( wdotN  , MN  * (2.0*(Rf2-Rb2)+2.0*(Rf3-Rb3)+2.0*(Rf4-Rb4)+2.0*(Rf5-Rb5)+1.0*(Rf6-Rb6)+1.0*(Rf7-Rb7)+1.0*(Rf8-Rb8)) )"
OwdotN2_GGS = "Eq( wdotN2 , MN2 * (-1.0*(Rf2-Rb2)-1.0*(Rf3-Rb3)-1.0*(Rf4-Rb4)-1.0*(Rf5-Rb5)-1.0*(Rf8-Rb8)) )"
OwdotNO_GGS = "Eq( wdotNO , MNO * (-1.0*(Rf6-Rb6)-1.0*(Rf7-Rb7)+1.0*(Rf8-Rb8)) )"

OZ_GGS  	= "Eq(Z, 10000.0/T)"
ORb00_GGS 	= "Eq(Rb00, Rf00/K0)"
ORb01_GGS 	= "Eq(Rb01, Rf01/K1)"
ORb02_GGS 	= "Eq(Rb02, Rf02/K2)"
ORb03_GGS 	= "Eq(Rb03, Rf03/K3)"
ORb04_GGS 	= "Eq(Rb04, Rf04/K4)"
ORb05_GGS 	= "Eq(Rb05, Rf05/K5)"
ORb06_GGS 	= "Eq(Rb06, Rf06/K6)"
ORb07_GGS 	= "Eq(Rb07, Rf07/K7)"
ORb08_GGS 	= "Eq(Rb08, Rf08/K8)"

ORf0_GGS 	= "Eq(Rf0 , Rf00*(0.001*rhoO2/MO2)*(0.001*rhoN/MN+0.001*rhoO/MO) )"
ORb0_GGS 	= "Eq(Rb0 , Rb00*(0.001*rhoO/MO)**2.0*(0.001*rhoN/MN+0.001*rhoO/MO) )"
ORf1_GGS 	= "Eq(Rf1 , Rf01*(0.001*rhoO2/MO2)*(0.001*rhoN2/MN2+0.001*rhoO2/MO2+0.001*rhoNO/MNO) )"
ORb1_GGS 	= "Eq(Rb1 , Rb01*(0.001*rhoO/MO)**2.0*(0.001*rhoN2/MN2+0.001*rhoO2/MO2+0.001*rhoNO/MNO) )"
ORf2_GGS 	= "Eq(Rf2 , Rf02*(0.001*rhoN2/MN2)*(0.001*rhoN/MN) )"
ORb2_GGS 	= "Eq(Rb2 , Rb02*(0.001*rhoN/MN)**3.0 )"
ORf3_GGS 	= "Eq(Rf3 , Rf03*(0.001*rhoN2/MN2)*(0.001*rhoO/MO) )"
ORb3_GGS 	= "Eq(Rb3 , Rb03*(0.001*rhoN/MN)**2.0*(0.001*rhoO/MO) )"
ORf4_GGS 	= "Eq(Rf4 , Rf04*(0.001*rhoN2/MN2)*(0.001*rhoN2/MN2+0.001*rhoO2/MO2) )"
ORb4_GGS 	= "Eq(Rb4 , Rb04*(0.001*rhoN/MN)**2.0*(0.001*rhoN2/MN2+0.001*rhoO2/MO2) )"
ORf5_GGS 	= "Eq(Rf5 , Rf05*(0.001*rhoN2/MN2)*(0.001*rhoNO/MNO) )"
ORb5_GGS 	= "Eq(Rb5 , Rb05*(0.001*rhoN/MN)**2.0*(0.001*rhoNO/MNO) )"
ORf6_GGS 	= "Eq(Rf6 , Rf06*(0.001*rhoNO/MNO)*(0.001*rhoO/MO+0.001*rhoO2/MO2+0.001*rhoN/MN+0.001*rhoN2/MN2+0.001*rhoNO/MNO) )"
ORb6_GGS	= "Eq(Rb6 , Rb06*(0.001*rhoN/MN)*(0.001*rhoO/MO)*(0.001*rhoO/MO+0.001*rhoO2/MO2+0.001*rhoN/MN+0.001*rhoN2/MN2+0.001*rhoNO/MNO) )"
ORf7_GGS 	= "Eq(Rf7 , Rf07*(0.001*rhoNO/MNO)*(0.001*rhoO/MO) )"
ORb7_GGS 	= "Eq(Rb7 , Rb07*(0.001*rhoO2/MO2)*(0.001*rhoN/MN) )"
ORf8_GGS 	= "Eq(Rf8 , Rf08*(0.001*rhoN2/MN2)*(0.001*rhoO/MO) )"
ORb8_GGS 	= "Eq(Rb8 , Rb08*(0.001*rhoNO/MNO)*(0.001*rhoN/MN) )"

OK0_GGS = "Eq(K0, exp(A00+A01*log(Z)+A02*Z+A03*Z**2.0+A04*Z**3.0) )"
OK1_GGS = "Eq(K1, exp(A10+A11*log(Z)+A12*Z+A13*Z**2.0+A14*Z**3.0) )"
OK2_GGS = "Eq(K2, exp(A20+A21*log(Z)+A22*Z+A23*Z**2.0+A24*Z**3.0) )"
OK3_GGS = "Eq(K3, exp(A30+A31*log(Z)+A32*Z+A33*Z**2.0+A34*Z**3.0) )"
OK4_GGS = "Eq(K4, exp(A40+A41*log(Z)+A42*Z+A43*Z**2.0+A44*Z**3.0) )"
OK5_GGS = "Eq(K5, exp(A50+A51*log(Z)+A52*Z+A53*Z**2.0+A54*Z**3.0) )"
OK6_GGS = "Eq(K6, exp(A60+A61*log(Z)+A62*Z+A63*Z**2.0+A64*Z**3.0) )"
OK7_GGS = "Eq(K7, exp(A70+A71*log(Z)+A72*Z+A73*Z**2.0+A74*Z**3.0) )"
OK8_GGS = "Eq(K8, exp(A80+A81*log(Z)+A82*Z+A83*Z**2.0+A84*Z**3.0) )"

ORf00_GGS = "Eq(Rf00, 1000.0*Cf0*(twoT**nf0)*exp(-thetaf0/twoT))"
ORf01_GGS = "Eq(Rf01, 1000.0*Cf1*(twoT**nf1)*exp(-thetaf1/twoT))"
ORf02_GGS = "Eq(Rf02, 1000.0*Cf2*(twoT**nf2)*exp(-thetaf2/twoT))"
ORf03_GGS = "Eq(Rf03, 1000.0*Cf3*(twoT**nf3)*exp(-thetaf3/twoT))"
ORf04_GGS = "Eq(Rf04, 1000.0*Cf4*(twoT**nf4)*exp(-thetaf4/twoT))"
ORf05_GGS = "Eq(Rf05, 1000.0*Cf5*(twoT**nf5)*exp(-thetaf5/twoT))"
ORf06_GGS = "Eq(Rf06, 1000.0*Cf6*(twoT**nf6)*exp(-thetaf6/twoT))"
ORf07_GGS = "Eq(Rf07, 1000.0*Cf7*(twoT**nf7)*exp(-thetaf7/twoT))"
ORf08_GGS = "Eq(Rf08, 1000.0*Cf8*(twoT**nf8)*exp(-thetaf8/twoT))"

# Park 1989 (GGS)
constants_GGS  = {
	'Cf0':'2.900e23', 	'nf0':'-2.0', 	'thetaf0':'5.975e4',	'A00':'2.855',	'A01':'0.988',	'A02':'-6.181',	'A03':'-0.023',	'A04':'-0.001',
	'Cf1':'9.680e22', 	'nf1':'-2.0', 	'thetaf1':'5.975e4',	'A10':'2.855',	'A11':'0.988',	'A12':'-6.181',	'A13':'-0.023',	'A14':'-0.001',
	'Cf2':'1.600e22', 	'nf2':'-1.6', 	'thetaf2':'1.132e5',	'A20':'1.858',	'A21':'-1.325',	'A22':'-9.856',	'A23':'-0.174',	'A24':'0.008',
	'Cf3':'4.980e22', 	'nf3':'-1.6', 	'thetaf3':'1.132e5',	'A30':'1.858',	'A31':'-1.325',	'A32':'-9.856',	'A33':'-0.174',	'A34':'0.008',
	'Cf4':'3.700e21', 	'nf4':'-1.6', 	'thetaf4':'1.132e5',	'A40':'1.858',	'A41':'-1.325',	'A42':'-9.856',	'A43':'-0.174',	'A44':'0.008',
	'Cf5':'4.980e21', 	'nf5':'-1.6', 	'thetaf5':'1.132e5',	'A50':'1.858',	'A51':'-1.325',	'A52':'-9.856',	'A53':'-0.174',	'A54':'0.008',
	'Cf6':'7.950e23', 	'nf6':'-2.0', 	'thetaf6':'7.550e4',	'A60':'0.792',	'A61':'-0.492',	'A62':'-6.761',	'A63':'-0.091',	'A64':'0.004',
	'Cf7':'8.370e12', 	'nf7':'0.0', 	'thetaf7':'1.945e4',	'A70':'-2.063',	'A71':'-1.480',	'A72':'-0.580',	'A73':'-0.114',	'A74':'0.005',
	'Cf8':'6.440e17', 	'nf8':'-1.0', 	'thetaf8':'3.837e4',	'A80':'1.066',	'A81':'-0.833',	'A82':'-3.095',	'A83':'-0.084',	'A84':'0.004',
}



##########################################################

# # Park 2001

#########################################################
constants_P01  = {
	'Cf0':'1.0e22', 	'nf0':'-1.5', 	'thetaf0':'5.936e4',	'A00':'1.578640',	'A01':'2.688744',	'A02':'4.215573',	'A03':'-8.091354',	'A04':'0.174260',
	'Cf1':'2.0e21', 	'nf1':'-1.5', 	'thetaf1':'5.936e4',	'A10':'1.578640',	'A11':'2.688744',	'A12':'4.215573',	'A13':'-8.091354',	'A14':'0.174260',
	'Cf2':'3.0e22', 	'nf2':'-1.6', 	'thetaf2':'1.132e5',	'A20':'-3.293682',	'A21':'0.998998',	'A22':'-8.237028',	'A23':'-5.526183',	'A24':'-0.582174',
	'Cf3':'7.0e21', 	'nf3':'-1.6', 	'thetaf3':'1.132e5',	'A30':'-3.293682',	'A31':'0.998998',	'A32':'-8.237028',	'A33':'-5.526183',	'A34':'-0.582174',
	'Cf4':'8.4e12', 	'nf4':'0.00', 	'thetaf4':'1.940e4',	'A40':'-1.840133',	'A41':'-1.768215',	'A42':'-4.759554',	'A43':'1.153872',	'A44':'-0.238985',
	'Cf5':'5.7e12', 	'nf5':'0.42', 	'thetaf5':'4.2938e4',	'A50':'-3.032189',	'A51':'0.0784648',	'A52':'-7.693047',	'A53':'1.411299',	'A54':'-0.517448',
}


# # Park 2001:
OwdotO_P01  = "Eq( wdotO  , MO * (2.0 *  (Rf0 - Rb0) + 2.0 * (Rf1 - Rb1) - 1.0 * (Rf4 - Rb4) - 1.0 * (Rf5 - Rb5)) )"
OwdotO2_P01 = "Eq( wdotO2 , MO2 * (-1.0 *(Rf0 - Rb0) - 1.0 * (Rf1 - Rb1) + 1.0 * (Rf4 - Rb4)) )"
OwdotN_P01  = "Eq( wdotN  , MN * (2.0 *  (Rf2 - Rb2) + 2.0 * (Rf3 - Rb3) + 1.0 * (Rf4 - Rb4) + 1.0 * (Rf5 - Rb5)) )"
OwdotN2_P01 = "Eq( wdotN2 , MN2 * (-1.0 *(Rf2 - Rb2) - 1.0 * (Rf3 - Rb3) - 1.0 * (Rf5 - Rb5)) )"
OwdotNO_P01 = "Eq( wdotNO , MNO * (-1.0 *(Rf4 - Rb4) + 1.0 * (Rf5 - Rb5)) )"


ORf00_P01 = "Eq(Rf00, 1000.0*Cf0*(twoT**nf0)*exp(-thetaf0/twoT))"
ORf01_P01 = "Eq(Rf01, 1000.0*Cf1*(twoT**nf1)*exp(-thetaf1/twoT))"
ORf02_P01 = "Eq(Rf02, 1000.0*Cf2*(twoT**nf2)*exp(-thetaf2/twoT))"
ORf03_P01 = "Eq(Rf03, 1000.0*Cf3*(twoT**nf3)*exp(-thetaf3/twoT))"
ORf04_P01 = "Eq(Rf04, 1000.0*Cf4*(twoT**nf4)*exp(-thetaf4/twoT))"
ORf05_P01 = "Eq(Rf05, 1000.0*Cf5*(twoT**nf5)*exp(-thetaf5/twoT))"

OZ_P01  	= "Eq(Z, 10000.0/T)"
ORb00_P01 	= "Eq(Rb00, Rf00/K0)"
ORb01_P01 	= "Eq(Rb01, Rf01/K1)"
ORb02_P01 	= "Eq(Rb02, Rf02/K2)"
ORb03_P01 	= "Eq(Rb03, Rf03/K3)"
ORb04_P01 	= "Eq(Rb04, Rf04/K4)"
ORb05_P01 	= "Eq(Rb05, Rf05/K5)"

ORf0_P01 	= "Eq(Rf0 , Rf00 * (0.001 * rhoO2 / MO2) * (0.001 * rhoN / MN + 0.001 * rhoO / MO) )"
ORb0_P01 	= "Eq(Rb0 , Rb00 * ((0.001 * rhoO / MO) ** 2) * (0.001 * rhoN / MN + 0.001 * rhoO / MO) )"
ORf1_P01 	= "Eq(Rf1 , Rf01 * (0.001 * rhoO2 / MO2) * (0.001 * rhoN2 / MN2 + 0.001 * rhoO2 / MO2 + 0.001 * rhoNO / MNO) )"
ORb1_P01 	= "Eq(Rb1 , Rb01 * ((0.001 * rhoO / MO) ** 2) * (0.001 * rhoN2 / MN2 + 0.001 * rhoO2 / MO2 + 0.001 * rhoNO / MNO) )"
ORf2_P01 	= "Eq(Rf2 , Rf02 * (0.001 * rhoN2 / MN2) * (0.001 * rhoN / MN + 0.001 * rhoO / MO) )"
ORb2_P01 	= "Eq(Rb2 , Rb02 * ((0.001 * rhoN / MN) ** 2) * (0.001 * rhoN / MN + 0.001 * rhoO / MO) )"
ORf3_P01 	= "Eq(Rf3 , Rf03 * (0.001 * rhoN2 / MN2) * (0.001 * rhoN2 / MN2 + 0.001 * rhoO2 / MO2 + (0.001 * rhoNO / MNO)) )"
ORb3_P01 	= "Eq(Rb3 , Rb03 * ((0.001 * rhoN / MN) ** 2) * (0.001 * rhoN2 / MN2 + 0.001 * rhoO2 / MO2 + (0.001 * rhoNO / MNO)) )"
ORf4_P01 	= "Eq(Rf4 , Rf04 * (0.001 * rhoNO / MNO) * (0.001 * rhoO / MO) )"
ORb4_P01 	= "Eq(Rb4 , Rb04 * (0.001 * rhoO2 / MO2) * (0.001 * rhoN / MN) )"
ORf5_P01 	= "Eq(Rf5 , Rf05 * (0.001 * rhoN2 / MN2) * (0.001 * rhoO / MO) )"
ORb5_P01 	= "Eq(Rb5 , Rb05 * (0.001 * rhoNO / MNO) * (0.001 * rhoN / MN) )"

OK0_P01 = "Eq(K0, exp(A00/Z+A01+A02*log(Z)+A03*Z+A04*(Z**2.0)))"
OK1_P01 = "Eq(K1, exp(A10/Z+A11+A12*log(Z)+A13*Z+A14*(Z**2.0)))"
OK2_P01 = "Eq(K2, exp(A20/Z+A21+A22*log(Z)+A23*Z+A24*(Z**2.0)))"
OK3_P01 = "Eq(K3, exp(A30/Z+A31+A32*log(Z)+A33*Z+A34*(Z**2.0)))"
OK4_P01 = "Eq(K4, exp(A40/Z+A41+A42*log(Z)+A43*Z+A44*(Z**2.0)))"
OK5_P01 = "Eq(K5, exp(A50/Z+A51+A52*log(Z)+A53*Z+A54*(Z**2.0)))"
