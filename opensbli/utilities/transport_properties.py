''' 
    transport_properties

    author. 
    last update 26.06.2023

    syntax: 
            RE = add_transport_properties(viscosity='viscosity', thermal_conductivity='viscosity') 
            
            constituent_eqns, substitutions = add_transport_properties(constituent_eqns, substitutions, viscosity='viscosity', thermal_conductivity='viscosity') 

            potential options: viscosity, ndsblottner
 
        update: 
        user doesn't have to do anything more. 

'''

from sympy import flatten, Idx, sqrt, Rational, pprint, factor, nsimplify, collect
from opensbli.core.opensbliobjects import ConstantObject, ConstantIndexed, Globalvariable
from opensbli.core.grid import GridVariable
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.kernel import Kernel
from opensbli.core.datatypes import Int
from opensbli.core.kernel import MultiConstantsToDeclare as MCTD
from opensbli.core.parsing import EinsteinEquation
# from opensbli.utilities.helperfunctions import substitute_simulation_parameters

class viscosity(object):
    '''
      viscosity viscosity models - 
        sutherland
    '''

    def __init__(self, constituent_eqns, substitutions, constants, model='sutherland'):
        # assign relations
        self.constituent_eqns = constituent_eqns
        self.substitutions = substitutions
        self.constants = constants
       
        self.model = model

    def __repr__(self):
        return self.contituent_eqns,self.substitutions, self.constants

    def add_viscosity(self, species):
        constituent_equations = []
        
        if self.model == 'sutherland':
            viscosity = 'Eq(mu  , (T**(1.5)) * (1.0 + SuthT/Tinf)/(T + SuthT/Tinf)'
        elif self.model == 'ndsblottner':
            viscosity = 'Eq(mu, ((yO2+yN2+yNO)*0.1*exp(-11.2202)*T**(0.021823*log(T)+0.34357)+(yO+yN)*0.1*exp(-11.7344)*T**(0.022652*log(T)+0.342509))*(1.0-exp(-0.010568*T)) )' # NDS simplication of Blottner/Sutherland
            # ((yO2+yN2+yNO)*0.1*exp(-11.2202)*T**(0.021823*log(T)+0.34357)+(yO+yN)*0.1*exp(-11.7344)*T**(0.022652*log(T)+0.342509))*(1.0-exp(-0.010568*T))
            viscosity = 'Eq(mu, (yO2 + yN2 + yNO)*0.1*exp(-11.2202)*T**(0.021823*log(T)+0.34357)           )'


        Tinf = ConstantObject('Tinf')
        constituent_equations += [viscosity]

        input_constants  = []
        input_constants += ['Cf0', 'nf0', 'thetaf0', 'A01', 'A02', 'A03', 'A04', 'A05']

        return constituent_equations, input_constants
 
class park01(object):
    '''
      park01 model (implementation of Ali Musawi's work) 

      Park, C., Jaffe, R. L.,Partridge, H. (2001). Chemical-Kinetic Parameters of Hyperbolic Earth Entry. 
                Journal of Thermophysics and Heat Transfer,15(1) , 76-90. https://doi.org/10.2514/2.6582

      Park, C. (1993). Review of Chemical-Kinetic Problems of Future NASA Missions, I: Earth Entries. 
                Journal of Thermophysics and Heat Transfer,(3) , 385–398. https://doi.org/10.2514/3.431

      (Park et al. 2001 pg 84, table 5, r = 1,2,11,12) 
        set of equations:

        Reactions

        R0 - O2 + M  ↔ 2O + M  (M = N, O); 					alphar1, betar1 
        R1 - O2 + M  ↔ 2O + M  (M = N2, O2, NO)				alphar2, betar2 

        R2 - N2 + M  ↔ 2N + M  (M = O, N)					alphar3, betar3
        R3 - N2 + M  ↔ 2N + M  (M = N2, O2, NO)				alphar4, betar4

        R4 - NO + O  ↔ 2N + M  (M = N2, O2)					alphar5, betar5
        R5 - N2 + O  ↔  N + NO 								alphar6, betar6

      stoichiometric coefficients are self explanatory, but are included for posterity 

    Rf_s and Rb_s equations are eqs.42 and 43 from viscosity paper. 
    '''

    def __init__(self, constituent_eqns, substitutions, constants, species, reaction_constants, two_temperature):
        # assign relations
        self.constituent_eqns = constituent_eqns
        self.substitutions = substitutions
        self.constants = constants
       
        self.two_temperature = two_temperature
        self.reaction_constants = reaction_constants
        # self.type = type

        self.constituent_eqns += self.rate_equations(species)
        self.substitutions    += self.reaction_dependencies()
        self.constants        += self.generate_constants()

    def __repr__(self):
        return self.contituent_eqns,self.substitutions, self.constants

    def rate_equations(self, species):
        constituent_equations = []

        if species == 5: 
            
            wdotO  = 'Eq( wdotO  , MO  * (2.0 *  (Rf0 - Rb0) + 2.0 * (Rf1 - Rb1) - 1.0 * (Rf4 - Rb4) - 1.0 * (Rf5 - Rb5)) )'
            wdotO2 = 'Eq( wdotO2 , MO2 * (-1.0 *(Rf0 - Rb0) - 1.0 * (Rf1 - Rb1) + 1.0 * (Rf4 - Rb4)) )'
            wdotN  = 'Eq( wdotN  , MN  * (2.0 *  (Rf2 - Rb2) + 2.0 * (Rf3 - Rb3) + 1.0 * (Rf4 - Rb4) + 1.0 * (Rf5 - Rb5)) )'
            wdotN2 = 'Eq( wdotN2 , MN2 * (-1.0 *(Rf2 - Rb2) - 1.0 * (Rf3 - Rb3) - 1.0 * (Rf5 - Rb5)) )'
            wdotNO = 'Eq( wdotNO , MNO * (-1.0 *(Rf4 - Rb4) + 1.0 * (Rf5 - Rb5)) )'

            constituent_equations += [wdotO2, wdotO, wdotN2, wdotN, wdotNO]

        return constituent_equations

    def reaction_dependencies(self):
        substitutions = []
        # forward and backward rate equations, Rfr, and Rbr respectively forming park01 framework
        Rf0 = 'Eq(Rf0, 1000.0 * kf0         * (0.001* rhoO2/MO2)**(1.0) * (0.001* rhoN/MN + 0.001* rhoO/MO)**(1.0))'
        Rb0 = 'Eq(Rb0, 1000.0 * (kf0 / Kc0) * (0.001*   rhoO/MO)**(2.0) * (0.001* rhoN/MN + 0.001* rhoO/MO)**(1.0))'

        Rf1 = 'Eq(Rf1, 1000.0 * kf1         * (0.001* rhoO2/MO2)**(1.0) * (0.001* rhoN2/MN2 + 0.001* rhoO2/MO2 + 0.001* rhoNO/MNO)**(1.0))'
        Rb1 = 'Eq(Rb1, 1000.0 * (kf1 / Kc1) * (0.001*   rhoO/MO)**(2.0) * (0.001* rhoN2/MN2 + 0.001* rhoO2/MO2 + 0.001* rhoNO/MNO)**(1.0))'

        Rf2 = 'Eq(Rf2, 1000.0 * kf2         * (0.001* rhoN2/MN2)**(1.0) * (0.001* rhoN/MN + 0.001* rhoO/MO)**(1.0))'
        Rb2 = 'Eq(Rb2, 1000.0 * (kf2 / Kc2) * (0.001*   rhoN/MN)**(2.0) * (0.001* rhoN/MN + 0.001* rhoO/MO)**(1.0))'

        Rf3 = 'Eq(Rf3, 1000.0 * kf3         * (0.001* rhoN2/MN2)**(1.0) * (0.001* rhoN2/MN2 + 0.001* rhoO2/MO2 + 0.001 * rhoNO/MNO)**(1.0))'
        Rb3 = 'Eq(Rb3, 1000.0 * (kf3 / Kc3) * (0.001*   rhoN/MN)**(2.0) * (0.001* rhoN2/MN2 + 0.001* rhoO2/MO2 + 0.001 * rhoNO/MNO)**(1.0))'

        Rf4 = 'Eq(Rf4, 1000.0 * kf4         * (0.001* rhoNO/MNO)**(1.0) * (0.001* rhoO/MO)**(1.0))'
        Rb4 = 'Eq(Rb4, 1000.0 * (kf4 / Kc4) * (0.001* rhoO2/MO2)**(1.0) * (0.001* rhoN/MN)**(1.0))'

        Rf5 = 'Eq(Rf5, 1000.0 * kf5         * (0.001* rhoN2/MN2)**(1.0) * (0.001* rhoO/MO)**(1.0))'
        Rb5 = 'Eq(Rb5, 1000.0 * (kf5 / Kc5) * (0.001* rhoNO/MNO)**(1.0) * (0.001* rhoN/MN)**(1.0))'
    
        substitutions += [Rf0,Rb0, Rf1,Rb1, Rf2,Rb2, Rf3,Rb3, Rf4,Rb4, Rf5,Rb5]

        # forward reaction rates, kfr, backward are already in the code since kbr = (kfr/ Kcr)
        kf0 = 'Eq(kf0, Cf0 * Tq**(nf0) * exp( - thetaf0/ Tq))' # here, Ef1/k values are provided in table IV, hence, theta1= Ef1/k
        kf1 = 'Eq(kf1, Cf1 * Tq**(nf1) * exp( - thetaf1/ Tq))' # 
        kf2 = 'Eq(kf2, Cf2 * Tq**(nf2) * exp( - thetaf2/ Tq))'
        kf3 = 'Eq(kf3, Cf3 * Tq**(nf3) * exp( - thetaf3/ Tq))'
        kf4 = 'Eq(kf4, Cf4 * Tq**(nf4) * exp( - thetaf4/ Tq))'
        kf5 = 'Eq(kf5, Cf5 * Tq**(nf5) * exp( - thetaf5/ Tq))'

        substitutions += [kf0, kf1, kf2, kf3, kf4, kf5]

        # equilibrium constant for reaction, r, employed by park01 curvefit
        Kc0 = 'Eq(Kc0, exp(A01/Z + A02 + A03*log(Z) + A04*Z + A05*(Z**2.0)))'
        Kc1 = 'Eq(Kc1, exp(A11/Z + A12 + A13*log(Z) + A14*Z + A15*(Z**2.0)))'
        Kc2 = 'Eq(Kc2, exp(A21/Z + A22 + A23*log(Z) + A24*Z + A25*(Z**2.0)))'
        Kc3 = 'Eq(Kc3, exp(A31/Z + A32 + A33*log(Z) + A34*Z + A35*(Z**2.0)))'
        Kc4 = 'Eq(Kc4, exp(A41/Z + A42 + A43*log(Z) + A44*Z + A45*(Z**2.0)))'
        Kc5 = 'Eq(Kc5, exp(A51/Z + A52 + A53*log(Z) + A54*Z + A55*(Z**2.0)))'
        substitutions += [Kc0, Kc1, Kc2, Kc3, Kc4, Kc5]

        # z-parameter for temperature and Tq value
        Z = 'Eq(Z, 10000/T)'
        Tq = 'Eq(Tq, Td)' # maybe get 2 temperature with Tq=Td? 

        if self.two_temperature:
            # two temperature model implementation, q = 0.5 
            Td = 'Eq(Td, (T*Tv)**0.5)'
            # or add in constant object q
        else:
            # two temperature model implementation, q = 0.5 
            Td = 'Eq(Td, T)'

        substitutions += [Z, Tq, Td]

        return substitutions
    
    def generate_constants(self):
        ''' constants extracted from table 5, 

        park01 constants
        '''
        
        Cf0, nf0, thetaf0, a01, a02, a03, a04, a05 = ConstantObject('Cf0'), ConstantObject('nf0'), ConstantObject('thetaf0'), ConstantObject('A01'), ConstantObject('A02'), ConstantObject('A03'), ConstantObject('A04'), ConstantObject('A05')
        Cf1, nf1, thetaf1, a11, a12, a13, a14, a15 = ConstantObject('Cf1'), ConstantObject('nf1'), ConstantObject('thetaf1'), ConstantObject('A11'), ConstantObject('A12'), ConstantObject('A13'), ConstantObject('A14'), ConstantObject('A15')
        Cf2, nf2, thetaf2, a21, a22, a23, a24, a25 = ConstantObject('Cf2'), ConstantObject('nf2'), ConstantObject('thetaf2'), ConstantObject('A21'), ConstantObject('A22'), ConstantObject('A23'), ConstantObject('A24'), ConstantObject('A25')
        Cf3, nf3, thetaf3, a31, a32, a33, a34, a35 = ConstantObject('Cf3'), ConstantObject('nf3'), ConstantObject('thetaf3'), ConstantObject('A31'), ConstantObject('A32'), ConstantObject('A33'), ConstantObject('A34'), ConstantObject('A35')
        Cf4, nf4, thetaf4, a41, a42, a43, a44, a45 = ConstantObject('Cf4'), ConstantObject('nf4'), ConstantObject('thetaf4'), ConstantObject('A41'), ConstantObject('A42'), ConstantObject('A43'), ConstantObject('A44'), ConstantObject('A45')
        Cf5, nf5, thetaf5, a51, a52, a53, a54, a55 = ConstantObject('Cf5'), ConstantObject('nf5'), ConstantObject('thetaf5'), ConstantObject('A51'), ConstantObject('A52'), ConstantObject('A53'), ConstantObject('A54'), ConstantObject('A55')
        

        input_constants  = []
        input_constants += ['Cf0', 'nf0', 'thetaf0', 'A01', 'A02', 'A03', 'A04', 'A05']
        input_constants += ['Cf1', 'nf1', 'thetaf1', 'A11', 'A12', 'A13', 'A14', 'A15']
        input_constants += ['Cf2', 'nf2', 'thetaf2', 'A21', 'A22', 'A23', 'A24', 'A25']
        input_constants += ['Cf3', 'nf3', 'thetaf3', 'A31', 'A32', 'A33', 'A34', 'A35']
        input_constants += ['Cf4', 'nf4', 'thetaf4', 'A41', 'A42', 'A43', 'A44', 'A45']
        input_constants += ['Cf5', 'nf5', 'thetaf5', 'A51', 'A52', 'A53', 'A54', 'A55']

        # kinetic model of park, park89 constants 
        Cf0.value, nf0.value, thetaf0.value, a01.value, a02.value, a03.value, a04.value, a05.value = 1.0e22, -1.5,  5.936e4,  1.578640, 2.688744 ,   4.215573, -8.091354,  0.174260
        Cf1.value, nf1.value, thetaf1.value, a11.value, a12.value, a13.value, a14.value, a15.value = 2.0e21, -1.5,  5.936e4,  1.578640, 2.688744 ,   4.215573, -8.091354,  0.174260
        Cf2.value, nf2.value, thetaf2.value, a21.value, a22.value, a23.value, a24.value, a25.value = 3.0e22, -1.6,  1.132e5, -3.293682, 0.998998 ,  -8.237028, -5.526183, -0.582174
        Cf3.value, nf3.value, thetaf3.value, a31.value, a32.value, a33.value, a34.value, a35.value = 7.0e21, -1.6,  1.132e5, -3.293682, 0.998998 ,  -8.237028, -5.526183, -0.582174
        Cf4.value, nf4.value, thetaf4.value, a41.value, a42.value, a43.value, a44.value, a45.value = 8.4e12, 0.00,  1.940e4, -1.840133, -1.768215,  -4.759554,  1.153872, -0.238985
        Cf5.value, nf5.value, thetaf5.value, a51.value, a52.value, a53.value, a54.value, a55.value = 5.7e12, 0.42, 4.2938e4, -3.032189, 0.0784648,  -7.693047,  1.411299, -0.517448

        # add constants to the .cpp file
        MCTD.add_constants([Cf0, nf0, thetaf0, a01, a02, a03, a04, a05])
        MCTD.add_constants([Cf1, nf1, thetaf1, a11, a12, a13, a14, a15])
        MCTD.add_constants([Cf2, nf2, thetaf2, a21, a22, a23, a24, a25])
        MCTD.add_constants([Cf3, nf3, thetaf3, a31, a32, a33, a34, a35])
        MCTD.add_constants([Cf4, nf4, thetaf4, a41, a42, a43, a44, a45])
        MCTD.add_constants([Cf5, nf5, thetaf5, a51, a52, a53, a54, a55])

        return input_constants
 

class park01am(object):
    '''
      park01 model (implementation of Ali Musawi's work with nasa 9 polynomial

      Park, C., Jaffe, R. L.,Partridge, H. (2001). Chemical-Kinetic Parameters of Hyperbolic Earth Entry. 
                Journal of Thermophysics and Heat Transfer,15(1) , 76-90. https://doi.org/10.2514/2.6582

      Park, C. (1993). Review of Chemical-Kinetic Problems of Future NASA Missions, I: Earth Entries. 
                Journal of Thermophysics and Heat Transfer,(3) , 385–398. https://doi.org/10.2514/3.431

      park model with nasa 9 polynomials, update of musawi. 

      (Park et al. 2001 pg 84, table 5, r = 1,2,11,12) 
        set of equations:

        Reactions


        R0 - N2 + M  ↔ 2N  + M  (M = N, O); 				alphar0, betar0
        R1 - N2 + M  ↔ 2N  + M  (M = N2, O2, NO)			alphar1, betar1

        R2 - O2 + M  ↔ 2O  + M  (M = O, N)					alphar2, betar2
        R3 - O2 + M  ↔ 2O  + M  (M = N2, O2, NO)			alphar3, betar3

        R4 - NO + O  ↔  N  + O + M (M = N, O, NO)			alphar4, betar4
        R5 - NO + M  ↔  N  + O + M (M = N2, O2) 			alphar5, betar5

        R6 - N2 + O  ↔  NO + N 								alphar6, betar6
        R7 - NO + O  ↔  O2 + N 								alphar7, betar7

      stoichiometric coefficients are self explanatory, but are included for posterity 

    Rf_s and Rb_s equations are eqs.42 and 43 from viscosity paper. 
    '''

    def __init__(self, constituent_eqns, substitutions, constants, species, reaction_constants, two_temperature):
        # assign relations
        self.constituent_eqns = constituent_eqns
        self.substitutions = substitutions
        self.constants = constants
       
        self.two_temperature = two_temperature
        self.reaction_constants = reaction_constants
        # self.type = type

        self.constituent_eqns += self.rate_equations(species)
        self.substitutions    += self.reaction_dependencies()
        self.constants        += self.generate_constants()

    def __repr__(self):
        return self.contituent_eqns,self.substitutions, self.constants

    def rate_equations(self, species):
        constituent_equations = []

        if species == 5: 
            
            wdotO  = "Eq( wdotO  , MO  * ( 2*(Rf2-Rb2) + 2*(Rf3-Rb3) + (Rf4-Rb4) + (Rf5-Rb5) - (Rf6-Rb6) - (Rf7-Rb7) ) )"
            wdotO2 = "Eq( wdotO2 , MO2 * ( - (Rf2-Rb2) - (Rf3-Rb3) + (Rf7-Rb7) ) )"
            wdotN  = "Eq( wdotN  , MN  * ( 2*(Rf0-Rb0) + 2*(Rf1-Rb1) + (Rf4-Rb4) + (Rf5-Rb5) + (Rf6-Rb6) + (Rf7-Rb7) ) )"
            wdotN2 = "Eq( wdotN2 , MN2 * ( - (Rf0-Rb0) - (Rf1-Rb1) - (Rf6-Rb6) ) )"
            wdotNO = "Eq( wdotNO , MNO * ( - (Rf4-Rb4) - (Rf5-Rb5) + (Rf6-Rb6) - (Rf7-Rb7) ) )"

            constituent_equations += [wdotO2, wdotO, wdotN2, wdotN, wdotNO]

        return constituent_equations

    def reaction_dependencies(self):
        substitutions = []
        # forward and backward rate equations, Rfr, and Rbr respectively forming park01 framework
        Rf0 = 'Eq(Rf0, 1000.0 * kf0  * (0.001* rhoN2/MN2)**(1.0) * (0.001* rhoO/MO + 0.001* rhoN/MN)**(1.0))'
        Rb0 = 'Eq(Rb0, 1000.0 * kb0 * (0.001*   rhoN/MN)**(2.0) * (0.001* rhoO/MO + 0.001* rhoN/MN)**(1.0))'

        Rf1 = 'Eq(Rf1, 1000.0 * kf1  * (0.001* rhoN2/MN2)**(1.0)  * (0.001* rhoO2/MO2 + 0.001* rhoN2/MN2 + 0.001* rhoNO/MNO)**(1.0))'
        Rb1 = 'Eq(Rb1, 1000.0 * kb1 * (0.001*   rhoN/MN)**(2.0)  * (0.001* rhoO2/MO2 + 0.001* rhoN2/MN2 + 0.001* rhoNO/MNO)**(1.0))'

        Rf2 = 'Eq(Rf2, 1000.0 * kf2  * (0.001* rhoO2/MO2)**(1.0) * (0.001* rhoO/MO + 0.001* rhoN/MN)**(1.0))'
        Rb2 = 'Eq(Rb2, 1000.0 * kb2 * (0.001*   rhoO/MO)**(2.0) * (0.001* rhoO/MO + 0.001* rhoN/MN)**(1.0))'

        Rf3 = 'Eq(Rf3, 1000.0 * kf3  * (0.001* rhoO2/MO2)**(1.0) * (0.001* rhoO2/MO2 + 0.001* rhoN2/MN2 + 0.001* rhoNO/MNO)**(1.0))'
        Rb3 = 'Eq(Rb3, 1000.0 * kb3 * (0.001*   rhoO/MO)**(2.0) * (0.001* rhoO2/MO2 + 0.001* rhoN2/MN2 + 0.001* rhoNO/MNO)**(1.0))'

        Rf4 = 'Eq(Rf4, 1000.0 * kf4  * (0.001* rhoNO/MNO)**(1.0) * (0.001* rhoO/MO + 0.001* rhoN/MN + 0.001* rhoNO/MNO)**(1.0))'
        Rb4 = 'Eq(Rb4, 1000.0 * kb4 * (0.001*   rhoN/MN)**(1.0) * (0.001* rhoO/MO)**(1.0) * (0.001* rhoO/MO + 0.001* rhoN/MN + 0.001* rhoNO/MNO)**(1.0))'

        Rf5 = 'Eq(Rf5, 1000.0 * kf5  * (0.001* rhoNO/MNO)**(1.0) * (0.001* rhoO2/MO2 + 0.001* rhoN2/MN2)**(1.0))'
        Rb5 = 'Eq(Rb5, 1000.0 * kb5 * (0.001*   rhoN/MN)**(1.0) * (0.001* rhoO/MO)**(1.0) * (0.001* rhoO2/MO2 + 0.001* rhoN2/MN2)**(1.0))'

        Rf6 = 'Eq(Rf6, 1000.0 * kf4  * (0.001* rhoN2/MN2)**(1.0) * (0.001* rhoO/MO)**(1.0))'
        Rb6 = 'Eq(Rb6, 1000.0 * kb4 * (0.001* rhoNO/MNO)**(1.0) * (0.001* rhoN/MN)**(1.0))'

        Rf7 = 'Eq(Rf7, 1000.0 * kf5  * (0.001* rhoNO/MNO)**(1.0) * (0.001* rhoO/MO)**(1.0))'
        Rb7 = 'Eq(Rb7, 1000.0 * kb5 * (0.001* rhoO2/MO2)**(1.0) * (0.001* rhoN/MN)**(1.0))'
    
        substitutions += [Rf0,Rb0, Rf1,Rb1, Rf2,Rb2, Rf3,Rb3, Rf4,Rb4, Rf5,Rb5, Rf6,Rb6, Rf7,Rb7]

        # forward reaction rates, kfr, backward are already in the code since kbr = (kfr/ Kcr)
        kf0 = 'Eq(kf0, Cf0 * Tq**(nf0) * exp( - thetaf0/ Tq))' # here, Ef1/k values are provided in table IV, hence, theta1= Ef1/k
        kf1 = 'Eq(kf1, Cf1 * Tq**(nf1) * exp( - thetaf1/ Tq))' # 
        kf2 = 'Eq(kf2, Cf2 * Tq**(nf2) * exp( - thetaf2/ Tq))'
        kf3 = 'Eq(kf3, Cf3 * Tq**(nf3) * exp( - thetaf3/ Tq))'
        kf4 = 'Eq(kf4, Cf4 * Tq**(nf4) * exp( - thetaf4/ Tq))'
        kf5 = 'Eq(kf5, Cf5 * Tq**(nf5) * exp( - thetaf5/ Tq))'
        kf6 = "Eq(kf6, Cf6 * Tq**(nf6) * exp( - thetaf6/ Tq))"
        kf7 = "Eq(kf7, Cf7 * Tq**(nf7) * exp( - thetaf7/ Tq))"

        substitutions += [kf0, kf1, kf2, kf3, kf4, kf5, kf6, kf7]

        # backward reaction rates, kbr, where kbr = kfr(T)/Kcr 
        kb0 = "Eq(kb0, (Cf0 * T **(nf0) * exp( - thetaf0/  T)) / Kc0)"
        kb1 = "Eq(kb1, (Cf1 * T **(nf1) * exp( - thetaf1/  T)) / Kc1)"
        kb2 = "Eq(kb2, (Cf2 * T **(nf2) * exp( - thetaf2/  T)) / Kc2)"
        kb3 = "Eq(kb3, (Cf3 * T **(nf3) * exp( - thetaf3/  T)) / Kc3)"
        kb4 = "Eq(kb4, (Cf4 * T **(nf4) * exp( - thetaf4/  T)) / Kc4)"
        kb5 = "Eq(kb5, (Cf5 * T **(nf5) * exp( - thetaf5/  T)) / Kc5)"
        kb6 = "Eq(kb6, (Cf6 * T **(nf6) * exp( - thetaf6/  T)) / Kc6)"
        kb7 = "Eq(kb7, (Cf7 * T **(nf7) * exp( - thetaf7/  T)) / Kc7)"

        substitutions += [kb0, kb1, kb2, kb3, kb4, kb5, kb6, kb7]

        # equilibrium constant for reaction, r, employed by (nasa9 in house am polynomials)
        Kc0 = "Eq(Kc0, exp(A00/(T**1.5) + A01/T + A02 + A03*log(T) + A04*T + A05*(T**2.0)))"
        Kc1 = "Eq(Kc1, exp(A10/(T**1.5) + A11/T + A12 + A13*log(T) + A14*T + A15*(T**2.0)))"
        Kc2 = "Eq(Kc2, exp(A20/(T**1.5) + A21/T + A22 + A23*log(T) + A24*T + A25*(T**2.0)))"
        Kc3 = "Eq(Kc3, exp(A30/(T**1.5) + A31/T + A32 + A33*log(T) + A34*T + A35*(T**2.0)))"
        Kc4 = "Eq(Kc4, exp(A40/(T**1.5) + A41/T + A42 + A43*log(T) + A44*T + A45*(T**2.0)))"
        Kc5 = "Eq(Kc5, exp(A50/(T**1.5) + A51/T + A52 + A53*log(T) + A54*T + A55*(T**2.0)))"
        Kc6 = "Eq(Kc6, exp(A60/(T**1.5) + A61/T + A62 + A63*log(T) + A64*T + A65*(T**2.0)))"
        Kc7 = "Eq(Kc7, exp(A70/(T**1.5) + A71/T + A72 + A73*log(T) + A74*T + A75*(T**2.0)))"

        substitutions += [Kc0, Kc1, Kc2, Kc3, Kc4, Kc5, Kc6, Kc7]

        # z-parameter for temperature and Tq value
        Z = 'Eq(Z, 10000/T)'
        Tq = 'Eq(Tq, Td)' # maybe get 2 temperature with Tq=Td? 

        if self.two_temperature:
            # two temperature model implementation, q = 0.5 
            Td = 'Eq(Td, (T*Tv)**0.5)'
            # or add in constant object q
        else:
            # two temperature model implementation, q = 0.5 
            Td = 'Eq(Td, T)'

        substitutions += [Z, Tq, Td]

        return substitutions
    
    def generate_constants(self):
        ''' constants extracted from table 5, 

        nasa9 constants 
        '''

        Cf0, nf0, thetaf0, a00, a01, a02, a03, a04, a05 = ConstantObject('Cf0'), ConstantObject('nf0'), ConstantObject('thetaf0'), ConstantObject('A00'), ConstantObject('A01'), ConstantObject('A02'), ConstantObject('A03'), ConstantObject('A04'), ConstantObject('A05')
        Cf1, nf1, thetaf1, a10, a11, a12, a13, a14, a15 = ConstantObject('Cf1'), ConstantObject('nf1'), ConstantObject('thetaf1'), ConstantObject('A10'), ConstantObject('A11'), ConstantObject('A12'), ConstantObject('A13'), ConstantObject('A14'), ConstantObject('A15')
        Cf2, nf2, thetaf2, a20, a21, a22, a23, a24, a25 = ConstantObject('Cf2'), ConstantObject('nf2'), ConstantObject('thetaf2'), ConstantObject('A20'), ConstantObject('A21'), ConstantObject('A22'), ConstantObject('A23'), ConstantObject('A24'), ConstantObject('A25')
        Cf3, nf3, thetaf3, a30, a31, a32, a33, a34, a35 = ConstantObject('Cf3'), ConstantObject('nf3'), ConstantObject('thetaf3'), ConstantObject('A30'), ConstantObject('A31'), ConstantObject('A32'), ConstantObject('A33'), ConstantObject('A34'), ConstantObject('A35')
        Cf4, nf4, thetaf4, a40, a41, a42, a43, a44, a45 = ConstantObject('Cf4'), ConstantObject('nf4'), ConstantObject('thetaf4'), ConstantObject('A40'), ConstantObject('A41'), ConstantObject('A42'), ConstantObject('A43'), ConstantObject('A44'), ConstantObject('A45')
        Cf5, nf5, thetaf5, a50, a51, a52, a53, a54, a55 = ConstantObject('Cf5'), ConstantObject('nf5'), ConstantObject('thetaf5'), ConstantObject('A50'), ConstantObject('A51'), ConstantObject('A52'), ConstantObject('A53'), ConstantObject('A54'), ConstantObject('A55')
        Cf6, nf6, thetaf6, a60, a61, a62, a63, a64, a65 = ConstantObject('Cf6'), ConstantObject('nf6'), ConstantObject('thetaf6'), ConstantObject('A60'), ConstantObject('A61'), ConstantObject('A62'), ConstantObject('A63'), ConstantObject('A64'), ConstantObject('A65')
        Cf7, nf7, thetaf7, a70, a71, a72, a73, a74, a75 = ConstantObject('Cf7'), ConstantObject('nf7'), ConstantObject('thetaf7'), ConstantObject('A70'), ConstantObject('A71'), ConstantObject('A72'), ConstantObject('A73'), ConstantObject('A74'), ConstantObject('A75')

        input_constants  = []
        input_constants += ['Cf0', 'nf0', 'thetaf0', 'A00', 'A01', 'A02', 'A03', 'A04', 'A05']
        input_constants += ['Cf1', 'nf1', 'thetaf1', 'A10', 'A11', 'A12', 'A13', 'A14', 'A15']
        input_constants += ['Cf2', 'nf2', 'thetaf2', 'A20', 'A21', 'A22', 'A23', 'A24', 'A25']
        input_constants += ['Cf3', 'nf3', 'thetaf3', 'A30', 'A31', 'A32', 'A33', 'A34', 'A35']
        input_constants += ['Cf4', 'nf4', 'thetaf4', 'A40', 'A41', 'A42', 'A43', 'A44', 'A45']
        input_constants += ['Cf5', 'nf5', 'thetaf5', 'A50', 'A51', 'A52', 'A53', 'A54', 'A55']
        input_constants += ['Cf6', 'nf6', 'thetaf6', 'A60', 'A61', 'A62', 'A63', 'A64', 'A65']
        input_constants += ['Cf7', 'nf7', 'thetaf7', 'A70', 'A71', 'A72', 'A73', 'A74', 'A75']

        # kinetic model of park, nasa9 constants, update ali 
        Cf0.value, nf0.value, thetaf0.value, a00.value, a01.value, a02.value, a03.value, a04.value, a05.value = 3.0e22, -1.6, 113200.0,  2.091870e+03,	-1.134452e+05, -6.427685e-01, 5.610081e-01, -3.300745e-04, 1.836761e-08
        Cf1.value, nf1.value, thetaf1.value, a10.value, a11.value, a12.value, a13.value, a14.value, a15.value = 7.0e21, -1.6, 113200.0,  2.091870e+03,	-1.134452e+05, -6.427685e-01, 5.610081e-01, -3.300745e-04, 1.836761e-08
        Cf2.value, nf2.value, thetaf2.value, a20.value, a21.value, a22.value, a23.value, a24.value, a25.value = 1.0e22, -1.5,  59360.0,  6.275151e+03,	-6.017939e+04,  3.082217e+00, 9.324624e-02, -2.147384e-04, 7.298221e-09
        Cf3.value, nf3.value, thetaf3.value, a30.value, a31.value, a32.value, a33.value, a34.value, a35.value = 2.0e21, -1.5,  59360.0,  6.275151e+03,	-6.017939e+04,  3.082217e+00, 9.324624e-02, -2.147384e-04, 7.298221e-09
        Cf4.value, nf4.value, thetaf4.value, a40.value, a41.value, a42.value, a43.value, a44.value, a45.value = 1.1e17, 0.00,  75500.0,  5.029763e+03,	-7.594795e+04,  7.563820e-01, 1.817658e-01, -6.685894e-06, 4.489270e-09
        Cf5.value, nf5.value, thetaf5.value, a50.value, a51.value, a52.value, a53.value, a54.value, a55.value = 5.0e15, 0.00,  75500.0,  5.029763e+03,	-7.594795e+04,  7.563820e-01, 1.817658e-01, -6.685894e-06, 4.489270e-09
        Cf6.value, nf6.value, thetaf6.value, a60.value, a61.value, a62.value, a63.value, a64.value, a65.value = 5.7e12, 0.42,  42938.0, -2.937898e+03,	-3.749727e+04, -1.399152e+00, 3.792425e-01, -1.086524e-04, 6.580220e-09
        Cf7.value, nf7.value, thetaf7.value, a70.value, a71.value, a72.value, a73.value, a74.value, a75.value = 8.4e12, 0.00,  19400.0, -1.245460e+03,	-1.576855e+04, -2.325883e+00, 8.852601e-02, -6.685894e-06, 4.489270e-09

        # add constants to the .cpp file
        MCTD.add_constants([Cf0, nf0, thetaf0, a00, a01, a02, a03, a04, a05])
        MCTD.add_constants([Cf1, nf1, thetaf1, a10, a11, a12, a13, a14, a15])
        MCTD.add_constants([Cf2, nf2, thetaf2, a20, a21, a22, a23, a24, a25])
        MCTD.add_constants([Cf3, nf3, thetaf3, a30, a31, a32, a33, a34, a35])
        MCTD.add_constants([Cf4, nf4, thetaf4, a40, a41, a42, a43, a44, a45])
        MCTD.add_constants([Cf5, nf5, thetaf5, a50, a51, a52, a53, a54, a55])
        MCTD.add_constants([Cf6, nf6, thetaf6, a60, a61, a62, a63, a64, a65])
        MCTD.add_constants([Cf7, nf7, thetaf7, a70, a71, a72, a73, a74, a75])

        return input_constants
    

class add_transport_properties(object):
    ''' adds rate equations into the opensbli python script'''
    # def __init__(self, split_type, ndim, constants, coordinate_symbol='x', conservative=True, viscosity=None, energy_formulation='none', debug=False):
    def __init__(self, constituent_eqns, substitutions, constants, simulation_type='perfectgas', viscosity='viscosity', thermal_conductivity='viscosity'):
        self.contituent_eqns = constituent_eqns
        self.substitutions   = substitutions
        # print(constants)
        self.constants = constants

        if simulation_type == 'perfectgas':
            
        else:
            print('generating source terms for equilibrium/nonequilibrum flow')
        
            if equation_type == 'viscosity':
                print('adding rate equation source terms of viscosity')
                viscosity(self.contituent_eqns, self.substitutions, self.constants, species, reaction_constants, two_temperature)
                # raise NotImplementedError('only frozen flow wdot terms are implemented')
            elif equation_type == 'park01':
                # raise NotImplementedError('only viscosity is currently implemented, updates will include park01')
                park01(self.contituent_eqns, self.substitutions, self.constants, species, reaction_constants, two_temperature)
            elif equation_type == 'park01am':
                # raise NotImplementedError('only viscosity is currently implemented, updates will include park01')
                park01am(self.contituent_eqns, self.substitutions, self.constants, species, reaction_constants, two_temperature)
            else:
                raise NotImplementedError('nothing included as of yet')

        return None

    def __repr__(self):
        return self.contituent_eqns,self.substitutions, self.constants