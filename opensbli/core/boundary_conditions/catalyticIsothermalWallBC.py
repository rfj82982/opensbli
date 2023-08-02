"""@brief This contains non catalytic isothermal wall boundary condition
   @authors Teja Ala
   @contributors 
   @details
"""

from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative, WallBC
from opensbli.physical_models.ns_physics import NSphysics, PhysicsVariable
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from sympy import sqrt, Matrix, pprint, exp, Mul, Pow
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.opensbliobjects import ConstantObject
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.core.grid import GridVariable


class catalyticIsothermalWallBC(ModifyCentralDerivative, BoundaryConditionBase, WallBC):
    """ catalytic wall condition, zero gradient dNO/dn = 0 over the boundary. (G.N. Coleman)
    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'catalyticIsothermalWall'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme

        self.deriv_indices = [[i for i in range(6)], [-i for i in range(6)]]
        self.deriv_coeffs = [[self.carp_coefficients()[0, i] for i in range(6)], [-self.carp_coefficients()[0, i] for i in range(6)]]
        return

    def carp_coefficients(self):
        """ Computes the finite-difference coefficients for the 1st order one-sided Carpenter wall boundary derivative.
        :returns: Matrix: bc4: Matrix of stencil coefficients."""
        R1 = -(2177.0*sqrt(295369.0)-1166427.0)/25488.0
        R2 = (66195.0*sqrt(53.0)*sqrt(5573.0)-35909375.0)/101952.0

        al4_0 = [-(216.0*R2+2160.0*R1-2125.0)/12960.0, (81.0*R2+675.0*R1+415.0)/540.0, -(72.0*R2+720.0*R1+445.0)/1440.0, -(108.0*R2+756.0*R1+421.0)/1296.0]
        al4_1 = [(81.0*R2+675.0*R1+415.0)/540.0, -(4104.0*R2+32400.0*R1+11225.0)/4320.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(216.0*R2+2160.0*R1+655.0)/4320.0]
        al4_2 = [-(72.0*R2+720.0*R1+445.0)/1440.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(4104.0*R2+32400.0*R1+12785.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0]
        al4_3 = [-(108.0*R2+756.0*R1+421.0)/1296.0, -(216.0*R2+2160.0*R1+655.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0, -(216.0*R2+2160.0*R1-12085.0)/12960.0]
        al4 = Matrix([al4_0, al4_1, al4_2, al4_3])

        ar4_0 = [(-1.0)/2.0, -(864.0*R2+6480.0*R1+305.0)/4320.0, (216.0*R2+1620.0*R1+725.0)/540.0, -(864.0*R2+6480.0*R1+3335.0)/4320.0, 0.0, 0.0]
        ar4_1 = [(864.0*R2+6480.0*R1+305.0)/4320.0, 0.0, -(864.0*R2+6480.0*R1+2315.0)/1440.0, (108.0*R2+810.0*R1+415.0)/270.0, 0.0, 0.0]
        ar4_2 = [-(216.0*R2+1620.0*R1+725.0)/540.0, (864.0*R2+6480.0*R1+2315.0)/1440.0, 0.0, -(864.0*R2+6480.0*R1+785.0)/4320.0, -1.0/12.0, 0.0]
        ar4_3 = [(864.0*R2+6480.0*R1+3335.0)/4320.0, -(108.0*R2+810.0*R1+415.0)/270.0, (864.0*R2+6480.0*R1+785.0)/4320.0, 0.0, 8.0/12.0, -1.0/12.0]
        ar4 = Matrix([ar4_0, ar4_1, ar4_2, ar4_3])
        # Form inverse and convert to rational
        al4_inv = al4.inv()
        bc4 = al4_inv*ar4
        return bc4

    def apply(self, arrays, block):

        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        wall_eqns = []
        direction, side = self.direction, self.side
        from_side_factor, to_side_factor = self.set_side_factor()
        NS = NSphysics(block)

        # append PhysicsVariable to the block
        PhysicsVariable.blocknumber = block.blocknumber
        PhysicsVariable.shape = block.shape
        PhysicsVariable.block = block

        for ar in arrays:
            if isinstance(ar, list):  # Set velocity components to zero on the wall
                rhs = [0 for i in range(len(ar))]
                wall_eqns += [OpenSBLIEq(x, y) for (x, y) in zip(ar, rhs)]

        # Define constants
        Rhat, Twall, pref = ConstantObject('Rhat'), ConstantObject('Twall'), ConstantObject('pref')
        thetavO2,thetavN2,thetavO,thetavN,thetavNO = ConstantObject('thetavO2'),ConstantObject('thetavN2'),ConstantObject('thetavO'),ConstantObject('thetavN'),ConstantObject('thetavNO')
        # constant variables for chemistry
        MNO,MN2,MN,MO,MO2 = ConstantObject('MNO'),ConstantObject('MN2'),ConstantObject('MN'),ConstantObject('MO'),ConstantObject('MO2') # setting rhoNO density for catalytic wall
        dhNO,dhN2,dhN,dhO,dhO2 = ConstantObject('dhNO'),ConstantObject('dhN2'),ConstantObject('dhN'),ConstantObject('dhO'),ConstantObject('dhO2') # setting rhoNO density for catalytic wall
        sigOtoNe = ConstantObject('sigOtoNe') # ratio of freestream O-N, currently definied in the main script

        # definite physical variables
        rN2v, rNv, rO2v, rOv, rNOv = PhysicsVariable('rhoN2'), PhysicsVariable('rhoN'), PhysicsVariable('rhoO2'), PhysicsVariable('rhoO'), PhysicsVariable('rhoNO')
        energy_storev, evib_storev = PhysicsVariable('rhoE'), PhysicsVariable('rhoev')

        # get variables
        rN2, rN, rO2, rO, rNO = rN2v.variable, rNv.variable, rO2v.variable, rOv.variable, rNOv.variable
        energy_store, evib_store = energy_storev.variable, evib_storev.variable

        cN2w, cNw, cO2w, cOw, cNOw, rhow = GridVariable('cN2w'), GridVariable('cNw'), GridVariable('cO2w'), GridVariable('cOw'), GridVariable('cNOw'), GridVariable('rhow')
        wall_eqns +=[OpenSBLIEq(GridVariable('cNw'), 0.0)]
        wall_eqns +=[OpenSBLIEq(GridVariable('cOw'), 0.0)]

        # Calculate cNO at the wall associated with  dcNO/dy = 0
        coeffs, grid_indices = self.deriv_coeffs[side], self.deriv_indices[side]
        # Evaluate yNO values at the wall and 5 points above/below it
        cNO_values = [GridVariable('cNO%d' % i) for i in range(len(grid_indices))]
        
        cNO_evaluations = [OpenSBLIEq(cNO_values[i], increment_dataset(rNO / (rN2 + rN + rO2 + rO + rNO), self.direction, i*to_side_factor))
                                   for i in range(1, len(cNO_values))]

        wall_eqns += cNO_evaluations
        deriv_equation = sum([c*cNO_values[i+1] for i, c in enumerate(coeffs[1:])])
        # Rearrange the equation to solve for the wall cNO
        deriv_equation = OpenSBLIEq(cNO_values[0], -1.0*deriv_equation/(coeffs[0]))  # set wall value of drhoNO/dy to 0...

        wall_eqns += [deriv_equation]
        wall_eqns += [OpenSBLIEq(GridVariable('cNOw'), cNO_values[0])]

        wall_eqns +=[OpenSBLIEq(GridVariable('cN2w'), (1.0-cNOw+cNOw/2.0*MO2/MNO*(1.0-sigOtoNe))/(1.0+sigOtoNe*MO2/MN2) )]
        wall_eqns +=[OpenSBLIEq(GridVariable('cO2w'), 1.0 - cN2w - cNOw)]

        # wall_eqns +=[OpenSBLIEq(GridVariable('rhow'), pref/(Rhat*Twall*(cN2w/MN2 + cNw/MN + cO2w/MO2 + cOw/MO + cNOw/MO)))]
        wall_eqns +=[OpenSBLIEq(GridVariable('rhow'), rN2 + rN + rO2 + rO + rNO )]

        wall_eqns += [OpenSBLIEq(rO, rhow * cOw)]
        wall_eqns += [OpenSBLIEq(rN, rhow * cNw)]
        wall_eqns += [OpenSBLIEq(rNO, rhow * cNOw)]
        wall_eqns += [OpenSBLIEq(rN2, rhow * cN2w)]
        wall_eqns += [OpenSBLIEq(rO2, rhow * cO2w)]

    
        # Calculate Temperature at the wall associated with  dT/dy = 0 <--------------------------------- only here to make all variables RW access, need to find a better way
        coeffs, grid_indices = self.deriv_coeffs[side], self.deriv_indices[side]
        # Evaluate temperature values at the wall and 5 points above/below it
        temperature_values = [GridVariable('T%d' % i) for i in range(len(grid_indices))]
        temperature_evaluations = [OpenSBLIEq(temperature_values[i], increment_dataset(NS.temperature(relation=True, conservative=True)+ evib_store, self.direction, i*to_side_factor))
                                   for i in range(1, len(temperature_values))]

        wall_eqns += temperature_evaluations
        deriv_equation = sum([c*temperature_values[i+1] for i, c in enumerate(coeffs[1:])])
        # Rearrange the equation to solve for the wall temperature
        deriv_equation = OpenSBLIEq(temperature_values[0], -1.0*deriv_equation/(coeffs[0] ))  # set wall value of dT/dy to 0...
        wall_eqns += [deriv_equation]

        # Set energy on the wall based on this calculated wall temperature and the density that comes from the continuity equation
        wall_eqns += [OpenSBLIEq(evib_store, Rhat*rN2*thetavN2/(MN2*exp(Mul(thetavN2, Pow(Twall, -1) )) -1.0) + Rhat*rO2*thetavO2/(MO2*exp(Mul(thetavO2, Pow(Twall,-1))) -1.0) + Rhat*rNO*thetavNO/(MNO*exp(Mul(thetavNO, Pow(Twall, -1))) -1.0) )]
        wall_eqns += [OpenSBLIEq(energy_store, evib_store + 4.1868e6*(rO*(dhO/MO)+rN*(dhN/MN)+rNO*(dhNO/MNO)) +  Mul(Twall, Rhat*((3.0/2.0)*((rO/MO)+(rN/MN))+(5.0/2.0)*((rO2/MO2)+(rN2/MN2)+(rNO/MNO)))) )]
        
        kernel.add_equation(wall_eqns)
        # Print out the current equations for the adiabatic wall conditio
        print("Printing equations for catalytic wall boundary condition (Carp-4):")  # <<<
        for eqn in kernel.equations:
            pprint(eqn)
        # exit()
        kernel.update_block_datasets(block)
        return kernel
