from sympy import flatten, Idx, sqrt, Rational, pprint, factor, nsimplify, collect
from opensbli.core.opensbliobjects import ConstantObject, ConstantIndexed, Globalvariable
from opensbli.core.grid import GridVariable
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.kernel import Kernel
from opensbli.core.datatypes import Int
from opensbli.core.parsing import EinsteinEquation

class NS_Split(object):
    """ Split forms for the convective parts of the Navier-Stokes equations with central/DRP schemes."""

    def __init__(self, split_type, ndim, constants, coordinate_symbol="x", conservative=True, viscosity=None):
        self.split_type = split_type
        self.conservative = conservative
        self.coordinate_symbol = coordinate_symbol
        self.constants = constants
        self.ndim = ndim
        self.viscosity = viscosity
        self.EE = EinsteinEquation()
        self.replace_factors = False
        # Storing either the conservative or primitive variables as the q vector to advance in time.
        if self.conservative:
            self.rhou = 'rhou'
            self.mom_lhs = 'rhou'
            self.energy_lhs = 'rhoE'
        else:
            self.rhou = 'rho*u'
            self.mom_lhs = 'u'
            self.energy_lhs = 'Et'
        # KGP coefficients
        if split_type == 'KGP':
            self.alpha = Rational(1,4)
            self.beta = Rational(1,4)
            self.delta = Rational(1,4)
            self.gamma = Rational(1,4)
            self.epsilon = 0
        # Diffusive terms
        self.substitutions = self.diffusive_terms()
        self.mass = self.continuity_eq()
        self.momentum = self.momentum_eq()
        self.energy = self.energy_eq()
        # self.diffusive = self.diffusive_eq()
        return

    def factor_replace(self, original_eqn):
        from sympy import symbols, count_ops, S
        print("Original operation count: {:}".format(original_eqn.count_ops()))
        a, b, c = ConstantObject('one_over_4'), ConstantObject('one_over_2'), ConstantObject('two_over_3')
        constant_dict = {a : Rational(1,4), b: Rational(1,2), c: Rational(2,3)}
        # Add the values
        for key, num in constant_dict.items():
            key.value = num
        reverse_dict = {v: k for k, v in constant_dict.items()}
        # Substitute the rational constants
        output = original_eqn.subs(reverse_dict)
        for key, value in constant_dict.items():
            output = collect(output, key)
        # Substitute simulation constants
        if ConstantObject('mu') in self.constants:
            output = collect(output, ConstantObject('mu')/ConstantObject('Re'))
        else: # variable viscosity
            output = collect(output, ConstantObject('Re'))
        # pprint(output)
        print("New operation count: {:}".format(output.count_ops()))
        return output

    def common_factors(self, eqn):
        """ Simplifies the equation by taking out common rational numbers."""
        lhs, rhs = eqn.lhs, eqn.rhs
        optimized = True
        if optimized:
            rhs = self.factor_replace(rhs)
        return OpenSBLIEq(lhs, rhs)

    def continuity_eq(self):
        if self.split_type == 'Feiereisen':
            out = "Eq(Der(rho, t), - Conservative(%s_j, x_j))" % self.rhou
        elif self.split_type == 'KGP':
            A, B, C, D = self.alpha, self.beta, self.gamma, self.delta
            if self.conservative:
                out = "Eq(Der(rho, t), - (%s*Conservative(rhou_j, x_j) + %s*Conservative(rhou_j, x_j) + %s*u_j*Der(rho, x_j) + %s*(rho*Der(u_j, x_j) + u_j*Der(rho, x_j))))" % (A, B, C, D)
            else:
                out = "Eq(Der(rho, t), - (%s*Conservative(rho*u_j, x_j) + %s*Conservative(rho*u_j, x_j) + %s*u_j*Der(rho, x_j) + %s*(rho*Der(u_j, x_j) + u_j*Der(rho, x_j))))" % (A, B, C, D)
        else:
            raise NotImplementedError("Only Feierisen and KGP splitting methods are implemented.")
        out = self.EE.expand(out, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)
        if self.replace_factors:
            out = self.common_factors(out)
        return out

    def momentum_eq(self):
        momentum = "Eq(Der(%s_i, t), - Der(p, x_i) + Der(tau_i_j, x_j))" % self.mom_lhs
        out = self.EE.expand(momentum, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)
        if self.split_type == 'Feiereisen':
            convective = "(1/2) * (Conservative(%s_i*u_j, x_j) + %s_j*Der(u_i,x_j) + u_i * Der(%s_j,x_j))" % (self.rhou, self.rhou, self.rhou)
        # Kennedy Gruber cubic split
        elif self.split_type == 'KGP':
            A, B, C, D = self.alpha, self.beta, self.gamma, self.delta
            if self.conservative:
                convective = "%s*Conservative(rhou_j*u_i, x_j) + %s*(u_i*Conservative(rhou_j, x_j) + rhou_j*Der(u_i, x_j)) + %s*(u_j*Conservative(rhou_i, x_j) + rhou_i*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*u_i, x_j) + u_i*u_j*Der(rho, x_j))" % (A, B, C, D)
            else:
                convective = "%s*Conservative(rho*u_j*u_i, x_j) + %s*(u_i*Conservative(rho*u_j, x_j) + rho*u_j*Der(u_i, x_j)) + %s*(u_j*Conservative(rho*u_i, x_j) + rho*u_i*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*u_i, x_j) + u_i*u_j*Der(rho, x_j))" % (A, B, C, D)
        else:
            raise NotImplementedError("Only Feierisen and KGP splitting methods are implemented.")
        # Add convective parts
        expanded_convective = self.EE.expand(convective, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)
        expanded_convective[0] = factor(expanded_convective[0])
        for no, value in enumerate(out):
            temp = OpenSBLIEq(out[no].lhs,  out[no].rhs - expanded_convective[no])
            if self.replace_factors:
                out[no] = self.common_factors(temp)
            else:
                out[no] = temp
        return out

    def energy_eq(self):
        if self.split_type == 'Feiereisen':
            if self.conservative:
                convective = "(1/2) * (Conservative(%s*u_j, x_j) + %s_j*Conservative(%s / rho, x_j) + (%s / rho) * Conservative(%s_j, x_j))" % (self.energy_lhs, self.rhou, self.energy_lhs, self.energy_lhs, self.rhou)
            else:
                convective = "(1/2) * (Conservative(rho*%s*u_j, x_j) + %s_j*Conservative(%s, x_j) + %s * Conservative(%s_j, x_j))" % (self.energy_lhs, self.rhou, self.energy_lhs, self.energy_lhs, self.rhou)
        elif self.split_type == 'KGP':
            # Split on phi = E
            A, B, C, D = self.alpha, self.beta, self.gamma, self.delta
            if self.conservative:
                convective = "(%s*Conservative(rhoE*u_j, x_j) + %s*((rhoE/rho)*Conservative(rhou_j, x_j) + rhou_j*Conservative(rhoE/rho, x_j)) + %s*(u_j*Der(rhoE, x_j) + rhoE*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*(rhoE/rho), x_j) + u_j*(rhoE/rho)*Der(rho, x_j)))" % (A, B, C, D)
            else:
                convective = "(%s*Conservative(rho*Et*u_j, x_j) + %s*(Et*Conservative(rho*u_j, x_j) + rho*u_j*Conservative(Et, x_j)) + %s*(u_j*Der(rho*Et, x_j) + rho*Et*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*Et, x_j) + u_j*Et*Der(rho, x_j)))" % (A, B, C, D)
        else:
            raise NotImplementedError("Only Feierisen and KGP splitting methods are implemented.")
        energy = "Eq(Der(%s, t), - %s - Conservative(p*u_j, x_j) + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j))" % (self.energy_lhs, convective)
        out = self.EE.expand(energy, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)
        if self.replace_factors:
            out = self.common_factors(out)
        return out

    def diffusive_terms(self):
        if self.viscosity == 'constant':
            stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)*Der(u_k,x_k)))" # *divV Der(u_k,x_k)
            heat_flux = "Eq(q_j, ((1.0/Re)/((gama-1)*Minf*Minf*Pr))*Der(T,x_j))"
        else:
            stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)*Der(u_k,x_k)))" # *divV Der(u_k,x_k)
            # stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)*divV))" # *divV Der(u_k,x_k)
            heat_flux = "Eq(q_j, ((mu/Re)/((gama-1)*Minf*Minf*Pr))*Der(T,x_j))"
        substitutions = [stress_tensor, heat_flux]
        return substitutions
