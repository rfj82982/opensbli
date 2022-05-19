""" David J. Lusher 09/21. Dispersion Relation Preserving (DRP) explicit filters."""

from opensbli import *
from sympy import pprint, Piecewise
from opensbli.core.opensbliobjects import DataObject, ConstantObject, GroupedPiecewise
from opensbli.equation_types.opensbliequations import OpenSBLIEquation
from opensbli.postprocess.post_process_eq import *
from opensbli.code_generation.algorithm.common import *
from opensbli.utilities.user_defined_kernels import UserDefinedEquations
from opensbli.core.kernel import ConstantsToDeclare as CTD

class DRPFilter(object):
    """ Selective filtering from Bogey & Bailly, A family of low dispersive and low dissipative explicit
    schemes for flow and noise computations, JoCP (2004) 194-214."""
    def __init__(self, block, width=11, q=None, optimized=False, sigma=0.1, wall_control=False):
        self.width, self.optimized = width, optimized
        print("Using a DRP filter with stencil width %d for block %d." % (self.width, block.blocknumber))
        self.depth = int(width/2.0)
        self.ndim = block.ndim
        self.block = block
        # Arrays to filter
        self.q_vector = [block.location_dataset(x) for x in flatten(q)]
        self.temp_arrays = [block.location_dataset('%s_RKold' % x.base.noblockname ) for x in self.q_vector]
        self.freq = ConstantObject('filter_frequency')
        self.freq.value = 25
        CTD.add_constant(self.freq)
        # Width and weightings of the filter
        self.generate_weights()
        self.sigma = ConstantObject('sigma_filt')
        self.sigma.value = sigma
        # Create the filter equations
        self.create_filter(block)
        return


    def detect_wall_boundaries(self):
        """ The shock-filter is turned off in the near-wall region. This function detects which directions, if any, have
        wall boundary conditions."""
        self.wall_boundaries = [[False, False] for _ in range(self.ndim)]
        try:
            for direction in range(self.ndim):
                for side in [0,1]:
                    if isinstance(self.block.boundary_types[direction][side], WallBC):
                        self.wall_boundaries[direction][side] = True
        except:
            raise ValueError("Please set boundary conditions on the block before calling the shock filter.")
        return

    def detect_interface_boundaries(self):
        """ The shock-filter is turned off close to block interfaces. This function detects which directions, if any, have
        interface boundary conditions."""
        self.interface_boundaries = [[False, False] for _ in range(self.ndim)]
        try:
            for direction in range(self.ndim):
                for side in [0,1]:
                    if isinstance(self.block.boundary_types[direction][side], InterfaceBC) or isinstance(self.block.boundary_types[direction][side], SharedInterfaceBC):
                        self.interface_boundaries[direction][side] = True
        except:
            raise ValueError("Please set boundary conditions on the block before calling the shock filter.")
        return

    def boundary_check(self, block):
        """ Checks if there are any non-periodic boundaries requirinig narrowing of the filter stencil."""
        self.non_periodic = [[False, False] for _ in range(block.ndim)]
        self.modify_directions = [False for _ in range(block.ndim)]
        for i, bc in enumerate(block.boundary_types):
            for side in [0,1]:
                if isinstance(bc[side], ModifyCentralDerivative):
                    self.non_periodic[i][side] = True
                    self.modify_directions[i] = True
        return

    def wall_control(self):
        """ Turns off the filter close to any of the walls or block interfaces in the problem."""
        buffer = 5
        wall_var = GridVariable('Wall')
        wall_conditions, wall_equations = [], []
        indexes = [OpenSBLIEq(GridVariable('Grid_%d' % direction), self.block.grid_indexes[direction]) for direction in range(self.ndim)]
        wall_equations += indexes
        # Disable the shock filter at any wall boundaries
        for direction in range(self.ndim):
            for side in [0,1]:
                wall = self.wall_boundaries[direction][side]
                if wall:
                    if side == 0:
                        wall_conditions += [ExprCondPair(0, indexes[direction].lhs <= buffer)]
                    else:
                        wall_conditions += [ExprCondPair(0, indexes[direction].lhs >= self.block.ranges[direction][side] - (buffer+1))]
        # Check for block interfaces
        for direction in range(self.ndim):
            for side in [0,1]:
                interface = self.interface_boundaries[direction][side]
                if interface:
                    if side == 0:
                        wall_conditions += [ExprCondPair(0, indexes[direction].lhs <= buffer)]
                    else:
                        wall_conditions += [ExprCondPair(0, indexes[direction].lhs >= self.block.ranges[direction][side] - (buffer+1))]
        # No wall or interface, default condition is the sensor is not turned off
        wall_conditions += [ExprCondPair(1, True)]
        wall_equations += [OpenSBLIEq(wall_var, Piecewise(*wall_conditions))]
        return wall_var, wall_equations

    def generate_weights(self):
        """ Weights are symmetric about the central point."""
        if self.width == 9:
            if self.optimized:
                self.weights = [0.008228661760,-0.045211119360,0.120007591680,-0.204788880640]
                self.weights += [0.243527493120] + self.weights[::-1]
            else:
                self.weights = [Rational(1,256),Rational(-1,32),Rational(7,64),Rational(-7,32)]
                self.weights += [Rational(35,128)] + self.weights[::-1]
        elif self.width == 11:
            if self.optimized:
                self.weights = [-0.002999540835,0.018721609157,-0.059227575576,0.123755948787,-0.187772883589]
                self.weights += [0.215044884112] + self.weights[::-1]
            else:
                self.weights = [Rational(-1,1024),Rational(5,512),Rational(-45,1024),Rational(15,128),Rational(-105,512)]
                self.weights += [Rational(63,256)] + self.weights[::-1]
        elif self.width == 13:
            if self.optimized:
                self.weights = [0.001254597714,-0.008520738659,0.029662754736,-0.069975429105,0.123632891797,-0.171503832236]
                self.weights += [0.190899511506] + self.weights[::-1]
            else:
                self.weights = [Rational(1,4096),Rational(-3,1024),Rational(33,2048),Rational(-55,1024),Rational(495,4096),Rational(-99,512)]
                self.weights += [Rational(231,1024)] + self.weights[::-1]
        self.locations = [i for i in range(-int(self.width/2.0), int(self.width/2.0)+1)]        
        return


    def create_stencil(self, direction):
        """ Indexes the datasets based on the width of the filter stencil."""
        output = []
        for dset_id, dset in enumerate(self.q_vector):
            stencil = []
            for i, location in enumerate(self.locations):
                stencil.append(self.weights[i]*increment_dataset(dset, direction, location))
            output += [OpenSBLIEq(self.temp_arrays[dset_id], sum(stencil))]
        return output

    def zero_temp_arrays(self):
        """ Ensure the temp arrays are zero everywhere."""
        zeroed = []
        for dset_id, dset in enumerate(self.temp_arrays):
            zeroed.append(OpenSBLIEq(dset, 0.0))
        return zeroed

    def create_equations(self, block, direction):
        # Create the indexed equations to calculate the filter
        application = self.create_stencil(direction)
        direction += 1
        # Update the q vector
        wall_var, wall_equations = self.wall_control()
        update = wall_equations[:]
        for dset_id, dset in enumerate(self.q_vector):
            update += [OpenSBLIEq(dset, dset - wall_var*self.sigma*self.temp_arrays[dset_id])]
        return application, update

    def create_UDF(self, block, equations, direction, order):
        UDF = UserDefinedEquations()
        UDF.algorithm_place = InTheSimulation(frequency=False)
        if order == 0:
            UDF.computation_name = 'Block %d: Zero the filter array' % block.blocknumber
        elif order == 1:
            UDF.computation_name = 'Block %d: DRP filter calculation direction %s' % (block.blocknumber, block.direction_labels[direction])
        else:
            UDF.computation_name = 'Block %d: DRP filter update direction %s' % (block.blocknumber, block.direction_labels[direction])
        # Place the filter at the very end
        UDF.order = 10000 + direction + order
        UDF.add_equations(equations)
        # # Attribute to modify ranges for non-periodic boundaries
        # if order is not 0 and self.modify_directions[direction]:
        #     UDF.modify_kernel_range = self.non_periodic
        #     UDF.depth = self.depth
        return UDF

    def create_filter(self, block):
        self.equation_classes = []
        # Zero the arrays
        zeroed = self.zero_temp_arrays()
        self.equation_classes += [self.create_UDF(block, zeroed, 0, 0)]
        # Check for non-periodic boundaries
        if self.wall_control:
            self.detect_wall_boundaries()
            self.detect_interface_boundaries()
        # Create a kernel at the end of the time loop, every iteration (no frequency)
        for direction in range(block.ndim):
            # Create the equations
            application, update = self.create_equations(block, direction)
            filter1 = self.create_UDF(block, application, direction, 1)
            filter2 = self.create_UDF(block, update, direction, 2)
            self.equation_classes += [filter1, filter2]
        return
