""" Max Walker & David J Lusher (09/2023): TVD non-linear filter for shock-capturing."""
from opensbli import *
from opensbli.equation_types.opensbliequations import NonSimulationEquations
from opensbli.postprocess.post_process_eq import *
from opensbli.schemes.spatial.TVD import *
from opensbli.filters.WENO_filter import NonLinearFilterBase

class TVDFilter(NonSimulationEquations, NonLinearFilterBase):
    """ Class to apply a TVD-based non-linear filter after a full time-step of a non-dissipative high order base scheme. The governing
    equations in the user script should be central derivatives in a skew-symmetric formulation to improve numerical stability."""
    def __init__(self, block, metrics=None, airfoil=False, passive_scalar=False):
        print("Using non-linear TVD filtering on block {:}.".format(block.blocknumber))
        self.passive_scalar = passive_scalar
        # Get the shared functionality between TVD/WENO non-linear filters
        NonLinearFilterBase.__init__(self, airfoil, block, metrics, optimize=False)
        # Main class to generate the filter
        self.main(block)
        return

    def TVD_filter_application(self, block):
        """ Applies the non-linear TVD filter by subtracting from the q vector after a full RK time-step."""
        resid_kernel = self.residual_kernels[0]
        filter_equations = []
        nvars = len(self.solution_vector)
        dt = ConstantObject('dt')
        for i, var in enumerate(self.solution_vector):
            filter_equations += [OpenSBLIEq(var, var + dt*resid_kernel.equations[i].rhs)]
        # Finish creating the kernel
        residual_kernel = self.create_kernel('Non-linear TVD Filter application', filter_equations, resid_kernel.halo_ranges, block)
        self.component_counter += 1
        self.add_kernel(residual_kernel)
        return

    def main(self, block):
        """ Main calling function to generate the kernels for the TVD filter."""
        # Counter to order the kernels. Put the TVD filtering kernels at the very end of the time loop
        self.component_counter = 1000 + block.blocknumber*1000
        # Create the equations for TVD
        if self.passive_scalar:
            eqn = self.Euler_equations_passive_scalar(block, 'TVD')
        else:
            eqn = self.Euler_equations(block, 'TVD')
        # Convert the equations to datasets on this block
        self.equations = self.convert_to_datasets(block, eqn)
        # Create a TVD scheme
        self.SF = TVDFlux(averaging=SimpleAverage([0, 1]), shock_filter=True, conservative=block.conservative, passive_scalar=self.passive_scalar)
        self.halo_type = set()
        self.halo_type.add(self.SF.halotype)
        # Start the discretisation and create residual arrays for the equations
        self.Kernels = []
        self.create_residual_arrays(block)
        CR, solution_vector, reductions = self.SF.discretise(self, block)
        # Q vector
        self.solution_vector = flatten(self.time_advance_arrays)
        # Swap over the TVD stencil if periodic boundaries
        # bc_kernels = self.update_periodic_boundary(block, self.SF.halotype)
        # Constituent relations evaluations on the Q vector at the end of the full RK time-step
        self.constituent_relations(block)
        # Zero the work arrays
        self.zero_work_arrays(block)
        # Create the TVD reconstruction kernels
        reconstruction_kernels = []
        for direction, ker in enumerate(self.reconstruction_kernels):
            halo_ranges = ker.halo_ranges
            reconstruction_kernels.append(self.create_kernel('TVD reconstruction direction %d' % direction, ker.equations, halo_ranges, block))
            self.component_counter += 1

        self.add_kernel(reconstruction_kernels)
        # Check if there any wall boundary conditions or interfaces defined on the block.
        self.detect_wall_boundaries()
        self.detect_interface_boundaries()
        # # Create the residual kernel
        self.TVD_filter_application(block)
        return
