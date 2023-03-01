#!/usr/bin/env python
from opensbli import *
from opensbli.utilities.user_defined_kernels import *
from opensbli.code_generation.algorithm.common import InTheSimulation, AfterSimulationEnds, BeforeSimulationStarts

def first_order_moments(ndim, q_vector, nsamples, conservative=True):
    # Create the statistics equations, this shows another way of writing the equations
    # create variables means required
    r_m = symbols('rhomean', **{'cls': DataObject})
    ru_m = symbols('rhou0:%dmean' % ndim, **{'cls': DataObject})
    E_mean = DataObject('E_mean')
    storage_arrays = [r_m] + list(ru_m) + [E_mean]
    # Step 1: Sum the quantities in time
    accumulation_equations = []
    # Accumulate density
    accumulation_equations += [Eq(r_m, r_m + q_vector[0])]
    # Accumulate_momentum
    for i in range(ndim):
        if conservative:
            accumulation_equations += [Eq(ru_m[i], ru_m[i] + q_vector[i+1])]
        else:
            accumulation_equations += [Eq(ru_m[i], ru_m[i] + q_vector[0]*q_vector[i+1])]
    # Accumulate total energy
    if conservative:
        accumulation_equations += [Eq(E_mean, E_mean + q_vector[-1]/q_vector[0])]
    else:
        accumulation_equations += [Eq(E_mean, E_mean + q_vector[-1])]
    # Step 2: Divide by the number of statistics samples
    normalise_equations = []
    # Density
    normalise_equations += [Eq(r_m, r_m/nsamples)]
    # Momentum
    for i in range(ndim):
        normalise_equations += [Eq(ru_m[i], ru_m[i]/nsamples)]
    # Energy
    normalise_equations += [Eq(E_mean, E_mean/nsamples)]
    return accumulation_equations, normalise_equations, storage_arrays

def primitive(ndim, q_vector):
    # Create the statistics equations, this shows another way of writing the equations
    # create variables means required
    u_m = symbols('u0:%dmean' % ndim, **{'cls': DataObject})
    niter = symbols("niter", **{'cls':ConstantObject})
    niter.datatype = Int()

    accumulation_equations, normalisation_equations = [], []

    grid_vars = [GridVariable('u%d' % i) for i in range(ndim)]

    a_ins = GridVariable("ains")
    
    p_ins = GridVariable("pins")

    t_ins = GridVariable("tins")

    M_ins = GridVariable("mins")

    mu_ins = GridVariable("muins")

    accumulation_equations += [Eq(grid_vars[i], q_vector[i+1]/q_vector[0]) for i in range(ndim)]

    pprint(accumulation_equations)
    # accumulate_momentum
    for i in range(ndim):
        accumulation_equations += [Eq(u_m[i], u_m[i] + grid_vars[i])]
    # accumulate_momentum
    for i in range(ndim):
        normalisation_equations += [Eq(u_m[i], u_m[i]/niter)]


    for i in range(ndim):
        for j in range(ndim):
            if i >= j:
                mean_var = DataObject("u%du%dmean" % (i, j))
                comp1 = grid_vars[i]
                comp2 = grid_vars[j]
                accumulation_equations += [Eq(mean_var, (comp1 * comp2) + mean_var)]
                normalisation_equations += [Eq(mean_var, mean_var/niter)]

   
    p_mean = DataObject('p_mean')
    pp_mean = DataObject('pp_mean')
    a_mean = DataObject('a_mean')
    T_mean = DataObject('T_mean')
    TT_mean = DataObject('TT_mean')
    M_mean = DataObject('M_mean')
    mu_mean = DataObject('mu_mean')
    
    p_ins = 0.4 * (q_vector[-1] - (0.5*q_vector[0]*(grid_vars[0]**2)+0.5*q_vector[0]*(grid_vars[1]**2)+0.5*q_vector[0]*(grid_vars[2]**2)))

    a_ins = (1.4 * p_ins / q_vector[0])**0.5 

    T_ins = p_ins*1.4*0.0955*0.0955/q_vector[0]

    mu_ins = T_ins**0.7

    M_ins = ((grid_vars[0]**2)+(grid_vars[1]**2)+(grid_vars[2]**2))**0.5 / a_ins

    accumulation_equations += [Eq(p_mean, p_mean + p_ins)]
    normalisation_equations += [Eq(p_mean, p_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(pp_mean, pp_mean + (p_ins*p_ins))]
    normalisation_equations += [Eq(pp_mean, pp_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(a_mean, a_mean + a_ins)]
    normalisation_equations += [Eq(a_mean, a_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(T_mean, T_mean + T_ins)]
    normalisation_equations += [Eq(T_mean, T_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(TT_mean, TT_mean + (T_ins*T_ins))]
    normalisation_equations += [Eq(TT_mean, TT_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(mu_mean, mu_mean + mu_ins)]
    normalisation_equations += [Eq(mu_mean, mu_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(M_mean, M_mean + M_ins)]
    normalisation_equations += [Eq(M_mean, M_mean/niter)]
    pprint(accumulation_equations)

    return accumulation_equations, normalisation_equations

def second_order_moments(ndim, q_vector, nsamples, conservative=True):
    # Matrix of favre reynolds stresses i.e rhou_i*u_j if i >= j, only lower triangular matrix is considered
    def RS(i,j):
        if i >= j:
            return DataObject("rhou%du%dmean" % (i, j))
        else:
            return 0
    stresses_reynolds = Matrix(ndim, ndim, RS)
    accumulation, normalisation, storage_arrays = [], [], []

    # Divide by density
    if conservative:
        rho_inv = GridVariable("rhoinv")
        accumulation += [Eq(rho_inv, 1.0/q_vector[0])]
        for i in range(ndim):
            for j in range(ndim):
                if i >= j:
                    mean_var = DataObject("rhou%du%dmean" % (j,i))
                    storage_arrays += [mean_var]
                    comp1 = q_vector[i+1]
                    comp2 = q_vector[j+1]*rho_inv
                    accumulation += [Eq(mean_var, (comp1 * comp2) + mean_var)]
                    normalisation += [Eq(mean_var, mean_var/nsamples)]
    else:
        for i in range(ndim):
            for j in range(ndim):
                if i >= j:
                    mean_var = DataObject("rhou%du%dmean" % (j,i))
                    storage_arrays += [mean_var]
                    accumulation += [Eq(mean_var, q_vector[0]*(q_vector[i+1] * q_vector[j+1]) + mean_var)]
                    normalisation += [Eq(mean_var, mean_var/nsamples)]
    return accumulation, normalisation, storage_arrays

def favre_averaged_stats(ndim, q_vector, conservative=True):
    # Create the UDF classes to hold the statistics
    accumulation = UserDefinedEquations()
    stat_frequency = symbols("stat_frequency", **{'cls':ConstantObject})
    stat_frequency.datatype = Int()
    accumulation.algorithm_place = InTheSimulation(frequency=stat_frequency)
    accumulation.order = 1e9
    nsamples = symbols("nsamples", **{'cls':ConstantObject})
    nsamples.value = 'niter/stat_frequency'
    nsamples.datatype = Int()
    # Divide at the end
    normalisation = UserDefinedEquations()
    normalisation.algorithm_place = AfterSimulationEnds()
    # Left hand side, names to be written by the HDF5 class
    storage_arrays = []
    # Calculate first moment statistics, one for calculation, one to divide by the number of samples
    first_moms_acc, first_mom_normalise, array_output  = first_order_moments(ndim, q_vector, nsamples, conservative=conservative)
    accumulation.add_equations(first_moms_acc)
    normalisation.add_equations(first_mom_normalise)
    storage_arrays += array_output
    # prim_acc, prim_norm = primitive(ndim, q_vector, conservative=conservative)
    # accumulation.add_equations(prim_acc)
    # normalisation.add_equations(prim_norm)
    # Calculate second moment statistics
    sec_moms_acc, sec_mom_normalise, array_output = second_order_moments(ndim, q_vector, nsamples, conservative=conservative)
    storage_arrays += array_output
    accumulation.add_equations(sec_moms_acc)
    normalisation.add_equations(sec_mom_normalise)
    print("Gathering {} statistics quantities: {}".format(len(storage_arrays), storage_arrays))
    return [accumulation, normalisation], storage_arrays
