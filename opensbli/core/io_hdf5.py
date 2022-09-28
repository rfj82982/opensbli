"""@brief
   @authors Satya Pramod Jammy
   @contributors David J Lusher
   @details
"""
from opensbli.code_generation.algorithm.common import InTheSimulation, AfterSimulationEnds, BeforeSimulationStarts
from opensbli.core.opensbliobjects import Globalvariable, ConstantObject
from opensbli.core.datatypes import Int
from opensbli.core.kernel import ConstantsToDeclare as CTD
from sympy import flatten, pprint


class opensbliIO(object):
    group_number = 0

    @staticmethod
    def increase_io_group_number():
        opensbliIO.group_number += 1
        return

class iohdf5(opensbliIO):
    def __new__(cls, arrays=None, save_every=None, **kwargs):
        ret = super(iohdf5, cls).__new__(cls)
        ret.order = 0
        ret.block_number = 0
        ret.group_number = cls.group_number
        cls.increase_io_group_number()
        if kwargs:
            ret.kwargs = {}
            for key in kwargs:
                if isinstance(key, str):
                    if isinstance(kwargs[key], str):
                        ret.kwargs[key.lower()] = kwargs[key].lower()
        else:
            # Default IO type is write to hdf5
            ret.kwargs = {'iotype': "write"}
            # Default write placement is the end of the simulation
        # Position of write calls in the output
        if 'position' not in ret.kwargs:
            ret.kwargs['position'] = 'end'
        ret.algorithm_place = []
        # Check if constants should be written to the HDF5 file
        if 'write_constants' in kwargs:
            cls.write_constants = kwargs['write_constants']
        else:
            cls.write_constants = True # by default always write the constants to HDF5 now
        # Constant for file write frequency
        if save_every:
            cls.save_every = ConstantObject('write_output_file', integer=True)
            cls.save_every._value = save_every
            cls.save_every.datatype = Int()
            CTD.add_constant(cls.save_every)
        else:
            cls.save_every = None
        ret.get_algorithm_location()
        ret.arrays = []
        if arrays:
            ret.add_arrays(arrays)
        return ret

    def get_algorithm_location(cls):
        if cls.save_every:
            cls.algorithm_place += [InTheSimulation(cls.save_every)]
        if cls.kwargs['iotype'] == "write":
            if cls.kwargs['position'] == "init":
                cls.algorithm_place = [BeforeSimulationStarts()]
            else:
                cls.algorithm_place += [AfterSimulationEnds()]
        elif cls.kwargs['iotype'] == "read":
            cls.algorithm_place = [BeforeSimulationStarts()]
        else:
            raise ValueError("")
        return

    def add_arrays(cls, arrays):
        cls.arrays += flatten(arrays)
        return

    def check_datasets(cls, block):
        """ Checks if the user has added any datasets to the IO class that are not defined within the simulation."""
        simulation_dsets = [str(ar) for ar in block.block_datasets.keys()]
        io_dsets = [str(ar) for ar in cls.arrays]
        missing_dsets = [x for x in io_dsets if x not in simulation_dsets]
        if len(missing_dsets) > 0:
            raise ValueError("The dataset(s): '%s' added to the HDF5 class are not defined in the simulation code. Please check the HDF5 add_arrays input in the problem script." % str(', '.join(missing_dsets)))
        return

    def set_read_from_hdf5_arrays(cls, block):
        if cls.kwargs['iotype'] == "read":
            if 'filename' in cls.kwargs:
                fname = cls.kwargs['filename']
            else:
                fname = 'data.h5'
            for ar in cls.arrays:
                if str(ar) in block.block_datasets.keys():
                    dset = block.block_datasets[str(ar)]
                    dset.read_from_hdf5 = True
                    dset.input_file_name = fname
                    block.block_datasets[str(ar)] = dset
                else:
                    block.block_datasets[str(ar)] = ar
                    block.block_datasets[str(ar)].read_from_hdf5 = True
                    dset.input_file_name = fname
            # Check if restarting from a previous solution or not
            if 'restart_simulation' in cls.kwargs.keys():
                if cls.kwargs['restart_simulation']:
                    for const in CTD.constants:
                        if str(const) == 'restart':
                            const._value = 1
        return

    def write_latex(cls, latex):
        string = ["HDF5 IO type %s on arrays" % (cls.kwargs['iotype'])]
        string += ["%s" % (d) for d in cls.arrays]
        latex.write_string(' '.join(string))
        return

    @property
    def opsc_code(cls):
        code = []
        if cls.kwargs['iotype'] == "write":
            code += cls.hdf5write_opsc_code()
        elif cls.kwargs['iotype'] == "read":
            code += cls.hdf5read_opsc_code()
        else:
            raise ValueError("Cant classify HDF5io")
        return code

    def set_output_constants(cls, constants):
        cls.constants_to_write += flatten([constants])
        return

    def hdf5write_opsc_code(cls):
        var_name = 'name%s' % cls.block_number
        code = []
        if "name" in cls.kwargs:
            if '.h5' in cls.kwargs["name"]:
                name = cls.kwargs["name"]
            elif '.' in cls.kwargs["name"]:
                raise ValueError("")
            else:
                name = cls.kwargs["name"] + '.h5'
            if cls.dynamic_fname:
                raise ValueError("dynamic fname not allowed ")
            filename = "\"%s\"" % name
        else:
            name = "opensbli_output"
            code += ['char %s[80];' % var_name]
            if cls.dynamic_fname:
                code += ['sprintf(%s, \"%s_%%06d.h5\", %s);' % (var_name, name, cls.control_parameter)]
            else:
                code += ['sprintf(%s, \"%s.h5\");' % (var_name, name)]
            filename = var_name
        dataset_write = []
        for ar in cls.arrays:
            block_name = ar.base.blockname
            dataset_write += ['ops_fetch_dat_hdf5_file(%s, %s);' % (ar, filename)]

        # generate the block name
        code += ['ops_fetch_block_hdf5_file(%s, %s);' % (block_name, filename)] + dataset_write
        # Write constants to the HDF5 output file, once per file (not per block)
        if cls.write_constants and cls.block_number == 0:
            # Generate the OPS API calls
            user_constants = [x for x in CTD.constants if isinstance(x, ConstantObject)]
            user_constants = [x for x in user_constants if not x.rational]
            for c in user_constants: # Write only user input constants, not rational factors and inverses
                code += ['ops_write_const_hdf5(\"%s\", 1, \"%s\", (char*)&%s, %s);' % (c.name, c.datatype.opsc(), c.name, filename)]
            # Constants to always write to HDF5
            code += ['ops_write_const_hdf5(\"iter\", 1, \"int\", (char*)&iter, %s);' % (filename)]
        return code

    def hdf5read_opsc_code(cls):
        """To keep the abstraction going return nothing for HDF5 OPSC code"""
        return []

    @property
    def evaluated_datasets(cls):
        evaluated = set()
        if cls.kwargs['iotype'] == "read":
            evaluated = evaluated.union(set([a.base for a in cls.arrays]))
        return evaluated


class IoGroup():
    pass
