""" Script to generate all of the current OpenSBLI test cases."""
import os, subprocess
# List of the current applications
cases = {
'/wave/'    : 'wave.py',
'/euler_wave/'  : 'euler_wave.py',
'/shu_osher/'   : 'shu_osher.py',
'/Sod_shock_tube/'  : 'Sod_shock_tube.py',
'/Lax_shock_tube/'  : 'Lax_shock_tube.py',
'/LeBlanc/' : 'LeBlanc.py',
'/taylor_green_vortex/' : 'taylor_green_vortex.py',
'/taylor_green_vortex/TGsym/'   : 'TG_IsoT.py',
'/taylor_green_vortex/TGsym/'   : 'TGsym.py',
'/viscous_shock_tube/'  : 'viscous_shock_tube.py',
'/kelvin_helmholtz/'    : 'kelvin_helmholtz.py',
'/inviscid_shock_reflection/'   : 'inviscid_shock.py',
'/katzer_SBLI/' : 'katzer_SBLI.py',
'/channel_flow/laminar_2D/' : 'laminar_channel.py',
'/channel_flow/turbulent_3D/'   : 'turbulent_channel.py',
'/channel_flow/compressible_TCF_Central/'   : 'turbulent_channel.py',
'/channel_flow/compressible_TCF_TENO/'  : 'turbulent_channel.py',
'/channel_flow/adiabatic_isothermal_channel'    : 'iso_adi_channel_heat_sink.py',
#'/Delery_bump/inviscid/'   : #'inviscid_shock_delery_aerofoil_forced.py',
#'/Delery_bump/viscous/'    : #'viscous_shock_delery_aerofoil.py',
'/transitional_SBLI/'   : 'transitional_SBLI.py',
'/cylinder/supersonic_cylinder/'    : 'supersonic_cylinder.py',
'/compressible_taylor_green_vortex/TGV_multi_block/'    : 'compressible_TGV_MB.py',
'/compressible_taylor_green_vortex/'    : 'compressible_TGV.py',
'/aerofoils/single_block/2D/'   : 'CRM_2D.py',
'/aerofoils/single_block/3D/'   : 'CRM_3D.py',
'/aerofoils/multi_block/2D/'    : 'transonic_MB.py',
'/aerofoils/multi_block/3D/'    : 'transonic_MB.py',
'/people/max/gaussian_bump/'    : 'gaussian_bump_3D_turbulent.py',
'/people/pushpender/'   : '2d_ramp_5deg_flat_top_invicid_Ly.py',
'/people/Ali/Mixlay_0D_001_TCNEQ_2T_TvNR_1ev_Clean/'    : 'Mixlay_2D_010.py',
'/people/Ali/Mixlay_2D_000_TNEQ_ML_1ev_Clean/'  : 'Mixlay_2D_000.py',
'/people/Ali/Mixlay_2D_000_TCNEQ_2T_TvNR_Transport_3ev_3rho/'   : 'Mixlay_0D_000.py',
'/people/teja/flatplate/'   : 'flatplate.py',
'/people/teja/transition/'  : 'transition17.py',
'/people/teja/mixtransition/'  : 'mixflat_transition.py',
'/people/teja/mixflat/'    : 'mixflat_N2.py',
'/people/teja/mixcylinder/' : 'mixcylinder.py',
'/vortex_core/2D/'  : 'vortex_core.py',
}

directories = [x for x in cases.keys()]
file_names = [cases[x] for x in directories]

assert len(directories) == len(file_names)
print('\33[4m' + "Found %d OpenSBLI applications." % len(file_names) + '\033[0m')
# Current working directory
owd = os.getcwd()
# Optional diff between the generated codes
check_diff = False
generate = True
compile_test = True
OPS_translator_path = '~/software/OPS/ops_translator/c/ops.py'
if check_diff:
    import difflib
    # Set a directory containing previously generated C codes
    old_code_dir = os.environ['two'] + 'apps/'


with open(os.devnull, 'w') as devnull:

    for fname, directory in zip(file_names, directories):
        print("Generating the %s application." % (directory+fname))
        output_code = subprocess.call(["python3.8 %s" % fname], shell=True, cwd=owd+directory, stdout=devnull)
        if output_code == 0:
            print('\33[92m' + "%s generated successfully." % fname + '\033[0m')
            # Compare the output code to a previously generated one
            if check_diff:
                file1, file2 = old_code_dir + directory + 'opensbli.cpp', owd + directory + 'opensbli.cpp'
                text1, text2 = open(file1).readlines(), open(file2).readlines()
                for line in difflib.unified_diff(text1, text2):
                    print(line)
            if generate:
                output_code = subprocess.call(["python3.8 {} opensbli.cpp".format(OPS_translator_path)], shell=True, cwd=owd+directory, stdout=devnull)
                if output_code == 0:
                    print('\33[92m' + "%s translated successfully." % fname + '\033[0m')
            if compile_test:
                if not os.path.isfile('Makefile'):
                    print('\33[91m' + "No Makefile for case: {}".format(directory) + '\033[0m')
                try:
                    proc = subprocess.Popen(["make -B opensbli_openmp"], shell=True, cwd=owd+directory, stdout=devnull, stderr=subprocess.PIPE)
                    for line in proc.stderr:
                        if "error" in str(line) and "linker command failed" not in str(line):
                            print('\33[91m' + "Compilation error: {}".format(line) + '\033[0m')
                except:
                        print("Compile failed.")
        else:
            print('\33[91m' + "Generation of %s%s has failed." % (fname, directory) + '\033[0m')
