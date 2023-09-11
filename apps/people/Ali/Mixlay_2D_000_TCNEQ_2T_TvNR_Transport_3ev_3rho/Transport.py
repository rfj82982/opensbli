# Coded by: Ali Musawi 01/08/2023
#       + Last edited: 13/08/2023

import re
import os

print("--------------------\n" + "Starting Transport Implementation")

####################################
# find the chosen viscosity and thermal conductivity model and the TC3 or TC1
# --------------------------------------------------------
# Open the file
fcpp        = open("Mixlay_0D_000.py",'r+')
filedatapy  = fcpp.read()
fcpp.close()
for line in filedatapy.split('\n'):
    if "Transport_model =" in line:
        dummy = line.split('\'')
        chosen_mu = dummy[1]
        chosen_TC = dummy[3]
        TC3 = dummy[5]
    if "Diffusion		= " in line:
        dummy = line.split('\'')
        chosen_Diff = dummy[1]
# --------------------------------------------------------
# Open the kernel file
f = open("opensbliblock00_kernels.h",'r+')
filedata  = f.read()
f.close()
Character = []
# Find where kappa and mu is in the kernel file
if TC3 == 'N':  # Single thermal conductivity
    for item in filedata.split("{"):
        for item2 in item.split("}"):
            if "kappatr_B0(0,0) =" in item2:
                Character_TC_tr = item2
            if "kappavib_B0(0,0) =" in item2:
                Character_TC_vib = item2
            elif "mu_B0(0,0) =" in item2:
                Character_mu = item2
elif TC3 == 'Y':       # multi thermal conductivity
    for item in filedata.split("{"):
        for item2 in item.split("}"):
            if "kappatr_B0(0,0) =" in item2:
                Character_TC_tr = item2
            elif "kappavibO2_B0(0,0) =" in item2:
                Character_TC_vibO2 = item2
            elif "kappavibN2_B0(0,0) =" in item2:
                Character_TC_vibN2 = item2
            elif "kappavibNO_B0(0,0) =" in item2:
                Character_TC_vibNO = item2
            elif "mu_B0(0,0) =" in item2:
                Character_mu = item2



# Allocate the equation for each of the chosen method
# -------------------------------------
if TC3 == 'N':
    temp3   = open("Transport_mu_1tc.txt",'r')
    mu_code = temp3.read()
    temp3.close()
    # Search the text file including the models and chose the correct one for implementation
    for item_2 in mu_code.split("// --break--"):
        if chosen_mu in item_2 and "viscosity" in item_2:
            mu = item_2.strip()
        elif chosen_TC in item_2 and "thermal conductivity tr" in item_2:
            TC_tr = item_2.strip()
        elif chosen_TC in item_2 and "thermal conductivity vib" in item_2:
            TC_vib = item_2.strip()
elif TC3 == 'Y':
    temp3   = open("Transport_mu_3tc.txt",'r')
    mu_code = temp3.read()
    temp3.close()
    # Search the text file including the models and chose the correct one for implementation
    for item_2 in mu_code.split("// --break--"):
        if chosen_mu in item_2 and "viscosity" in item_2:
            mu = item_2.strip()
        elif chosen_TC in item_2 and "thermal conductivity tr" in item_2:
            TC_tr = item_2.strip()
        elif chosen_TC in item_2 and "thermal conductivity vibO2" in item_2:
            TC_vibO2 = item_2.strip()
        elif chosen_TC in item_2 and "thermal conductivity vibN2" in item_2:
            TC_vibN2 = item_2.strip()
        elif chosen_TC in item_2 and "thermal conductivity vibNO" in item_2:
            TC_vibNO = item_2.strip()






temp4   = open("Transport_Diffusion.txt",'r')
diffusion_code = temp4.read()
temp4.close()
# Search the text file including the models and chose the correct one for implementation
if chosen_Diff == 'Lee85':
    for item in filedata.split("{"):
        for item2 in item.split("}"):
            if "DO_B0(0,0) =" in item2:
                Character_diffO = item2
            elif "DO2_B0(0,0) =" in item2:
                Character_diffO2 = item2
            elif "DN_B0(0,0) =" in item2:
                Character_diffN = item2
            elif "DN2_B0(0,0) =" in item2:
                Character_diffN2 = item2
            elif "DNO_B0(0,0) =" in item2:
                Character_diffNO = item2
    # -----------------------------------------------
    for item_2 in diffusion_code.split("// --break--"):
        if 'diffusionO-' in item_2:
            diffO = item_2.strip()
        elif 'diffusionO2' in item_2:
            diffO2 = item_2.strip()
        elif 'diffusionN-' in item_2:
            diffN = item_2.strip()
        elif 'diffusionN2' in item_2:
            diffN2 = item_2.strip()
        elif 'diffusionNO' in item_2:
            diffNO = item_2.strip()
    # -----------------------------------------------
    filedata = filedata.replace(Character_diffO, diffO)
    filedata = filedata.replace(Character_diffO2, diffO2)
    filedata = filedata.replace(Character_diffN, diffN)
    filedata = filedata.replace(Character_diffN2, diffN2)
    filedata = filedata.replace(Character_diffNO, diffNO)
    print('Diffusion changed')




# make a backup of the kernel file
filedata_backup = filedata

# Implementing the models in the kernel code
# -------------------------------------
# Replace the viscosity
if TC3 == 'N':
    if chosen_mu == 'none':
        print('Viscosity not changed')
    elif not chosen_mu == 'none':
        filedata = filedata.replace(Character_mu, mu)
        print('Viscosity changed')
    # Replace the thermal conductivity
    if chosen_TC == 'none':
        print('Thermal Conductivity not changed')
    elif not chosen_TC == 'none' :
        filedata = filedata.replace(Character_TC_tr, TC_tr)
        filedata = filedata.replace(Character_TC_vib, TC_vib)
        print('Thermal Conductivity changed')
elif TC3 == 'Y':
    if chosen_mu == 'none':
        print('Viscosity not changed')
    elif not chosen_mu == 'none':
        filedata = filedata.replace(Character_mu, mu)
        print('Viscosity changed')
    # Replace the thermal conductivity
    if chosen_TC == 'none':
        print('Thermal Conductivity not changed')
    elif not chosen_TC == 'none' :
        filedata = filedata.replace(Character_TC_tr, TC_tr)
        filedata = filedata.replace(Character_TC_vibO2, TC_vibO2)
        filedata = filedata.replace(Character_TC_vibN2, TC_vibN2)
        filedata = filedata.replace(Character_TC_vibNO, TC_vibNO)
        print('3x Thermal Conductivity changed')

print('----------------------------------------- \n-------------------------------- ')

# Implementing the constants in the cpp code
# -------------------------------------
# Reading the cpp file
temp4       = open("opensbli.cpp",'r+')
filedatacpp = temp4.read()
filedatacpp_backup = filedatacpp
temp4.close()

# Extract the constants from the .txt
constants = mu_code.split("// --break-- CONSTANTS")
# add the constants to the cpp file if they are both NOT "none"
if not (chosen_TC == 'none' and chosen_mu == 'none'):
    for linecpp in filedatacpp.split('\n'):
        if 'double' in linecpp:                                     # finds where the first "souble" is
            linecpp2 = linecpp  + constants[1]
            filedatacpp = filedatacpp.replace(linecpp, linecpp2)
            break
# declare the constants for OPS
if not (chosen_TC == 'none' and chosen_mu == 'none'):
    for linecpp in filedatacpp.split('\n'):
        if 'ops_decl_const(' in linecpp:                            # finds where the first "ops_decl_const" is
            linecpp2 = linecpp  + constants[2]
            filedatacpp = filedatacpp.replace(linecpp, linecpp2)
            break

# Write out the newly created files: .cpp and kernel
# -------------------------------------
# Write out the file
f = open("opensbliblock00_kernels.h",'w')
f.write(filedata)
f.close()
f = open("opensbli.cpp",'w')
f.write(filedatacpp)
f.close()

