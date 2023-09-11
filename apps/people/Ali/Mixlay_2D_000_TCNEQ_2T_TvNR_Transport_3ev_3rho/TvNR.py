# Coded by: Ali Musawi
#   + Last edited: 24/08/2023

import re
import os
import sys

print("--------------------\n" + "Starting Tv replacement")

####################################
# Open the file
f = open("opensbliblock00_kernels.h",'r+')
filedata = f.read()                 # set "filedata" as input
f.close()


# See if replacement is enabled in the file
# --------------------------------------------------------
# Open the file
fcpp        = open("Mixlay_0D_000.py",'r+')
filedatapy  = fcpp.read()
fcpp.close()
for line in filedatapy.split('\n'):
    if "TvNR_py" in line and "\'Y\'" in line:
        answer = 'Y'
    elif "TvNR_py" in line and "\'N\'" in line:
        answer = 'N'
        print("Tv NOT Implemented" + "\n--------------------")
        sys.exit()

Character = []
for item in filedata.split("{"):
    for item2 in item.split("}"):
        if "Tv_B0(0,0) = Tvref_B0(0,0)" in item2:
            Character.append(item2.strip())


# Replace Tvref by Tv
CharacterNew = Character[0].replace('Tvref_B0(0,0)', 'Tv_B0(0,0)')

# CharacterNewer = "unsigned short ppp = 1; \n Tv_B0(0,0) = Tvref_B0(0,0) + 100.0 ;\n while(ppp <= 2) { \n"+ CharacterNew +"    \n    ppp++; \n }"
CharacterNewer = "    double oldTv_B0,error, errtolTvNR = 1.0e-8; \n" \
                 "    error = 1.0;\n" \
                 "    Tv_B0(0,0) = Tvref_B0(0,0) + 200 ;\n" \
                 "    while(error > errtolTvNR ) { \n" \
                 "      oldTv_B0 = Tv_B0(0,0); \n"\
                 +"      " + CharacterNew +"\n" \
                 "      error = fabs(Tv_B0(0,0)-oldTv_B0); \n}"


# Replace
newdata = filedata.replace(         # Replaces first arguemnt with second
    Character[0], CharacterNewer )

# Write out the file
if answer == 'Y':
    f = open("opensbliblock00_kernels.h",'w')
    f.write(newdata)
    f.close()
    print("Tv Implemented" + "\n--------------------")


###################################
# print(Character[0])
# print(CharacterNewer)
