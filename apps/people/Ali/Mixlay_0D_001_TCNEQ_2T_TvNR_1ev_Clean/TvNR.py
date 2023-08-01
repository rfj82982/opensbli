import re
import os

print("--------------------\n" + "Starting Tv replacement")

####################################
# Open the file
f = open("opensbliblock00_kernels.h",'r+')
filedata = f.read()                 # set "filedata" as input
f.close()


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
f = open("opensbliblock00_kernels.h",'w')
f.write(newdata)
f.close()
###################################
print(Character[0])
print("-----------------------------------")
print(CharacterNewer)
print("Tv Implemented"+ "\n--------------------")
