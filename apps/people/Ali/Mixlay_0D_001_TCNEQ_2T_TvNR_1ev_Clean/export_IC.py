import re
import os

print("--------------------\n" + "Starting replacement")

####################################
# Open the file
f = open("opensbli.cpp",'r+')
filedata = f.read()                 # set "filedata" as input
f.close()

# -------------------------------------
# Find the specified number of iteration:

# Find the name of python file where the number of iteration is read from
for file in os.listdir('.'):
    if file.startswith("Mixlay_"):
        fname = file
# -------------------------------------

L = open(fname,'r+')
filedata2 = L.read()                 # set "filedata" as input
L.close()

for item in filedata2.split("\n"):
    if "save_every=" in item:
        Character= item.strip()
# Strip the letters from the line + strip the first two digits
iter = re.sub('\D', '', Character)
iter = iter[2:len(iter)]

print('Print Every  :' + iter)

# Replace the line in Cpp file
newdata = filedata.replace(         # Replaces first arguemnt with second
    '(fmod(iter + 1,'+iter+') == 0)', '(fmod(iter + 1,'+iter+') == 0 || (iter + 1) == 1)')
    # "|| (iter + 1) == 1)" adds the argument of when iteration number is 1

# Write out the file
f = open("opensbli.cpp",'w')
f.write(newdata)
f.close()
###################################

print("Replaced"+ "\n--------------------")
