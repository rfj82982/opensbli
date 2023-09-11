# Coded by: Ali Musawi - 12/08/2023
#   + Last edited: 13/08/2023
#   + The purpose of this is to change inf values to zero for 3rho simulations
#       and it's done by adding a single if statement

# -----------------------------------
# Reading the cpp file
temp4       = open("opensbliblock00_kernels.h",'r+')
filedatacpp = temp4.read()
filedatacpp_backup = filedatacpp
temp4.close()


# -----------------------------------
# Define the line to find and the if statement to add
evNO = "   evNO_B0(0,0) ="
if_state_evNO = "\n   if (evNO_B0(0,0) != evNO_B0(0,0) or evNO_B0(0,0) == 1.0/0.0 or evNO_B0(0,0) == -1.0/0.0){\n" \
            +"       evNO_B0(0,0) = 0.0;\n" \
            +"   }"
evN2 = "   evN2_B0(0,0) ="
if_state_evN2 = "\n   if (evN2_B0(0,0) != evN2_B0(0,0) or evN2_B0(0,0) == 1.0/0.0 or evN2_B0(0,0) == -1.0/0.0){\n" \
            +"       evN2_B0(0,0) = 0.0;\n" \
            +"   }"
evO2 = "   evO2_B0(0,0) ="
if_state_evO2 = "\n   if (evO2_B0(0,0) != evO2_B0(0,0) or evO2_B0(0,0) == 1.0/0.0 or evO2_B0(0,0) == -1.0/0.0){\n" \
            +"       evO2_B0(0,0) = 0.0;\n" \
            +"   }"

# -----------------------------------
# Add the if statements
for item_2 in filedatacpp.split('\n'):
    if evNO in item_2:
        filedatacpp = filedatacpp.replace(item_2, item_2+if_state_evNO)
    elif evN2 in item_2:
        filedatacpp = filedatacpp.replace(item_2, item_2+if_state_evN2)
    elif evO2 in item_2:
        filedatacpp = filedatacpp.replace(item_2, item_2+if_state_evO2)


# See if replacement is enabled in the file
# --------------------------------------------------------
# Open the file
fcpp        = open("Mixlay_0D_000.py",'r+')
filedatapy  = fcpp.read()
fcpp.close()
for line in filedatapy.split('\n'):
    if "inf_to_zero" in line and "\'Y\'" in line:
        answer = 'Y'
    elif "inf_to_zero" in line and "\'N\'" in line:
        answer = 'N'
# -----------------------------------

# write the file out
if answer == 'Y':
    f = open("opensbliblock00_kernels.h",'w')
    f.write(filedatacpp)
    f.close()
    print("inf_to_zero_3rho APPLIED")
else:
    print("inf_to_zero_3rho NOT APPLIED")

print("------------------------")
print("------------------------")