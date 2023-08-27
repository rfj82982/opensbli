# import numpy 
import numpy as np
# import matplotlib for plotting figures
import matplotlib.pyplot as plt
# import scipy for the PSD calculation
from scipy import signal

#----------------------------------------------------------------------#
# read the probe file, made of float values
probe = np.loadtxt("cylinder_probes.log", dtype='f', delimiter=',', skiprows=1)

#----------------------------------------------------------------------#
# construct the iterations array
iterations = probe[:,0] # iterations
time       = probe[:,1] # iteration x dt_code => nondimensional

for i in range(2,9):

  q1 = probe[:,i]
  
  #----------------------------------------------------------------------#
  # get the sampling frequency
  fs = 1/(time[1]-time[0])
  
  #calculate the PSD for the lift coefficients
  f, PSD = signal.periodogram(q1, fs)
  
  #----------------------------------------------------------------------#
  # plot the drag and the lift coefficients
  plt.subplot(2,7,i-1)
  plt.plot(time,q1)
  plt.xlabel("Time")
  plt.grid(True)
  
  # plot the lift PSD distributions
  plt.subplot(2,7,7+i-1)
  plt.loglog(f,PSD)
  plt.xlabel("Frequency")
  plt.ylabel("PSD")
  plt.grid(True)
  
# show the plot
plt.savefig('./simulation_plots/probe_data.pdf', bbox_inches='tight')
