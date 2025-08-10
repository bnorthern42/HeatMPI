import matplotlib.pyplot as plt
import numpy as np
import csv

# add your file handling code instead of the next command
b = np.genfromtxt("output.csv", delimiter=",")
#b = np.loadtxt("output.csv", delimiter=',')

#---
plt.imshow(b, cmap='hot', interpolation='nearest')
plt.show()