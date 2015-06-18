


# Numpy, Matplotplib
import numpy as np
import matplotlib.pyplot as plt



a = 3.4
b = -1.5
x = np.arange(-10,10,0.5)
y = x**2 + a*x + b
##plt.plot(x,y)
##plt.show()

nx = x**2
nnx = np.sqrt(nx)

z = nx + a*np.sqrt(nx) + b

z1 = np.interp(nnx,x,y)

plt.plot(nx,z)
plt.plot(nx,z1)
plt.show()

