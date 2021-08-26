#!/usr/bin/env python
# coding: utf-8

# In[2]:


# disc code
import numpy as np
import matplotlib.pyplot as plt

#initialize arrays for the positions and velocities
x = []
y = []
z = []
Vx = []
Vy = []
Vz = []
R =[]
t = 0
energy = []
time = []

#initial conditions
x_val = 0
y_val = 2300 * 3.086e+16
z_val = 0
Vx_val = 238000
Vy_val = 0
Vz_val = 0
h = 1000000000000

# set the constants
mu_e = 1.2*10**11 * 1.989e+30 * 6.67e-11
N = 2*10**9*365.25*24*60*60 #orbit time

# height and length of MW
a = 3.086e+16 * 4500
b = 3.086e+16 * 260

while t < N:
    
    # the MN potentials
    pot = mu_e * (x_val**2 + y_val**2 + (a + np.sqrt(b**2+z_val**2) )**2 )**-1.5
    pot_z = mu_e * (x_val**2 + y_val**2 + (a+np.sqrt(b**2+z_val**2))**2)**-1.5 * (2*z_val + 2*z_val*a*(b**2+z_val**2)**-0.5)
    
    Vx_val = Vx_val - h * pot * x_val    
    x_val = x_val + h * Vx_val
    x.append(x_val)
    Vx.append(Vx_val)

    Vy_val = Vy_val - h * pot * y_val
    y_val = y_val + h * Vy_val
    y.append(y_val)
    Vy.append(Vy_val)
    
    z_val = z_val + h * Vz_val
    Vz_val = Vz_val - h * pot_z
    Vz_val = Vz_val - h * pot_z
    z.append(z_val)
    Vz.append(Vz_val)

    t = t + h
    
plt.plot(x,y)


# In[ ]:





# In[ ]:




