#!/usr/bin/env python
# coding: utf-8

# In[9]:


import numpy as np
import matplotlib.pyplot as plt

m_sun = 2e+30
r = 5.3e+12
G = 6.67e-11
m_1 = 2e+30
m_2 = 4e+30
grav_c = G*m_1*m_2
period = 365*24*60*60

# positions of the two stars
x_1 = 0
y_1 = 0
z_1 = 0
x_2 = 1.5e+11
y_2 = 0
z_2 = 0

# velocities of the two stars
Vx_1 = 0
Vy_1 = 0
Vz_1 = 0
Vx_2 = 0
Vy_2 = 10000
Vz_2 = 0

t = 0
daysec = 24*60*60
h = 0.01*daysec

x1list = []
y1list = []
z1list = []

x2list = []
y2list = []
z2list = []

while t < 1000*daysec:
    
    rx = x_2 - x_1
    ry = y_2 - y_1
    rz = z_2 - z_1
    
    r3 = (rx**2 + ry**2 + rz**2)**1.5
    
    fx = -grav_c*rx/r3
    fy = -grav_c*ry/r3
    fz = -grav_c*rz/r3
    
    Vx_1 = Vx_1 - fx*h/m_1
    Vy_1 = Vy_1 - fy*h/m_1
    Vz_1 = Vz_1 - fz*h/m_1
    
    x_1 = Vx_1*h
    y_1 = Vy_1*h
    z_1 = Vz_1*h
    
    x1list.append(x_1)
    y1list.append(y_1)
    z1list.append(z_1)
    
    Vx_2 = Vx_2 + fx*h/m_2
    Vy_2 = Vy_2 + fy*h/m_2
    Vz_2 = Vz_2 + fz*h/m_2
    
    x_2 = Vx_2*h
    y_2 = Vy_2*h
    z_2 = Vz_2*h
    
    x2list.append(x_2)
    y2list.append(y_2)
    z2list.append(z_2)
    
    t = t + h
    
plt.plot(x1list,y1list)
plt.plot(x2list,y2list)


# In[ ]:




