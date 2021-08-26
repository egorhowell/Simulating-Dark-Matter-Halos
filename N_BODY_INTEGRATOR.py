#!/usr/bin/env python
# coding: utf-8

# In[16]:


import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d

def getAcc(pos, mass, G, softening,rs,v0):

    N = pos.shape[0]
    a = np.zeros((N,3))

    # solar system gravity
    for i in range(N):
        for j in range(N-1):
            dx = pos[j,0] - pos[i,0]
            dy = pos[j,1] - pos[i,1]
            dz = pos[j,2] - pos[i,2]
            inv_r3 = (dx**2 + dy**2 + dz**2 + softening**2)**(-1.5)
            a[i,0] +=  G * (dx * inv_r3) * mass[j]
            a[i,1] +=  G * (dy * inv_r3) * mass[j]
            a[i,2] +=  G * (dz * inv_r3) * mass[j]
     
        #kepler   
        for i in range(N):
            dx = pos[-1,0] - pos[i,0]
            dy = pos[-1,1] - pos[i,1]
            dz = pos[-1,2] - pos[i,2]
            inv_r3 = (dx**2 + dy**2 + dz**2 + softening**2)**(-1.5)
            a[i,0] +=  G * (dx * inv_r3) * mass[-1]
            a[i,1] +=  G * (dy * inv_r3) * mass[-1]
            a[i,2] +=  G * (dz * inv_r3) * mass[-1]
       
        #plummer    
        #for i in range(N):
         #   dx = pos[-1,0] - pos[i,0]
          #  dy = pos[-1,1] - pos[i,1]
           # dz = pos[-1,2] - pos[i,2]
            #inv_r3 = (dx**2 + dy**2 + dz**2 + rs**2 + softening**2)**(-1.5)
            #a[i,0] +=  G * (dx * inv_r3) * mass[-1]
            #a[i,1] +=  G * (dy * inv_r3) * mass[-1]
            #a[i,2] +=  G * (dz * inv_r3) * mass[-1]
            
        #logarithmic
        #for i in range(N):
         #   dx = pos[-1,0] - pos[i,0]
          #  dy = pos[-1,1] - pos[i,1]
           # dz = pos[-1,2] - pos[i,2]
            #inv_r3 = (dx**2 + dy**2 + (dz/0.8)**2 + rs**2 + softening**2)**(-1)
            #a[i,0] +=  (dx * inv_r3 * v0**2)
            #a[i,1] +=  (dy * inv_r3 * v0**2) 
            #a[i,2] +=  (dz * inv_r3 * v0**2)
            
        
    return a

G = 6.67e-11
m_sun = 1.989e+30
mass = np.array([1,1.7e-7,2.44e-6,3e-6,3.2e-7,9.5e-4,2.84e-4,4.34e-5,5.1e-5,0]) * m_sun 
softening = 0.00000000000001

pos = np.array([[0,0,0],
               [0.46e+11,0,0], # merc
               [1.08e+11,0,0], # venus
               [1.49e+11,0,0], # earth  
               [2.06e+11,0,0], # mars
               [7.4e+11,0,0],  # jupiter
               [13.52e+11,0,0],# saturn
               [27.4e+11,0,0], # uranus
               [44.44e+11,0,0],# neptune
               [0,0,0]]) # dark matter

vel = np.array([[0,0,0],
               [0,59000,0], # merc
               [0,35000,0], # venus
               [0,30000,0], # earth
               [0,26500,0], # mars
               [0,13720,0], # jupiter
               [0,10180,0], # saturn
               [0,7110,0],  # uranus
               [0,5500,0],  # neptune
               [0,0,0]]) # dark matter


t = 0
dt = 10000
period  = 10e+7
rs = 3e+11
v0 = 10000


x1 = []
y1 = []
z1 = []
x2 = []
y2 = []
z2 = []
x3 = []
y3 = []
z3 = []
x4 = []
y4 = []
z4 = []
x5 = []
y5 = []
z5 = []
x6 = []
y6 = []
z6 = []
x7 = []
y7 = []
z7 = []
x8 = []
y8 = []
z8 = []

v1 = []
v2 = []
v3 = []
v4 = []
v5 = []
v6 = []
v7 = []
v8 = []

a = getAcc(pos, mass, G, softening, rs,v0)

while t < period:  

    x1.append(pos[1][0])
    y1.append(pos[1][1])
    z1.append(pos[1][2])
    x2.append(pos[2][0])
    y2.append(pos[2][1])
    z2.append(pos[2][2])
    x3.append(pos[3][0])
    y3.append(pos[3][1])
    z3.append(pos[3][2])
    x4.append(pos[4][0])
    y4.append(pos[4][1])
    z4.append(pos[4][2])
    x5.append(pos[5][0])
    y5.append(pos[5][1])
    z5.append(pos[5][2])
    x6.append(pos[6][0])
    y6.append(pos[6][1])
    z6.append(pos[6][2])
    x7.append(pos[7][0])
    y7.append(pos[7][1])
    z7.append(pos[7][2])
    x8.append(pos[8][0])
    y8.append(pos[8][1])
    z8.append(pos[8][2])
    
    v1.append( np.sqrt((vel[1][0]**2) + (vel[1][1]**2)) )
    v2.append( np.sqrt((vel[2][0]**2) + (vel[2][1]**2)) )
    v3.append( np.sqrt((vel[3][0]**2) + (vel[3][1]**2)) )
    v4.append( np.sqrt((vel[4][0]**2) + (vel[4][1]**2)) )
    v5.append( np.sqrt((vel[5][0]**2) + (vel[5][1]**2)) )
    v6.append( np.sqrt((vel[6][0]**2) + (vel[6][1]**2)) )
    v7.append( np.sqrt((vel[7][0]**2) + (vel[7][1]**2)) )
    v8.append( np.sqrt((vel[8][0]**2) + (vel[8][1]**2)) )
    
    
    vel = vel + a * dt/2

    pos = pos + vel * dt
    
    a = getAcc(pos, mass, G, softening, rs,v0)
    
    vel = vel + a * dt/2
       
    t = t + dt
    
x1 = np.array(x1) * 1e-10
y1 = np.array(y1) * 1e-10
z1 = np.array(z1) * 1e-10
x2 = np.array(x2) * 1e-11
y2 = np.array(y2) * 1e-11
x3 = np.array(x3) * 1e-11
y3 = np.array(y3) * 1e-11
x4 = np.array(x4) * 1e-11
y4 = np.array(y4) * 1e-11
    
x5 = np.array(x5) * 1e-11
y5 = np.array(y5) * 1e-11
x6 = np.array(x6) * 1e-11
y6 = np.array(y6) * 1e-11
x7 = np.array(x7) * 1e-11
y7 = np.array(y7) * 1e-11
x8 = np.array(x8) * 1e-11
y8 = np.array(y8) * 1e-11


plt.figure(2)
fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(1, 1, 1)
plt.plot(x1,y1, label = 'Mercury')
#plt.plot(x2,y2, label = 'Venus')
#plt.plot(x3,y3, label = 'Earth')
#plt.plot(x4,y4, label = 'Mars')
#plt.plot(x5,y5, label = 'Jupiter')
#plt.plot(x6,y6, label = 'Saturn')
#plt.plot(x7,y7, label = 'Uranus')
#plt.plot(x8,y8, label = 'Neptune')
plt.legend(loc = 'upper left')
plt.xlabel('X-Position ($10^{11}$m)',fontsize=14)
plt.ylabel('Y-Position ($10^{11}$m)',fontsize=14)
plt.rc('font', family='serif')
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.plot(0,0,'ko',label='Sun and Halo')
plt.legend(loc="upper left")
plt.savefig('merc.png',dpi=300, bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection='3d')
ax.plot3D(x1, y1, z1)
ax.set_xlabel('X-Position ($10^{10}$m)', fontsize=20,labelpad=15)
ax.set_ylabel('Y-Position ($10^{10}$m)', fontsize=20,labelpad=15)
ax.set_zlabel('Z-Position ($10^{10}$m)', fontsize=20,labelpad=15)
ax.xaxis.set_tick_params(labelsize=20)
ax.yaxis.set_tick_params(labelsize=20)
ax.zaxis.set_tick_params(labelsize=20)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

# Now set color to white (or whatever is "invisible")
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')
plt.savefig('updare.png')



mean_vel = np.array([np.mean(v1),np.mean(v2), np.mean(v3), np.mean(v4), np.mean(v5), np.mean(v6), np.mean(v7), np.mean(v8)])*1e-3
radius = np.array([0.3871, 0.7233, 1, 1.5273, 5.2, 9.5388, 19.19, 30])

plt.plot(radius,mean_vel,color='black')
plt.xlabel('Semi-Major Axis (AU)',fontsize=14)
plt.ylabel('Velocity (km/s)', fontsize=14)
plt.rc('font', family='serif')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('nodmvel.png',dpi=300, bbox_inches='tight')


# In[ ]:





# In[ ]:




