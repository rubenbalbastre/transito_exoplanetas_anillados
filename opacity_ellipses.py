# This script plots the opacity profiles of the different models


# IMPORT LIBRARYS needed to compute the models

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# 2 rings. Powlaw model
rad_central = 0.8588570681698539
gamma_decay = 1.5895456911756274
r_powlaw = np.linspace(rad_central,6,1000)
e_planet = 0.3
radius_planet = 0.29
planet_separation = 0
r0 = radius_planet+planet_separation
e_rings = 0.997913661355117

theta_angle = 0.05706727918143067 * 180/np.pi

# 2 rings. Powlawsum model
radius_central_1 = 0.8788602780343863
radius_central_2 =  0.7289851845216327  
gamma_decay_1 =  1.5614611129237697  
radius_limit =  0.999935719291165 
r_powlawsum = np.linspace(radius_limit,6,1000)

# 3 rings. Powlaw and powlaw

fig,ax = plt.subplots(figsize=(6,6))

planet = mpatches.Ellipse((0,0),r0,r0*np.sqrt(1-e_planet**2.),theta_angle)
ring1 = mpatches.Ellipse((0,0),rad_central,rad_central*np.sqrt(1-e_rings**2.),theta_angle,color = "k",alpha=0.9,)
ring2 = mpatches.Ellipse((0,0),2,2*np.sqrt(1-e_rings**2.),theta_angle,color = "k",alpha=0.5)

ax.add_patch(ring1)
ax.add_patch(ring2)
ax.add_patch(planet)
ax.set_xlabel(r"x $[R_{\odot}]$")
ax.set_ylabel(r"y $[R_{\odot}]$")
ax.set_axis_off()
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
plt.show()
"""
import opacity_models as om

x = np.linspace(-1, 1, 21)
r=15
y1 = np.sqrt((r**2 - x**2)*(1-e_rings**2))
y2 = -y1
X1,Y1 = np.meshgrid(x,y1)
X2,Y2 = np.meshgrid(x,y2)
R1 = np.sqrt(x**2+y1**2)
R2 = np.sqrt(x**2+y2**2)
Z1 = om.opacity_powlaw(R1,rad_central,gamma_decay_1)

Z2 = om.opacity_powlaw(R2,rad_central,gamma_decay_1)

plt.contourf(X1,Y1,Z1)
plt.show()
"""