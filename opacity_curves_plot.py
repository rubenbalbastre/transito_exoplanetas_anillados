# This script plots the opacity profiles of the different models


# IMPORT LIBRARYS needed to compute the models

import numpy as np
import matplotlib.pyplot as plt
import opacity_models as om 

rmin = 0.01
rmax = 24
o_lim = 0.005
r = np.linspace(rmin,rmax,10000)
fig = plt.figure(figsize=(6,6))

# CASO 1.

# 2 rings. Powlaw model
rad_central = 0.8618565116529561
gamma_decay = 1.5945382368765897
print("Extensión anillo 4.2.1",rad_central*pow(o_lim, -1./gamma_decay))
plt.plot(r,om.opacity_powlaw(r,rad_central,gamma_decay),label="4.2.1")

# 2 rings. Powlawsum model
radius_central_1 = 0.8792040877776828
radius_central_2 =  0.7291666067340421
gamma_decay_1 =  1.561093780448864
radius_limit =  1.000268335095424
plt.plot(r,om.opacity_powlawsum(r,radius_central_1,radius_central_2,gamma_decay_1,radius_limit),label="4.2.2")

# 3 rings. Powlaw and powlaw
radius_central_1 =  0.9182743667960531 
radius_limit =  5.980580742532575 
gamma_decay_1 =  1.6181189941921197 
gamma_decay_2 =  1.4215463854327959
plt.plot(r,om.opacity_powlaw_powlaw(r,radius_central_1,radius_limit,gamma_decay_1,gamma_decay_2),label="4.2.3")



# 2 rings. Powlaw model
rad_central = 0.8502320238645875
gamma_decay = 1.5929497977878695
print("Extensión anillo 4.2.1 (separado)",rad_central*pow(o_lim, -1./gamma_decay))
plt.plot(r,om.opacity_powlaw(r,rad_central,gamma_decay),label="4.2.1 (separado)")

# 2 rings. Powlawsum model
radius_central_1 = 0.8788522056488338
radius_central_2 =  0.7289846195667176
gamma_decay_1 =  1.5614661403425687
radius_limit =  0.9999762479306954
plt.plot(r,om.opacity_powlawsum(r,radius_central_1,radius_central_2,gamma_decay_1,radius_limit),label="4.2.2 (separado)")

# 3 rings. powlaw & powlaw
radius_central_1 =  0.9150883679005799 
radius_limit =  6.000910219136926 
gamma_decay_1 =  1.619195972467201
gamma_decay_2 =  1.4272263549623854
plt.plot(r,om.opacity_powlaw_powlaw(r,radius_central_1,radius_limit,gamma_decay_1,gamma_decay_2),label="4.2.3 (separado)")



#plt.plot(r,om.opacity_powlaw_powlaw(r,radius_central_1,radius_limit,gamma_decay_1,gamma_decay_2),label="3 anillos. Ley de potencias y ley de potencias")
plt.legend()
plt.xlim([0.01,5])
plt.ylim([0.05,1.05])
plt.xlabel(r"r $[R_{\star}]$")
plt.ylabel("Opacidad")
plt.show()
fig.savefig("/Users/rubenbalbastre/Desktop/Física/CUARTO AÑO/2do cuatrimestre/TFG/Presentación/perfiles_opacidad_ampl.pdf")

