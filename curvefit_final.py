
#################################################################################
####### THIS SCRIPT MINIMIZES THE CHI^2 TEST OF A MODEL COMPARED WITH THE DATA ##
#################################################################################

from scipy.optimize import curve_fit
from pandas import read_table
import modelos_planet_ring
from libtransits import generate_randpoints_in_circle_limbdark
from time import time
import numpy as np 
import numpy.random as rd
rd.seed(314) # to repeat the same calculus

start_time = time()


# DATA 
data_of_interest = read_table('/Users/rubenbalbastre/Desktop/Física/CUARTO AÑO/2do cuatrimestre/TFG/Curvas de luz/data_of_interest.txt',delim_whitespace=True).to_numpy()
Nrp_limb = 200000 # Number of points
rand_points_in_circle = generate_randpoints_in_circle_limbdark(Nrp_limb) 
t_steps=data_of_interest[:,0] 
x_data = data_of_interest[:,0] # kepler date
y_data = data_of_interest[:,1] # relative flux

# MODELS & INITIAL CONDITIONS

# 2 anillos. Ley de potencias

def f(t_steps,T0, velocity, yp, radius_planet,radius_empty_ellipse, rad_central, gamma_decay, theta_angle, e_ellipse,e_planet, initial_opacity):
    
    alpha = [T0, velocity, yp, radius_planet,radius_empty_ellipse, rad_central, gamma_decay, theta_angle, e_ellipse,initial_opacity,e_planet]
    relflux =  modelos_planet_ring.relflux_model_allt_planet_ellipse_powerlaw(t_steps,rand_points_in_circle,alpha)
    return relflux

T0 =  792.7353157372115
velocity = 5.755936043130253
yp =  -0.6641242337473513
radius_planet =  0.2892221874081098
planet_separation =  9.907112659448818e-11
rad_central =  0.8618565116529561
gamma_decay =  1.5945382368765897
theta_angle =  0.05882554462900485
e_ellipse =  0.9979368067666516
e_planet =  0.2522929758800985
initial_opacity =  0.9993418216478868

initial_conditions = [T0,velocity,yp,radius_planet,planet_separation,rad_central,gamma_decay,theta_angle,e_ellipse,e_planet,initial_opacity]
lista = ["T0","velocity","yp","radius_planet","planet_separation","rad_central","gamma_decay","theta_angle","e_ellipse","e_planet","initial_opacity"]
n = len(initial_conditions)

bounds = ((792,1,-1,0.1,0,0,0,0,0,0,0.9),(793,10,0,0.5,0.5,5,2,1,1,0.5,1))


# 2 anillo. Suma de 2 leyes de potencias
"""
def f(t_steps,T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity,e_planet):
    
    alpha_powlawsum = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity,e_planet]
    relflux = modelos_planet_ring.relflux_model_planet_powlaw_sum(t_steps,rand_points_in_circle,alpha_powlawsum)

    return relflux

T0 =  792.7358476160401 
velocity =  6.058794124422742 
yp =  -0.6235537413816092 
radius_planet =  0.2759797061204283 
planet_separation =  1e-10  
radius_central_1 =  0.8792193098593314 
radius_central_2 =  0.7291647151869433 
gamma_decay_1 =  1.5611697192298082 
radius_limit =  0.9997928162832527 
theta_angle =  0.05689651407788253 
e_ellipse =  0.9983864934019999  
initial_opacity =  0.9924713789696227 
e_planet =  0.21

lista = ["T0", "velocity", "yp","radius_planet","planet_separation","radius_central_1","radius_central_2","gamma_decay_1", "radius_limit","theta_angle",  "e_ellipse", "initial_opacity","e_planet"]
initial_conditions = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity,e_planet]
n = len(initial_conditions)
bounds = ((792,1,-1,0,0,0,0,0,0,0,0.9,0.9,0),(793,10,-0.1,0.5,0.5,5,0.99,2,5,0.5,1,1,0.5))
"""
# 3 anillos. Ley de potencias y ley de potencias.
"""
def f(t_steps,T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity,e_planet):
    alpha_powlaw_powlaw = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity,e_planet]
    relflux = modelos_planet_ring.relflux_powlaw_powlaw(t_steps,rand_points_in_circle,alpha_powlaw_powlaw)
    return relflux

T0 =  792.736801060297
velocity =   5.770398790535912
yp =   -0.556889563007003  
radius_planet =  0.276
planet_separation = 0
radius_central_1 =  0.9190885385025974
radius_limit =  6.00026831997475
gamma_decay_1 =  1.6182487170629745
gamma_decay_2 =  1.4206411918789255 
theta_angle = 0.07067671931673186
e_ellipse =  0.9988367203091169
initial_opacity =  0.99993933102073
e_planet =  0.25

initial_conditions =   [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity,e_planet]
lista = ["T0", "velocity", "yp","radius_planet", "planet_separation","radius_central_1","radius_limit", "gamma_decay_1","gamma_decay_2", "theta_angle", "e_ellipse", "initial_opacity","e_planet"]
n = len(initial_conditions)
bounds = ((792.7,1,-1,0,0,0.2,1,1,1,0,0.9,0.9,0),(793,10,0,1,1,5,15,3,3,0.3,1,1,0.5))
"""
"""
# MINIMIZATION
optimal_values, cov = curve_fit(f,x_data,y_data,initial_conditions,bounds=bounds,diff_step = 1e-2)
print("\n Datos de la minimización de chi2 \n")
errores = np.sqrt(np.diag(cov))


# RESULTS
k = -1      
for n in initial_conditions:
    k = k+1
    name = lista[k]
    print(name,"= ",optimal_values[k]," Error = ", errores[k],"\n")
"""
# Chi^2
expected_value = f(x_data,*initial_conditions)
chi2 = (y_data-expected_value)**2 / expected_value
print("Initial contidions Chi^2: ",sum(chi2))
"""
expected_value = f(x_data,*optimal_values)
chi2 = (y_data-expected_value)**2 / expected_value
print("Optimal values Chi^2: ",sum(chi2))
"""
# Coeficiente de determinación R^2
# Suma de los cuadrados de los residuos
ss_res = np.sum( (y_data - expected_value)**2  )

# Suma total de cuadrados
ss_tot = np.sum( (y_data - np.mean(y_data) )**2  )

R2 = 1 - (ss_res / ss_tot)
print('R2 = {:10.8f}'.format(R2) )



# show the running time
print("\n Code runned for : ",(time()-start_time)," s \n")
"""
# PLOTTING
import matplotlib.pyplot as plt
plt.plot(x_data,y_data,'.',label="Data")
plt.plot(x_data,f(x_data,*initial_conditions),label="Model")
plt.legend()
plt.show()
"""


#SPECIAL PLOT

std = np.empty(len(x_data)) # standard deviation between data and our model
for step in np.arange(len(x_data)):
    std[step] = y_data[step] - expected_value[step]

import matplotlib.pyplot as plt 
fig2,ax=plt.subplots(2,1,figsize=(6,6),sharex=True,gridspec_kw={'height_ratios': [10, 2]})
ax[0].plot(x_data,y_data,'.',label="Datos")
ax[0].plot(x_data,expected_value,label="Mejor ajuste del modelo")
ax[0].legend()
ax[0].set_ylabel("Flujo relativo")
ax[1].plot(x_data,std[:])
ax[1].set_ylim([-0.008,0.008])
ax[1].set_ylabel("Residuos")
ax[1].set_xlabel("Tiempo [Fecha Kepler]")
fig2.tight_layout()
plt.show()
fig2.savefig("/Users/rubenbalbastre/Desktop/Física/CUARTO AÑO/2do cuatrimestre/TFG/Presentación/mejor_ajuste_protoplaneta.pdf")

