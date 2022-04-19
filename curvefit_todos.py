
#################################################################################
####### THIS SCRIPT MINIMIZES THE CHI^2 TEST OF A MODEL COMPARED WITH THE DATA ##
#################################################################################

# SELECT A MODEL AND ITS CONDITIONS TO EXECUTE

from scipy.optimize import curve_fit
from pandas import read_table
import modelos_planet_ring
from libtransits import generate_randpoints_in_circle_limbdark
from time import time
import numpy as np 

start_time = time()


# DATA
data_of_interest = read_table('/Users/rubenbalbastre/Desktop/Física/CUARTO AÑO/2do cuatrimestre/TFG/Curvas de luz/data_of_interest.txt',delim_whitespace=True).to_numpy()
Nrp_limb = 200000 # Number of points
rand_points_in_circle = generate_randpoints_in_circle_limbdark(Nrp_limb)
t_steps=data_of_interest[:,0] 
x_data = data_of_interest[:,0]
y_data = data_of_interest[:,1]

# INITIAL CONDITIONS

# VALUES OF THE MODEL PARAMETERS. THEY ARE ABLE TO BE CHANGED WHEN TESTING DIFFERENT MODELS
"""
T0 =  792.72  
velocity = 2.5
yp = -0.7
radius_planet = 0.387
planet_separation = 0.09
radius1 = 3.5
radius2 = 3.8
radius3 = 7
radius4 = 13
val1 = 0.7
val2 = 0
val3 = 0.13
val4 = 0.01
theta_angle = 2.8*np.pi /180
e_ellipse = 0.9999
e_planet = 0
"""

T0 =  792.7524126791726
velocity =  6.310054468250201
yp =  -0.6729189254203921
planet_separation =  0.049996378212313994 
radius_central =  0.9476180688040844
radius_limit =  0.9456470784547976
gamma_decay =  1.5796377789108569 
theta_angle =  0.057680977281105536
e_ellipse =  0.9981398755339095 
initial_opacity =  0.9997279400633751 
radius_planet =  0.2766418729581342
e_planet =  0.09999860675476276
# MODELS

# 2 anillos. Ley de potencias
"""
def f(t_steps,T0, velocity, yp, radius_planet,planet_separation, rad_central, gamma_decay, theta_angle, e_ellipse, initial_opacity, e_planet):
    
    alpha = [T0, velocity, yp, radius_planet,planet_separation, rad_central, gamma_decay, theta_angle, e_ellipse, initial_opacity,e_planet]
    relflux =  modelos_planet_ring.relflux_model_allt_planet_ellipse_powerlaw(t_steps,rand_points_in_circle,alpha)
    return relflux

initial_conditions = [T0, velocity, yp, radius_planet,planet_separation, rad_central, gamma_decay, theta_angle, e_ellipse, initial_opacity,e_planet]
lista = ["T0", "velocity", "yp", "radius_planet","planet_separation", "rad_central", "gamma_decay", "theta_angle","e_ellipse", "initial_opacity","e_planet"]

"""
# 3 anillos. Ley de potencias y ley de potencias * exp
"""
def f(t_steps,T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,gamma_decay_2, theta_angle,  e_ellipse, initial_opacity, e_planet):
    
    alpha = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,gamma_decay_2, theta_angle,  e_ellipse, initial_opacity, e_planet]
    relflux =  modelos_planet_ring.relflux_model_planet_powlaw_and_powlawexp(t_steps,rand_points_in_circle,alpha)
    return relflux

initial_conditions=[T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,gamma_decay_2, theta_angle,  e_ellipse, initial_opacity, e_planet]
lista = ["T0", "velocity", "yp","radius_planet","planet_separation","radius_central_1","radius_central_2","gamma_decay_1","gamma_decay_2", "theta_angle", "e_ellipse", "initial_opacity", "e_planet"]

"""
# 2 anillos. Ley de potencias por exponencial

def f(t_steps,T0, velocity, yp,planet_separation,radius_central, radius_limit,gamma_decay, theta_angle,  e_ellipse, initial_opacity,radius_planet,e_planet):
    
    alpha_powlawexp = [T0, velocity, yp,planet_separation,radius_central, radius_limit,gamma_decay, theta_angle,  e_ellipse, initial_opacity,radius_planet,e_planet]
    relflux = modelos_planet_ring.rel_flux_model_allt_planet_powerlawexp(t_steps,rand_points_in_circle,alpha_powlawexp)
    return relflux

initial_conditions=[T0, velocity, yp,planet_separation,radius_central, radius_limit,gamma_decay, theta_angle,  e_ellipse, initial_opacity,radius_planet,e_planet]
lista = ["T0","velocity", "yp","planet_separation","radius_central", "radius_limit","gamma_decay", "theta_angle", "e_ellipse", "initial_opacity","radius_planet","e_planet"]
bounds =((792,1,-1,0,0,0,1,0,0.9,0.2,0,0),(793,10,0,1,10,10,3,0.3,1,1,1,0.8))

# 2 anillo. Suma de 2 leyes de potencias
"""
def f(t_steps,T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity,e_planet):
    
    alpha_powlawsum = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity,e_planet]
    relflux = modelos_planet_ring.relflux_model_planet_powlaw_sum(t_steps,rand_points_in_circle,alpha_powlawsum)

    return relflux

initial_conditions=[T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity,e_planet]
lista = ["T0", "velocity", "yp","radius_planet","planet_separation","radius_central_1","radius_central_2","gamma_decay_1", "radius_limit","theta_angle",  "e_ellipse", "initial_opacity","e_planet"]

"""
# 3 anillos. Ley de potencias y exponencial
"""
def f(t_steps,T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1, theta_angle,  e_ellipse, initial_opacity,e_planet):
   
    alpha_powlaw_exp = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1, theta_angle,  e_ellipse, initial_opacity,e_planet]
    relflux = modelos_planet_ring.relflux_model_planet_powlaw_and_exp(t_steps,rand_points_in_circle,alpha_powlaw_exp)
    return relflux

initial_conditions=[T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1, theta_angle,  e_ellipse, initial_opacity,e_planet]
lista = ["T0", "velocity", "yp","radius_planet","planet_separation","radius_central_1","radius_central_2","gamma_decay_1","theta_angle", "e_ellipse", "initial_opacity","e_planet"]

"""
# 3 anillos. Suma de 2 leyes de potencias y ley de potencias
"""
def f(t_steps,T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1, gamma_decay_2,radius_central_3,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity,e_planet ):
    
    alpha_powlawsum_powlaw = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1, gamma_decay_2,radius_central_3,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity,e_planet ]
    relflux = modelos_planet_ring.relflux_model_planet_powlawsum_powlaw(t_steps,rand_points_in_circle,alpha_powlawsum_powlaw)
    return relflux

initial_conditions=[T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1, gamma_decay_2,radius_central_3,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity,e_planet ]
lista = ["T0", "velocity", "yp","radius_planet", "planet_separation","radius_central_1","radius_central_2","gamma_decay_1","gamma_decay_2","radius_central_3","radius_limit_1","radius_limit_2", "theta_angle", "e_ellipse", "initial_opacity","e_planet" ]

"""
# 3 anillos. Suma de dos leyes de potencias y suma de dos leyes de potencias
"""
def f(t_steps,T0, velocity, yp,radius_planet, planet_separation,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity,e_planet):
    alpha_powlawsum_powlawsum = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity,e_planet]
    relflux_powlawsum_powlawsum = modelos_planet_ring.relflux_model_planet_powlawsum_powlawsum(t_steps,rand_points_in_circle,alpha_powlawsum_powlawsum)
    return relflux_powlawsum_powlawsum

initial_conditions=[T0, velocity, yp,radius_planet, planet_separation,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity,e_planet]
lista = ["T0", "velocity", "yp","radius_planet", "planet_separation","radius_central_1","gamma_decay_1","radius_central_2","radius_central_3","gamma_decay_3","radius_central_4","radius_limit_1","radius_limit_2", "theta_angle", "e_ellipse", "initial_opacity","e_planet"]

"""
# 3 anillos. Suma de 2 leyes de potencias y ley de potencias por exponencial
"""
def f(t_steps,T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1, radius_central_3,gamma_decay_3,radius_limit, theta_angle, e_ellipse, initial_opacity, e_planet):
    alpha_powlawsum_powlawexp = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1, radius_central_3,gamma_decay_3,radius_limit, theta_angle, e_ellipse, initial_opacity, e_planet]
    relflux_powlawsum_powlawexp = modelos_planet_ring.relflux_model_planet_powlawsum_powlawexp(t_steps,rand_points_in_circle,alpha_powlawsum_powlawexp)
    return relflux_powlawsum_powlawexp

initial_conditions=[T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1, radius_central_3,gamma_decay_3,radius_limit, theta_angle, e_ellipse, initial_opacity, e_planet]
lista =["T0", "velocity", "yp","radius_planet", "planet_separation","radius_central_1","radius_central_2","gamma_decay_1", "radius_central_3","gamma_decay_3","radius_limit", "theta_angle", "e_ellipse", "initial_opacity", "e_planet"]

"""
# 3 anillos. Ley de potencias y ley de potencias.
"""
def f(t_steps,T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity,e_planet):
    alpha_powlaw_powlaw = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity,e_planet]
    relflux = modelos_planet_ring.relflux_powlaw_powlaw(t_steps,rand_points_in_circle,alpha_powlaw_powlaw)
    return relflux


initial_conditions =   [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity,e_planet]
lista = ["T0", "velocity", "yp","radius_planet", "planet_separation","radius_central_1","radius_limit", "gamma_decay_1","gamma_decay_2", "theta_angle", "e_ellipse", "initial_opacity","e_planet"]

"""
# 4 anillos de opacidad constante
"""
def f(t_steps,T0, velocity, yp,radius_planet, planet_separation,radius1,radius2,radius3,radius4,val1,val2,val3,val4, theta_angle, e_ellipse,e_planet):
    alpha_cst = [T0, velocity, yp,radius_planet, planet_separation,radius1,radius2,radius3,radius4,val1,val2,val3,val4, theta_angle, e_ellipse,e_planet]
    relflux = modelos_planet_ring.relflux_model_planet_cst(t_steps,rand_points_in_circle,alpha_cst)
    return relflux

initial_conditions = [T0, velocity, yp,radius_planet, planet_separation,radius1,radius2,radius3,radius4,val1,val2,val3,val4, theta_angle, e_ellipse,e_planet] 
bounds = ((788,2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0),(796,10,0,1,1,10,20,20,30,1,1,1,1,1,1,0.35))
lista = ["T0", "velocity", "yp","radius_planet", "planet_separation","radius1","radius2","radius3","radius4","val1","val2","val3","val4", "theta_angle", "e_ellipse","e_planet"] 

"""








"""
# MINIMIZATION

optimal_values, cov = curve_fit(f,x_data,y_data,initial_conditions,bounds=bounds,diff_step=1e-2)

print("\n Datos de la minimización de chi2 \n")
errores = np.sqrt(np.diag(cov))


# RESULTS
k = -1      
for n in initial_conditions:
    k = k+1
    name = lista[k]
    print(name,"= ",optimal_values[k]," Error = ",errores[k],"\n")
print("\n Code runned for : ",(time()-start_time)," s \n")
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
ax[1].set_ylabel("Residuos")
ax[1].set_ylim([-0.012,0.012])
ax[1].set_xlabel("Tiempo [Fecha Kepler]")
fig2.tight_layout()
plt.show()
fig2.savefig("/Users/rubenbalbastre/Desktop/Física/CUARTO AÑO/2do cuatrimestre/TFG/Presentación/figura4_ampl.pdf")