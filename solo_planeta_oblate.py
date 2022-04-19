
#################################################################################
####### THIS SCRIPT MINIMIZES THE CHI^2 TEST OF AN OBLATE PLANET MODEL ##
#################################################################################


from scipy.optimize import curve_fit
from pandas import read_table
import modelos_planet_ring
from libtransits import generate_randpoints_in_circle_limbdark
from time import time
import numpy as np 
import pandas as pd 

start_time = time() 

# Data 
data_of_interest = read_table('/Users/rubenbalbastre/Desktop/Física/CUARTO AÑO/2do cuatrimestre/TFG/Curvas de luz/data_of_interest.txt',delim_whitespace=True).to_numpy()
Nrp_limb = 2000 # Number of points
rand_points_in_circle = generate_randpoints_in_circle_limbdark(Nrp_limb)
t_steps=data_of_interest[:,0] 
x_data = data_of_interest[:,0]
y_data = data_of_interest[:,1]

# Model to minimize

# 2 anillos. Ley de potencias

def f(t_steps,T0, velocity, yp, radius_planet,e_planet,theta_angle):

    alpha = [T0, velocity, yp, radius_planet,e_planet,theta_angle]
    relflux =  modelos_planet_ring.relflux_oblate_planet(t_steps,rand_points_in_circle,alpha)
    return relflux

# Initial conditions

T0 =  792.7329244132397 
velocity =  2
yp =  -0.6234316496989329
radius_planet =  0.4
separation_planet =  0
radius_central_1 =  0.8788602780343863 
radius_central_2 =  0.7289851845216327
gamma_decay_1 =  1.5614611129237697 
radius_limit =  0.999935719291165 
theta_angle =  0.056920249775347154 
e_ellipse =  0.9984006654603454  
initial_opacity =  0.992665399035574 
e_planet =  0.15

initial_conditions = [T0, velocity, yp,radius_planet,e_planet,theta_angle]
# Name of the minimized parameters
lista = ["T0", "velocity", "yp","radius_planet","e_planet","theta_angle"]

# Bounds for the minimization

# Minimize
optimal_values, cov = curve_fit(f,x_data,y_data,initial_conditions)

print(cov) # show covariance matrix
print("\n Datos de la minimización de chi2 \n")
errores = np.sqrt(np.diag(cov)) # show the uncertanty

# Show the minimization results
k = -1      
for n in initial_conditions:
    k = k+1
    name = lista[k]
    print(name,"= ",optimal_values[k]," Error = ",errores[k],"\n")

# Chi^2
expected_value = f(x_data,*initial_conditions)
a = (y_data-expected_value)**2 / expected_value
print("Chi^2: ",sum(a))

# show the running time
print("\n Code runned for : ",(time()-start_time)," s \n")

# plot the results
import matplotlib.pyplot as plt
plt.plot(x_data,y_data,'.',label="Data")
plt.plot(x_data,expected_value,label="Model")
plt.legend()
plt.show()