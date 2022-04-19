
#################################################################################
####### THIS SCRIPT MINIMIZES THE CHI^2 TEST OF A MODEL COMPARED WITH THE DATA ##
# in this script we use the least_squares function, which makes us redefine our functions
#################################################################################

from scipy.optimize import least_squares
from pandas import read_table
import modelos_planet_ring
from libtransits import generate_randpoints_in_circle_limbdark
from time import time
import numpy as np 

start_time = time()

# DATA 
data_of_interest = read_table('/Users/rubenbalbastre/Desktop/Física/CUARTO AÑO/2do cuatrimestre/TFG/Curvas de luz/data_of_interest.txt',delim_whitespace=True).to_numpy()
Nrp_limb = 20000 # Number of points
rand_points_in_circle = generate_randpoints_in_circle_limbdark(Nrp_limb)
t_steps=data_of_interest[:,0] # time to integrate: the time column of the observed data # np.linspace(780,810,m)
x_data = data_of_interest[:,0]
y_data = data_of_interest[:,1]

# MODELS

# 2 anillos. Ley de potencias

def f(T0, velocity, yp, radius_planet,radius_empty_ellipse, rad_central, gamma_decay, theta_angle, e_ellipse,e_planet, initial_opacity,t_steps):
    
    alpha = [T0, velocity, yp, radius_planet,radius_empty_ellipse, rad_central, gamma_decay, theta_angle, e_ellipse,initial_opacity,e_planet]
    relflux =  modelos_planet_ring.relflux_model_allt_planet_ellipse_powerlaw(t_steps,rand_points_in_circle,alpha)
    return relflux

T0 =  792.7330149230421 
velocity =  5.892690452128485 
yp =  -0.6449754594452981  
radius_planet =  0.276 
planet_separation =  0
rad_central =  0.8588642775999662  
gamma_decay =  1.5894790896643616  
theta_angle =  0.056737515258601504  
e_ellipse =  0.9979695522067488  
initial_opacity =  0.9993692921997381
e_planet =  0.15

initial_conditions = [T0,velocity,yp,radius_planet,planet_separation,rad_central,gamma_decay,theta_angle,e_ellipse,e_planet,initial_opacity]
bounds = ((788,1,-1,0,0,0,0,0,0,0,0),(796,10,1,1,10,30,10,10,1,1,1))
lista = ["T0","velocity","yp","radius_planet","planet_separation","rad_central","gamma_decay","theta_angle","e_ellipse","e_planet","initial_opacity"]


# 2 anillo. Suma de 2 leyes de potencias
"""
def f(t_steps,T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity,e_planet):
    
    alpha_powlawsum = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity,e_planet]
    relflux = modelos_planet_ring.relflux_model_planet_powlaw_sum(t_steps,rand_points_in_circle,alpha_powlawsum)

    return relflux

T0 =  792.7329244132397 
velocity =  6.058795884931286 
yp =  -0.6234316496989329
radius_planet =  0.27600000000000013 
radius_empty_ellipse =  0.27599999999999997  
radius_central_1 =  0.8788602780343863 
radius_central_2 =  0.7289851845216327
gamma_decay_1 =  1.5614611129237697 
radius_limit =  0.999935719291165 
theta_angle =  0.056920249775347154 
e_ellipse =  0.9984006654603454  
initial_opacity =  0.992665399035574 
e_planet =  0.15

initial_conditions = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity,e_planet]
bounds = ((788,1,-1,0,0,0,0,0,0,0,0,0,0),(796,10,1,1,5,30,30,10,10,3,1,1,1))

"""
# 3 anillos. Ley de potencias y ley de potencias.
"""
def f(t_steps,T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity,e_planet):
    alpha_powlaw_powlaw = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity,e_planet]
    relflux = modelos_planet_ring.relflux_powlaw_powlaw(t_steps,rand_points_in_circle,alpha_powlaw_powlaw)
    return relflux
"""










# MINIMIZATION

def fun(param):
    T0, velocity, yp, radius_planet,radius_empty_ellipse, rad_central, gamma_decay, theta_angle, e_ellipse,e_planet, initial_opacity=param
    return f(T0, velocity, yp, radius_planet,radius_empty_ellipse, rad_central, gamma_decay, theta_angle, e_ellipse,e_planet, initial_opacity,x_data)-y_data

res = least_squares(fun,initial_conditions,bounds=bounds)

# RESULTS
k = -1      
for n in initial_conditions:
    k = k+1
    name = lista[k]
    print(name,"= ",res.x[k]," Error = ",res.fun[k],"\n")

cov = - np.linalg.inv(np.dot(np.matrix.transpose(res.jac),res.jac)) # calculate covariance matrix
print(cov) # show the covariance matrix
print(np.sqrt(np.diag(cov))) # show the "error" of the parameters

# Chi^2
expected_value = f(*res.x,x_data)
chi2 = (y_data-expected_value)**2 / expected_value
print("Chi^2: ",sum(chi2))

# PLOTTING
import matplotlib.pyplot as plt
plt.plot(x_data,y_data,'.',label="Data")
plt.plot(x_data,f(*res.x,x_data),label="Model")
plt.legend()
plt.show()
#https://stats.stackexchange.com/questions/231868/relation-between-covariance-matrix-and-jacobian-in-nonlinear-least-squares