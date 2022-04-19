#################################################################################
####### THIS SCRIPT CONTAINS A VARIETY OF TRANSIT MODELS ########################
#################################################################################



# Import the needed librarys

import modelos 

# This scripts implements both the planet model and the rings system model. Several models are proposed.

""" Model 1.

2 rings

Opacity following power-law * exp opacity.

Non-spherical planet

"""

def rel_flux_model_allt_planet_powerlawexp(time_steps,rand_points_in_circle,alpha):

    # Data for the function
    T0, velocity, yp,radius_empty_ellipse,radius_central, radius_limit,gamma_decay, theta_angle,  e_ellipse, initial_opacity,radius_planet,e_planet = alpha

    # Functions data to calculate the planet and ring fluxes
    alpha_planet = [T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_empty_ellipse,radius_planet,radius_central,radius_limit,gamma_decay, theta_angle,  e_ellipse, initial_opacity ]
                
    # Calculation of the fluxes                             
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps,rand_points_in_circle,alpha_planet)
    relflux_ellipse = modelos.relflux_model_allt_powerlawexp(time_steps, rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux


""" Model 2.

 2 rings.

Opacity following powlaw opacity.

Non-spherical planet
    
"""

def relflux_model_allt_planet_ellipse_powerlaw (time_steps,rand_points_in_circle,alpha):

    # Data for the function                              
    T0, velocity, yp, radius_planet,radius_empty_ellipse, rad_central, gamma_decay, theta_angle, e_ellipse, initial_opacity,e_planet = alpha
    
    # Data to compute the planet and ring fluxes
    alpha_planet = [T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet,radius_empty_ellipse,rad_central, gamma_decay, theta_angle, e_ellipse, initial_opacity ]

    # Calculation of the fluxes                             
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps,rand_points_in_circle,alpha_planet)
    relflux_ellipse = modelos.relflux_model_allt_powerlaw(time_steps, rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux
        
""" Model 3.

2 rings.

Opacity following powlaw opacity. 

Non-spherical planet
    
"""

def relflux_model_e_planet_ellipse_powlaw(time_steps,rand_points_in_circle,alpha):

    T0, velocity, yp,radius_planet, rad_central,radius_empty_ellipse, gamma_decay, theta_angle, e_ellipse, radius_planet, e_planet ,initial_opacity= alpha

    # Data to compute the planet and ring fluxes
    alpha_planet = [T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet,radius_empty_ellipse,rad_central, gamma_decay, theta_angle, e_ellipse, initial_opacity ]
    
    # Calculation of the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet)                         
    relflux_ellipse = modelos.relflux_model_allt_powerlaw(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1 - (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux
    
""" Model 4.

2 rings

Opacity following the power-law sum. 

Non-spherical planet

"""

def relflux_model_planet_powlaw_sum(time_steps,rand_points_in_circle,alpha):

    # Data for the function   
    T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse, initial_opacity, e_planet = alpha

    # Data to compute the planet and ring fluxes
    alpha_planet = [T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse,  initial_opacity]

    # Calculation of the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet)                    
    relflux_ellipse = modelos.relflux_model_allt_powerlaw_sum(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux
    
""" Model 5.

3 rings

Opacity following the power-law and exp. 

Non-spherical planet

"""

def relflux_model_planet_powlaw_and_exp(time_steps,rand_points_in_circle,alpha):

    # Data for the function
    T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1, theta_angle,  e_ellipse, initial_opacity,e_planet=alpha
        
    # Data to compute the planet and ring fluxes
    alpha_planet = [T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1, theta_angle,  e_ellipse, initial_opacity]
    
    # Calculation of the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet) 
    relflux_ellipse = modelos.relflux_model_allt_powerlaw_and_exp(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux
    
""" Model 6.

3 rings

Opacity following the power-law and pow-law*exp

Non-spherical planet

"""

def relflux_model_planet_powlaw_and_powlawexp(time_steps,rand_points_in_circle,alpha):

    # Data for the function
    T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,gamma_decay_2, theta_angle,  e_ellipse, initial_opacity, e_planet=alpha
        
    # Data to compute the fluxes
    alpha_planet = [T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,gamma_decay_2, theta_angle,  e_ellipse, initial_opacity]
    
    # Calculate the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet) 
    relflux_ellipse = modelos.relflux_model_allt_powerlaw_and_powlawexp(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux
    
"""

Model 7.

Opacity powlaw for nrings

"""
def relflux_model_planet_powlaw_nrings(time_steps,rand_points_in_circle,alpha):

    # Data for the function
    T0, velocity, yp,radius_planet,planet_separation, theta_angle,  e_ellipse, initial_opacity,rvector,gvector,e_planet=alpha
        
    # Data to compute the fluxes
    alpha_planet = [T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet,planet_separation, theta_angle,  e_ellipse, initial_opacity,rvector,gvector]

    # Calculate the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet)  
    relflux_ellipse = modelos.relflux_model_allt_powerlaw_nrings(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux

""" 

Model 8. 

Opacity powlaw_sum & powlaw*exp for 3 rings

"""

def relflux_model_planet_powlawsum_powlawexp(time_steps,rand_points_in_circle,alpha):

    # Data for the function
    T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_central_3,gamma_decay_3,radius_limit, theta_angle, e_ellipse, initial_opacity, e_planet = alpha
    # Data to compute the fluxes
    alpha_planet = [T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_central_3,gamma_decay_3,radius_limit, theta_angle, e_ellipse, initial_opacity]

    # Calculate the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet)  
    relflux_ellipse = modelos.relflux_model_allt_powerlawsum_and_powlawexp(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux

"""

Model 9. 

Opacity powlaw_sum & powlaw for 3 rings

"""

def relflux_model_planet_powlawsum_powlaw(time_steps,rand_points_in_circle,alpha):

    # Data for the function
    T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1, gamma_decay_2,radius_central_3,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity,e_planet = alpha

    # Data to compute the fluxes
    alpha_planet = [T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_central_3,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity ]

    # Calculate the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet)  
    relflux_ellipse = modelos.relflux_model_allt_powerlawsum_and_powlaw(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux

""" 

Model 10.

Opacity powlawsum & powlawsum for 3 rings

"""

def relflux_model_planet_powlawsum_powlawsum(time_steps,rand_points_in_circle,alpha):

    # Data
    T0, velocity, yp,radius_planet, planet_separation,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity,e_planet= alpha

    # Data to fluxes
    alpha_planet =[T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity]

    # Calculate the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet)  
    relflux_ellipse = modelos.relflux_model_allt_powerlawsum_and_powlawsum(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux

def relflux_powlaw_powlaw(time_steps,rand_points_in_circle,alpha):

    # Data
    T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity,e_planet= alpha
    # Data to fluxes
    alpha_planet =[T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity]
    # Calculate the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet)  
    relflux_ellipse = modelos.relflux_model_allt_powerlaw_and_powlaw(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux


""" Model 11. 

Constant opacity 4 rings

"""

def relflux_model_planet_cst(time_steps,rand_points_in_circle,alpha):

    # Data
    T0, velocity, yp,radius_planet, planet_separation,radius1,radius2,radius3,radius4,val1,val2,val3,val4, theta_angle, e_ellipse,e_planet= alpha

    # Data to fluxes
    alpha_planet =[T0, velocity, yp, radius_planet, e_planet, theta_angle]
    alpha_ellipse = [T0, velocity, yp,radius_planet, planet_separation,radius1,radius2,radius3,radius4,val1,val2,val3,val4, theta_angle, e_ellipse]

    # Calculate the fluxes
    relflux_planet = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha_planet)  
    relflux_ellipse = modelos.relflux_model_allt_cst(time_steps,rand_points_in_circle,alpha_ellipse)
    
    relflux = 1- (( 1 - relflux_planet ) + (1 -relflux_ellipse))
    
    return relflux

""" Model 12.

Only oblate planet
"""

def relflux_oblate_planet(time_steps,rand_points_in_circle,alpha):

    T0, velocity, yp, radius_planet,e_planet,theta_angle_planet = alpha

    alpha = [T0, velocity, yp, radius_planet,e_planet,theta_angle_planet]
    relflux = modelos.relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,alpha)  

    return relflux


