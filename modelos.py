#################################################################################
####### THIS SCRIPT CONTAINS A VARIETY OF RING MODELS ########################
#################################################################################

# IMPORT LIBRARYS

import numpy as np
import opacity_models
import libtransits as transits
import numpy.random as rd
rd.seed(314)                

opac_lim = 0.001 # If we change it 0.001, running time increases considerably. The same occurs when we decrease the number of points

""" PLANET MODELS """

""" Model 1.

Calculates the flux that passes through a non-spherical planet whose projection has an eccentricity e_planet and inclination angle theta_angle_planet.
    
"""

def relflux_model_singlet_nonsph_planet(time,rand_points_in_circle, T0, velocity, yp,radius_planet,e_planet,theta_angle_planet):

    Nptotal = len(rand_points_in_circle)

    xp = velocity*(time - T0)

    # Define relative coordinates to the planet
    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle_planet)
    st = np.sin(theta_angle_planet)
    
    # Define the rotated coordinates of the ellipse (planet)
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # radius of the rotated ellipse
    radrot_planet = np.sqrt( xrot*xrot + yrot*yrot/(1-e_planet**2.) )

    #In/Out conditions
    points_in_planet = radrot_planet < radius_planet
    
    counts_in_planet = np.array(points_in_planet,float)

    # Flux 
    relflux = 1 - counts_in_planet[points_in_planet==True].sum()/Nptotal

    return relflux

def relflux_model_allt_nonsph_planet(time_steps, rand_points_in_circle,
                                       alpha):
                                       
    T0, velocity, yp, radius_planet,e_planet,theta_angle_planet = alpha # parameters of the model

    Nt = len(time_steps)

    rflux = np.empty(Nt, float)

    for i in np.arange(Nt): # integrate for each time step
        rflux[i] = relflux_model_singlet_nonsph_planet(time_steps[i],rand_points_in_circle, T0,velocity, yp,radius_planet, e_planet,theta_angle_planet)
    return rflux

    

""" RING MODELS """

""" Model 1.

    2 rings. 
    
    Pow-law opacity with separation distance between the planet and the ring included.
    
"""

def relflux_model_singlet_powerlaw(time, rand_points_in_circle, T0, velocity, yp, radius_planet,planet_separation, rad_central, gamma_decay,  theta_angle, e_ellipse, initial_opacity):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    #  Rotated coordinates
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    a_max = rad_central*pow(opac_lim, -1./gamma_decay) # Define the maximum ellipse we will consider 

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse<radius_empty_ellipse)
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)
    points_in_ellipse = (radrot_ellipse <= a_max)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse,float)
    
    # Evaluate the opacity of the points inside the ellipse (ring)
    counts_in_ellipse[counts_in_ellipse == 1] =initial_opacity * transits.opacity_decay_powlaw(radrot_ellipse[counts_in_ellipse == 1],rad_central, gamma_decay)

    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==False].sum()
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==False].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlaw (time_steps, rand_points_in_circle, alpha):

    T0, velocity, yp,radius_planet,planet_separation,rad_central, gamma_decay, theta_angle, e_ellipse, initial_opacity = alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlaw(time_steps[i],rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,
                                                         rad_central,
                                                         gamma_decay,
                                                         theta_angle,
                                                         e_ellipse,
                                                         initial_opacity)

    return rflux

""" Model 2.

    2 rings

    Power-law * exp opacity function with the separation distance between the planet and the ring included.
   
"""

def relflux_model_singlet_powerlawexp (time,  rand_points_in_circle, T0, velocity, yp, radius_planet,planet_separation,radius_central,radius_limit,gamma_decay,  theta_angle, e_ellipse, initial_opacity):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    # Relative coordinates to the planet

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    # Rotated coordinates

    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points'
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    # Define the maximum ellipse we will consider. Two options:
    
    a_max = radius_central*pow(opac_lim, -1./gamma_decay)
        
    # a_max = - np.log(opac_lim)/gamma_decay_2
    # a_max = rad_central*pow(opac_lim, -1./gamma_decay)

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2)
    
    points_in_ellipse = (radrot_ellipse <= a_max)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)

    counts_in_ellipse[counts_in_ellipse == 1] =initial_opacity * opacity_models.opacity_powlawexp(radrot_ellipse[counts_in_ellipse == 1],radius_central,radius_limit,gamma_decay)

    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlawexp(time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,planet_separation,radius_planet,radius_central,radius_limit,gamma_decay, theta_angle,  e_ellipse, initial_opacity = alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlawexp(time_steps[i],  rand_points_in_circle, T0, velocity, yp, radius_planet,planet_separation,radius_central,radius_limit,gamma_decay,  theta_angle, e_ellipse, initial_opacity)
    
    return rflux
    
    
""" Model 3.

    2 rings

    Pow-law-sum opacity function with the separation distance between the planet and the ring included.
  
"""

def relflux_model_singlet_powerlaw_sum (time, rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit,  theta_angle, e_ellipse, initial_opacity):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    a_max = radius_central_1 * pow(opac_lim, -1./gamma_decay_1)

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_ellipse = (radrot_ellipse <= a_max)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)

    counts_in_ellipse[counts_in_ellipse == 1] = initial_opacity * opacity_models.opacity_powlawsum(radrot_ellipse[counts_in_ellipse == 1],radius_central_1,radius_central_2,gamma_decay_1,radius_limit)

    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlaw_sum (time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit, theta_angle,  e_ellipse,  initial_opacity= alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlaw_sum(time_steps[i],rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_limit,  theta_angle, e_ellipse, initial_opacity)

    return rflux
    
    
""" Model 4.

    3 rings

    Pow-law and exp opacity function with the separation distance between the planet and the ring included.
    
"""

def relflux_model_singlet_powerlaw_and_exp (time,  rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,  theta_angle, e_ellipse, initial_opacity):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    a_max = radius_central_1*pow(opac_lim, -1./gamma_decay_1) # Define the maximum ellipse we will consider. Con esto ya estaría bien porque dos potencias caen más rápido que una.

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_ellipse = (radrot_ellipse <= a_max)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)
    
    counts_in_ellipse[counts_in_ellipse == 1] = initial_opacity* opacity_models.opacity_powlaw_and_exp(radrot_ellipse[counts_in_ellipse == 1],radius_central_1,gamma_decay_1,radius_central_2)

    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlaw_and_exp (time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1, theta_angle,  e_ellipse, initial_opacity= alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlaw_and_exp(time_steps[i],rand_points_in_circle,
                                                         T0, velocity, yp,radius_planet, planet_separation,
                                                             radius_central_1,radius_central_2,gamma_decay_1,
                                                         theta_angle,
                                                         e_ellipse, initial_opacity)

    return rflux


""" Model 5.

    Pow-law and pow-law-exp opacity function with the separation distance between the planet and the ring included.
    
"""

def relflux_model_singlet_powerlaw_and_powlawexp (time,  rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1, gamma_decay_2, theta_angle, e_ellipse, initial_opacity):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    a_max = radius_central_1*pow(opac_lim, -1./gamma_decay_1) # Define the maximum ellipse we will consider. Con esto ya estaría bien porque dos potencias caen más rápido que una.

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_ellipse = (radrot_ellipse <= a_max)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)
    
    counts_in_ellipse[counts_in_ellipse == 1] = initial_opacity * opacity_models.opacity_powlaw_and_powlawexp(radrot_ellipse[counts_in_ellipse == 1],radius_central_1,gamma_decay_1,radius_central_2,gamma_decay_2)

    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlaw_and_powlawexp (time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,radius_planet,planet_separation,radius_central_1,radius_central_2,gamma_decay_1,gamma_decay_2, theta_angle,  e_ellipse, initial_opacity = alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlaw_and_powlawexp(time_steps[i],rand_points_in_circle,
                                                         T0, velocity, yp,radius_planet, planet_separation,
                                                             radius_central_1,radius_central_2,gamma_decay_1,gamma_decay_2,
                                                         theta_angle,
                                                         e_ellipse,initial_opacity)

    return rflux

""" Model 6.

    Pow-law n rings opacity function with the separation distance between the planet and the ring included.
    
"""

def relflux_model_singlet_powerlaw_nrings (time,  rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation, theta_angle, e_ellipse, initial_opacity,rvector,gvector):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    l=len(rvector)
    a_max = rvector[l-1]*pow(opac_lim, -1./gvector[l-1]) # Define the maximum ellipse we will consider.

    # In/Out conditions
    radius_empty_ellipse = planet_separation+radius_planet
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_ellipse = (radrot_ellipse <= a_max)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)

    counts_in_ellipse[counts_in_ellipse == 1] = initial_opacity * opacity_models.opacity_piecewise(radrot_ellipse[counts_in_ellipse == 1],rvector,gvector)

    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlaw_nrings (time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,radius_planet,planet_separation, theta_angle,  e_ellipse, initial_opacity,rvector,gvector = alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlaw_nrings(time_steps[i],rand_points_in_circle,
                                                         T0, velocity, yp,radius_planet, planet_separation,
                                                         theta_angle,
                                                         e_ellipse,initial_opacity,rvector,gvector)

    return rflux
    

""" 
Model 7.

    Pow-law-sum and powlaw*exp opacity for 3 rings opacity. Separation distance between the planet and the ring included.
"""

def relflux_model_singlet_powerlawsum_and_powlawexp (time,  rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_central_3,gamma_decay_3,radius_limit, theta_angle, e_ellipse, initial_opacity):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    a_max = radius_central_1*pow(opac_lim, -1./gamma_decay_1) # Define the maximum ellipse we will consider. Con esto ya estaría bien porque dos potencias caen más rápido que una.

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_ellipse = (radrot_ellipse <= a_max)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)
    
    counts_in_ellipse[counts_in_ellipse == 1] = initial_opacity * opacity_models.opacity_powlawsum_and_powlawexp(radrot_ellipse[counts_in_ellipse==1],radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_limit)

    # Flux
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlawsum_and_powlawexp (time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_central_3,gamma_decay_3,radius_limit, theta_angle, e_ellipse, initial_opacity = alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlawsum_and_powlawexp(time_steps[i],rand_points_in_circle,T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_central_3,gamma_decay_3,radius_limit, theta_angle, e_ellipse, initial_opacity)
    return rflux

""" 
Model 8.

    Pow-law-sum and powlaw opacity for 3 rings opacity. Separation distance between the planet and the ring included.
"""

def relflux_model_singlet_powerlawsum_and_powlaw(time,  rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_central_3,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    a_max = radius_central_1*pow(opac_lim, -1./gamma_decay_1) # Define the maximum ellipse we will consider. Con esto ya estaría bien porque dos potencias caen más rápido que una.

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_ellipse = (radrot_ellipse <= a_max)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)
    
    counts_in_ellipse[counts_in_ellipse == 1] = initial_opacity * opacity_models.opacity_powlawsum_and_powlaw(radrot_ellipse[counts_in_ellipse==1],radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,radius_limit_1,radius_limit_2)

    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlawsum_and_powlaw (time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_central_3,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity = alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlawsum_and_powlaw(time_steps[i],rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_central_2,gamma_decay_1,radius_central_3,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity)
    return rflux


""" 
Model 9.

    Pow-law-sum and powlaw-sum opacity for 3 rings opacity. Separation distance between the planet and the ring included.
"""

def relflux_model_singlet_powerlawsum_and_powlawsum(time,  rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    a_max = radius_central_1*pow(opac_lim, -1./gamma_decay_1) # Define the maximum ellipse we will consider. Con esto ya estaría bien porque dos potencias caen más rápido que una.

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_ellipse = (radrot_ellipse <= a_max)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)
    
    counts_in_ellipse[counts_in_ellipse == 1] = initial_opacity * opacity_models.opacity_powlawsum_and_powlawsum(radrot_ellipse[counts_in_ellipse==1],radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2)
    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlawsum_and_powlawsum (time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,radius_planet, planet_separation,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity = alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlawsum_and_powlawsum(time_steps[i],  rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2, theta_angle, e_ellipse, initial_opacity)
    
    return rflux 


""" Model 10.

    Pow-law and pow-law opacity function with the separation distance between the planet and the ring included.
    
"""

def relflux_model_singlet_powerlaw_and_powlaw (time,  rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit,gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    a_max = radius_central_1*pow(opac_lim, -1./gamma_decay_1) # Define the maximum ellipse we will consider. Con esto ya estaría bien porque dos potencias caen más rápido que una.

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_ellipse = (radrot_ellipse <= a_max)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)
    
    counts_in_ellipse[counts_in_ellipse == 1] =  initial_opacity * opacity_models.opacity_powlaw_powlaw(radrot_ellipse[counts_in_ellipse == 1],radius_central_1,radius_limit,gamma_decay_1,gamma_decay_2)

    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_powerlaw_and_powlaw (time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit, gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity = alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_powerlaw_and_powlaw(time_steps[i],rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius_central_1,radius_limit,gamma_decay_1,gamma_decay_2, theta_angle, e_ellipse, initial_opacity)

    return rflux



""" Model 11.

    4 constant opacity rings separated
    
"""

def relflux_model_singlet_cst (time,  rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius1,radius2,radius3,radius4,val1,val2,val3,val4,theta_angle, e_ellipse):

    Nptotal = len(rand_points_in_circle) # number of total points

    # Relative coordinates
    
    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st

    # 'Elliptical radius for the points', which correspond to the value
    # of the semi-major axis of the elliptical annulus to which this point
    # belongs
    
    radrot_ellipse = np.sqrt(xrot**2. + (yrot**2./(1 - e_ellipse**2.)))
    
    a_max = 35 # Define the maximum ellipse we will consider. Con esto ya estaría bien porque dos potencias caen más rápido que una.

    # In/Out conditions
    radius_empty_ellipse = radius_planet+planet_separation
    points_in_empty_ellipse = (radrot_ellipse < radius_empty_ellipse)
    
    points_in_ellipse = (radrot_ellipse <= a_max)
    
    points_in_circle = (xrel*xrel+yrel*yrel < radius_planet**2.)

    # Opacity
    
    counts_in_ellipse = np.array(points_in_ellipse, float)
    
    counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)
    
    counts_in_ellipse[counts_in_ellipse == 1] =  opacity_models.opacity_cst(radrot_ellipse[counts_in_ellipse == 1],radius_empty_ellipse,radius1,radius2,radius3,radius4,val1,val2,val3,val4)

    # Flux
    
    N_in_ellipse = counts_in_ellipse[points_in_circle==0].sum()
    
    N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==0].sum()

    relflux = 1. - ((N_in_ellipse-N_in_empty_ellipse)/Nptotal)

    return relflux


def relflux_model_allt_cst (time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp,radius_planet, planet_separation, radius1,radius2,radius3,radius4,val1,val2,val3,val4,theta_angle, e_ellipse = alpha

    Nt = len(time_steps)
    
    rflux = np.empty(Nt, float)
    
    for i in range(Nt):
    
        rflux[i] = relflux_model_singlet_cst(time_steps[i],rand_points_in_circle, T0, velocity, yp,radius_planet, planet_separation,radius1,radius2,radius3,radius4,val1,val2,val3,val4, theta_angle, e_ellipse)

    return rflux








# Arreglar modelo ¿?

"""

PROBABLY TO DELETE

    Model XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX.
    
    Short description of the general idea:
    
    Here it is presented a model that calculates the flux across a system of n rings. The opacity of the rings follow a pow-law function.
    Some of their properties are listed as input vectors:
    
    * The internal radius
    
        r_int: int_ring_radius_vector
        
    * The external radius
    
        r_ext: ext_ring_radius_vector
        
    * The gamma exponents of the pow-law
    
        G: gamma_decay_vector       where gamma_decay_vector[0] = 0   indicates a first ring with opacity = initial_opacity
    
    The opacity constants R are derived from the opacity continuity equations. The opacity for the first ring takes the form:
    
    O_0(r) = initial_opacity      r < r_ext[0]     where R[0] = int_ring_radius[0]      and A[0] = initial_opacity * (r_ext[0]/R[0])**(-G[0])
            (r/R[0])**(-G[0])     r > r_ext[0]
    
    The expression of the m ring's opacity is, for m > 0.
    
    O_m(r) =  0               r < r_ext[m]          where R[m] = (A[m-1]**(1/G[m]) r_int[m-1]      and A[m] = (r_ext[m]/R[m])**(-G[m])
            (r/R[m])**(-G[m]) r > r_ext[m]
                
    So we have two constant vectors for opacity, R and A. Both the ring's system and the planet share the same orientation angle but not the eccentricity. Therefore, the needed parameters of the model are:

    """
"""
def relflux_model_singlet_nrings_powlaw(time,rand_points_in_circle,T0,yp,radius_planet,velocity,theta_angle,initial_opacity,ext_ring_radius_vector,int_ring_radius_vector,number_rings,gamma_decay_vector,e_ellipse):

    # Create the ring separation vector
        
    ring_separation = ext_ring_radius_vector - int_ring_radius_vector
    
    # Create the R and A opacity's constants vectors
    
    R = np.empty(len(ext_ring_radius_vector))
    A = np.empty(len(ext_ring_radius_vector))
    
    # set the initial constants for the first ring
    
    R[0] = int_ring_radius_vector[0]
    A[0] =  (ext_ring_radius_vector[0]/R[0])**(-gamma_decay_vector[0])
    
    for number in np.arange(number_rings-1):
    
        ring=number+1
    
        R[ring] = A[ring-1] ** (1/gamma_decay_vector[ring]) * int_ring_radius_vector[ring]
        
        A[ring] = (ext_ring_radius_vector[ring] / R[ring])**(-gamma_decay_vector[ring])
        
    # Probe that input data has sense
    
    for element in np.arange(len(ring_separation)):
    
        if ring_separation[element]<0:
        
            import sys
            sys.exit("Input error. int_ring_radius must be greater than radius_planet")
    
    if len(gamma_decay_vector) == len(ring_separation):
    
        gamma_decay_vector=gamma_decay_vector
        
    else:
        print("Gamma_decay_vector is not the same length as ext_ring_radius and int_ring_radius")
        
    # General conditions of the problem
    
    Nptotal=len(rand_points_in_circle)
    
    xp=velocity*(time-T0)
    
    xrel = rand_points_in_circle[:,0] - xp
    yrel = rand_points_in_circle[:,1] - yp
    
    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st
        
    # Flux calculus for the n rings
    
    relflux = 1
    
    for ring in np.arange(number_rings): # Calculate the flux for each ring
    
        # Select the vector element's
                
        a_ellipse = ext_ring_radius_vector[ring]
        a_empty_ellipse = int_ring_radius_vector[ring]
        g_decay = gamma_decay_vector[ring]
    
        b_ellipse = a_ellipse * np.sqrt(1. - e_ellipse*e_ellipse)
        b_empty_ellipse = a_empty_ellipse * np.sqrt(1. - e_ellipse*e_ellipse)
    
        radrot_ellipse = np.sqrt( pow(xrot/a_ellipse,2) + pow(yrot/b_ellipse,2) )
        radrot_empty_ellipse = np.sqrt( pow(xrot/a_empty_ellipse,2) + pow(yrot/b_empty_ellipse,2) )
    
        # In/Out conditions
    
        points_in_circle = (xrel*xrel + yrel*yrel < radius_planet*radius_planet)

        points_in_ellipse = ( radrot_ellipse <= 1)
        counts_in_ellipse = np.array(points_in_ellipse, float)
    
        points_in_empty_ellipse = ( radrot_empty_ellipse < 1)
        counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)
        
        # Opacity function
        
        counts_in_ellipse[counts_in_ellipse==1] = opacity_models.opacity_powlaw(radrot_ellipse[counts_in_ellipse==1],R[ring],g_decay,a_empty_ellipse)

        # flux

        N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle=Flase].sum()

        N_in_ellipse = counts_in_ellipse[points_in_circle==False].sum()
        
        relflux = relflux - (N_in_ellipse-N_in_empty_ellipse)/Nptotal
        
        # Calculate the opacity vector that we will need in the next step ARREGLAR
    
        #constant_opacity_condition = (radrot_ellipse = a_ellipse )
    
        #opacity_vector[ring] = opacity_decay_powlaw_modelos(a_ellipse,radius_central,g_decay) # verify if this is correct
        
    opacity_vector = 0.3*np.ones(number_rings-1) # ?¿?¿?¿?
    
    # Flux calculus for the n rings-separation
    
    for separation in np.arange(number_rings-1):
    
        # Select the vector's element
    
        a_empty_ellipse = ext_ring_radius_vector[separation]
        a_ellipse = int_ring_radius_vector[separation + 1]
        
        b_ellipse = a_ellipse * np.sqrt(1. - e_ellipse*e_ellipse)
        b_empty_ellipse = a_empty_ellipse * np.sqrt(1. - e_ellipse*e_ellipse)
        
        # Radius is ellipse rotated system
    
        radrot_ellipse = np.sqrt( pow((xrot)/a_ellipse,2) + pow((yrot)/b_ellipse,2) )
        radrot_empty_ellipse = np.sqrt( pow((xrot)/a_empty_ellipse,2) + pow((yrot)/b_empty_ellipse,2) )
        
        # In/Out conditions
        
        points_in_circle = xrel*xrel + yrel*yrel < radius_planet*radius_planet

        points_in_ellipse = (radrot_ellipse <= 1)
        points_in_empty_ellipse = radrot_empty_ellipse < 1
        
        counts_in_ellipse = np.array(points_in_ellipse,float)
        counts_in_empty_ellipse = np.array(points_in_empty_ellipse,float)
        
        # Opacity ARREGLAR 
        
        cst = opacity_vector[separation]
        
        counts_in_ellipse[counts_in_ellipse==1] = opacity_constant( radrot_ellipse[counts_in_ellipse==1], a_empty_ellipse, cst) ## Arreglar
        
        # Flux
        
        N_in_ellipse = counts_in_ellipse[points_in_circle==False].sum()
        N_in_empty_ellipse = counts_in_empty_ellipse[points_in_circle==False].sum()

        relflux = relflux - (N_in_ellipse-N_in_empty_ellipse)/Nptotal
        
    
    return relflux
    
    
def relflux_model_allt_nrings_powlaw(time_steps,rand_points_in_circle,alpha):


    T0,yp,velocity,theta_angle,initial_opacity,ext_ring_radius_vector,int_ring_radius_vector,number_rings,gamma_decay_vector,e_ellipse = alpha

    Np=len(time_steps)
    
    relflux = np.empty(Np)
    
    for step in np.arange(Np):
        
        relflux[step] = relflux_model_singlet_nrings_powlaw(time_steps[step],rand_points_in_circle,T0,yp,velocity,theta_angle,initial_opacity,ext_ring_radius_vector,int_ring_radius_vector,number_rings,gamma_decay_vector,e_ellipse)
        
    return relflux
"""     

# Arreglar ¿?

""" model of n rings equi-spaciated """

"""
def relflux_model_singlet_nrings_powlaw(time,rand_points_in_circle,T0,yp,velocity,theta_angle,radius_planet,initial_opacity,ext_ring_radius,int_ring_radius,number_rings,gamma_decay_vector,rad_central_vector):

    # Creamos los vectores con los radios internos y externos equiespaciados
        
    rad = ext_ring_radius-int_ring_radius
    ring_separation = int_ring_radius-radius_planet
    
    # Probe that input data has sense
    
    if ring_separation<0:
        
        import sys
        sys.exit("Input error. int_ring_radius must be greater than radius_planet")
        
    int_ring_radius_vector=np.empty(number_rings)
    int_ring_radius_vector[0] = int_ring_radius
    
    ext_ring_radius_vector = np.empty(number_rings)
    ext_ring_radius_vector[0] = ext_ring_radius
    
    for number in np.arange(number_rings-1):
        int_ring_radius_vector[number+1] = ring_separation + ext_ring_radius_vector[number]
        ext_ring_radius_vector[number+1] = int_ring_radius_vector[number+1]+ rad
        
    # definimos las condiciones del problema como hemos hecho anteriormente
    
    Nptotal=len(rand_points_in_circle)
    
    xp=velocity*(time-T0)
    
    xrel = rand_points_in_circle[:,0] - xp
    yrel = rand_points_in_circle[:,1] - yp
    
    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st
    
    # Flux calculus for the n rings
    
    relflux = 1
    
    for ring in np.arange(number_rings): # HACER EL LOOP. Calcular el flujo para cada anillo
    
        a_ellipse = ext_ring_radius_vector[ring]
        a_empty_ellipse = int_ring_radius_vector[ring]
    
        b_ellipse = a_ellipse * np.sqrt(1. - e_ellipse*e_ellipse)
        b_empty_ellipse = a_empty_ellipse * np.sqrt(1. - e_ellipse*e_ellipse)
    
        radrot_ellipse = np.sqrt( pow((xrot)/a_ellipse,2) + pow((yrot)/b_ellipse,2) )
        radrot_empty_ellipse = np.sqrt( pow((xrot)/a_empty_ellipse,2) + pow((yrot)/b_empty_ellipse,2) )
    
        # In/Out conditions
    
        points_in_ellipse = ( radrot_ellipse <= 1)
        counts_in_ellipse = np.array(points_in_ellipse, float)
    
        points_in_empty_ellipse = ( radrot_empty_ellipse < 1)
        counts_in_empty_ellipse = np.array(points_in_empty_ellipse, float)
        
        # Opacity
        
        radius_central = rad_central_vector[ring]
        g_decay = gamma_decay_vector[ring]
    
        counts_in_ellipse[counts_in_ellipse==1] = initial_opacity* transits.opacity_decay_powlaw(radrot_ellipse[counts_in_ellipse==1],radius_central,g_decay)
        
        # flux
    
        N_in_ellipse = counts_in_ellipse[points_in_empty_ellipse==False].sum()
        
        relflux = relflux - N_in_ellipse/Nptotal
        
        
    return relflux
    
    """
    
"""
def relflux_model_allt_nrings_powlaw(time_steps,rand_points_in_circle,alpha):

    T0,yp,velocity,theta_angle,radius_planet,e_planet,e_ellipse,initial_opacity,ext_ring_radius,int_ring_radius,number_rings,gamma_decay = alpha

    Np = len(time_steps)
    
    relflux=np.empty(Np)
    
    for step in np.arange(Np):
        
        relflux[step] = relflux_model_singlet_nrings_powlaw(time_steps[step],rand_points_in_circle,T0,yp,velocity,theta_angle,radius_planet,e_planet,e_ellipse,initial_opacity,ext_ring_radius,int_ring_radius,number_rings,gamma_decay)
    
    return relflux
        
"""



    
