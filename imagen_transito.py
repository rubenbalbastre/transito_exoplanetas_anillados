"""
import h5py
import numpy as np

archivo = h5py.File("/Users/rubenbalbastre/Desktop/mcmc_prueba.h5","r")
cadena = np.array(archivo['mcmc']['chain'])
log_prob = np.argmax(np.exp(np.array(archivo['mcmc']['log_prob'])))


import emcee
reader = emcee.backends.HDFBackend("/Users/rubenbalbastre/Desktop/mcmc_prueba.h5")
# Mostrar mejor ajuste
nwalkers, ndim = reader.shape
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

mpl.rcParams["text.usetex"] = False
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'




T0 =  792.7322779935188 
velocity =  5.722180392461998  
yp =  -0.6501855195751064 
radius_planet =  0.290880718974143 
planet_separation =  0.10643083452085439 
rad_central =  0.8502320238645875 
gamma_decay =  1.5929497977878695 
theta_angle =  0.05817662612225648 
e_ellipse =  0.9979368067666516 
e_planet =  0.28277742729832883 
initial_opacity =  0.9952465189487167
rlim = 10

def transit_diagram_planet_ellipse(yp, rad_planet, theta_angle, e_ellipse,
                                   opacity_law, opacity_param_1,
                                   opacity_param_2,
                                   opacity_step=0.025, star_color='orange',
                                   planet_color='skyblue',
                                   rings_color='black',
                                   ax=None):
    """
    Function to draw a diagram showing the configuration of the system at the
    time of inferior conjunction.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(10,10))

    # First, draw the circles corresponding to the star
    circle_star = mpatches.Circle([0, 0], 1, ec="black", fc=star_color,
                                  alpha=1.0)
    ax.add_patch(circle_star)

    # For model 1, draw a single ellipse with alpha set to its opacity
    if opacity_law == 1:
        a_ellipse = opacity_param_2
        opac = opacity_param_1
        b_ellipse = a_ellipse*np.sqrt(1. - e_ellipse*e_ellipse)

        ellipse = mpatches.Ellipse([0, yp], 2*a_ellipse, 2*b_ellipse,
                                   angle=(theta_angle*180./np.pi),
                                   ec='black', fc=rings_color, alpha=opac)
        ax.add_patch(ellipse)

    # In the other cases, draw the ellipses corresponding to the
    # required levels
    else:
        opacity_levels = np.arange(opacity_step, 1.0 + opacity_step, opacity_step)
        for this_level in opacity_levels:

            # First, obtain semi-major axis of ellipse
            if opacity_law == 2:
                a_level = opacity_param_1 - opacity_param_2*np.log(this_level)
            elif opacity_law == 3:
                a_level = opacity_param_1*pow(this_level, -1./opacity_param_2)
            elif opacity_law == 4:
                a_level = np.sqrt(-2*(opacity_param_2**2) *
                                  np.log(this_level/opacity_param_1))
            elif opacity_law == 5:
                a_level = opacity_param_1*pow(this_level, -1./opacity_param_2)

            # Now, get corresponding semi-minor axis and draw ellipse
            b_level = a_level*np.sqrt(1. - e_ellipse*e_ellipse)

            ellipse = mpatches.Ellipse([0, yp], 2*a_level, 2*b_level,
                                       angle=(theta_angle*180./np.pi),
                                       ec='black', fc=rings_color,
                                       alpha=opacity_step)
            ax.add_patch(ellipse)
        
        ellipse = mpatches.Ellipse([0, yp], 2*rad_central, 2*rad_central*np.sqrt(1. - e_ellipse*e_ellipse),
                                       angle=(theta_angle*180./np.pi),
                                       ec='black', fc=rings_color,
                                       alpha=1)
        ax.add_patch(ellipse)

        # add empty ellipse
        radius_empty_ellipse = radius_planet +  planet_separation
        ellipse = mpatches.Ellipse([0, yp], 2*radius_empty_ellipse, 2*radius_empty_ellipse*np.sqrt(1. - e_ellipse*e_ellipse),
                                       angle=(theta_angle*180./np.pi),ec='black', fc=star_color)
        ax.add_patch(ellipse)

    e_planet =  0.2522929758800985 
    # Finally, draw the circle corresponding to the planet
    ellipse_planet = mpatches.Ellipse([0,yp], 2*rad_planet,2*rad_planet*np.sqrt(1-e_planet**2),
                                   angle=(theta_angle*180./np.pi), ec="black", fc=planet_color,
                                  alpha=1.0)
    ax.add_patch(ellipse_planet)

    # And set some global properties of the plot
    ax.set_xlim(-rlim, rlim)
    ax.set_ylim(-2, 2)
    ax.set_aspect('equal')
    ax.set_xlabel(r" x [$R_{\star}$]")
    ax.set_ylabel(r" y [$R_{\star}$]")
    ax.set_axis_off()

    fig.savefig("/Users/rubenbalbastre/Desktop/Física/CUARTO AÑO/2do cuatrimestre/TFG/Presentación/protoplanet_representation.pdf")

    return ax


transit_diagram_planet_ellipse(yp, radius_planet, theta_angle, e_ellipse, opacity_law=3, opacity_param_1=rad_central, opacity_param_2=gamma_decay)
plt.show()
