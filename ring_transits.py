#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:30:37 2017

@author: pablo

Main library of functions for the ring-transits repository

Most code comes from my original libtransits.py file, using the 'Simpons-2'
approach to integration.
"""


import numpy as np

import scipy.integrate as integrate

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# FUNCTIONS TO DESCRIBE THE LIMB DARKENING
def _unnorm_radial_intensity_limb(r_squared, limb_params=[0.371,  0.8674,
                                                          -0.7499,  0.2507]):
    """
    Un-normalised limb darkening function for the non-linear [Claret (2000)]
    model.
    The default values of the parameters are those corresponding to Tabby star,
    from the ATLAS tables in Claret (2000).
    Takes as input value r_squared=r^2
    """

    mu_param = np.sqrt(1. - r_squared)

    c1, c2, c3, c4 = limb_params

    intensity = \
        1 - c1*(1 - mu_param**(1./2.)) - c2*(1 - mu_param) \
        - c3*(1 - mu_param**(3./2.)) - c4*(1 - mu_param**2)

    return intensity


def radial_intensity_limb(r_squared, limb_params=[0.371,  0.8674,
                                                  -0.7499,  0.2507]):
    """
    Function that returns the normalised (total flux of star = 1) radial
    intensity function for the non-linear limb darkening model of
    Claret (2000).
    Simpler models (e.g. linear, quadratic) can be recovered from here by
    setting some of the model parameters to 0.

    The default values of the parameters are those corresponding to Tabby star,
    from the ATLAS tables in Claret (2000).
    Takes as input value r_squared=r^2.
    """

    r_squared = np.atleast_1d(r_squared)
    do_calc = (r_squared <= 1)

    c1, c2, c3, c4 = limb_params

    # Normalisation obtained from the Mathematica calculations
    I0 = 210./(np.pi*(210. - 42.*c1 - 70.*c2 - 90.*c3 - 105.*c4))

    un_norm_flux = _unnorm_radial_intensity_limb(r_squared[do_calc],
                                                 limb_params)

    I_out = np.zeros(r_squared.shape)
    I_out[do_calc] = I0 * un_norm_flux

    return I_out


def relflux_onlyplanet(time, T0, velocity, yp, radius_planet,
                       limb_params=[0.371,  0.8674, -0.7499,  0.2507],
                       Nx=51, Ny=51):
    """
    Function to get the relative flux at a given time for a simple
    model with a single component: a planet (i.e. disc).
    We use here a Simpson's integration over the disc of the *planet*.
    Nx and Ny define the grid to use.
    limb_params define the parameters to use to model the limb darkening
    of the star. The default values correspond to Tabby star following
    Claret (2000)

    Parameters:
    * T0: central time of transit [days]
    * velocity: velocity of transit [stellar radii/day]
    * yp: impact parameter [stellar radii]
    * radius_planet: radius of planet/disc [stellar radii]

    """

    time = np.atleast_1d(time)
    xp = velocity*(time - T0)

    # Create grid of xrel, yrel values
    xrel = np.linspace(-radius_planet, radius_planet, Nx)
    yrel = np.linspace(-radius_planet, radius_planet, Ny)
    xstep = 2*radius_planet/float(Nx - 1)
    ystep = 2*radius_planet/float(Ny - 1)

    # Assign opacities to points in the grid
    # (Shape: Nx,Ny)
    rel_radius_2 = np.add.outer(xrel**2, yrel**2)
    opacity_rel = np.zeros(rel_radius_2.shape)
    opacity_rel[rel_radius_2 <= radius_planet**2] = 1.0

    # Now, calculate the corresponding points in the axes centered at the star
    xvals = np.add.outer(xp, xrel)  # Shape (Nt, Nx)
    yvals = yp + yrel               # Shape (Ny)

    # And compute the corresponding radial distances to the star
    # and relative fluxes (default values for limb darkening)
    dist_2_to_star_centre = np.add.outer(xvals**2, yvals**2)  # Shape (Nt, Nx, Ny)
    rel_intensity = radial_intensity_limb(r_squared=dist_2_to_star_centre,
                                          limb_params=limb_params)   # Shape (Nt, Nx, Ny)

    # Finally, compute the relative intensity blocked by each pixel
    # at each time step, and integrate to obtain the total blocked line
    # at each time
    blocked_flux_pixels = rel_intensity*opacity_rel

    total_blocked_flux = \
        integrate.simps(integrate.simps(blocked_flux_pixels, dx=ystep),
                        dx=xstep)

    return 1. - total_blocked_flux


def relflux_planet_ellipse_model(time, T0, velocity, yp, rad_planet,
                                 theta_angle, e_ellipse,
                                 opacity_law, opacity_param_1,
                                 opacity_param_2,
                                 limb_params=[0.371, 0.8674, -0.7499, 0.2507],
                                 opacity_cut=0.001, Nx=51, Ny=51):
    """
    Function to get the relative flux at a given time for a model which
    includes a planet (disc) and the corresponding ring system (ellipse).
    The opacity for the ellipse can follow different laws, depending on the
    parameter opacity_law:
        * opacity_law=1: constant opacity inside an ellipse of fixed
                         size
        * opacity_law=2: opacity=1 inside a 'central' ellipse, and exponential
                         decay outside of it
        * opacity_law=3: opacity=1 inside a 'central' ellipse, and power-law
                         decay outside of it
        * opacity_law=4: opacity decays from centre following a Gaussian
                         function.

    We use here a Simpson's integration over the blocking object, up to the point
    where the opacity==opacity_cut (except for opacity_law=1, where the limit
    is set by the size of the ellipse). Nx and Ny define the grid to use.
    If Ny=0, it will be set so that the spacing of the grid is the same
    in the x and y directions, or Ny=51 (whichever sets a larger Ny).
    limb_params define the parameters to use to model the limb darkening
    of the star. The default values correspond to Tabby star following
    Claret (2000)


    Parameters:
    * T0: central time of transit [days]
    * velocity: velocity of transit [stellar radii/day]
    * yp: impact parameter [stellar radii]
    * rad_planet: radius of the planet [stellar radii]
    * theta_angle: orientation angle of the ellipses, measured from the
        direction of transit [radians]
    * e_ellipse: ellipticity of both ellipses

    * opacity_param_1, opacity_param_2: have different meanings depending on
      the opacity_law used:
          - opacity_law=1:
              + opacity_param_1: opacity of the ellipse
              + opacity_param_2: semi-major axis of the ellipse [stellar radii]
          - opacity_law=2:
              + opacity_param_1: semi-major axis of the central ellipse
                                 [stellar radii]
              + opacity_param_2: characteristic scale of the exponential decay
                                 along the semi-major axis [stellar radii]
          - opacity_law=3:
              + opacity_param_1: semi-major axis of the central ellipse
                                 [stellar radii]
              + opacity_param_2: exponent of the power law decay
          - opacity_law=4:
              + opacity_param_1: opacity at the centre of the ellipse (r=0),
                                 which sets the normalisation of the Gaussian
                                 function
              + opacity_param_2: standard deviation of the Gaussian function
                                 [stellar radii]

    NOTES:
        Actual radial dependence of the different opacity laws (always along
        the semi-major axis):

        * opacity_law = 1:

            op0 = opacity_param_1
            r0 = opacity_param_2

            opac(r) = {op0   if r <= r0;
                       0     if r >  r0}

        * opacity_law = 2:

            r_cent = opacity_param_1
            r_exp = opacity_param_2

            opac(r) = {1                          if  r <= r_cent;
                       Exp[-(r - r_cent)/r_exp]   if  r >  r_cent}

        * opacity_law = 3:

            r_cent = opacity_param_1
            gamma = opacity_param_2

            opac(r) = {1                     if  r <= r_cent;
                       (r/r_cent)^(-gamma)   if r >  r_cent}

        * opacity_law = 4:

            op_cent = opacity_param_1
            sigma = opacity_param_2

            opac(r) = op_cent*exp(- r^2 / (2*sigma^2))
    """

    assert opacity_law in [1, 2, 3, 4]

    time = np.atleast_1d(time)

    xp = velocity*(time - T0)

    # Define the limits over which to generate our grid of points
    # for integration
    #
    # First, obtain, depending on the opacity_law, the maximum extension
    # of the ellipse in both directions
    if opacity_law == 1:
        a_max = opacity_param_2
    elif opacity_law == 2:
        a_max = opacity_param_1 - opacity_param_2*np.log(opacity_cut)
    elif opacity_law == 3:
        a_max = opacity_param_1*pow(opacity_cut, -1./opacity_param_2)
    elif opacity_law == 4:
        a_max = np.sqrt(-2*(opacity_param_2**2) *
                        np.log(opacity_cut/opacity_param_1))

    b_max = a_max*np.sqrt(1 - e_ellipse**2)

    # Now, set the region to integrate taking into account both the ellipse
    # and the planet, and generate the grid

    xlim = max(rad_planet, a_max)
    ylim = max(rad_planet, b_max)

    # For the case Ny=0, need to re-define it
    if Ny == 0:
        Ny_rel = int(Nx*ylim/xlim)
        Ny = max(Ny_rel, 51)

    xrel = np.linspace(-xlim, xlim, Nx)
    yrel = np.linspace(-ylim, ylim, Ny)
    xstep = 2*xlim/float(Nx - 1)
    ystep = 2*ylim/float(Ny - 1)

    # Now, for the grid of points, compute the opacities
    # *due to the ellipse (rings)*, following the given opacity law
    # (Shape: Nx,Ny)
    rad_ell_values = np.sqrt(np.add.outer(xrel**2, (yrel**2)/(1-e_ellipse**2)))

    opacity_rel = np.zeros(rad_ell_values.shape)
    sel_for_calc = (rad_ell_values <= a_max)

    if opacity_law == 1:
        opacity_rel[sel_for_calc] = opacity_param_1
    elif opacity_law == 2:
        opacity_rel[sel_for_calc] =  \
            np.exp(-(rad_ell_values[sel_for_calc] - opacity_param_1) /
                   opacity_param_2)
        if opacity_param_1 > 0:
            opacity_rel[rad_ell_values <= opacity_param_1] = 1.0
    elif opacity_law == 3:
        opacity_rel[sel_for_calc] = \
            (rad_ell_values[sel_for_calc]/opacity_param_1)**(-opacity_param_2)
        opacity_rel[rad_ell_values <= opacity_param_1] = 1.0
    elif opacity_law == 4:
        opacity_rel[sel_for_calc] = \
            opacity_param_1*np.exp(-(rad_ell_values[sel_for_calc]**2) /
                                   (2*(opacity_param_2**2)))


    # Now, if needed, compute the opacities
    # *due to the planet*
    if rad_planet > 0:
        rel_radius_2 = np.add.outer(xrel**2, yrel**2)
        opacity_rel[rel_radius_2 <= rad_planet**2] = 1.0

    # Now, calculate the corresponding points in the axes centered
    # at the star
    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)

    # First, 'undo' the rotation
    x_i = np.add.outer(xrel*ct, -yrel*st)   # (Nx, Ny)
    y_i = np.add.outer(xrel*st, yrel*ct)  # (Nx, Ny)

    # And now, get the final values
    xvals = np.add.outer(xp, x_i)     # (Nt, Nx, Ny)
    yvals = yp + y_i                  # (Nx, Ny)

    # And compute the corresponding radial distances to the star
    # and relative fluxes  (default values for limb darkening)
    dist_2_to_star_centre = (xvals**2 + yvals**2)               # Shape (Nt, Nx, Ny)
    rel_intensity = radial_intensity_limb(r_squared=dist_2_to_star_centre,
                                          limb_params=limb_params)   # Shape (Nt, Nx, Ny)

    # Finally, compute the relative intensity blocked by each pixel
    # at each time step, and integrate to obtain the total blocked line
    # at each time
    blocked_flux_pixels = rel_intensity*opacity_rel

    total_blocked_flux = \
        integrate.simps(integrate.simps(blocked_flux_pixels, dx=ystep),
                        dx=xstep)

    return 1. - total_blocked_flux


def transit_diagram_planet_ellipse(yp, rad_planet, theta_angle, e_ellipse,
                                   opacity_law, opacity_param_1,
                                   opacity_param_2,
                                   opacity_step=0.1, star_color='lightcoral',
                                   planet_color='skyblue',
                                   rings_color='black',
                                   ax=None):
    """
    Function to draw a diagram showing the configuration of the system at the
    time of inferior conjunction.
    """

    if ax is None:
        fig, ax = plt.subplots()

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
        opacity_levels = np.arange(opacity_step, 1.0 + opacity_step,
                                   opacity_step)
        for this_level in opacity_levels:

            # First, obtain semi-major axis of ellipse
            if opacity_law == 2:
                a_level = opacity_param_1 - opacity_param_2*np.log(this_level)
            elif opacity_law == 3:
                a_level = opacity_param_1*pow(this_level, -1./opacity_param_2)
            elif opacity_law == 4:
                a_level = np.sqrt(-2*(opacity_param_2**2) *
                                  np.log(this_level/opacity_param_1))

            # Now, get corresponding semi-minor axis and draw ellipse
            b_level = a_level*np.sqrt(1. - e_ellipse*e_ellipse)

            ellipse = mpatches.Ellipse([0, yp], 2*a_level, 2*b_level,
                                       angle=(theta_angle*180./np.pi),
                                       ec='black', fc=rings_color,
                                       alpha=opacity_step)
            ax.add_patch(ellipse)

    # Finally, draw the circle corresponding to the planet
    circle_planet = mpatches.Circle([0, yp], rad_planet, ec='black',
                                    fc=planet_color, alpha=1)
    ax.add_patch(circle_planet)

    # And set some global properties of the plot
    ax.set_xlim(-4, 4)
    ax.set_ylim(-2, 2)
    ax.set_aspect('equal')
    ax.set_xlabel("X (stellar radii)")
    ax.set_ylabel("Y (stellar radii)")

    return ax
