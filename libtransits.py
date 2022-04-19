# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 13:05:12 2015

Just a library with some functtions I will use for the modelling of transit
light curves


@author: pablo
"""


import numpy as np
import numpy.random as rd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

rd.seed(314)

def radial_bright_func_limb(r, radius=1.0, limb_constant_u=0.56):

    return (1. - limb_constant_u) + \
        limb_constant_u*np.sqrt(1. - ((r/radius)**2.))

def generate_randpoints_in_circle_limbdark(N, radius=1.0,
                                           limb_constant_u=0.56):
    """
    Follow http://astrowww.phys.uvic.ca/~tatum/stellatm/atm6.pdf to
    introduce the limb darkening.
    """
    xi1, xi2 = rd.rand(2, N)
    # Get r-coordinate from xi1
    r = radius*np.sqrt(xi1)
    # Get phi coordinate (uniform in 0-2pi) from xi2
    phi = 2*np.pi*xi2

    xi3 = rd.rand(N)

    selection = (xi3 <
                 radial_bright_func_limb(r=r, radius=radius,
                                         limb_constant_u=limb_constant_u))

    #print("Selected a fraction of %g points from total") % \
    #   (selection.sum()/float(N))

    xout = r[selection]*np.sin(phi[selection])
    yout = r[selection]*np.cos(phi[selection])

    return np.array((xout, yout)).T


def relflux_model_singlet_onlyplanet(time, T0, velocity, yp,
                                           radius_planet,
                                           rand_points_in_circle):
    """
    Function to get the relative flux at a given time for a simple
    model with a single component: a planet (i.e. disc)

    Parameters:
    * T0: central time of transit [days]
    * velocity: velocity of transit [stellar radii/day]
    * yp: impact parameter [stellar radii]
    * radius_planet: radius of planet/disc [stellar radii]
    """
    Nptotal = len(rand_points_in_circle)

    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp
    

    condition_circle = ((xrel*xrel) + (yrel*yrel) > (radius_planet*radius_planet))

    relflux = condition_circle.sum()/float(Nptotal)
    return relflux


def relflux_model_pred_allt_onlyplanet(time_steps, rand_points_in_circle,
                                       alpha):
    T0, velocity, yp, radius_planet = alpha

    Nt = len(time_steps)
    rflux = np.empty(Nt, float)
    for i in range(Nt):
        rflux[i] = relflux_model_singlet_onlyplanet(time_steps[i], T0,
                                                          velocity, yp,
                                                          radius_planet,
                                                          rand_points_in_circle)
    return rflux


def relflux_model_singlet_opacdisc(time, T0, velocity, yp,
                                         radius_planet, theta_angle, a_ellipse,
                                         e_ellipse, opacity_ellipse,
                                         rand_points_in_circle):
    """
    Function to get the relative flux at a given time for a model with two
    components: a planet (i.e. opaque disc), and an ellipse with fixed
    opacity.

    Parameters:
    * T0: central time of transit [days]
    * velocity: velocity of transit [stellar radii/day]
    * yp: impact parameter [stellar radii]
    * radius_planet: radius of planet/disc [stellar radii]
    * theta_angle: orientation angle of the ellipse [radians]
    * a_ellipse: semi-major axis of ellipse [stellar radii]
    * opacity ellipse: opacity of the ellipse, which is fixed for all of it.
        Must be between 0 and 1 (the code does not check for this!)
    """

    if 0 <= opacity_ellipse <= 1:
        opacity_ellipse=opacity_ellipse 
    else: 
        print("Please select an ellipse opacity between 0 and 1")

    Nptotal = len(rand_points_in_circle)

    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    condition_circle = ((xrel*xrel) + (yrel*yrel) > (radius_planet*radius_planet))

    b_ellipse = a_ellipse*np.sqrt(1. - e_ellipse*e_ellipse)
    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)

    condition_ellipse = ((pow((xrel*ct + yrel*st)/a_ellipse, 2) +
                         pow((yrel*ct - xrel*st)/b_ellipse, 2)) > 1)

    cond_both = condition_circle*condition_ellipse

    N_out_both = cond_both.sum()

    N_in_circle_and_ellipse = ((condition_ellipse == False)*
                               (condition_circle == False)).sum()

    N_in_ellipse = (condition_ellipse == False).sum()

    Nout = N_out_both + (1. - opacity_ellipse) * \
        (N_in_ellipse - N_in_circle_and_ellipse)

    relflux = Nout/float(Nptotal)
    return relflux


def relflux_model_pred_allt_opacdisc(time_steps, rand_points_in_circle, alpha):

    T0, velocity, yp, radius_planet, theta_angle, a_ellipse, e_ellipse, opacity_ellipse = alpha

    Nt = len(time_steps)
    rflux = np.empty(Nt, float)
    for i in range(Nt):
        rflux[i] = relflux_model_singlet_opacdisc(time_steps[i],
                                                        T0, velocity, yp,
                                                        radius_planet,
                                                        theta_angle, a_ellipse,
                                                        e_ellipse,
                                                        opacity_ellipse,
                                                        rand_points_in_circle)
    return rflux


def relflux_model_singlet_opacdisc_exp_fix(time, T0, velocity, yp,
                                                 radius_planet, theta_angle,
                                                 a_ellipse, e_ellipse,
                                                 rand_points_in_circle,
                                                 n_efolds=5.0):

    # In this case, opacity decreases exponentially linearly along the
    # major semi-axis of the ellipse,
    # It is =opacity_ellipse/e at xrot=a/n_efolds (by default, n_efolds=5),
    # and =opacity_ellipse at xrot=0
    # (the center)
    # In this version, opacity_ellipse is set to 1
    opacity_ellipse = 1.0

    Nptotal = len(rand_points_in_circle)

    xp = velocity*(time - T0)

    xrel = rand_points_in_circle[:, 0] - xp
    yrel = rand_points_in_circle[:, 1] - yp

    b_ellipse = a_ellipse*np.sqrt(1. - e_ellipse*e_ellipse)
    ct = np.cos(theta_angle)
    st = np.sin(theta_angle)
    xrot = xrel*ct + yrel*st
    yrot = yrel*ct - xrel*st
    radrot_ellipse = np.sqrt((xrot/a_ellipse)**2. + (yrot/b_ellipse)**2.)

    points_in_circle = ((xrel*xrel) + (yrel*yrel) <= (radius_planet*radius_planet))

    points_in_ellipse = (radrot_ellipse <= 1)

    # Make points in ellipse count only as 'opacity (function of xrot)'
    # opacity_points = opacity_ellipse*np.exp(-5.*radrot_ellipse)
    # counts_in_ellipse = np.array(points_in_ellipse, float)*opacity_points

    counts_in_ellipse = np.array(points_in_ellipse, float)
    counts_in_ellipse[counts_in_ellipse == 1] = opacity_ellipse*np.exp(-n_efolds * radrot_ellipse[counts_in_ellipse == 1])

    # Count no. of points taken out by the circle (planet)
    N_in_circle = points_in_circle.sum()

    # Count no. of points taken out by the ellipse, but do not double-count
    # those already taken out by circle
    N_in_ellipse = counts_in_ellipse[points_in_circle == False].sum()

    # Count total no. of points taken out by the model
    N_out_total = N_in_circle + N_in_ellipse

    # And derive the relative flux
    relflux = 1. - (N_out_total/Nptotal)

    return relflux
    

# Very inefficient!!
def relflux_model_pred_allt_opacdisc_exp_fix(time_steps, rand_points_in_circle,
                                             alpha, n_efolds=5.0):
    T0, velocity, yp, radius_planet, theta_angle, a_ellipse, e_ellipse = alpha

    Nt = len(time_steps)
    rflux = np.empty(Nt, float)
    for i in range(Nt):
        rflux[i] = relflux_model_singlet_opacdisc_exp_fix(time_steps[i],T0, velocity,yp, radius_planet,
        theta_angle,a_ellipse, e_ellipse,rand_points_in_circle,n_efolds=n_efolds)
    return rflux


def opacity_decay_exp(r, radius_central, radius_decay):
    opac = np.exp(-(r-radius_central)/radius_decay)
    opac[r < radius_central] = 1
    return opac


def relflux_model_singlet_opac_exp_central(time, T0, velocity, yp,
                                                 rad_central, rad_exp_decay,
                                                 theta_angle, e_ellipse,
                                                 rand_points_in_circle):
    """
    Function to get the relative flux for a given time and for a set of
    model parameters.

    In this case the model describes a elliptical disc with two 'components':
    * A central ellipse with opacity=1
    * A disc surrounding this ellipse with opacity decreasing exponentially,
      from the edge of the central ellipse

    Both components share the same orientation angle and ellipticity.

    Parameters of the model:

    * T0: central time of transit [days]
    * velocity: velocity of transit [stellar radii/day]
    * yp: impact parameter [stellar radii]
    * rad_central: semi-major axis of the central ellipse [stellar radii]
    * rad_exp_decay: characteristic scale of the exponential decay
        along the semi-major axis. For r = rad_central+rad_exp_decay,
        the opacity will be equal to 1/e. We'll do our integration up to
        r=rad_central + 5*rad_exp_decay. [stellar radii]
    * theta_angle: orientation angle of the ellipses, measured from the
        direction of transit [radians]
    * e_ellipse: ellipticity of both ellipses

    Other parameters:
    * time: time at which we calculate the flux [days]
    * rad_points_in_circle: set of random points distributed following the
        brightness profile of the star [x/y positions in stellar radii]
    """

    Nptotal = len(rand_points_in_circle)

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

    # Define the maximum ellipse we will consider
    a_max = rad_central + (5.*rad_exp_decay)

    # First, define the points which are inside the maximum ellipse considered
    # los que cumplen la condición, se les pone la etiqueta True. Los que no, False.
    points_in_ellipse = (radrot_ellipse <= a_max)

    # And now, get the 'counts' taking into account the changing opacity
    counts_in_ellipse = np.array(points_in_ellipse, float)

    # Cuando seleccionamos ==1, solo cogemos los True.
    counts_in_ellipse[counts_in_ellipse == 1] = opacity_decay_exp(radrot_ellipse[counts_in_ellipse == 1], rad_central, rad_exp_decay)


    # Get the total count and calculate the flux
    N_in_ellipse = counts_in_ellipse.sum()

    relflux = 1. - (N_in_ellipse/Nptotal)

    return relflux


def relflux_model_pred_allt_opac_exp_central(time_steps, rand_points_in_circle,
                                             alpha):

    T0, velocity, yp, rad_central, rad_exp_decay,theta_angle, e_ellipse = alpha

    Nt = len(time_steps)
    rflux = np.empty(Nt, float)
    for i in range(Nt):
        rflux[i] = \
            relative_flux_model_singlet_opac_exp_central(time_steps[i],
                                                         T0, velocity, yp,
                                                         rad_central,
                                                         rad_exp_decay,
                                                         theta_angle,
                                                         e_ellipse,
                                                         rand_points_in_circle)
    return rflux


def opacity_decay_powlaw(r, radius_central, gamma_decay):
    opac = pow(r/radius_central, -gamma_decay)
    opac[r < radius_central] = 1
    return opac


def relflux_model_singlet_powerlaw_central(time, T0, velocity, yp,
                                                 rad_central, gamma_decay,
                                                 theta_angle, e_ellipse,
                                                 rand_points_in_circle,
                                                 opac_lim=0.001):
    """
    Function to get the relative flux for a given time and for a set of
    model parameters.

    In this case the model describes a elliptical disc with two 'components':
    * A central ellipse with opacity=1
    * A disc surrounding this ellipse with opacity decreasing following a
      power law.
    Mathematically, we define this model using two parameters, as

    opac(x) = {1                            if  r <= rad_central,
               (x/rad_central)^-gamma_decay if  r>rad_centreal }

    We do the integration up to the point where opac=opac_lim (default 0.001).

    Both components share the same orientation angle and ellipticity.

    Parameters of the model:

    * T0: central time of transit [days]
    * velocity: velocity of transit [stellar radii/day]
    * yp: impact parameter [stellar radii]
    * rad_central: semi-major axis of the central ellipse [stellar radii]
    * gamma_decay: exponent of the power law decay, as described above
    * theta_angle: orientation angle of the ellipses, measured from the
        direction of transit [radians]
    * e_ellipse: ellipticity of both ellipses

    Other parameters:
    * time: time at which we calculate the flux [days]
    * rad_points_in_circle: set of random points distributed following the
        brightness profile of the star [x/y positions in stellar radii]
    """

    Nptotal = len(rand_points_in_circle)

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

    # Define the maximum ellipse we will consider
    # a_max = rad_central + (5.*rad_exp_decay)
    a_max = rad_central*pow(opac_lim, -1./gamma_decay)

    # First, define the points which are inside the maximum ellipse considered
    points_in_ellipse = (radrot_ellipse <= a_max)

    # And now, get the 'counts' taking into account the changing opacity
    counts_in_ellipse = np.array(points_in_ellipse, float)

    counts_in_ellipse[counts_in_ellipse == 1] = opacity_decay_powlaw(radrot_ellipse[counts_in_ellipse == 1],rad_central, gamma_decay)

    # Get the total count and calculate the flux
    N_in_ellipse = counts_in_ellipse.sum()

    relflux = 1. - (N_in_ellipse/Nptotal)

    return relflux


def relflux_model_pred_allt_powerlaw_central(time_steps, rand_points_in_circle,
                                             alpha, opac_lim=0.001):

    T0, velocity, yp, rad_central, gamma_decay, theta_angle, e_ellipse = alpha

    Nt = len(time_steps)
    rflux = np.empty(Nt, float)
    for i in range(Nt):
        rflux[i] =   relflux_model_singlet_powerlaw_central(time_steps[i],
                                                         T0, velocity, yp,
                                                         rad_central,
                                                         gamma_decay,
                                                         theta_angle,
                                                         e_ellipse,
                                                         rand_points_in_circle,
                                                         opac_lim)

    return rflux


def geometry_plot(yp=0.0, rad_disc=0.1, theta_ellipse=0.0, a_ellipse=0.2,
                  e_ellipse=0.5):
    """
    Function to make a plot representing the geometry of the transit for
    the case of two components: a disc and an ellipse with the same center.
    We represent the geometry at the central point of the transit.

    Parameters:
    * yp: impact factor
    * rad_disc: radius of the central disc
    * theta_ellipse: orientation angle of the ellipse [radians]
    * a_ellipse: semi-major axis of the ellipse
    * e_ellipse: ellipticity of the ellipse
    """

    fig, ax = plt.subplots(figsize=(8, 8))

    b_ellipse = a_ellipse*np.sqrt(1. - e_ellipse*e_ellipse)

    circle_star = mpatches.Circle([0, 0], 1, ec="black", fc='red', alpha=0.5)
    ax.add_patch(circle_star)

    circle_planet = mpatches.Circle([0, yp], rad_disc, ec='black', fc='blue',
                                    alpha=0.5)
    ax.add_patch(circle_planet)

    ellipse_disc = mpatches.Ellipse([0, yp], 2*a_ellipse, 2*b_ellipse,
                                    angle=(theta_ellipse*180./np.pi),
                                    ec='black', fc='green', alpha=0.2)
    ax.add_patch(ellipse_disc)

    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.grid(b=True)
    ax.set_aspect('equal')
    ax.set_xlabel("X (stellar radii)")
    ax.set_ylabel("Y (stellar radii)")
    plt.show()
