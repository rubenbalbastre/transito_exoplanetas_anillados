#################################################################################
####### THIS SCRIPT CONTAINS SEVERAL MODELS FOR THE OPACITY FUNCTION ######
#################################################################################


# IMPORT LIBRARYS needed to compute the models

import numpy as np

""" OPACITY MODELS """

"""Gaussian distribution of the ring opacity NON-USED"""

def gauss(x,a_ellipse,a_empty_ellipse):

    mu = (a_ellipse - a_empty_ellipse)/2
    
    sigma = mu  # 30% of the mean value \mu

    return 1./np.sqrt(2*np.pi) / sigma * np.exp(- 1/2 *((x-mu)/sigma)**2 )

""" constant opacity model. 4 rings """

def opacity_cst(xdata,radius_empty_ellipse,radius1,radius2,radius3,radius4,val1,val2,val3,val4):
    opac = np.empty(len(xdata))
    for n in np.arange(len(xdata)):
        r = xdata[n]
        if r < radius_empty_ellipse:
            opac[n] = 0
        elif radius_empty_ellipse < r < radius1:
            opac[n] = val1
        elif radius1 < r < radius2:
            opac[n] = val2
        elif radius2 < r < radius3:
            opac[n] = val3
        elif radius3 < r < radius4:
            opac[n] = val4
        elif r > radius4:
            opac[n] = 0
    return opac

""" Power-law of libtransits.py 

 This model consists on a piecewise function defined as:
 O(r) = (r/r0)^(-gamma_decay)  if r > r0
      = 1                      if r < r0
"""
def opacity_powlaw(xdata,r0,gammaDecay):

    l=len(xdata)
    opac = np.empty(l)

    for n in np.arange(l):
        r=xdata[n]
        opac[n] = pow(r/r0,-gammaDecay)
        if r < r0:
            opac[n] = 1
    return opac
    

""" Opacity power-law * exp 

 This model consists on a piecewise function defined as:
 O(r) = (r/radius_central)^(-gamma_decay) * exp(-r*gamma_decay_2)  if r > radius_limit
      = 1                                                          if r < radius_limit
 The coefficent gamma_decay_2 is obtanid through the continuity conditions of the opacity
"""

def opacity_powlawexp (r,radius_central,radius_limit,gamma_decay):


    gamma_decay_2 = gamma_decay / radius_limit * np.log(radius_central/radius_limit) # continuity condition
    opac = pow(r/radius_central,-gamma_decay) * np.exp(-gamma_decay_2*r)
    opac[r<radius_limit] = 1
    return opac


""" Opacity Pow-law-sum 

"""

def opacity_powlawsum(r,radius_central_1,radius_central_2,gamma_decay_1,radius_limit):

    # Continuity & sense condition 
    if radius_central_2>=radius_limit:
        raise TypeError("radius_central_2 > radius_limit has no sense")
    elif radius_limit<=radius_central_1:
        raise TypeError("Non-sense")

    gamma_decay_2 = np.log(1 - (radius_limit/radius_central_1)**(-gamma_decay_1)) / np.log(radius_central_2/radius_limit)

    opac = (r/radius_central_1)**(-gamma_decay_1) + (r/radius_central_2)**(-gamma_decay_2)
        
    opac[r<radius_limit] = 1
    
    if gamma_decay_2<0:
        raise TypeError("gamma_decay_2 <0 ")

    for n in np.arange(len(r)):
        if  opac[n] < 0:
            raise TypeError("opacity<0 has no sense")
        elif opac[n] > 1:
            raise TypeError("opacity>1 has no sense")
    return opac

""" Opacity power-law and exp. 3 rings

"""


def opacity_powlaw_and_exp(r,radius_central_1,gamma_decay_1,radius_central_2):

    opac = pow(r/radius_central_1,-gamma_decay_1)
    opac[r<radius_central_1] = 1
    
    # continuity condition. It is easier to get gamma_decay_2 from radius_central_2 than otherwise.
    gamma_decay_2 = gamma_decay_1/radius_central_2 * np.log(radius_central_2/radius_central_1)
    
    for n in np.arange(len(r)): 
        if r[n]>radius_central_2:
            opac[n] = np.exp(-gamma_decay_2*r[n])

    return opac

# Modelo cambiado para ver si va más rápido. Va igual de lento. aprox 260s. aprox 330s
"""
def opacity_powlaw_and_exp(r,radius_central_1,gamma_decay_1,radius_central_2):

    # continuity condition. It is easier to get gamma_decay_2 from radius_central_2 than otherwise.
    gamma_decay_2 = gamma_decay_1/radius_central_2 * np.log(radius_central_2/radius_central_1)

    opac = np.empty(len(r))

    for n in np.arange(len(r)):
        rad = r[n]
        if rad < radius_central_1:
            opac[n] = 1
        elif radius_central_1 < rad < radius_central_2:
            opac[n] = pow(rad/radius_central_1,-gamma_decay_1)
        elif rad > radius_central_2:
            opac[n] = np.exp(-gamma_decay_2*rad)
    
    return opac

"""


"""
 Opacity power-law & power-law * exp
"""

def opacity_powlaw_and_powlawexp(r,radius_central_1,gamma_decay_1,radius_central_2,gamma_decay_2):

    opac = pow(r/radius_central_1,-gamma_decay_1)
    opac[r<radius_central_1] = 1
    
    # continuity condition
    gamma_decay_3 = gamma_decay_1/radius_central_2 * np.log(radius_central_2/radius_central_1)
    
    for n in np.arange(len(r)):
        if r[n] > radius_central_2:
            opac[n] = pow(r[n]/radius_central_2,-gamma_decay_2) * np.exp(-gamma_decay_3*r[n])
    return opac

""" 
Opacity power-law sum & power-law * exp

"""

def opacity_powlawsum_and_powlawexp(r,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_limit_1):

    # Continuity conditions
    gamma_decay_2 = np.log( 1 - pow(radius_limit_1/radius_central_1,-gamma_decay_1) ) / np.log(radius_central_2/radius_limit_1)
    gamma_decay_4 = - 1/radius_central_3 * (np.log(pow(radius_central_3/radius_central_1,-gamma_decay_1)  + pow(radius_central_3/radius_central_2,-gamma_decay_2)))

    opac = pow(r/radius_central_1,-gamma_decay_1) + pow(r/radius_central_2,-gamma_decay_2)

    opac[r < radius_central_1] = 1

    for n in np.arange(len(r)):
        if r[n] > radius_central_3:
            opac[n] = pow(r[n]/radius_central_3,-gamma_decay_3) * np.exp(-gamma_decay_4*r[n])

    return opac 

""" 
Opacity power-law sum & power-law

"""

def opacity_powlawsum_and_powlaw(r,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,radius_limit_1,radius_limit_2):

    if not radius_limit_1 < radius_limit_2:
        raise TypeError("Sorry, radius_limit_1 is not lower than radius_limit_2")
    elif radius_central_3>radius_limit_2:
        raise TypeError("Sorry, radius_central_3 has to be greater than radius_limit_2")

    # Continuity conditions
    gamma_decay_2 = np.log( 1 - pow(radius_limit_1/radius_central_1,-gamma_decay_1) ) / np.log(radius_central_2/radius_limit_1)
    gamma_decay_3 = np.log( pow(radius_limit_2/radius_central_1,-gamma_decay_1) + pow(radius_limit_2/radius_central_2,-gamma_decay_2) )/np.log(radius_central_3/radius_limit_2)

    if gamma_decay_3<0:
        raise TypeError("G3<0")

    opac = pow(r/radius_central_1,-gamma_decay_1) + pow(r/radius_central_2,-gamma_decay_2)

    opac[r < radius_limit_1] = 1

    for n in np.arange(len(r)):
        x = r[n]
        if x > radius_limit_2:
            opac[n] = pow(x/radius_central_3,-gamma_decay_3)
    return opac 

""" 
Opacity power-law sum & power-law sum 

"""

def opacity_powlawsum_and_powlawsum(r,radius_central_1,gamma_decay_1,radius_central_2,radius_central_3,gamma_decay_3,radius_central_4,radius_limit_1,radius_limit_2):

    # Continuity conditions
    gamma_decay_2 = np.log( 1 - pow(radius_limit_1/radius_central_1,-gamma_decay_1) ) / np.log(radius_central_2/radius_limit_1)
    gamma_decay_4 = np.log( pow(radius_limit_2/radius_central_1,-gamma_decay_1) + pow(radius_limit_2/radius_central_2,-gamma_decay_2) - pow(radius_limit_2/radius_central_3,-gamma_decay_3) ) / np.log(radius_central_4/radius_limit_2)

    opac = pow(r/radius_central_1,-gamma_decay_1) + pow(r/radius_central_2,-gamma_decay_2)

    opac[r < radius_limit_1] = 1

    for n in np.arange(len(r)):
        x = r[n]
        if x > radius_limit_2:
            opac[n] = pow(x/radius_central_3,-gamma_decay_3) + pow(x/radius_central_4,-gamma_decay_4)
    return opac 

def opacity_powlaw_powlaw(r,radius_central_1,radius_limit,gamma_decay_1,gamma_decay_2):

    radius_central_2 = radius_limit * pow(radius_limit/radius_central_1,-gamma_decay_1/gamma_decay_2)

    opac = pow(r/radius_central_1,-gamma_decay_1)
    opac[r<radius_central_1] = 1

    for n in np.arange(len(r)):
        x=r[n]
        if x > radius_limit:
            opac[n] = pow(x/radius_central_2,-gamma_decay_2)
    return opac




















""" 
THE REST OF THE SCRIPT WILL BE PROBABLY TO DELETED AS IT IS NOT USED IN THE TFG

Opacity power-law for n rings 

Input:

* r vector. (r0,r1,r2,...). Contains the radius of the ring parts
* g vector. (g0,g1,g2,...). Contains the gamma decay coefficient of the ring parts

We calculate the R vector of constants, through the application of the continuity conditions.

"""

def opacity_piecewise(xdata,rvector,gvector):

    if len(rvector) != len(gvector): # probe that input data is correct
        
        import sys
        sys.exit("\n Error in the opacity_piecewise() input. \n The rvector and gvector do not have the same length. \n")
    
    # Create the Rvector 
    l = len(rvector)
    Rvector = np.zeros(l)
    Rvector[0] = rvector[0]
    
    for k in np.arange(l-1): # Calculate the Rvector
        
        k=k+1
        a = rvector[k]
        b = gvector[k]
        b_ = gvector[k-1]
        
        Rvector[k] = a * pow(Rvector[k-1]/a,b_/b)
        
    opac = np.empty(len(xdata))
    
    # Create the limit of the last ring and add it to rvector. This is computed in a simple way but can be modified
    
    rmax = pow(0.01,-1/gvector[l-1])*Rvector[l-1]
   
    rvector = np.concatenate( (rvector, np.array([rmax])) ,axis=None)
    
    for n in np.arange(len(xdata)): # Loop to calculate the opacity of each point depending on where it is placed
    
        r = xdata[n]
        
        if r < rvector[0]:
        
            opac[n] = 1
            
        for ring in np.arange(l):
        
            if rvector[ring+1] > r > rvector[ring]:
        
                opac[n] =pow(r/Rvector[ring],-gvector[ring])
            
    return opac
    
"""
Piecewise function 2 POWLAW SUM with intervals

This function is prepared for n rings.

Input:

* rvector. (r0,r1,r2,...). Contains the radius of the ring parts, including the constant opacity rings
* gvector. (g0,g1,g2,...). Contains the gamma decay coefficient of the ring parts. For constant opacity vectors, g will be equal to zero.
* cstvector. (cst0,cst1,cst2,...). When it refers to a ring with variable opacity, its value is one. When it refers to a constant opacity ring, its value is the opactiy.

We calculate the R vector of constants, through the application of the continuity conditions

"""
"""
def opacity_piecewise_sep(xdata,rvector,gvector,cstvector):

    # Just to check data input it is ok

    if len(rvector) != len(gvector):
        
        import sys
        sys.exit("\n Error in the opacity_piecewise_sep() input. \n The rvector and gvector do not have the same length. \n")
    
    if (len(rvector)) != len(cstvector):
    
        import sys
        sys.exit("\n Error in the opacity_piecewise_sep() input. \n The rvector and cstvector do not have the same length. \n")
    
    # Create the Rvector

    l = len(rvector)
    Rvector = np.zeros(l)
    Rvector[0] = rvector[0] # initial condition
    
    for k in np.arange(l-1): # Loop to calculate the Rvector
        
        k=k+1
        a = rvector[k]
        b = gvector[k]
        b_ = gvector[k-1]
        Rvector[k] = a * pow(Rvector[k-1]/a,b_/b)
        
    opac = np.empty(len(xdata)) # Create opacity vector that we will get as an output
    
    # Create the limit of the last ring and add it to rvector. This is computed in a simple way but can be modified
    
    rmax = pow(0.01,-1/gvector[l-1])*Rvector[l-1]
   
    rvector = np.concatenate( (rvector, np.array([rmax])) ,axis=None)
   
    for n in np.arange(len(xdata)): # Loop to calculate the opacity of each point depending on where it is placed
    
        r = xdata[n]
        
        if r < rvector[0]:
        
            opac[n] = 1
            
        for ring in np.arange(l):
        
            if rvector[ring+1] > r > rvector[ring]:
        
                opac[n] = cstvector[ring] *pow(r/Rvector[ring],-gvector[ring])
            
    #print(opac)
    
    return opac
    
"""
# This shows how the continuity conditions are satisfied in the n rings models.
"""
plt.figure()

x = np.linspace(0,10,100)
rvector = np.array([2,6,9])
gvector = np.array([1.2,0.9,1])
plt.plot(x,opacity_piecewise(x,rvector,gvector),'.',label="No separation")
rvector = np.array([2,4,5,7])
gvector = np.array([1.2,0.7,1e-6,1.5])
cstvector = np.array([1,1,0.6,1.2])
plt.plot(x,opacity_piecewise_sep(x,rvector,gvector,cstvector),'.',label="Separation")

plt.ylabel("Opacity")
plt.xlabel("r [$R_{\odot}$]")
plt.legend()
plt.show()
"""

# This shows the continuity of the opacity function calculated with modelos.relflux_model_allt_nrings_powlaw()
"""
x1 = np.linspace(0.3,5,30)
x2=np.linspace(5.2,7,30)
x3=np.linspace(7.3,10,30)

r1 = 0.3
r2 = 2.96233107
r3 = 3.41818377

g1=0.3
g2=1.5
g3=1.7

#plt.plot(x1,Opacity_1(x1,r1,g1),label="power law")
#plt.plot(x2,Opacity_1(x2,r2,g2),label="power law")
#plt.plot(x3,Opacity_1(x3,r3,g3),label="power law")
"""

"""
        
x=np.linspace(0,8,100)
r0=3
gammaDecay=1.56

radius_central_1= r0/2
radius_central_2 = r0
gamma_decay_1=0.5
gamma_decay_2=gammaDecay


"""


