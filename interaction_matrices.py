import numpy as np
from math import factorial
#%Interaction matrices
#%Numerical diagonalization of fluxonium and antennae
#%
#%This function calculates the full interaction matrix as well as
#%the symmetric one ( Re[I] = cos() ) and
#%the anti-symmetric one ( Im[I] = sin() ).
#%
#%It depends on the energy offset "m",
#%the relative energy difference "a" between the levels,
#%and the strenth of the Flux ZPF "z_x" 
#
#%Strength of the flux ZPF
#%zr = 0.0540;
#%zq = 2.8903;

#%z = zq;
#
#%Hilbert space dimension
#%M_max = 20;
#
#% Interaction matrix
M_int = np.zeros((M_max, M_max))
#%Symmetric interaction matrix "cos()"
S_int = np.zeros((M_max, M_max))
#%Anti-symmetric interaction matrix "cos()"
A_int = zeros((M_max, M_max))

#%Interaction function 1
C1_az = np.zeros(M_max)


def interaction_function_1(a, z):
    """
    Interaction function 1
    Numerical diagonalization of fluxonium and antennae
    
    This function is independent of the "base" energy quantum
    (energy offset), it only takes into account the relative 
    energy difference between the levels.
    
    Order parameter "a" (difference of energy quanta)
    z coefficient "z" (relative strength of the "zero point fluctiations")
    """
    return ( 1/( sqrt( factorial(a)*(2**a) ) ) )*(z**a)*np.exp(-(1/4)*z**2)

def interaction_function_2(m, a, z):
    """
    Numerical diagonalization of fluxonium and antennae
    This function is depends on both the "base" energy quanta
    (energy offset) and the the relative energy difference 
    between the levels.

    Offset energy quanta "m"
    Order parameter "a" (difference of energy quanta)
    z coefficient "z" (relative strength of the "zero point fluctiations")
    """
    if m==0:
        F_maz = 1
    else:
        Gmma = np.sqrt(factorial(m+a) / ( factorial(a) *factorial(m) ) )
        Pmma = 0
        for vv in range(m+1):
            Pmma = Pmma + ((-1)**vv)* ( factorial(a) *factorial(m) / \
                ( (2**vv) * factorial(vv)*factorial(vv+a)*factorial(m-vv) ) )\
                *( z**(2*vv) )
        F_maz = Gmma*Pmma
    return F_maz

Csg_a = 0
for a in range(M_max):
    F_az = interaction_function_1(a, z)
    if mod(a, 2) == 0:
        Csg_a = (-1)**(a/2)
    else:
        Csg_a = (-1)**((a-1)/2)
    C1_az[a] = F_az*Csg_a

for m in range(M_max):
    for a in range(M_max-m):
        F_maz = interaction_function_2(m, a, z)
        
        #Full interaction matrix
        M_int[m,m+a] = C1_az[a]*F_maz
        M_int[m+a,m] = C1_az[a]*F_maz
        
        if mod(a,2) == 0:
            #%Symmetric
            S_int[m,m+a] = C1_az[a]*F_maz
            S_int[m+a,m] = C1_az[a]*F_maz
            #%Anti-symmetric
            A_int[m,m+a] = 0
            A_int[m+a,m] = 0
        else:
            #%Symmetric
            S_int[m,m+a] = 0
            S_int[m+a,m] = 0
            #%Anti-symmetric
            A_int[m,m+a] = C1_az[a]*F_maz
            A_int[m+a,m] = C1_az[a]*F_maz