# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 18:34:08 2021

@author: CodingCow
"""
# ------------------------------------------------------------  
# ----- Import standard packages -----
import math
import numpy as np




# ------------------------------------------------------------  
# ---- Constants and parameters -----
CKM = {
       "s12": {"val": 0.22650, "sdev": 0.00048, 
               "sdev_plus": 0.00048, "sdev_minus": 0.00048},
       "s13": {"val": 0.00361, "sdev": 0.00010, 
               "sdev_plus": 0.00011, "sdev_minus": 0.00009},
       "s23": {"val": 0.04053, "sdev": 0.00072, 
               "sdev_plus": 0.00083, "sdev_minus": 0.00061},
       "delta": {"val": 1.196, "sdev": 0.0044, 
                 "sdev_plus": 0.045, "sdev_minus": 0.043}
       }

PMNS = { 
        None
        }

Yukawa = {
    "u": {"val": 3.29456 * 10**(-6), "sdev": 1.41 * 10**(-9)},
    "c": {"val": 0.00165737, "sdev": 7.4 * 10**(-7)},
    "t": {"val": 0.497757, "sdev": 2. * 10**(-6)},
    "d": {"val": 0.0000244146, "sdev": 5. * 10**(-10)},
    "s": {"val": 0.000486184, "sdev": 2.21 * 10**(-7)},
    "b": {"val": 0.0237974, "sdev": 4. * 10**(-7)},
    "e": {"val": 2.96535 * 10**(-6), "sdev": 3. * 10**(-11)},
    "mu": {"val": 0.000624451, "sdev": 3. * 10**(-9)},
    "tau": {"val": 0.0106073, "sdev": 1. * 10**(-7)}
    }



# ------------------------------------------------------------  
# ----- Matrix utils -----
def randomMatrix(dim, abs_min=0.01, abs_max=4.0):
    """
    Random complex square matrix of dimension dim * dim
        with entry norms uniformly choosen between [abs_min, abs_max)
    """
    return ( abs_min + np.random.rand(dim, dim)*(abs_max-abs_min) ) + \
                    np.random.rand(dim, dim)*np.pi*2j



def randomUnitaryMatrix(dim, abs_min=0.01, abs_max=4.0):
    """
    Random square unitary matrix of dimension dim * dim
        with entry norms between [abs_min, abs_max)
    """
    return np.linalg.qr( randomMatrix(dim, abs_min, abs_max) )[0]



# ------------------------------------------------------------  
# ----- 
def getCKM(unitaryMatrix): 
    """

    Parameters
    ----------
    unitaryMatrix : TYPE
        DESCRIPTION.

    Returns
    -------
    s12 : TYPE
        DESCRIPTION.
    s13 : TYPE
        DESCRIPTION.
    s23 : TYPE
        DESCRIPTION.
    delta : TYPE
        DESCRIPTION.
    phi_left : TYPE
        DESCRIPTION.
    phi_right : TYPE
        DESCRIPTION.

    """
    # separate norms and phases
    abss = np.absolute(unitaryMatrix)
    args = np.angle(unitaryMatrix)
    
    # CKM angles, s for sin, c for cos
    s13 = abss[0][2]; c13 = np.sqrt(1 - s13*s13)
    c12 = abss[0][0] / c13; s12 = np.sqrt(1 - c12*c12)
    c23 = abss[2][2] / c13; s23 = np.sqrt(1 - c23*c23)
    
    # calculating delta_CKM require calculating the quark rotation phases
    u1 = 0.5*np.angle(
        - np.exp([-(1j/3)*(
            4*args[0][0] + args[0][1] + 2*args[1][2] + 2*args[2][2]
            )]) * 
        (
            np.exp([ 1j*(args[1][0] + args[2][2]) ]) * c23 * abss[1][0] - 
            np.exp([ 1j*(args[1][2] + args[2][0]) ]) * s23 * abss[2][0] 
        ) 
    )
    
    phi_left = np.array([
        u1,
        -u1 + (1/3)*(- args[0][0] - args[0][1] - 2*args[1][2] + args[2][2]),
        -u1 + (1/3)*(- args[0][0] - args[0][1] + args[1][2] - 2*args[2][2])
        ]).T[0] # to convert to numpy 1d array
    
    phi_right = np.array([
        u1 + args[0][0],                  
        u1 + args[0][1],                 
        -u1 + (1/3)*(-args[0][0] - args[0][1] + args[1][2] + args[2][2])
        ]).T[0]
    
    # delta_CKM
    delta = (np.remainder(np.pi - (1/3)*(6*u1 + args[0][0] + args[0][1] + 
                          3*args[0][2] - args[1][2] - args[2][2]), 2 * np.pi) 
             - np.pi)[0]    
    
    return s12, s13, s23, delta, phi_left, phi_right
# ------------------------------------------------------------  



# ------------------------------------------------------------  
def getUnitaryMatrix(s12, s13, s23, delta):
    mat23 = np.array(
        [[1, 0, 0], 
         [0, np.sqrt(1-s23*s23), s23], 
         [0, -s23, np.sqrt(1-s23*s23)]]
        )
    mat13 = np.array(    
        [[np.sqrt(1-s13*s13), 0, s13*np.exp(-1j*delta)],               
         [0, 1, 0],               
         [-s13*np.exp(1j*delta), 0, np.sqrt(1-s13*s13)]]    
        )
    mat12 = np.array(            
        [[np.sqrt(1-s12*s12), s12, 0],                            
         [-s12, np.sqrt(1-s12*s12), 0],                
         [0, 0, 1]]
        )
    return mat23.dot(mat13).dot(mat12)
# ------------------------------------------------------------  



# ------------------------------------------------------------  
def fermionProfile(cL, z, zuv=1.01, zir=10.):
    return np.piecewise(cL, [cL != 0.5, cL == 0.5], [
        lambda cL: np.sqrt(
            ( 0.5 - cL ) / ( np.power(zir, 1-2*cL) - np.power(zuv, 1-2*cL) )
            ) * np.power(z, 2-cL),
        lambda cL: 1 / np.sqrt( 2 * (np.log(zir) - np.log(zuv)) )
        ])
# ------------------------------------------------------------  



# ------------------------------------------------------------  
def ytilde(y, cL_array, cR_minus_array, z=1.01, zuv=1.01, zir=10.):
    cL_mat = np.diag(fermionProfile(cL_array, z, zuv, zir))
    cR_minus_mat = np.diag(fermionProfile(cR_minus_array, z, zuv, zir))
    return cL_mat.dot(y).dot(cR_minus_mat)
# ------------------------------------------------------------  



# ------------------------------------------------------------
def predictCKM(x, c_array):
    cL = c_array[0:3]
    cR1_minus = c_array[3:6]
    cR2_minus = c_array[6:]
    
    y1 = x[0][0]
    y2 = x[0][1]
    
    y1_tilde = ytilde(y1, cL, cR1_minus)
    y2_tilde = ytilde(y2, cL, cR2_minus)
    
    u1, s1, _ = np.linalg.svd(y1_tilde)
    _, s2, vh2 = np.linalg.svd(y2_tilde)
    
    ckm_unitary = u1.dot(vh2)
    s12, s13, s23, delta, phi_l, phi_r = getCKM(ckm_unitary)
    
    return np.concatenate((s1, s2, np.array([s12, s13, s23, delta])))
# ------------------------------------------------------------