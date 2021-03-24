# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 18:34:08 2021

@author: CodingCow
"""
# ------------------------------------------------------------  
# ----- Import standard packages -----
import math
import numpy as np
from scipy import linalg





# ------------------------------------------------------------  
# ----- Matrix utils -----
def generateRandomMatrix(
        dim:int, abs_min:float=0.01, abs_max:float=4.0
        ) -> np.array:
    """
    Generate random square COMPLEX matrix of dimension (dim, dim)
        with matrix element norms UNIFORMLY generated between [abs_min, abs_max)
        and matrix element phases UNIFORMLY generated between [0, 2 Pi).

    Parameters
    ----------
    dim : int
        Dimension of the desired output matrix.
    abs_min : float, optional
        Min of the matrix element norms. The default is 0.01.
    abs_max : float, optional
        Max of the matrix element norms. The default is 4.0.

    Returns
    -------
    np.array
        A single matrix of dimension (dim, dim).

    """
    return ( abs_min + (abs_max-abs_min) * np.random.rand(dim, dim) ) * \
                       np.exp( np.random.rand(dim, dim)*np.pi*2j )



def generateRandomUnitaryMatrix(
        dim:int, abs_min:float=0.01, abs_max:float=4.0
        ) -> np.array:
    """
    Generate random square UNITARY matrix of dimension (dim, dim)
        with matrix element norms UNIFORMLY generated between [abs_min, abs_max)
        and matrix element phases UNIFORMLY generated between [0, 2 Pi).

    Method: Perform QR decomposition of a random complex matrix.

    Dependent functions:
        generateRandomMatrix()

    Parameters
    ----------
    dim : int
        Dimension of the desired output matrix.
    abs_min : float, optional
        Min of the matrix element norms. The default is 0.01.
    abs_max : float, optional
        Max of the matrix element norms. The default is 4.0.

    Returns
    -------
    np.array
        A single matrix of dimension (dim, dim).

    """
    return np.linalg.qr( generateRandomMatrix(dim, abs_min, abs_max) )[0]



def generateRandomParameters(dim:int, abs_min=0.01, abs_max=4.0) -> np.array:
    """
    Generate flattened (dim, dim) matrix with REAL parameters.
    
    Note: For optimization purpose, prefer flattened real arguments. 

    Parameters
    ----------
    dim : int
        Dimension of the desired output matrix.
    abs_min : TYPE, optional
        Min of the matrix element norms. The default is 0.01.
    abs_max : TYPE, optional
        Max of the matrix element norms. The default is 4.0.

    Returns
    -------
    np.array
        A single array of dimension (dim * dim, ).

    """
    return np.concatenate((
        abs_min + np.random.rand(dim*dim)*(abs_max-abs_min), 
        np.random.rand(dim*dim)
        ))



def generateRandomPhases(dim:int) -> np.array:
    """
    Generate flattened (dim, dim) matrix with elements being phases only. 
    
    Used in conjunction with generateRandomParameters() 
        to generate flattened COMPLEX matrices. 
        
    Note: For optimization purpose, prefer flattened real arguments. 
    
    Parameters
    ----------
    dim : int
         Dimension of the desired output matrix.

    Returns
    -------
    np.array
        A single array of dimension (dim * dim, ).

    """
    return np.random.rand(dim*dim)



def convertToRandomMatrix(params:np.array) -> np.array:
    """
    Generate random square COMPLEX matrix of dimension (dim, dim)
        with matrix element norms UNIFORMLY generated between [abs_min, abs_max)
        and matrix element phases UNIFORMLY generated between [0, 2 Pi).
        FROM a flattened array of real parameters of size (2 dim * dim, )

    generateRandomMatrix(dim) serves the same purpose as
    convertToRandomMatrix( np.concatenate(( generateRandomParameters(dim),
                                            generateRandomPhases(dim) )) )

    Parameters
    ----------
    params : np.array
        Real parameter input.
        [0, dim * dim) : real parts
         [dim * dim, -1] : phases 

    Returns
    -------
    np.array
        A single matrix of dimension (dim, dim).

    """
    real_part = params[:9]
    real_part = real_part.reshape((3, 3))
    img_part = params[9:]
    img_part = 2j * np.pi * img_part.reshape(3, 3)
    return real_part * np.exp( img_part )



# ------------------------------------------------------------  
# ----- 
def getCKM(unitaryMatrix:np.array) -> tuple: 
    """
    Obtaining CKM parameter from random unitary CKM matrix. 

    Parameters
    ----------
    unitaryMatrix : np.array
        numpy array of dimension (3, 3).

    Returns
    -------
    s12 : float
        sin(theta_{12}).
    s13 : float
        sin(theta_{13}).
    s23 : float
        sin(theta_{23}).
    delta : float
        delta^{CP} between -pi and pi.
    phi_left : numpy array of dimension (3, )
        Phases on V_L.
    phi_right : numpy array of dimension (3, )
        Phases on V_R.

    """
    # separate norms and phases
    abss = np.absolute(unitaryMatrix)
    args = np.angle(unitaryMatrix)
    
    # CKM angles, s for sin, c for cos
    s13 = abss[0][2] / math.sqrt(abss[0][0]**2 + abss[0][1]**2 + abss[0][2]**2)
    #c13 = math.sqrt(abss[0][0]**2 + abss[0][1]**2 ) / math.sqrt(abss[0][0]**2 + abss[0][1]**2 + abss[0][2]**2)
    #c12 = abss[0][0] / math.sqrt(abss[0][0]**2 + abss[0][1]**2 )
    s12 = abss[0][1] / math.sqrt(abss[0][0]**2 + abss[0][1]**2 )
    c23 = abss[2][2] / math.sqrt(abss[1][2]**2 + abss[2][2]**2 )
    s23 = abss[1][2] / math.sqrt(abss[1][2]**2 + abss[2][2]**2 )
    
    # calculating delta_CKM require calculating the quark rotation phases
    u1 = 0.5*np.angle(
        - np.exp([ - (1j/3) * (    
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
    
    # delta_CP
    delta = (np.remainder(np.pi - (1/3)*(6*u1 + args[0][0] + args[0][1] + 
                          3*args[0][2] - args[1][2] - args[2][2]), 2 * np.pi) 
             - np.pi)[0]    
    
    return s12, s13, s23, delta, phi_left, phi_right
# ------------------------------------------------------------  



# ------------------------------------------------------------  
def getUnitaryMatrix(s12:float, s13:float, s23:float, delta:float) -> np.array:
    """
    From CKM parameter, convert to a specific unitary CKM matrix.
    Assuming standard parametrization. 

    Parameters
    ----------
    s12 : float
        sin(theta_{12}).
    s13 : float
        sin(theta_{13}).
    s23 : float
        sin(theta_{23}).
    delta : float
        delta^{CP} between -pi and pi.

    Returns
    -------
    np.array
        Unitary matrix of dimension (3, 3).

    """
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
def fermionProfile(
        cL:np.array, z:float, zuv:float=1.01, zir:float=10.
        ) -> np.array:
    """
    Obtain fermion profiles from fermion profile c-parameters.

    Parameters
    ----------
    cL : np.array
        Numpy array of cL or -cR.
    z : float
        5D location.
    zuv : float, optional
        The default is 1.01.
    zir : float, optional
        The default is 10..

    Returns
    -------
    np.array
        Numpy array of f(cL) or f(-cR).

    """
    return np.piecewise(cL.astype('float'), [cL != 0.5, cL == 0.5], [
        lambda cL: np.sqrt(
            ( 0.5 - cL ) / ( np.power(zir, 1.-2.*cL) - np.power(zuv, 1.-2.*cL) )
            ) * np.power(z, 2. - cL),
        lambda cL: 1. / np.sqrt( 2. * (np.log(zir) - np.log(zuv)) )
        ])
# ------------------------------------------------------------  



# ------------------------------------------------------------  
def ytilde(y, cL_array, cR_minus_array, z=1.01, zuv=1.01, zir=10.):
    """
    Scale Yukawa by 5D profile parameter

    Parameters
    ----------
    y : np.array
        Matrix of dimension (3, 3).
    cL_array : np.array
        Left c-parameter array of dimension (3, ).
    cR_minus_array : np.array
        Right - c-parameter array of dimension (3, ).
    z : float, optional
        5D location. The default is 1.01.
    zuv : float, optional
        The default is 1.01.
    zir : float, optional
        The default is 10..

    Returns
    -------
    np.array
        Matrix of dimension (3, 3).

    """
    cL_mat = np.diag(fermionProfile(cL_array, z, zuv, zir))
    cR_minus_mat = np.diag(fermionProfile(cR_minus_array, z, zuv, zir))
    return cL_mat.dot(y).dot(cR_minus_mat)
# ------------------------------------------------------------  






# ------------------------------------------------------------
def predictCKMfromParams(params:np.array, c_array:np.array):
    """
    Function used for optimization.
    Predict observable quantities from given params and c_array.

    Parameters
    ----------
    params : np.array
        All Yukawa parameters, shape (36, ).
        (0:9, ) real part of y1
        (9:18, ) phases of y1
        (18:27, ) real part of y2
        (27:36, ) phases of y2
    c_array : np.array
        All c-paramters.
        (0:3, ) cL
        (3:6, ) -cU
        (6:9, ) -cD

    Returns
    -------
    np.array
        (3 yu diagonal Yukawas, 
         3 yd diagonal Yukawas, 
         s12, s13, s23, delta).

    """
    # extract cL and cR
    cL = c_array[0:3]
    cR1_minus = c_array[3:6]
    cR2_minus = c_array[6:]
    
    # extract complex Yukawa matrices
    p1 = params[0, :18]
    y1 = convertToRandomMatrix(p1)
    p2 = params[0, 18:]
    y2 = convertToRandomMatrix(p2)
    
    y1_tilde = ytilde(y1, cL, cR1_minus)
    y2_tilde = ytilde(y2, cL, cR2_minus)
    
    u1, s1, _ = linalg.svd(y1_tilde, check_finite=True)
    _, s2, vh2 = linalg.svd(y2_tilde, check_finite=True)
    
    ckm_unitary = u1.dot(vh2)
    s12, s13, s23, delta, phi_l, phi_r = getCKM(ckm_unitary)
    
    return np.concatenate((s1, s2, np.array([s12, s13, s23, delta])))
# ------------------------------------------------------------



# ------------------------------------------------------------
def predictCKMfromParams2(phases, x_array):
    """
    Function used for optimization.
    Predict observable quantities from given params and c_array.

    Parameters
    ----------
    phases : np.array
        All Yukawa phases, shape (18, ).
        (0:9, ) phases of y1
        (9:18, ) phases of y2
    x_array : np.array
        All Yukawa real parameters
        (0:9, ) real part of y1
        (9:18, ) real part of y2
        All c-paramters.
        (18:21, ) cL
        (21:24, ) -cU
        (24:27, ) -cD

    Returns
    -------
    np.array
        (3 yu diagonal Yukawas, 
         3 yd diagonal Yukawas, 
         s12, s13, s23, delta).

    """
    p1 = x_array[:9]
    p2 = x_array[9:18]
    
    if len(phases.shape) == 1:
        phase1 = phases[:9]
        phase2 = phases[9:]
    else: 
        phase1 = phases[0, :9]
        phase2 = phases[0, 9:]
    
    cL = x_array[18:21]
    cR1_minus = x_array[21:24]
    cR2_minus = x_array[24:27]
    
    y1 = convertToRandomMatrix( np.concatenate((p1, phase1)) )
    y2 = convertToRandomMatrix( np.concatenate((p2, phase2)) )
    
    y1_tilde = ytilde(y1, cL, cR1_minus)
    y2_tilde = ytilde(y2, cL, cR2_minus)
    
    u1, s1, _ = linalg.svd(y1_tilde, check_finite=True)
    _, s2, vh2 = linalg.svd(y2_tilde, check_finite=True)
    
    ckm_unitary = u1.dot(vh2)
    s12, s13, s23, delta, phi_l, phi_r = getCKM(ckm_unitary)
    
    return np.concatenate((s1, s2, np.array([s12, s13, s23, delta])))


