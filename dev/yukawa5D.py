#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 18:25:21 2021

@author: Peggy
"""
# Packages
# ------------------------------------------------------------  
from math import sqrt
import numpy as np
from numpy.random import rand
# ------------------------------------------------------------  





# Methods for random matrices
# ------------------------------------------------------------  
def generate_random_matrix(
        dim:int, abs_min:float=0.001, abs_max:float=4.0,
        output_type:str="complex", flattened:bool=False
        ) -> np.array:
    """
    Generate random square COMPLEX matrix of dimension (dim, dim)
    with matrix element norms UNIFORMLY generated between [abs_min, abs_max)
    and matrix element phases UNIFORMLY generated between [- Pi, Pi).

    Parameters
    ----------
    dim : int
        Dimension of the desired output square matrix.
    abs_min : float, optional
        Min of the matrix element norms. The default is 0.001.
    abs_max : float, optional
        Max of the matrix element norms. The default is 4.0.

    Returns
    -------
    np.array
        A single matrix of dimension (dim, dim).
    """

    abss = abs_min + (abs_max - abs_min) * rand(dim, dim)
    args = - np.pi + 2 * np.pi * rand(dim, dim)

    if output_type == "real": 
        return abss.reshape((dim * dim)) if flattened else abss

    elif output_type == "phase":
        return args.reshape((dim * dim)) if flattened else args

    elif output_type == "complex":
        if flattened: 
            return np.concatenate((abss.reshape((dim * dim)), 
                                   args.reshape((dim * dim))))
        return abss * np.exp(args * 1j)

    elif output_type == "unitary":
        res = np.linalg.qr(abss * np.exp(args * 1j))[0]
        if flattened: 
            abss = np.absolute(res)
            args = np.angle(res)
            return np.concatenate((abss.reshape((dim * dim)), 
                                   args.reshape((dim * dim))))
        return res

    else: 
        raise ValueError(
            "Variable `output` must be either \"complex\", \"real\", " + 
            "or \"phase\" or \"unitary\"."
            )
# ------------------------------------------------------------  

# ------------------------------------------------------------
def unflatten_random_matrix(elms:np.array, is_complex:bool=True) -> np.array:
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
    elms : np.array
        Real parameter input.
        [0, dim * dim) : real parts
        [dim * dim, -1] : phases 

    Returns
    -------
    np.array
        A single matrix of dimension (dim, dim).
    """

    if is_complex:
        dim = int(sqrt(len(elms) / 2))
        assert dim ** 2 == len(elms) / 2, "Incorrect input length."
        
        real_part = elms[:(dim*dim)].reshape((dim, dim))
        img_part = elms[(dim*dim):].reshape((dim, dim))

        return real_part * np.exp(1j * img_part)

    else:
        dim = int(sqrt(len(elms)))
        assert dim ** 2 == len(elms), "Incorrect input length."
    
        return elms.reshape((dim, dim))
# ------------------------------------------------------------

# ------------------------------------------------------------
def convert_phases(phases:np.array, start:float=-np.pi):
    return np.remainder(phases - start, 2*np.pi) + start
# ------------------------------------------------------------

# ------------------------------------------------------------
def swap(arr, start_index, last_index, axis=0):
    if axis == 0:
        arr[[start_index, last_index], :] = arr[[last_index, start_index], :] 
    elif axis == 1:
        arr[:, [start_index, last_index]] = arr[:, [last_index, start_index]] 
    else:
        raise ValueError
# ------------------------------------------------------------



# Methods for profile
# ------------------------------------------------------------
def fermionProfile(
        cL:np.array, z:float, zuv:float=1.001, zir:float=1.e8
        ) -> np.array:
    """
    Obtain fermion profiles from fermion c-parameters.

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
    return np.piecewise(cL.astype('float'), [cL!=0.5, cL==0.5], [
        lambda cL: np.sqrt( 
            (0.5 - cL) * z**(2. - cL) / (zir**(1.-2.*cL) - zuv**(1.-2.*cL)) 
            ),
        lambda cL: 1. / np.sqrt(2. * (np.log(zir) - np.log(zuv)))
        ])
# ------------------------------------------------------------  

# ------------------------------------------------------------  
def ytilde(y, cL_array, cR_minus_array, z=1.01, zuv=1.01, zir=1.e8):
    """
    Scale Yukawa by left and right fermion profiles.

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





# Methods for optimization purpose
# ------------------------------------------------------------
def getChiSquared(y, ym, err):
    return np.sum(((y - ym)/err)**2)
# ------------------------------------------------------------
