#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 18:25:21 2021

Python module for 5D Yukawa calculation. 
These methods are universal, independent of which the model used 
(quark CKM, or lepton PMNS matrix).

Methods for dealing with 5D Yukawa matrices include:
    
    - 5D Matrix random generator (output could be in flattened form). 
    Corresponding unflattener for generated array-form matrix.
    
        generate_random_matrix()
        unflatten_random_matrix()
        
    - Phase converter to deal with different phase conventions.
    Some function output phases between [0, 2 pi), 
    while the CKM matrix is specified between [-pi, pi).
                                               
        convert_phases()
        
    - Element swap function to deal 
    with the ordering of (scipy) SVD calculator.
        
        swap()
        
Methods for dealing with profile parameters include:
    
    - Calculate profiles given c-inputs.
    Calculate the scaled 5D Yukawa from the profiles.
    
        fermionProfile()
        ytilde()
        
Methods for dealing with optimization include: 
    
    - Standard Chi Squared function (for testing purpose, 
    this function is included in Iminuit standard package)

        get_chi_squared()
    
        

@author: Minh D. Nguyen
"""

# Packages
# ------------------------------------------------------------  
from math import sqrt
import numpy as np
from numpy.random import rand
# ------------------------------------------------------------  





# Methods for matrices
# ------------------------------------------------------------  
def generate_random_matrix(
        dim:int, abs_min:float=0.001, abs_max:float=4.0,
        output_type:str="complex", flattened:bool=False
        ) -> np.array:
    """
    Generate random square REAL matrix of dimension (dim, dim)
        with matrix elements UNIFORMLY generated between [abs_min, abs_max)
    or
    generate random square COMPLEX matrix of dimension (dim, dim)
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
    output_type : str, optional
        Select output type among ["complex", "real", "phase", "unitary"].
        The default value is "complex".
    flattened : bool, optional
        Select whether to flatten output.
        If False output standard matrix. 
        If True output 1D real array. 
        The default value is False.

    Returns
    -------
    np.array
        A singla matrix of dimension of dimension (dim, dim)
        If flattened: 
            If output_type is "complex" : (2 * dim * dim, )
                (First norms, then phases)
            If output_type is "real" :  (dim * dim, ) 
            If output_type is "phase" : (dim * dim, )
            If output_type is "unitary" : (2 * dim * dim, ).
    """
    # generate norms and phases separately
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
    Generate random square REAL matrix of dimension (dim, dim)
        with matrix elements UNIFORMLY generated between [abs_min, abs_max)
        FROM a flattened array of real parameters of size (dim * dim, ),
    or,
    generate random square COMPLEX matrix of dimension (dim, dim)
        with matrix element norms UNIFORMLY generated between [abs_min, abs_max)
        and matrix element phases UNIFORMLY generated between [0, 2 Pi).
        FROM a flattened array of real parameters of size (2 dim * dim, ).

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
def convert_phases(phases:np.array, start:float=-np.pi) -> np.array:
    """
    Convert phases into standardized range [start, start + 2 pi)

    Parameters
    ----------
    phases : np.array
    start : float, optional. The default is -np.pi.

    Returns
    -------
    np.array
        Converted phases.

    """
    return np.remainder(phases - start, 2 * np.pi) + start
# ------------------------------------------------------------

# ------------------------------------------------------------
def swap(matrix:np.array, i_0:int, i_f:int, axis:int=0):
    """
    Swap two particular rows (axis=0) or columns (axis=1) of given matrix.

    Parameters
    ----------
    matrix : np.array
        Input matrix.
    i_0 : int
        Index of the first row/column to swap.
    i_f : int
        Index of the second row/column to swap.
    axis : int, optional
        0 for row, 1 for column. The default is 0.

    Raises
    ------
    ValueError
        `axis` value must be 0 or 1.

    Returns
    -------
    None (in-place).

    """
    if axis not in [0, 1]:
        raise ValueError("`axis` value must be 0 (for row) and 1 (for column).")

    if axis == 0:
        matrix[[i_0, i_f], :] = matrix[[i_f, i_0], :] 
        
    matrix[:, [i_0, i_f]] = matrix[:, [i_f, i_0]] 
# ------------------------------------------------------------





# Methods for fermion profile
# ------------------------------------------------------------
def fermionProfile(
        cL:np.array, z:float, zuv:float=1.001, zir:float=1.e8
        ) -> np.array:
    """
    Obtain fermion profiles given fermion c-parameters.

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
    return np.piecewise(
        cL.astype('float'), 
        [cL!=0.5, cL==0.5], 
        [lambda cL: np.sqrt((0.5 - cL)/(zir**(1.-2.*cL)-zuv**(1.-2.*cL))) * z**(2. - cL),
        lambda cL: z**1.5 / np.sqrt(2. * (np.log(zir) - np.log(zuv)))
        ])
# ------------------------------------------------------------  

# ------------------------------------------------------------  
def ytilde(y, cL_array, cR_minus_array, 
           zuv=1., zir=1.e8, 
           localization:str='uv') -> np.array:
    """
    Scale Yukawa by left and right fermion profiles for UV or bulk case.

    Dependencies
    ----------
    fermionProfile()

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
        The default is 10.
    localization : str, optional    
        Type of fermion localization Must be in ['uv', 'bulk']. 
        The default is 'uv'.

    Returns
    -------
    np.array
        Matrix of dimension (3, 3).
    """
    localization_cases = ['uv', 'bulk']
    if localization not in localization_cases:
        raise ValueError(f"`localization` must be among {localization_cases}")

    cL_mat = np.diag(fermionProfile(cL_array, zuv, zuv, zir))
    cR_minus_mat = np.diag(fermionProfile(cR_minus_array, zuv, zuv, zir))

    if localization == 'uv':
        return cL_mat @ y @ cR_minus_mat

    integral_mat = np.ones((3, 3))
    for i in range(3):
        for j in range(3):
            cL, cRm = cL_array[i], cR_minus_array[j]
            cm = cL + cRm
            integral_mat[i, j] = np.piecewise(
                cm, 
                [cm!=0, cm==0], 
                [lambda cm: (zuv**(-cm) - zir**(-cm)) / cm, 
                 lambda cm: np.log(zir) - np.log(zuv)]
                )
    return 2 * cL_mat @ (y * integral_mat) @ cR_minus_mat
# ------------------------------------------------------------





# Methods for optimization purpose
# ------------------------------------------------------------
def get_chi_squared(y:np.array, ym:np.array, yerr:np.array) -> float:
    """
    Standard Chi Squared function.

    Parameters
    ----------
    y : np.array
        Prediction.
    ym : np.array
        Experimental values.
    yerr : np.array
        Experimental uncertainties.

    Returns
    -------
    float
        Chi Squared.

    """
    return np.sum(((y - ym) / yerr) ** 2)
# ------------------------------------------------------------
