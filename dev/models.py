#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 16:49:09 2021

@author: Peggy
"""

from math import sqrt, log
import numpy as np
from numpy.random import rand
from scipy.linalg import svd
#from numpy.linalg import svd

from constants import Const
from yukawa5D import generate_random_matrix, unflatten_random_matrix, convert_phases, swap
from yukawa5D import ytilde

class Model():
    """
    Super-class encompassing all models
    """
    @classmethod
    def __init_subclass__(self, input_names, output_names):
        self.input_names = input_names
        self.output_names = output_names


class Quarks(
        Model,
        input_names = [
            # Real parts of 5D Yu matrix entries
            "yu11", "yu12", "yu13",
            "yu21", "yu22", "yu23",
            "yu31", "yu32", "yu33",
            # Phases of 5D Yu matrix entries
            "pu11", "pu12", "pu13",
            "pu21", "pu22", "pu23",
            "pu31", "pu32", "pu33",
            # Real parts of 5D Yd matrix entries
            "yd11", "yd12", "yd13",
            "yd21", "yd22", "yd23",
            "yd31", "yd32", "yd33",
            # Phases of 5D Yd matrix entries
            "pd11", "pd12", "pd13",
            "pd21", "pd22", "pd23",
            "pd31", "pd32", "pd33",
            # Left and right profile parameters
            "cL1", "cL2", "cL3",
            "cu1", "cu2", "cu3",
            "cd1", "cd2", "cd3"
            ],
        output_names = [
            # 4D quark (diagonal) couplings at 10^8 - 10^10 GeV
            "log_y4D_u", "log_y4D_c", "log_y4D_t", 
            "log_y4D_d", "log_y4D_s", "log_y4D_b", 
            # CKM parameters (do not run significantly, taken at expm. value)
            "s12", "s13", "s23", "delta"
            ]
        ):
    
    def __init__(self):
        self.input_limits = np.array(
            [(0.01, 4.)] * 9 + 
            [(-np.pi, np.pi)] * 9 +  
            [(0.01, 4.)] * 9 + 
            [(-np.pi, np.pi)] * 9 + 
            [(-2., 5.)] * 9
            )
        self.output_values = np.array([
            [log(Const.Yukawa[y]["val"]) for y in ["u", "c", "t", "d", "s", "b"]] +
            [Const.CKM[y]["val"] for y in ["s12", "s13", "s23", "delta"]]
            ])
        self.output_errors = np.array([
            [Const.Yukawa[y]["std"] / Const.Yukawa[y]["val"] for y in ["u", "c", "t", "d", "s", "b"]] +
            [Const.CKM[y]["std"] for y in ["s12", "s13", "s23", "delta"]]
            #[log(10) * 0.1 * Const.CKM[y]["val"] for y in ["s12", "s13", "s23", "delta"]]
            ])
        
        # det version
        self.output_det_values = np.array([
            [sum([log(Const.Yukawa[y]["val"]) for y in ["u", "c", "t"]])] +
            [sum([log(Const.Yukawa[y]["val"]) for y in ["d", "s", "b"]])] + 
            [Const.CKM[y]["val"] for y in ["s12", "s13", "s23", "delta"]]
            ])
        self.output_det_errors = np.array([
            [sum([Const.Yukawa[y]["std"] / Const.Yukawa[y]["val"] for y in ["u", "c", "t"]])] + 
            [sum([Const.Yukawa[y]["std"] / Const.Yukawa[y]["val"] for y in ["d", "s", "b"]])] + 
            [Const.CKM[y]["std"] for y in ["s12", "s13", "s23", "delta"]]
            ])

    def generate_random_inputs(
            self, dim:int=3, abs_min:float=0.01, abs_max:float=4., 
            c_min:float=0.1, c_max:float=1.
            ):
    	return np.concatenate((
    		generate_random_matrix(dim, abs_min, abs_max, output_type="complex", flattened=True), 
    		generate_random_matrix(dim, abs_min, abs_max, output_type="complex", flattened=True), 
    		c_min + (c_max - c_min) * rand(3 * dim)
    	))

    def _calculate_CKM(self, unitaryMatrix:np.array) -> tuple: 
        """
        Calculate CKM parameter from random unitary CKM matrix. 
        
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
        s13 = abss[0][2] / sqrt(abss[0][0]**2 + abss[0][1]**2 + abss[0][2]**2)
        s12 = abss[0][1] / sqrt(abss[0][0]**2 + abss[0][1]**2)
        c23 = abss[2][2] / sqrt(abss[1][2]**2 + abss[2][2]**2)
        s23 = abss[1][2] / sqrt(abss[1][2]**2 + abss[2][2]**2)
        
        # calculating delta_CKM require calculating the quark rotation phases
        u1 = 0.5 * np.angle(
            - np.exp([-(1j/3) * (    
                4*args[0][0] + args[0][1] + 2*args[1][2] + 2*args[2][2]
                )]) * 
            (
                np.exp([1j*(args[1][0] + args[2][2])]) * c23 * abss[1][0] - 
                np.exp([1j*(args[1][2] + args[2][0])]) * s23 * abss[2][0] 
            ) 
        )
        
        phi_left = np.array([
            u1,
            -u1 + (1/3)*(-args[0][0] - args[0][1] - 2*args[1][2] + args[2][2]),
            -u1 + (1/3)*(-args[0][0] - args[0][1] + args[1][2] - 2*args[2][2])
            ]).T[0] # to convert to numpy 1d array
        
        phi_right = np.array([
            u1 + args[0][0],                  
            u1 + args[0][1],                 
            -u1 + (-args[0][0] - args[0][1] + args[1][2] + args[2][2]) / 3
            ]).T[0]
        
        # delta_CP
        # defined in range [-pi, pi)
        delta = -(6*u1 + args[0][0] + args[0][1] + 3*args[0][2] - args[1][2] - args[2][2]) / 3
        
        return s12, s13, s23, convert_phases(delta)[0], phi_left, phi_right
    
    def get_output_full(self, x:np.array, inputs:np.array):
        """
        Main function used for optimization.
        Get observable quantities from input parameters.
    
        Parameters
        ----------

    
        Returns
        -------
        np.array
            (3 yu diagonal Yukawas, 
             3 yd diagonal Yukawas, 
             s12, s13, s23, delta).
        """     
        # extract complex Yukawa matrices
        y1 = unflatten_random_matrix(inputs[0:18], is_complex=True)
        y2 = unflatten_random_matrix(inputs[18:36], is_complex=True)
        # extract cL and cR
        cL, cR1_minus, cR2_minus = inputs[18:21], inputs[21:24], inputs[24:27]   
        
        y1_tilde = ytilde(y1, cL, cR1_minus)
        y2_tilde = ytilde(y2, cL, cR2_minus)
        
        # # using numpy
        # u1, s1, _ = svd(y1_tilde,)
        # _, s2, vh2 = svd(y2_tilde)        

        # using scipy        
        u1, s1, _ = svd(y1_tilde, check_finite=True)
        _, s2, vh2 = svd(y2_tilde, check_finite=True)
        
        s1 = np.sort(s1)
        s2 = np.sort(s2)
            
        swap(u1, 0, 2, axis=1)
        swap(vh2, 0, 2)
        
        ckm_unitary = u1.dot(vh2)
        s12, s13, s23, delta, phi_l, phi_r = self._calculate_CKM(ckm_unitary)
        
        return np.concatenate((np.log(s1), np.log(s2), 
                                   np.array([s12, s13, s23, delta])))

    def get_output_masses(self, x:np.array, inputs:np.array):
        # extract complex Yukawa matrices
        y1 = unflatten_random_matrix(inputs[0:18], is_complex=True)
        y2 = unflatten_random_matrix(inputs[18:36], is_complex=True)
        # extract cL and cR
        cL, cR1_minus, cR2_minus = inputs[18:21], inputs[21:24], inputs[24:27]   
        
        y1_tilde = ytilde(y1, cL, cR1_minus)
        y2_tilde = ytilde(y2, cL, cR2_minus)
        
        # # using numpy
        # _, s1, _ = svd(y1_tilde,)
        # _, s2, _ = svd(y2_tilde)        

        # using scipy        
        _, s1, _ = svd(y1_tilde, check_finite=True)
        _, s2, _ = svd(y2_tilde, check_finite=True)

        s1 = np.sort(s1)
        s2 = np.sort(s2)
        
        return np.concatenate((np.log(s1), np.log(s2)))

    def get_output_det(self, x:np.array, inputs:np.array):
        y1 = unflatten_random_matrix(inputs[0:18], is_complex=True)
        y2 = unflatten_random_matrix(inputs[18:36], is_complex=True)
        # extract cL and cR
        cL, cR1_minus, cR2_minus = inputs[18:21], inputs[21:24], inputs[24:27]   
        
        y1_tilde = ytilde(y1, cL, cR1_minus)
        y2_tilde = ytilde(y2, cL, cR2_minus)
        
        u1, s1, _ = svd(y1_tilde, check_finite=True)
        _, s2, vh2 = svd(y2_tilde, check_finite=True)
        
        ckm_unitary = u1.dot(vh2)
        s12, s13, s23, delta, phi_l, phi_r = self._calculate_CKM(ckm_unitary)
        
        return np.array([
            sum(np.log(s1)), sum(np.log(s2)), s12, s13, s23, delta
            ])  

    def get_output_CKM(self, x:np.array, inputs:np.array):
        # extract complex Yukawa matrices
        y1 = unflatten_random_matrix(inputs[0:18], is_complex=True)
        y2 = unflatten_random_matrix(inputs[18:36], is_complex=True)
        # extract cL and cR
        cL, cR1_minus, cR2_minus = inputs[18:21], inputs[21:24], inputs[24:27]   
        
        y1_tilde = ytilde(y1, cL, cR1_minus)
        y2_tilde = ytilde(y2, cL, cR2_minus)
        
        u1, s1, _ = svd(y1_tilde, check_finite=True)
        _, s2, vh2 = svd(y2_tilde, check_finite=True)
        
        ckm_unitary = u1.dot(vh2)
        s12, s13, s23, delta, phi_l, phi_r = self._calculate_CKM(ckm_unitary)
        
        return np.array([s12, s13, s23, delta]) 

    def check_mass_order(self, inputs:np.array):
        pass








