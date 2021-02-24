# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 16:54:05 2021

@author: Minh D. Nguyen
"""
# ----- Import standard packages -----
import numpy as np
# from matplotlib import pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares


# ----- Update path to script folder -----
import os
path = os.path.dirname(os.path.abspath(__file__)) # prefer ('')
os.chdir(path)

# ----- Import home-made packages -----
import yukawa5d

# ----- To reproduce the result -----
np.random.seed(3)


# ----- Workspace -----
y1 = yukawa5d.randomMatrix(3)
y2 = yukawa5d.randomMatrix(3)

cL = np.array([0.2, 0.3, 0.4])
cR1_minus = np.array([2.1, 2.2, 3.2])
y1_tilde = yukawa5d.ytilde(y1, cL, cR1_minus)

cR2_minus = np.array([5.1, 7.2, 10.2])
y2_tilde = yukawa5d.ytilde(y2, cL, cR2_minus)

u1, s1, vh1 = np.linalg.svd(y1)

u2, s2, vh2 = np.linalg.svd(y2)

ckm_unitary = u1.dot(vh2)

s12, s13, s23, delta, phi_l, phi_r = yukawa5d.getCKM(ckm_unitary)

c_array = np.array([0.2, 0.3, 0.4, 2.1, 2.2, 3.2, 5.1, 7.2, 10.2])
data_x = np.array([[y1, y2]])

res = yukawa5d.predictCKM(data_x, c_array)


data_y = np.array([[yukawa5d.CKM["s12"]["val"], 
                   yukawa5d.CKM["s13"]["val"],
                   yukawa5d.CKM["s23"]["val"],
                   yukawa5d.CKM["delta"]["val"],
                   yukawa5d.Yukawa["u"]["val"],
                   yukawa5d.Yukawa["c"]["val"],
                   yukawa5d.Yukawa["t"]["val"],
                   yukawa5d.Yukawa["d"]["val"],
                   yukawa5d.Yukawa["s"]["val"],
                   yukawa5d.Yukawa["b"]["val"]
                   ]])

data_yerr = np.array([[yukawa5d.CKM["s12"]["sdev"], 
                     yukawa5d.CKM["s13"]["sdev"],
                     yukawa5d.CKM["s23"]["sdev"],
                     yukawa5d.CKM["delta"]["sdev"],
                     yukawa5d.Yukawa["u"]["sdev"],
                     yukawa5d.Yukawa["c"]["sdev"],
                     yukawa5d.Yukawa["t"]["sdev"],
                     yukawa5d.Yukawa["d"]["sdev"],
                     yukawa5d.Yukawa["s"]["sdev"],
                     yukawa5d.Yukawa["b"]["sdev"]
                     ]])

least_squares = LeastSquares(data_x, data_y, data_yerr, yukawa5d.predictCKM)

m = Minuit(least_squares, (1., 2., 20., 3.1, 3.2, 10.3, 4.1, 4.2, 10.3) )

m.migrad()  # finds minimum of least_squares function
m.hesse()   # accurately computes uncertainties

# ----- Check if the CKM parameters are correct -----
"""
m2 = yukawa5d.randomUnitaryMatrix(3)

s12, s13, s23, delta, phi_l, phi_r = yukawa5d.getCKM(m2)

m2_rotated = yukawa5d.getUnitaryMatrix(s12, s13, s23, delta)

m2_rotated_check = np.diag(np.exp(1j*phi_l)).dot(m2).dot(np.diag(np.exp(-1j*phi_r)))
"""

# print(yukawa5d.fermionProfile(np.array([0.2, 0.3, 0.4]), 2.))


# u1, s1, vh1 = np.linalg.svd(m1)
# print( np.allclose(m1, np.dot( np.dot( u1, np.diag(s1) ), vh1 )) )


# check if unitary
# print( np.allclose(np.identity(3), m2.dot(m2.conj().T) ) )


#foo = yukawa5d.CKM["s12"]["val"]
#print(foo)