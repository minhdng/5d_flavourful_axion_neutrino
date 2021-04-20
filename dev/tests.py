#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 17:54:57 2021

@author: Peggy
"""
import numpy as np
import pandas as pd
import unittest


def is_unitary(m:np.array) -> bool:
    return np.allclose(np.eye(len(m)), m.dot(m.T.conj()))


# import saved data
data = pd.read_csv("../data/quark_08.csv")

test_case = np.array(data.iloc[1, 1:])

# testing matrixUtil methods
# ------------------------------------------------------------  #
from yukawa5D import generate_random_matrix
class test_random_matrix(unittest.TestCase):
    def test_unitary(self):
        self.assertTrue(
            is_unitary(generate_random_matrix(3, output_type="unitary"))
            )
        print("Unitary tests passed.")


# testing models
# Quark
from yukawa5D import unflatten_random_matrix, getChiSquared
yu = unflatten_random_matrix(test_case[:18], is_complex=True)
yd = unflatten_random_matrix(test_case[18:36], is_complex=True)
cL, cR1_minus, cR2_minus = test_case[36:39], test_case[39:42], test_case[42:]

from models import Quarks
model = Quarks()


x = np.array([None])
y_temp = model.get_output_full(x, test_case)
chi_squared_temp = getChiSquared(y_temp, model.output_values, model.output_errors)

print(y_temp[-4:])



# from scipy.linalg import svd
# svd(np.array( [[ 4.87822635e-03+1.88803452e-03, 4.96658382e-06+4.62887475e-06j,
#   -1.40903534e-03-2.49001362e-04j],
#   [-2.67510967e-01+4.26055398e-01j, -5.29183784e-05+9.28918649e-05j,
#     3.19391722e-04+2.65564278e-02j],
#   [ 9.80459976e-04+4.18661688e-05j,  2.86459201e-07-1.33207817e-06j,
#     7.66456399e-05+2.40972598e-05j]]))


# ckm_unitary = [[ 7.53287870e-01+0.61747274j, -2.07652458e-01-0.09028878j,
#    2.89976957e-03+0.00218838j],
#  [-1.17660474e-01+0.19330571j, -2.59261785e-01+0.93803267j,
#    4.05046355e-02+0.00513076j],
#  [ 6.77378580e-04-0.00855928j, 8.50874442e-03-0.03916668j,
#    9.97850930e-01+0.05112116j]]

# model._calculate_CKM(ckm_unitary)
