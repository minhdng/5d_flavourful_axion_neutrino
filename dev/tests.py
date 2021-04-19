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
data = pd.read_csv("../data/quark_7.csv")

test_case = np.array(data.iloc[0, 1:])

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


model._calculate_CKM(
    [
      [-0.968095 - 0.123777j, -0.133887 + 0.170621j, 0.0206427 - 0.00296574j], 
      [-0.0355603 - 0.214707j, 0.907812 + 0.355236j, -0.0100886 - 0.0471135j], 
      [-0.00539309 - 0.00889722j, -0.0466364 + 0.0217527j, -0.861182 - 0.505578j]]
    )

x = np.array([None])
y_temp = model.get_output_full(x, test_case)
chi_squared_temp = getChiSquared(y_temp, model.output_values, model.output_errors)

print(y_temp[-4:])



# svd(np.array( [[ 4.87822635e-03+1.88803452e-03, 4.96658382e-06+4.62887475e-06j,
#   -1.40903534e-03-2.49001362e-04j],
#  [-2.67510967e-01+4.26055398e-01j, -5.29183784e-05+9.28918649e-05j,
#    3.19391722e-04+2.65564278e-02j],
#  [ 9.80459976e-04+4.18661688e-05j,  2.86459201e-07-1.33207817e-06j,
#    7.66456399e-05+2.40972598e-05j]]))

