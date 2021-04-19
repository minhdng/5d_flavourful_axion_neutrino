#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 17:54:57 2021

@author: Peggy
"""
import numpy as np

from model import Quarks
from matrixUtil import * 
import unittest


def is_unitary(m:np.array) -> bool:
    return np.allclose(np.eye(len(m)), m.dot(m.T.conj()))


# testing matrixUtil methods
# ------------------------------------------------------------  #
class test_random_matrix(unittest.TestCase):
    def test_unitary(self):
        self.assertTrue(
            is_unitary(generate_random_matrix(3, output_type="unitary"))
            )
        print("Unitary tests passed.")
