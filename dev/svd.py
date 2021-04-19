#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 11:15:46 2021

@author: Peggy
"""

import numpy as np
from scipy.linalg import svd

from yukawa5D import generate_random_matrix


def my_svd(m:np.array):
    m_squared_l = m.T.conj().dot(m)
    m_squared_r = m.dot(m.T.conj())
    # np.linalg.eig(m_squared_r)
    return (
        np.real(np.sqrt(np.linalg.eigvals(m_squared_l))), 
        np.real(np.sqrt(np.linalg.eigvals(m_squared_r)))
        )


# %% test

def swap(arr, start_index, last_index, axis=0):
    if axis == 0:
        arr[[start_index, last_index], :] = arr[[last_index, start_index], :] 
    elif axis == 1:
        arr[:, [start_index, last_index]] = arr[:, [last_index, start_index]] 
    else:
        raise ValueError


m = generate_random_matrix(3)

u, s, vh = svd(m)

swap(u, 0, 2, axis=1)
swap(vh, 0, 2)

print(m)
print(u @ np.diag(np.sort(s)) @ vh)



