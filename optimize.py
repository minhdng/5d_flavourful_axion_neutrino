# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:01:59 2021

@author: CodingCow
"""
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares

from matplotlib import pyplot as plt



def func(x, para):
    return x.dot(para)


np.random.seed(2)
data_x = np.mgrid[ 0:5:0.5, 0:5:0.5 ].reshape(2, -1).T


data_yerr = 0.1
data_y = func(data_x, (2, 3))

least_squares = LeastSquares(data_x, data_y, data_yerr, func)

m = Minuit(least_squares, (2, 3))

m.migrad()
m.hesse()

"""
def fcn(x, y, z):
    return (x - 2) ** 2 + (y - 3) ** 2 + (z - 4) ** 2

fcn.errordef = Minuit.LEAST_SQUARES

m = Minuit(fcn, x=0, y=0, z=0)

m.migrad()  # run optimiser
print(m.values)  # x: 2, y: 3, z: 4

m.hesse()   # run covariance estimator
print(m.errors)  # x: 1, y: 1, z: 1
"""
