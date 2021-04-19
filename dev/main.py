#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 16:23:41 2021

@author: Peggy
"""

import numpy as np
import pandas as pd

from tqdm import tqdm

from yukawa5D import generate_random_matrix, unflatten_random_matrix, getChiSquared
from models import Quarks

from iminuit import Minuit
from iminuit.cost import LeastSquares

np.random.seed(2021)


# Testing random matrix generators
# m = generate_random_matrix(3, flattened=True)
# m = unflatten_random_matrix(m)


# %% Basic examle optimizing parameters simultaneously

# Initialize model
model = Quarks()

# Define chi squared using output values and least square function
x = np.array([None])
params = model.generate_random_inputs()
y = model.get_output_full(x, params)
chi_squared = getChiSquared(y, model.output_values, model.output_errors)
print("Initial chi squared: \n", chi_squared, "\n")

least_squares = LeastSquares(x, model.output_values, model.output_errors,
                              model.get_output_full)

# Define optimizer
optimizer = Minuit(least_squares, params, 
                    name=model.input_names)
optimizer.limits = model.input_limits
# optimizer.fixed[:36] = True

# Run optimizer
optimizer.migrad()

print("Chi squared after round 1: \n", optimizer.fval)
print("Log of Yukawa couplings: \n")
print(model.output_values)
print("At generated minimum:")
print(model.get_output_full(x, np.array(optimizer.values)))




# %% Basic example in loop to check efficiency
# Not very efficient, sub-percent efficiency

# number of 
chi_squared_threshold = 20.
N = 2000

# initialize output
res_0 = []

for i in tqdm(range(N)):
    model = Quarks()
    x = np.array([None])
 
    params = np.array(optimizer.values) # feed forward
    least_squares = LeastSquares(x, 
                                   model.output_values, 
                                   model.output_errors,
                                   model.get_output_full)
    optimizer = Minuit(least_squares, params, name=model.input_names)
    optimizer.limits = model.input_limits
    try:
        optimizer.migrad()
    except ValueError:
        pass
    
    if optimizer.fval and optimizer.fval < chi_squared_threshold:
        print("Match found!")
        res_0.append([optimizer.fval] + list(optimizer.values))

# Log the result

# res_pd = pd.read_csv("../data/quark_01.csv")

# res_names = ["chi_squared"] + model.input_names
# res_pd_temp = pd.DataFrame(res, columns=res_names)
# res_pd_temp = pd.concat([res_pd, res_pd_temp], axis=0).reset_index(drop=True)
# res_pd_temp.to_csv("../data/quark_02.csv", index=False)

# # Fix index
# res_pd = res_pd.drop("index", axis=1)
# res_pd.to_csv("../data/quark_01.csv", index=False)




# %% 2-step, masses then full
# Initialize model
model = Quarks()

# Define chi squared using output values and least square function
x = np.array([None])
params = model.generate_random_inputs()
y = model.get_output_masses(x, params)
chi_squared = getChiSquared(y, model.output_values[:, :6], model.output_errors[:, :6])
print("Initial chi squared: ", chi_squared)

least_squares = LeastSquares(x, model.output_values[:, :6], model.output_errors[:, :6],
                              model.get_output_masses)

# Define optimizer
optimizer = Minuit(least_squares, params, 
                    name=model.input_names)
optimizer.limits = model.input_limits

# fix all the phases
# optimizer.fixed[9:18] = True
# optimizer.fixed[27:36] = True

# Run optimizer
optimizer.migrad()

print("Chi squared after round 1:", optimizer.fval)
# print("Log of Yukawa couplings:")
# print(model.output_values[:, :6], "\n")
# print("At generated minimum")
# print(model.get_output_masses(x, np.array(optimizer.values)), "\n")
# print("The c-parameters are:")
# print(optimizer.values[-9:], "\n")


print("Step 1 done. Move on to step 2.\n")
params_2 = np.array(optimizer.values)
y_2 = model.get_output_full(x, params_2)
chi_squared_2 = getChiSquared(y_2, model.output_values, model.output_errors)
print("Initial chi squared round 2: ", chi_squared_2, "\n")

least_squares_2 = LeastSquares(x, model.output_values, model.output_errors,
                               model.get_output_full)

# Define optimizer
optimizer_2 = Minuit(least_squares_2, params_2, 
                     name=model.input_names)
optimizer_2.limits = model.input_limits

# fix all c-values
# optimizer_2.fixed[36:45] = True

# fix all the phases
#optimizer_2.fixed[9:18] = True
#optimizer_2.fixed[27:36] = True


# Run optimizer
optimizer_2.migrad()

print("Chi squared after round 2: ", optimizer_2.fval, "\n")
print("Log of Yukawa couplings: \n")
print(model.output_values, "\n")
print("At generated minimum:\n")
print(model.get_output_full(x, np.array(optimizer_2.values)), "\n")
# print("The c-parameters are:\n")
# print(optimizer_2.values[-9:], "\n")




# %% 2-step, masses then full, looped

chi_squared_threshold = 20.
N = 5

# initialize output
res = []

for i in tqdm(range(N)):
    model = Quarks()
    x = np.array([None])
    
    # Step 1
    params = model.generate_random_inputs()   
    least_squares = LeastSquares(x, 
                                 model.output_values[:, :6], 
                                 model.output_errors[:, :6], 
                                 model.get_output_masses)
    optimizer = Minuit(least_squares, params, name=model.input_names)
    optimizer.limits = model.input_limits
    optimizer.migrad()
    
    # Step 2
    params_2 = np.array(optimizer.values) # feed forward
    least_squares_2 = LeastSquares(x, 
                                   model.output_values, 
                                   model.output_errors,
                                   model.get_output_full)
    optimizer_2 = Minuit(least_squares_2, params_2, name=model.input_names)
    optimizer_2.limits = model.input_limits
    # optimizer_2.fixed[36:45] = True
    try:
        optimizer_2.migrad()
    except ValueError:
        pass
    
    if optimizer_2.fval and optimizer_2.fval < chi_squared_threshold:
        print("Match found!")
        res.append([optimizer_2.fval] + list(optimizer_2.values))

# Log the result

# res_pd = pd.read_csv("../data/quark_01.csv")

# res_names = ["chi_squared"] + model.input_names
# res_pd_temp = pd.DataFrame(res, columns=res_names)
# #res_pd_temp = pd.concat([res_pd, res_pd_temp], axis=0).reset_index(drop=True)
# res_pd_temp.to_csv("../data/quark_04.csv", index=False)

# # Fix index
# res_pd = res_pd.drop("index", axis=1)
# res_pd.to_csv("../data/quark_01.csv", index=False)



# %% 2-step, masses then full, looped with cache for learning

chi_squared_threshold = 20.
N = 5

# initialize output
res = []
cache = []

for i in tqdm(range(N)):
    model = Quarks()
    x = np.array([None])
    
    # Step 1
    params = model.generate_random_inputs()   
    least_squares = LeastSquares(x, 
                                 model.output_values[:, :6], 
                                 model.output_errors[:, :6], 
                                 model.get_output_masses)
    optimizer = Minuit(least_squares, params, name=model.input_names)
    optimizer.limits = model.input_limits
    optimizer.migrad()
    
    # Step 2
    params_2 = np.array(optimizer.values) # feed forward
    y_2 = model.get_output_full(x, params_2)
    chi_squared_2 = getChiSquared(y_2, model.output_values, model.output_errors)
    least_squares_2 = LeastSquares(x, 
                                   model.output_values, 
                                   model.output_errors,
                                   model.get_output_full)
    optimizer_2 = Minuit(least_squares_2, params_2, name=model.input_names)
    optimizer_2.limits = model.input_limits
    # optimizer_2.fixed[36:45] = True
    
    try:
        optimizer_2.migrad()
    except ValueError:
        pass
    
    if optimizer_2.fval:
        if optimizer_2.fval < chi_squared_threshold:
            print("Match found!")
            res.append([optimizer_2.fval] + list(optimizer_2.values))
            outcome = 1 # optimization successful
        else:
            outcome = 0 # optimization unsuccessful
    else:
        outcome = 2 # optimization encountered error

    cache.append([outcome, chi_squared_2] + list(optimizer.values))




# # Log the result

# appending = True

# res_names = ["chi_squared"] + model.input_names
# res_pd_temp = pd.DataFrame(res, columns=res_names)

# if appending: 
#     res_pd = pd.read_csv("../data/quark_05.csv")
#     res_pd_temp = pd.concat([res_pd, res_pd_temp], axis=0).reset_index(drop=True)

# res_pd_temp.to_csv("../data/quark_all.csv", index=False)

# # # Fix index if needed
# # res_pd = res_pd.drop("index", axis=1)
# # res_pd.to_csv("../data/quark_01.csv", index=False)

# # Log the cache
# cache_names = ["outcome", "chi_squared"] + model.input_names
# cache_pd_temp = pd.DataFrame(cache, columns=cache_names)

# if appending:
#     cache_pd = pd.read_csv("../data/quark_cache_05.csv")    
#     cache_pd_temp = pd.concat([cache_pd, cache_pd_temp], axis=0).reset_index(drop=True)

# cache_pd_temp.to_csv("../data/quark_cache_all.csv", index=False)





# %% 2-steps, masses then det
# Initialize model
model = Quarks()

# Define chi squared using output values and least square function
x = np.array([None])
params = model.generate_random_inputs()
y = model.get_output_masses(x, params)
chi_squared = getChiSquared(y, model.output_values[:, :6], model.output_errors[:, :6])
print("Initial chi squared: ", chi_squared)

least_squares = LeastSquares(x, model.output_values[:, :6], model.output_errors[:, :6],
                              model.get_output_masses)

# Define optimizer
optimizer = Minuit(least_squares, params, 
                    name=model.input_names)
optimizer.limits = model.input_limits

# fix all the phases
# optimizer.fixed[9:18] = True
# optimizer.fixed[27:36] = True

# Run optimizer
optimizer.migrad()

print("Chi squared after round 1:", optimizer.fval)
# print("Log of Yukawa couplings:")
# print(model.output_values[:, :6], "\n")
# print("At generated minimum")
# print(model.get_output_masses(x, np.array(optimizer.values)), "\n")
# print("The c-parameters are:")
# print(optimizer.values[-9:], "\n")


print("Step 1 done. Move on to step 3.\n")
params_3 = np.array(optimizer.values)
y_3 = model.get_output_det(x, params_3)
chi_squared_3 = getChiSquared(y_3, model.output_det_values, model.output_det_errors)
print("Initial chi squared round 3: ", chi_squared_3, "\n")

least_squares_3 = LeastSquares(x, model.output_det_values, model.output_det_errors,
                               model.get_output_det)

# Define optimizer
optimizer_3 = Minuit(least_squares_3, params_3, 
                     name=model.input_names)
optimizer_3.limits = model.input_limits

# fix all c-values
optimizer_3.fixed[36:45] = True

# fix all the phases
# optimizer_3.fixed[9:18] = True
# optimizer_3.fixed[27:36] = True


# Run optimizer
optimizer_3.simplex().migrad()

print("Chi squared after round 3: ", optimizer_3.fval, "\n")
global_chi_squared = getChiSquared(
    model.get_output_full(x, np.array(optimizer_3.values)),
    model.output_values, 
    model.output_errors)
print("Global chi squared after round 3: ", global_chi_squared, "\n")
# print("Log of Yukawa couplings: \n")
# print(model.output_det_values, "\n")
# print("At generated minimum:\n")
# print(model.get_output_det(x, np.array(optimizer_3.values)), "\n")
# print("The c-parameters are:\n")
# print(optimizer_2.values[-9:], "\n")


# %% Puttting the det algorithm in a loop
# for testing efficiency

# number of 
chi_squared_threshold = 20.
N = 100

# initialize output
res_3 = []

for i in tqdm(range(N)):
    model = Quarks()
    x = np.array([None])

    # step 1    
    params = model.generate_random_inputs()
    least_squares = LeastSquares(x, model.output_values[:, :6], model.output_errors[:, :6],
                                  model.get_output_masses)
    optimizer = Minuit(least_squares, params, 
                        name=model.input_names)
    optimizer.limits = model.input_limits
    optimizer.migrad()
    
    # step 3
    params_3 = np.array(optimizer.values)
    least_squares_3 = LeastSquares(x, model.output_det_values, model.output_det_errors,
                                   model.get_output_det)
    
    optimizer_3 = Minuit(least_squares_3, params_3, 
                         name=model.input_names)
    optimizer_3.limits = model.input_limits    
    optimizer_3.fixed[36:45] = True
    
    try:
        optimizer_3.simplex().migrad()
    except ValueError:
        pass
    
    if optimizer_3.fval and optimizer_3.fval < chi_squared_threshold:
        print("Match found!")
        res_3.append([optimizer_3.fval] + list(optimizer_3.values))


# %%
# Initialize model
model = Quarks()

# Define chi squared using output values and least square function
x = np.array([None])
params = model.generate_random_inputs()
y = model.get_output_masses(x, params)
chi_squared = getChiSquared(y, model.output_values[:, :6], model.output_errors[:, :6])
print("Initial chi squared: ", chi_squared)

least_squares = LeastSquares(x, model.output_values[:, :6], model.output_errors[:, :6],
                              model.get_output_masses)

# Define optimizer
optimizer = Minuit(least_squares, params, 
                    name=model.input_names)
optimizer.limits = model.input_limits

# fix all the phases
# optimizer.fixed[9:18] = True
# optimizer.fixed[27:36] = True

# Run optimizer
optimizer.migrad()

print("Chi squared after round 1:", optimizer.fval)
# print("Log of Yukawa couplings:")
# print(model.output_values[:, :6], "\n")
# print("At generated minimum")
# print(model.get_output_masses(x, np.array(optimizer.values)), "\n")
# print("The c-parameters are:")
# print(optimizer.values[-9:], "\n")


print("Step 1 done. Move on to step 2.\n")
params_4 = np.array(optimizer.values)
y_4 = model.get_output_CKM(x, params_4)
chi_squared_4 = getChiSquared(y_4, model.output_values[:, -4:], model.output_errors[:, -4:])
print("Initial chi squared round 2: ", chi_squared_4, "\n")

least_squares_4 = LeastSquares(x, model.output_values[:, -4:], model.output_errors[:, -4:],
                               model.get_output_CKM)

# Define optimizer
optimizer_4 = Minuit(least_squares_4, params_4, 
                     name=model.input_names)
optimizer_4.limits = model.input_limits

# fix all c-values
optimizer_4.fixed[36:45] = True

# fix all the phases
#optimizer_2.fixed[9:18] = True
#optimizer_2.fixed[27:36] = True


# Run optimizer
optimizer_4.migrad()

print("Chi squared after round 2: ", optimizer_4.fval, "\n")
global_chi_squared = getChiSquared(
    model.get_output_full(x, np.array(optimizer_4.values)),
    model.output_values, 
    model.output_errors)
print("Global chi squared after round 2: ", global_chi_squared, "\n")
print("Log of Yukawa couplings: \n")
print(model.output_values, "\n")
print("At generated minimum:\n")
print(model.get_output_full(x, np.array(optimizer_4.values)), "\n")
# print("The c-parameters are:\n")
# print(optimizer_2.values[-9:], "\n")

    
    
