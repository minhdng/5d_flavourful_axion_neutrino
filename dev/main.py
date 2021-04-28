#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 16:23:41 2021

@author: Peggy
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
from yukawa5D import generate_random_matrix, unflatten_random_matrix, get_chi_squared
from models import Quarks
from iminuit import Minuit
from iminuit.cost import LeastSquares



np.random.seed(2021)


# Global varibales
TRIAL = 6
model_name = "quark"
boundary_type = "uv"
logging = True
appending = False

# initialize output
res = []
cache = []

# initilize model
model = Quarks(boundary_type=boundary_type)
chi_squared_threshold = 10.
N = 100  # number of trials
x = np.array([None])    

# 2-step optimization process
# fit masses then fit full
for i in tqdm(range(N)):    
    
    # Step 1: fitting the masses only
    params = model.generate_random_inputs(c_min=1., c_max=5.)   
    least_squares = LeastSquares(x, 
                                 model.output_values[:, :6], 
                                 model.output_errors[:, :6], 
                                 model.get_output_masses)
    optimizer = Minuit(least_squares, params, name=model.input_names)
    optimizer.limits = model.input_limits
    optimizer.fixed[:36] = True
    try:
        optimizer.migrad()
        optimizer.hesse()
    except ValueError:
        pass

    # Step 2: fitting all at once  
    # Only run if step 1 is error-free
    if optimizer.fval:
        params_2 = np.array(optimizer.values) # feed forward
        y_2 = model.get_output_full(x, params_2)
        chi_squared_2 = get_chi_squared(y_2, model.output_values, model.output_errors)
        least_squares_2 = LeastSquares(x, 
                                       model.output_values, 
                                       model.output_errors,
                                       model.get_output_full)
        optimizer_2 = Minuit(least_squares_2, params_2, name=model.input_names)
        optimizer_2.limits = model.input_limits
        
        try:
            optimizer_2.migrad()
            optimizer_2.hesse()
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

# Logging
if logging: 
    # Log the result
    res_names = ["chi_squared"] + model.input_names
    res_url = f"../data/{model_name}_{boundary_type}_{TRIAL:02d}.csv"
    res_pd_temp = pd.DataFrame(res, columns=res_names)
    res_pd_temp.to_csv(res_url, index=False)
    
    if appending: 
        res_url_all = f"../data/{model_name}_{boundary_type}_all.csv"
        res_pd_all = pd.read_csv(res_url_all)
        res_pd_all = pd.concat([res_pd_all, res_pd_temp], axis=0).reset_index(drop=True)
        res_pd_all.to_csv(res_url_all, index=False)
    
    # Log the cache
    cache_names = ["outcome", "chi_squared"] + model.input_names
    cache_url = f"../data/{model_name}_{boundary_type}_cache_{TRIAL:02d}.csv"
    cache_pd_temp = pd.DataFrame(cache, columns=cache_names)
    cache_pd_temp.to_csv(cache_url, index=False)
    
    if appending:
        cache_url_all = f"../data/{model_name}_{boundary_type}_cache_all.csv"
        cache_pd_all = pd.read_csv(cache_url_all)    
        cache_pd_all = pd.concat([cache_pd_all, cache_pd_temp], axis=0).reset_index(drop=True)
        cache_pd_all.to_csv(cache_url_all, index=False)
    



# # Fix index if needed
# res_pd = res_pd.drop("index", axis=1)
# res_pd.to_csv("../data/quark_01.csv", index=False)

    
