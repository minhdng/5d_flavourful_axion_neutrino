# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 17:03:19 2021

@author: CodingCow
"""

class Const: 
    # ------------------------------------------------------------  
    # ---- CKM angles and phases -----
    CKM = {
           "s12": {"val": 0.22650, "std": 0.00048, 
                   "std_plus": 0.00048, "std_minus": 0.00048},
           "s13": {"val": 0.00361, "std": 0.00010, 
                   "std_plus": 0.00011, "std_minus": 0.00009},
           "s23": {"val": 0.04053, "std": 0.00072, 
                   "std_plus": 0.00083, "std_minus": 0.00061},
           "delta": {"val": 1.196, "std": 0.0044, 
                     "std_plus": 0.045, "std_minus": 0.043}
           }
    
    # ------------------------------------------------------------
    # PMNS phases
    PMNS = {
        None
        }

    
    # Neutrino masses difference, run to ... GeV
    neutrino = {
        None
        }
    
    # ------------------------------------------------------------     
    # Quark and Lepton 4D Yukawa parameters run to 10^8 GeV 
    Yukawa = {
        "u": {"val": 3.29456e-6, "std": 1.7129e-6},
        "c": {"val": 0.00165737, "std": 0.0000324975},
        "t": {"val": 0.497757, "std": 0.00350594},
        "d": {"val": 0.0000244146, "std": 4.59907e-6},
        "s": {"val": 0.000486184, "std": 0.0000287442},
        "b": {"val": 0.0237974, "std": 0.000172465},
        "e": {"val": 2.96535e-6, "std": 8.70457e-16},
        "mu": {"val": 0.000624451, "std": 1.35932e-11},
        "tau": {"val": 0.0106073, "std": 7.16363e-7}
        }

    
    