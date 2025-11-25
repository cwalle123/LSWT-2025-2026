"""This file is used to find the ambient conditions in the wind tunnel for further calculations in other files"""

####################################################################################################

# External imports
import numpy as np

# Internal imports

####################################################################################################
"""Functions"""

v_air= 4 # m/s
a_air = 4 # m/s
M_air = v_air / a_air # Mach number
p_air = 101325 # Pa

def ideal_gas_law(p: float, rho: float, R: float , T: float):
    if p == None:
        p = 1
    return p

####################################################################################################
"""Run this file"""

def main():
    print("Hello World")

if __name__ == "__main__":
    main() # makes sure this only runs if you run *this* file, not if this file is imported somewhere else