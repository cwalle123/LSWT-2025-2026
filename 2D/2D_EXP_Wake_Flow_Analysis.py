"""This file is used to compute the wake flow velocity profile at any angle of attack
   Written by: """

####################################################################################################

# External imports
import numpy as np
import matplotlib.pyplot as plt

# Internal imports
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from Ambient_Conditions import q_inf_2d, p_air_2d, two_dimensional_data

####################################################################################################
"""Functions"""

# ADD FUCTIONS HERE

# TIP: THE WAKE IS FROM P050 TO P097 (Im not sure what happens in P098 to P113 tho, maybe thats another parts of the pressure wake)

####################################################################################################
"""Run this file"""

def main():
    print("Hello World")

if __name__ == "__main__":
    main() # makes sure this only runs if you run *this* file, not if this file is imported somewhere else