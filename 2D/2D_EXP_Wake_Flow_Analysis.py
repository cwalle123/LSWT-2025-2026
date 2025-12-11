"""This file is used to compute the wake flow velocity profile at any angle of attack
   Written by: Clifton-John Walle"""

####################################################################################################

# External imports
import numpy as np
import matplotlib.pyplot as plt

# Internal imports
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from Ambient_Conditions import q_inf_2d, v_air_2d, two_dimensional_data
from Data_Handling_and_Processing import angle_of_attack_to_row

####################################################################################################
"""Functions"""

def compute_wake_velocity(aoa):

    # Wake probe positions (mm)
    wake_positions_mm = np.array([
        0, 12, 21, 27, 33, 39, 45, 51, 57, 63, 69, 72, 75, 78,
        81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114,
        117, 120, 123, 126, 129, 132, 135, 138, 141, 144,
        147, 150, 156, 162, 168, 174, 180, 186, 195, 207, 219])

    # Shift so it starts at 0 mm
    # wake_positions_mm = wake_positions_mm - wake_positions_mm[0]

    # Continue as before...
    row_index = angle_of_attack_to_row(aoa)
    row = two_dimensional_data.iloc[row_index]

    wake_cols = [f"P{str(i).zfill(3)} (Pa)" for i in range(50, 97)]
    wake_cols = [c for c in wake_cols if c in two_dimensional_data.columns]

    cp_wake = np.array([row[c] / q_inf_2d for c in wake_cols])

    if np.any(cp_wake > 1):
        print("\nSome values of Cp in the wake are more than 1. \nHere it is not possible to compute the velocity, which indicates the limits of neglecting viscosity\n")

    v_wake = v_air_2d * np.sqrt(1 - cp_wake)

    # Trim list lengths if needed
    n = min(len(cp_wake), len(wake_positions_mm))
    cp_wake = cp_wake[:n]
    v_wake = v_wake[:n]
    wake_positions_mm = wake_positions_mm[:n]

    return wake_positions_mm, v_wake, cp_wake

def plot_wake_velocities(aoa):
    x, v, cp = compute_wake_velocity(aoa)

    plt.plot(x, v)
    plt.xlabel("X (wake)")
    plt.ylabel("Velocity (m/s)")
    plt.title(f"Wake Velocity Profile at AoA = {aoa}°")
    plt.xlim(0, 220)
    plt.grid(True)
    plt.show()

def plot_multiple_wake_velocities(aoa_list):
    """
    Plot wake velocity profiles for multiple angles of attack with different node styles.
    
    Inputs:
        aoa_list : list of AoA values (e.g. [-2, 0, 2, 4, 6])
    """

    plt.figure(figsize=(10, 5))

    # Define a list of marker styles
    markers = ['s', 'o', '^', 'D', 'v', 'P', '*', 'X', '<', '>']  # circles, squares, triangles, diamonds, etc.

    for i, aoa in enumerate(aoa_list):
        # Compute wake velocity for this AoA
        x_wake, v_wake, cp_wake = compute_wake_velocity(aoa)

        # Cycle through markers if more AoAs than markers
        marker_style = markers[i % len(markers)]

        # Plot line with different marker
        plt.plot(x_wake, v_wake, marker=marker_style, linestyle='-', label=f"AoA {aoa}°")

    # Formatting
    plt.xlabel("X (wake) [mm]")
    plt.ylabel("Velocity (m/s)")
    plt.title("Wake Velocity Profiles for Multiple AoA Values")
    plt.grid(True)
    plt.xlim(0, 220)
    plt.legend()
    plt.tight_layout()
    plt.show()

####################################################################################################
"""Run this file"""

def main():
    # aoa = 18
    # plot_wake_velocities(aoa)
    plot_multiple_wake_velocities([0, 7, 18])
if __name__ == "__main__":
    main() # makes sure this only runs if you run *this* file, not if this file is imported somewhere else