"""This file is used to compute the wake flow velocity profile at any angle of attack
   Written by: Clifton-John Walle"""

####################################################################################################

# External imports
import numpy as np
import matplotlib.pyplot as plt

# Internal imports
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from Ambient_Conditions import q_inf_2d, v_air_2d, p_air_2d, two_dimensional_data
from Data_Handling_and_Processing import angle_of_attack_to_row

####################################################################################################
"""Functions"""

def compute_wake_velocity(aoa):
    """
    Computes wake velocity using:
        Cp(x) = (P_wake_meas - P_static_meas_interp) / q_inf_2d
        U(x)  = U_inf * sqrt(1 - Cp)
    """

    # --------------------------
    # 1. Probe positions
    # --------------------------
    total_wake_positions_mm = np.array([
        0, 12, 21, 27, 33, 39, 45, 51, 57, 63, 69, 72, 75, 78,
        81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114,
        117, 120, 123, 126, 129, 132, 135, 138, 141, 144,
        147, 150, 156, 162, 168, 174, 180, 186, 195, 207, 219
    ])

    static_wake_positions_mm = np.array([
        43.5, 55.5, 67.5, 79.5, 91.5, 103.5,
        115.5, 127.5, 139.5, 151.5, 163.5, 175.5
    ])

    # --------------------------
    # 2. Select AoA row
    # --------------------------
    row = two_dimensional_data.iloc[angle_of_attack_to_row(aoa)]

    # --------------------------
    # 3. Wake total-pressure taps
    # --------------------------
    wake_cols = [f"P{str(i).zfill(3)} (Pa)" for i in range(50, 97)]
    wake_cols = [c for c in wake_cols if c in two_dimensional_data.columns]

    P_wake_meas = np.array([row[c] for c in wake_cols])  # This is Pt - P_off

    # --------------------------
    # 4. Static pressure taps P098–P109
    # --------------------------
    static_cols = [f"P{str(i).zfill(3)} (Pa)" for i in range(98, 110)]
    P_static_meas = np.array([row[c] for c in static_cols])  # This is Ps - P_off

    # --------------------------
    # 5. Interpolate static pressure to wake positions
    #    No extrapolation → outside range = NaN
    # --------------------------
    P_static_interp = np.interp(
        total_wake_positions_mm,
        static_wake_positions_mm,
        P_static_meas,
        left=np.nan,
        right=np.nan)

    # --------------------------
    # 6. True pressure difference: Pt - Ps
    # --------------------------
    dP_true = P_wake_meas - P_static_interp

    # --------------------------
    # 7. Cp and velocity
    # --------------------------
    cp_wake = dP_true / q_inf_2d

    if np.any(cp_wake > 1):
        print("\nSome values of Cp in the wake are more than 1. \nHere it is not possible to compute the velocity, which indicates the limits of neglecting viscosity and not having tunnel corrections\n")

    v_wake = v_air_2d * np.sqrt(1 - cp_wake)

    return total_wake_positions_mm, v_wake, cp_wake

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