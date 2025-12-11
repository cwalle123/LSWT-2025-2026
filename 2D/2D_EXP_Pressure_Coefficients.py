"""This file is used to compute the 2D Pressure Coefficients along the chord as a function of angle of attack
   Written by: Clifton-John Walle"""

####################################################################################################

# External imports
import numpy as np
import matplotlib.pyplot as plt

# Internal imports
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from Ambient_Conditions import q_inf_2d, p_air_2d, two_dimensional_data
from Data_Handling_and_Processing import angle_of_attack_to_row, load_airfoil_pressure_tap_geometry

####################################################################################################
"""Functions"""

csv_path = os.path.join(os.path.dirname(__file__), "..", "Data", "Processed_Data", "airfoil_pressure_tap_coordinates.csv")
csv_path = os.path.abspath(csv_path)
tap_data = load_airfoil_pressure_tap_geometry(csv_path)

tap_names        = tap_data["tap_names"]
upper_mask       = tap_data["upper_mask"]
lower_mask       = tap_data["lower_mask"]
upper_sorted_idx = tap_data["upper_sorted_idx"]
lower_sorted_idx = tap_data["lower_sorted_idx"]
tap_x_upper      = tap_data["tap_x_upper"]
tap_x_lower      = tap_data["tap_x_lower"]

def compute_cp_for_an_aoa(aoa):
    """
    Returns Cp values for P001–P049 in the correct tap order
    as defined in the CSV.
    """
    row = two_dimensional_data.iloc[angle_of_attack_to_row(aoa)]

    cp_values = np.array([row[f"{tap} (Pa)"] / q_inf_2d for tap in tap_names]) # pressure is already Δp

    return cp_values

def plot_all_cp(aoa):
    cp = compute_cp_for_an_aoa(aoa)

    # Use linear x distribution only for this debug plot
    x = np.linspace(0, 1, len(cp))

    plt.figure(figsize=(8, 4))
    plt.plot(x, cp, marker='o')
    plt.gca().invert_yaxis()
    plt.xlabel("x / c")
    plt.ylabel("Cp")
    plt.title(f"Cp Distribution (AoA = {aoa}°)")
    plt.grid()
    plt.show()

def plot_cp_distribution(aoa):
    cp = compute_cp_for_an_aoa(aoa)

    # Separate upper/lower Cp using CSV masks and sort indices
    cp_upper = cp[upper_mask][upper_sorted_idx]
    cp_lower = cp[lower_mask][lower_sorted_idx]

    x_upper = tap_x_upper
    x_lower = tap_x_lower

    plt.figure(figsize=(10, 5))

    plt.plot(x_upper, cp_upper, 'o-', label="Upper surface")
    plt.plot(x_lower, cp_lower, 'o-', label="Lower surface")

    plt.gca().invert_yaxis()
    plt.xlabel("x / c")
    plt.ylabel("Cp")
    plt.title(f"Cp Distribution at AoA = {aoa}°")
    plt.grid()
    plt.legend()
    plt.show()

def plot_multiple_cp_distributions(aoa_list):
    from matplotlib.ticker import MultipleLocator

    plt.figure(figsize=(12, 6))

    for aoa in aoa_list:
        cp = compute_cp_for_an_aoa(aoa)

        cp_upper = cp[upper_mask][upper_sorted_idx]
        cp_lower = cp[lower_mask][lower_sorted_idx]

        x_u = tap_x_upper
        x_l = tap_x_lower

        x_total  = np.concatenate([x_u, x_l[::-1]])
        cp_total = np.concatenate([cp_upper, cp_lower[::-1]])

        plt.plot(x_total, cp_total, marker="o", label=f"AoA {aoa}°")

    plt.gca().invert_yaxis()
    plt.xlabel("x / c")
    plt.ylabel("Cp")
    plt.title("Cp Distributions for Multiple AoA Values")
    plt.xlim(0,1)

    ax = plt.gca()
    plt.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(0.1))  # major grid every 0.1
    ax.xaxis.set_minor_locator(MultipleLocator(0.05)) # minor grid every 0.05
    ax.yaxis.set_major_locator(MultipleLocator(0.5))  # adjust based on Cp range
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    plt.grid(which='major', linestyle='--', color='gray', alpha=1)
    plt.grid(which='minor', linestyle=':', color='lightgray', alpha=1)
    plt.legend()
    plt.tight_layout()
    plt.show()

####################################################################################################
"""Run this file"""

def main():
    # aoa = 7
    # plot_all_cp(aoa)
    # plot_cp_distribution(aoa)
    plot_multiple_cp_distributions([-1, 3, 7])

if __name__ == "__main__":
    main() # makes sure this only runs if you run *this* file, not if this file is imported somewhere else