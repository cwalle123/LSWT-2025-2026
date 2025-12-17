"""This file is used to compute the 2D Pressure Coefficients along the chord as a function of angle of attack
   Written by: Clifton-John Walle"""

####################################################################################################

# External imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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

def plot_cp_distribution(aoa, make_csv=False):
    cp = compute_cp_for_an_aoa(aoa)

    # Separate upper/lower Cp using masks and sort indices
    cp_upper = cp[upper_mask][upper_sorted_idx]
    cp_lower = cp[lower_mask][lower_sorted_idx]

    x_upper = tap_x_upper
    x_lower = tap_x_lower

    # Plotting
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

    # Export to CSV if requested
    if make_csv:
        # Find max length
        max_len = max(len(x_upper), len(x_lower))

        # Pad shorter arrays with NaN
        x_upper_pad = np.pad(x_upper, (0, max_len - len(x_upper)), constant_values=np.nan)
        cp_upper_pad = np.pad(cp_upper, (0, max_len - len(cp_upper)), constant_values=np.nan)
        x_lower_pad = np.pad(x_lower, (0, max_len - len(x_lower)), constant_values=np.nan)
        cp_lower_pad = np.pad(cp_lower, (0, max_len - len(cp_lower)), constant_values=np.nan)

        data = {
            "x_upper": x_upper_pad,
            "Cp_upper": cp_upper_pad,
            "x_lower": x_lower_pad,
            "Cp_lower": cp_lower_pad
        }

        df = pd.DataFrame(data)
        filename = f"cp_distribution_aoa_{aoa}.csv"
        df.to_csv(filename, index=False)
        print(f"Cp values saved to {filename}")

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
    # plot_cp_distribution(1, True)
    plot_multiple_cp_distributions([-1, 3, 7])

if __name__ == "__main__":
    main() # makes sure this only runs if you run *this* file, not if this file is imported somewhere else