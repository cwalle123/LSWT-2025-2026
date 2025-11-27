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

####################################################################################################
"""Functions"""

def compute_cp_for_row(row_index):
    """
    Compute pressure coefficients (Cp) for all P### (Pa) columns in a specific row.
    
    Returns:
        cp_values: numpy array of Cp values, sorted by tap number
    """
    # Extract row
    row = two_dimensional_data.iloc[row_index]

    # Identify all pressure tap columns (P001 (Pa) ... P113 (Pa))
    p_columns = [col for col in two_dimensional_data.columns if col.startswith("P") and col[1:4].isdigit()]

    # Sort taps numerically: P001 (Pa), P002 (Pa), ...
    sorted_columns = sorted(p_columns, key=lambda x: int(x[1:4]))

    # Compute Cp for each tap
    cp_values = np.array([(row[col]) / q_inf_2d for col in sorted_columns]) # The pressure measured is already the difference in the local pressure and the refrence pressure
    cp_values = cp_values[0:49]

    return cp_values

def plot_all_cp(row_index):
    """
    Plot all Cp values against a normalized chord axis (0 → 1)
    without splitting into upper/lower surfaces.
    """

    # Extract row
    row = two_dimensional_data.iloc[row_index]

    # Identify all pressure tap columns (P001 (Pa) ... P113 (Pa))
    p_columns = [col for col in two_dimensional_data.columns if col.startswith("P") and col[1:4].isdigit()]

    # Sort taps numerically: P001 (Pa), P002 (Pa), ...
    sorted_columns = sorted(p_columns, key=lambda x: int(x[1:4]))

    # Compute Cp for each tap
    cp_values = np.array([(row[col]) / q_inf_2d for col in sorted_columns])

    # Create evenly spaced x/c positions
    num_points = len(cp_values)
    x = np.linspace(0, 1, num_points)

    # Plot
    plt.figure(figsize=(8, 4))
    plt.plot(x, cp_values, marker='o', linestyle='-')
    plt.gca().invert_yaxis()  # aerodynamic Cp convention
    plt.xlabel("x / c")
    plt.ylabel("Cp")
    plt.title(f"Cp Distribution (Row {row_index})")
    plt.grid(True)
    plt.xlim(0, 1)
    plt.show()

def plot_cp_distribution(row_index):
    """
    Plot Cp distribution along the chord for a given row.
    Assumes Cp array is ordered such that:
      - first half = upper surface taps (LE → TE)
      - second half = lower surface taps (TE → LE)
    """

    # Compute Cp array
    cp_values = compute_cp_for_row(row_index)

    # Split into upper and lower arrays
    cp_upper = cp_values[0:25]
    cp_lower = cp_values[25:49]

    # Number of points excluding the first
    n_upper = len(cp_upper) - 1
    n_lower = len(cp_lower) - 1

    # Create x arrays for points 2 → end (scaled 0 → 1)
    x_upper = np.linspace(0, 1, n_upper)
    x_lower = np.linspace(0, 1, n_lower)

    # Prepend first point at the same location as the second point (x = 0)
    x_upper_adjusted = np.concatenate(([x_upper[0]], x_upper))
    x_lower_adjusted = np.concatenate(([x_lower[0]], x_lower))

    # Adjust Cp arrays to match lengths
    cp_upper_adjusted = np.concatenate(([cp_upper[0]], cp_upper[1:]))
    cp_lower_adjusted = np.concatenate(([cp_lower[0]], cp_lower[1:]))

    plt.figure(figsize=(10, 5))

    # Plot upper surface
    plt.plot(x_upper_adjusted, cp_upper_adjusted, marker='o', linestyle='-', color="tab:blue")

    # Plot lower surface
    plt.plot(x_lower_adjusted, cp_lower_adjusted, marker='o', linestyle='-', color="tab:blue")

    # Aerodynamic Cp convention (inverted y-axis)
    plt.gca().invert_yaxis()

    plt.xlabel("x / c")
    plt.ylabel("Cp")
    plt.title(f"Pressure Coefficient Distribution (Row {row_index})")
    plt.xlim(-0.001, 1)
    plt.ylim(1, -1.6)
    plt.xticks(np.arange(0, 1.01, 0.05))
    plt.yticks(np.arange(-1.6, 1, 0.1))
    plt.grid(True)
    # plt.legend()
    plt.show()

####################################################################################################
"""Constants and Calculations"""

####################################################################################################
"""Run this file"""

def main():
    aoa = 4

    if -3 <= aoa < 12: # WITH STEP SIZE 1.0 DEG
        row = aoa + 3
    elif 12 <= aoa <= 18: # WITH STEP SIZE 0.5 DEG
        row = 2 * aoa - 9

    plot_cp_distribution(row)

if __name__ == "__main__":
    main() # makes sure this only runs if you run *this* file, not if this file is imported somewhere else