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

def pressure_tap_locations(n_points, upper=False):
    """
    Generate the required custom x-distribution.
    If upper=True, inserts 0.375 immediately after passing 0.35.
    """

    base_x = [0, 0.005, 0.015, 0.04, 0.075, 0.12, 0.16, 0.20]

    x = base_x.copy()

    # Fill remaining points
    while len(x) < n_points - 1:

        # If upper curve and last point passed 0.35, insert 0.375 once
        if upper and x[-1] < 0.35:
            # Compute the next normal point
            nxt = x[-1] + 0.045
            if nxt >= 0.35:
                # Force next to be exactly 0.375
                x.append(0.375)
                continue

        # Normal behavior
        nxt = x[-1] + 0.045

        if nxt >= 0.93:
            break

        x.append(nxt)

    # If we still need more points and haven't hit 0.93 yet
    while len(x) < n_points - 1 and x[-1] < 0.93:
        x.append(0.93)

    # Final point always = 1.0
    x.append(1.0)

    # Trim or pad (safety)
    return np.array(x[:n_points])

def angle_of_attack_to_row(aoa):
    """
    This function uses an angle of attack and finds the row that contains the data of that angle of attack.
    
    Returns:
        The row that contains the data for the given angle of attack
    """

    if -3 <= aoa < 12: # WITH STEP SIZE 1.0 DEG
        if aoa % 1 != 0:
            raise ValueError('If you are providing an angle of attack in the range of -3 and 12, please provide an integer')
        row = aoa + 3
    elif 12 <= aoa <= 18: # WITH STEP SIZE 0.5 DEG
        if aoa % 0.5 != 0:
            raise ValueError('If you are providing an angle of attack in the range of 12 and 18, please provide an integer followed by .0 or .5')
        row = int(2 * aoa - 9)
    else:
        raise ValueError('Please provide an angle of attack in the range of -3 and 18')

    # print(row)
    return row

def compute_cp_for_an_aoa(aoa):
    """
    Compute pressure coefficients (Cp) for all P### (Pa) columns for a specific angle of attack.
    
    Returns:
        cp_values: numpy array of Cp values
    """

    # Extract row
    row_index = angle_of_attack_to_row(aoa)
    row = two_dimensional_data.iloc[row_index]

    # Identify all pressure tap columns (P001 (Pa) ... P113 (Pa))
    p_columns = [col for col in two_dimensional_data.columns if col.startswith("P") and col[1:4].isdigit()]

    # Sort taps numerically: P001 (Pa), P002 (Pa), ...
    sorted_columns = sorted(p_columns, key=lambda x: int(x[1:4]))

    # Compute Cp for each tap
    cp_values = np.array([(row[col]) / q_inf_2d for col in sorted_columns]) # The pressure measured is already the difference in the local pressure and the refrence pressure
    cp_values = cp_values[0:49]

    return cp_values

def plot_all_cp(aoa):
    """
    Plot all Cp values against a normalized chord axis (0 → 1)
    """

    # Extract row
    row_index = angle_of_attack_to_row(aoa)
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

def plot_cp_distribution(aoa):
    """
    Plot Cp distribution along the chord for a given angle of attack.
    """

    cp_values = compute_cp_for_an_aoa(aoa)

    cp_upper = cp_values[0:25]
    cp_lower = cp_values[25:49]

    # Upper uses special 0.375 rule
    x_upper = pressure_tap_locations(len(cp_upper), upper=True)

    # Lower uses normal rule
    x_lower = pressure_tap_locations(len(cp_lower), upper=False)

    plt.figure(figsize=(10, 5))

    plt.plot(x_upper, cp_upper, marker='o', linestyle='-', color="tab:blue", label="Upper")
    plt.plot(x_lower, cp_lower, marker='o', linestyle='-', color="tab:orange", label="Lower")

    plt.gca().invert_yaxis()
    plt.xlabel("x / c")
    plt.ylabel("Cp")
    plt.title(f"Cp Distribution (AoA = {aoa}°)")
    plt.xticks(np.arange(0, 1.01, 0.05))
    plt.grid(True)
    plt.legend()
    plt.show()

def plot_multiple_cp_distributions(aoa_list):
    """
    Plot Cp distributions for multiple angles of attack on one figure.
    """

    plt.figure(figsize=(12, 6))

    for aoa in aoa_list:

        # Compute Cp
        cp_values = compute_cp_for_an_aoa(aoa)

        # Split upper/lower
        cp_upper = cp_values[0:25]
        cp_lower = cp_values[25:49]

        # Apply spacing rules
        x_upper = pressure_tap_locations(len(cp_upper), upper=True)
        x_lower = pressure_tap_locations(len(cp_lower), upper=False)

        # Plot upper & lower as one combined curve per AoA
        # (Airfoil convention: plot upper first, then lower)
        x_total = np.concatenate((x_upper, x_lower[::-1]))
        cp_total = np.concatenate((cp_upper, cp_lower[::-1]))

        plt.plot(x_total, cp_total, marker='o', linestyle='-',
                 label=f"AoA {aoa}°")

    # Formatting
    plt.gca().invert_yaxis()
    plt.xlabel("x / c")
    plt.ylabel("Cp")
    plt.title("Pressure Coefficient Distribution for Multiple AoA Values")
    plt.xticks(np.arange(0, 1.01, 0.05))
    plt.grid(True)
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