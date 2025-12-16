#3.2.3 for github

"""This part contains the python code for the low speed wind tunnel test part 3.2.3.
(Aerodynamic Forces). It contains code with calculations, simulations and plots. It 
also contains conversions for data collected during the wind tunnel test and 
other functions imported from other files. Written by: Aram Kchadourian"""

"""Imports-------------------------------------------------------------------"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

#Import modules
import numpy as np
import matplotlib.pyplot as plt

#Import needed variables from other files
from Ambient_Conditions import(
    q_inf_2d,
    v_air_2d,
    rho_air_2d,
    p_air_2d
)
#Import functions from other files
from Data_Handling_and_Processing import (
    load_airfoil_pressure_tap_geometry,
    angle_of_attack_to_row
)
#imports functions for cp and wake velocity
from EXP_2D_Pressure_Coefficients import compute_cp_for_an_aoa
from EXP_2D_Wake_Flow_Analysis import compute_wake_velocity
from Ambient_Conditions import ref_static_p_inf_2d

#Import pressure tap location and cp_values from previous files
csv_path = "../Data/Processed_Data/airfoil_pressure_tap_coordinates.csv"
tap_data = load_airfoil_pressure_tap_geometry(csv_path)
#x/c position of each tap 
x_upper = tap_data["tap_x_upper"]
x_lower = tap_data["tap_x_lower"]
#Boolean to verify if tap is up or down
upper_mask = tap_data["upper_mask"]
lower_mask = tap_data["lower_mask"]
#Ordering for taps from leading- to trailing edge
upper_sorted_idx = tap_data["upper_sorted_idx"]
lower_sorted_idx = tap_data["lower_sorted_idx"]

#Import aoa values and convert to list
from Ambient_Conditions import two_dimensional_data
AOA_LIST = two_dimensional_data["Alpha (degrees)"].to_numpy()

"""--------------------------------------------------------------------------"""

"""Creates a common chordwise grid using upper tap locations and manually 
interpolates lower Cp values onto this grid. Returns:
        x_grid          : x/c positions (upper taps)
        cp_upper_sorted : Cp on upper surface
        cp_lower_interp : Cp on lower surface interpolated to x_grid
        """
def get_common_cp_grid(aoa):
    # 1) Get Cp values for this AoA (one value per pressure tap)
    cp_all = compute_cp_for_an_aoa(aoa)
    # 2) Split Cp into upper and lower using masks
    cp_upper = cp_all[upper_mask]
    cp_lower = cp_all[lower_mask]
    # 3) Sort Cp values from leading edge to trailing edge
    cp_upper_sorted = cp_upper[upper_sorted_idx]
    cp_lower_sorted = cp_lower[lower_sorted_idx]
    # 4) Use upper x/c positions as the common grid
    x_grid = x_upper
    # 5) Manually interpolate lower Cp onto upper x/c grid
    cp_lower_interp = np.zeros(len(x_grid))
    for i, x in enumerate(x_grid):
        # Loop over lower tap segments
        for j in range(len(x_lower) - 1):
            x1 = x_lower[j]
            x2 = x_lower[j + 1]
            # Check if x lies between two lower taps
            if x1 <= x <= x2:
                cp1 = cp_lower_sorted[j]
                cp2 = cp_lower_sorted[j + 1]
                # Linear interpolation
                cp_lower_interp[i] = cp1 + (cp2 - cp1) * (x - x1) / (x2 - x1)
                break
    return x_grid, cp_upper_sorted, cp_lower_interp

"""--------------------------------------------------------------------------"""

"""The following part of the code determines all the required coefficients, then
transforms them into aerodynamic forces and finally plots the most importatn ones."""

#Geometry & reference values
c = 0.16          # chord length [m] (from manual)
S = c * 1.0       # reference area for 2D (unit span)

#Storage arrays
Cn_list     = []
Cl_list     = []
Cd_list     = []
Cm_LE_list  = []
Cm_c4_list  = []

x_cp_c_list = []

L_list      = []
D_list      = []
M_LE_list   = []
M_c4_list   = []

#Loop over angles of attack
for aoa in AOA_LIST:
    aoa = int(aoa)
    # ---- 1. Common Cp grid ----
    x, cp_upper, cp_lower = get_common_cp_grid(aoa)
    # ---- 2. Normal force coefficient Cn ----
    delta_cp = cp_lower - cp_upper
    Cn = 0.0
    for i in range(len(x) - 1):
        dx = x[i+1] - x[i]
        avg_delta_cp = 0.5 * (delta_cp[i] + delta_cp[i+1])
        Cn += avg_delta_cp * dx
    # ---- 3. Moment coefficient about leading edge ----
    Cm_LE = 0.0
    for i in range(len(x) - 1):
        dx = x[i+1] - x[i]
        x_mid = 0.5 * (x[i] + x[i+1])
        avg_delta_cp = 0.5 * (delta_cp[i] + delta_cp[i+1])
        Cm_LE -= avg_delta_cp * x_mid * dx
    # ---- 4. Moment coefficient about quarter chord ----
    Cm_c4 = Cm_LE + 0.25 * Cn
    # ---- 5. Drag coefficient from wake (manual Eq. 2.3) ----
    wake_x, wake_v, _ = compute_wake_velocity(aoa)
    Cd = 0.0
    for i in range(len(wake_x) - 1):
        # skip invalid wake points
        if np.isnan(wake_v[i]) or np.isnan(wake_v[i+1]):
            continue
        dy = (wake_x[i+1] - wake_x[i]) / 1000  # mm → m
        uU_i  = wake_v[i]   / v_air_2d
        uU_ip = wake_v[i+1] / v_air_2d
        integrand_i  = uU_i  * (1 - uU_i)
        integrand_ip = uU_ip * (1 - uU_ip)
        # trapezoidal integration
        Cd += 0.5 * (integrand_i + integrand_ip) * dy
    Cd *= 2 / c
    # ---- 6. Lift coefficient (Eq. 2.5, manual) ----
    alpha = np.deg2rad(aoa)
    Cl = Cn * (np.cos(alpha) + (np.sin(alpha)**2)/np.cos(alpha)) - Cd * np.tan(alpha)
    # ---- 7. Calculate centre of pressure
    x_cp_c = - Cm_LE / Cn
    # ---- 8. Forces and moments ----
    L = Cl * q_inf_2d * S
    D = Cd * q_inf_2d * S
    M_LE = Cm_LE * q_inf_2d * S * c
    M_c4 = Cm_c4 * q_inf_2d * S * c
    # ---- 9. Store results ----
    Cn_list.append(Cn)
    Cl_list.append(Cl)
    Cd_list.append(Cd)
    Cm_LE_list.append(Cm_LE)
    Cm_c4_list.append(Cm_c4)

    x_cp_c_list.append(x_cp_c)

    L_list.append(L)
    D_list.append(D)
    M_LE_list.append(M_LE)
    M_c4_list.append(M_c4)

plt.figure()
plt.plot(AOA_LIST, Cl_list, marker='o')
plt.xlabel("Angle of attack α [deg]")
plt.ylabel("Lift coefficient $C_l$")
plt.title("$C_l$ vs Angle of Attack")
plt.grid(True)
plt.show()

plt.figure()
plt.plot(AOA_LIST, Cd_list, marker='o')
plt.xlabel("Angle of attack α [deg]")
plt.ylabel("Drag coefficient $C_d$")
plt.title("$C_d$ vs Angle of Attack")
plt.grid(True)
plt.show()

plt.figure()
plt.plot(AOA_LIST, Cm_c4_list, marker='o')
plt.xlabel("Angle of attack α [deg]")
plt.ylabel("Moment coefficient $C_{m,c/4}$")
plt.title("$C_{m,c/4}$ vs Angle of Attack")
plt.grid(True)
plt.show()

plt.figure()
plt.plot(Cl_list, Cd_list, 'o-', label=r"$C_d$ vs $C_l$")
plt.plot(Cl_list, Cm_c4_list, 's--', label=r"$C_{m,c/4}$ vs $C_l$")
plt.xlabel(r"Lift coefficient $C_l$")
plt.ylabel(r"Coefficient value")
plt.title(r"$C_d$ and $C_{m,c/4}$ as functions of $C_l$")
plt.grid(True)
plt.legend()
plt.show()

plt.figure()
plt.plot(AOA_LIST, x_cp_c_list, marker='o')
plt.xlabel("Angle of attack α [deg]")
plt.ylabel("Center of pressure x_cp / c [-]")
plt.title("Center of Pressure Location vs Angle of Attack")
plt.grid(True)
plt.tight_layout()
plt.show()

"""--------------------------------------------------------------------------"""

"""This next part focusses on the code to compute drag from wake rake data. This
includes plots to compare wake rake data with pressure tap data. """

# Storage lists
D_wake_list  = []   # Drag force per unit span [N/m]
Cd_wake_list = []   # Drag coefficient [-]
# Loop over angles of attack
for aoa in AOA_LIST:
    aoa = int(aoa)
    # 1. Get wake data for this AoA
    wake_y_mm, wake_v, cp_wake = compute_wake_velocity(aoa)
    # Convert y to meters
    wake_y = wake_y_mm / 1000.0
    # 2. Compute drag force D from wake equation
    #     D = rho ∫ (U∞ − U)U dy  +  ∫ (p∞ − p(y)) dy
    #     (turbulence term neglected)
    D = 0.0   # drag force per unit span [N/m]
    for i in range(len(wake_y) - 1):
        # Skip invalid wake points
        if np.isnan(wake_v[i]) or np.isnan(wake_v[i+1]):
            continue
        dy = wake_y[i+1] - wake_y[i]
        U_mid  = 0.5 * (wake_v[i] + wake_v[i+1])
        Cp_mid = 0.5 * (cp_wake[i] + cp_wake[i+1])
        # Term 1: momentum deficit (positive drag)
        D += rho_air_2d * (v_air_2d - U_mid) * U_mid * dy
        # Term 3: pressure recovery (reduces drag)
        D += Cp_mid * q_inf_2d * dy
    # 3. Convert drag force to drag coefficient
    Cd = D / (q_inf_2d * c)
    # 4. Store results
    D_wake_list.append(D)
    Cd_wake_list.append(Cd)

plt.figure(figsize=(7,5))
plt.plot(AOA_LIST, Cd_list, 'o-', label='Cd (surface pressure)')
plt.plot(AOA_LIST, Cd_wake_list, 's--', label='Cd (wake rake)')
plt.xlabel('Angle of attack α [deg]')
plt.ylabel('Drag coefficient Cd [-]')
plt.title('Drag coefficient vs angle of attack')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(7,5))
plt.plot(Cl_list, Cd_list, 'o-', label='Cd vs Cl (surface pressure)')
plt.plot(Cl_list, Cd_wake_list, 's--', label='Cd vs Cl (wake rake)')
plt.xlabel('Lift coefficient Cl [-]')
plt.ylabel('Drag coefficient Cd [-]')
plt.title('Drag polar comparison')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

"""--------------------------------------------------------------------------"""

"""The code below saves all lists into a txt file."""

# ============================================================
# Save all result lists to a text file
# ============================================================

output_path = "LSWT_3_2_3_results.txt"   # saved in current working folder

with open(output_path, "w") as f:

    f.write("==== 3.2.3 Aerodynamic Forces Results ====\n\n")

    f.write("Angle of attack [deg]:\n")
    f.write(str(list(AOA_LIST)) + "\n\n")

    f.write("Cn list:\n")
    f.write(str(Cn_list) + "\n\n")

    f.write("Cl list:\n")
    f.write(str(Cl_list) + "\n\n")

    f.write("Cd (surface pressure) list:\n")
    f.write(str(Cd_list) + "\n\n")

    f.write("Cd (wake rake) list:\n")
    f.write(str(Cd_wake_list) + "\n\n")

    f.write("Cm_LE list:\n")
    f.write(str(Cm_LE_list) + "\n\n")

    f.write("Cm_c/4 list:\n")
    f.write(str(Cm_c4_list) + "\n\n")

    f.write("Center of pressure x_cp/c list:\n")
    f.write(str(x_cp_c_list) + "\n\n")

    f.write("Lift force L [N/m]:\n")
    f.write(str(L_list) + "\n\n")

    f.write("Drag force D (surface) [N/m]:\n")
    f.write(str(D_list) + "\n\n")

    f.write("Drag force D (wake) [N/m]:\n")
    f.write(str(D_wake_list) + "\n\n")

    f.write("Moment about LE [Nm/m]:\n")
    f.write(str(M_LE_list) + "\n\n")

    f.write("Moment about c/4 [Nm/m]:\n")
    f.write(str(M_c4_list) + "\n\n")

print(f"Results successfully written to {output_path}")






