"""This file is used to find the ambient conditions in the wind tunnel for further calculations in other files
   Written by: Clifton-John Walle"""

####################################################################################################

# External imports
import numpy as np

####################################################################################################
"""Functions"""

def ideal_gas_law(p: float = None, rho: float = None, R: float = 287.057, T: float = None):
    """
    Compute one missing variable of the ideal gas law: p = rho * R * T

    Parameters
    ----------
    p : float or None
        Pressure in Pascals (Pa). Set to None if unknown.
    rho : float or None
        Density in kilograms per cubic meter (kg/m^3). Set to None if unknown.
    R : float or None
        Specific gas constant in Joules per kilogram per Kelvin (J/(kg·K)). Set to None if unknown.
    T : float or None
        Temperature in Kelvin (K). Set to None if unknown.

    Returns
    -------
    The missing variable as a float.

    Notes
    -----
    Exactly one variable must be None. If more than one or none are None,
    the function raises a ValueError.
    """
    vars_list = [p, rho, R, T]
    missing = vars_list.count(None)

    if missing != 1:
        raise ValueError("Exactly one variable must be None to compute it.")

    # Compute missing variable
    if p is None:
        p = rho * R * T
        return p
    elif rho is None:
        rho = p / (R * T)
        return rho
    elif R is None:
        R = p / (rho * T)
        return R
    elif T is None:
        T = p / (rho * R)
        return T

def dynamic_viscosity(T: float):
    """
    Compute the dynamic viscosity of air using Sutherland's law.

    Sutherland's Law:
        μ = μ₀ * (T / T₀)^(3/2) * (T₀ + S) / (T + S)

    Parameters
    ----------
    T : float
        Temperature in Kelvin (K).

    Returns
    -------
    float
        Dynamic viscosity of air in kilograms per meter-second (kg/(m·s)).

    Constants
    ---------
    μ₀ : float
        Reference dynamic viscosity of air at T₀ = 273.15 K.
        μ₀ = 1.716 × 10⁻⁵ kg/(m·s)

    T₀ : float
        Reference temperature (273.15 K).

    S : float
        Sutherland's constant for air (110.4 K).
    """
    mu_0 = 1.716 * (10**(-5)) # kg/(m·s)
    T_0 = 273.15 # K
    S = 110.4 # K
    mu = mu_0 * (T / T_0)**(1.5) * (T_0 + S) / (T + S)
    return mu

def dynamic_pressure(Delta_P_b: float):
    """
    Compute the dynamic pressure from a measured pressure difference
    using a quadratic calibration curve.

    The equation used is:
        q_inf = 0.21180 4 + 1.928442 · ΔP_b + 1.879374 * 10⁻⁴ · (ΔP_b)²

    Parameters
    ----------
    Delta_P_b : float
        Measured pressure difference (ΔP_b), in Pascals (Pa).

    Returns
    -------
    float
        Dynamic pressure q_inf in Pascals (Pa).
    """
    q_inf = 0.211804 + 1.928442 * Delta_P_b + 1.879374 * 10**(-4) * Delta_P_b**2
    return q_inf

def bernoulli(p1=None, v1=None, p2=None, v2=None, p_total_1=None, p_total_2=None, rho: float = None):
    """
    Solve one unknown in the horizontal, incompressible Bernoulli equation
    along a streamline, using optional total pressures.

    Behavior:
    1. If a total pressure (p_total_1 or p_total_2) is given AND one variable 
       from the other side (p2, v2, or p1, v1), compute the missing variable
       on the opposite side directly and return it.
    2. Otherwise, use the standard Bernoulli equation with individual pressures
       and velocities to compute the missing variable.

    Parameters
    ----------
    p1, p2 : float or None
        Static pressure at point 1 or 2, in Pascals (Pa).
    v1, v2 : float or None
        Flow velocity at point 1 or 2, in meters per second (m/s).
    p_total_1, p_total_2 : float or None
        Total (stagnation) pressure at point 1 or 2, in Pascals (Pa).
    rho : float
        Fluid density in kilograms per cubic meter (kg/m^3).

    Returns
    -------
    float
        The computed unknown variable.

    Notes
    -----
    - Exactly one unknown variable must exist on the side being solved.
    """

    if rho is None: 
        raise ValueError("Please provide a density.")

    # --- 1. Total pressure shortcut ---
    # If p_total_1 is given and one variable from side 2
    if p_total_1 is not None:
        if p2 is not None and v2 is None:  # solve v2
            return ((2 / rho) * (p_total_1 - p2)) ** 0.5
        if v2 is not None and p2 is None:  # solve p2
            return p_total_1 - 0.5 * rho * v2**2

    # If p_total_2 is given and one variable from side 1
    if p_total_2 is not None:
        if p1 is not None and v1 is None:  # solve v1
            return ((2 / rho) * (p_total_2 - p1)) ** 0.5
        if v1 is not None and p1 is None:  # solve p1
            return p_total_2 - 0.5 * rho * v1**2

    # --- 2. Standard Bernoulli calculation ---
    unknowns = [p1, v1, p2, v2]
    if unknowns.count(None) != 1:
        raise ValueError("Exactly one unknown variable must be provided.")

    # Compute total from whichever side is fully known
    if p1 is not None and v1 is not None:
        total = p1 + 0.5 * rho * v1**2
    elif p2 is not None and v2 is not None:
        total = p2 + 0.5 * rho * v2**2
    else:
        raise ValueError("Not enough information to compute Bernoulli relation.")

    # Solve the missing variable
    if p1 is None:
        return total - 0.5 * rho * v1**2
    elif v1 is None:
        return ((2 / rho) * (total - p1)) ** 0.5
    elif p2 is None:
        return total - 0.5 * rho * v2**2
    elif v2 is None:
        return ((2 / rho) * (total - p2)) ** 0.5

def reynolds_number(rho: float, V: float, L: float, mu: float) -> float:
    """
    Compute the Reynolds number for flow over a surface or in a pipe.

    Parameters
    ----------
    rho : float
        Fluid density in kilograms per cubic meter (kg/m^3).
    V : float
        Flow velocity in meters per second (m/s).
    L : float
        Characteristic length in meters (m).
    mu : float
        Dynamic viscosity of the fluid in kilograms per meter-second (kg/(m·s)).

    Returns
    -------
    float
        Reynolds number (dimensionless).
    """
    return (rho * V * L) / mu

####################################################################################################
"""Constants"""

# Calculate air density using ideal gas law (2 marks)
R_air = 287.057 # J/(kg·K)
p_air = 101389 # Pa                                    
T_air = 273.15 + 23.1 # K                                  
rho_air = ideal_gas_law(p=p_air, T=T_air) # kg/m^3

# Calculate dynamic viscosity of air μ using Sutherland’s law [2] (2 marks)
mu_air = dynamic_viscosity(T_air) # kg/(m·s)

# Compute reference dynamic pressure q∞ (1 marks)
Delta_P_b = 190 # Pa                                        TODO: FIND EXPERIMENTALLY   
q_inf = dynamic_pressure(Delta_P_b) # Pa

# Compute reference static pressure p∞ and compare it with the measurement in the tunnel (1 marks)
p_total_inf = 101375 # Pa                                   TODO: FIND EXPERIMENTALLY 
ref_static_p_inf = p_total_inf - q_inf # Pa

# Compute reference free-stream velocity U∞ and the Reynolds number (1 marks)
gamma_air = 1.4
v_air= bernoulli(p_total_1=p_total_inf, p2=ref_static_p_inf, rho=rho_air) # m/s
a_air = np.sqrt(R_air * gamma_air * T_air) # m/s                        
M_air = v_air / a_air # Mach number
chord_length = 0.16 # m 
Re = reynolds_number(rho_air, v_air, chord_length, mu_air)

####################################################################################################
"""Run this file"""

def main():
    # --- Print results ---
    print("\n=== Airflow Calculations ===")
    print(f"Air density (rho): {rho_air:.4f} kg/m³")
    print(f"Dynamic viscosity (mu): {mu_air:.6e} kg/(m·s)")
    print(f"Reference dynamic pressure (q∞): {q_inf:.4f} Pa")
    print(f"Free-stream velocity (U∞): {v_air:.4f} m/s")
    print(f"Reynolds number (Re): {Re:.2e}")
    print("============================\n")

if __name__ == "__main__":
    main() # makes sure this only runs if you run *this* file, not if this file is imported somewhere else