import numpy as np
import math
import matplotlib.pyplot as plt


#requirement. The payload and service module shall withstand a temperature difference of +-65°C for a duration of 5 minutes, starting from a temperature of 15°C. 

delta_temp=50 #C

# Given data
L_payload = 0.102  # m (Length of payload)
D = 0.29  # m (Diameter of the cylindrical payload)
t_skin = 0.001  # m (Thickness of aluminum skin)
R = D / 2  # m (Outer radius of payload)

# Material properties (thermal conductivity in W/m-K)
k_alu = 167  # Aluminum 6061-T6
k_insulation_materials = {
    "silica": 0.024,  # Silica Aerogel
    "fiberglass": 0.036  # Fiberglass (corrected value)
}

# Thermal properties of air inside payload
rho_air = 1.225  # kg/m^3 (Density of air at sea level)
C_p_payload = 1005  # J/kg-K (Specific heat capacity of payload)

# Temperature conditions
T_initial = 15 # C (Initial temperature of payload)
T_outside_hot = T_initial+delta_temp # C (External high temperature)
T_outside_cold = T_initial-delta_temp  # C (External low temperature)
  
T_max = 50  # C (Maximum functional temperature of payload)
T_min = -10  # C (Minimum functional temperature of payload)

time_flight = 5 * 60  # s (Total launch duration in seconds)

# Function to calculate the final temperature based on insulation thickness
def calculate_temperature(t_insulator, k_insulator, T_outside, condition):
    r_aluminium = R  # Outer radius
    r_insulator = R - t_skin  # Radius after aluminum skin
    r_air = R - t_skin - t_insulator  # Inner radius after insulation
    
    # Ensure valid radius values
    if r_air <= 0:
        return float('-inf')  # Return an invalid temperature to indicate failure
    
    # Thermal resistances (cylindrical conduction model)
    Resistance_alu = np.log(r_aluminium / r_insulator) / (2 * np.pi * k_alu * L_payload)
    Resistance_insulator = np.log(r_insulator / r_air) / (2 * np.pi * k_insulator * L_payload)
    
    # Total radial thermal resistance
    R_cylindrical = Resistance_alu + Resistance_insulator
    
    # Mass of payload module 
    mass_payload = 6.4  # kg (Assumed constant payload mass)
    
    # Lumped capacitance method
    tau = (mass_payload * C_p_payload) * R_cylindrical  # Thermal time constant
    
    # Transient heat balance (cooling or heating equation based on condition)
    if condition == "heating":
        T_final = T_outside - (T_outside - T_initial) * np.exp(-time_flight / tau)
    else:
        T_final = T_outside + (T_initial - T_outside) * np.exp(-time_flight / tau)
    
    return T_final

# Select insulation material
insulation_material = "fiberglass"  # Change this to "fiberglass" to test different materials
k_insulator = k_insulation_materials[insulation_material]  # Get thermal conductivity

# Iterate to find required insulation thickness for both conditions
t_insulator_values = np.linspace(0.0001, 0.01, 1000)
optimal_t_hot = None
optimal_t_cold = None
T_final_values_hot = []
T_final_values_cold = []

for t_insulator in t_insulator_values:
    T_final_hot = calculate_temperature(t_insulator, k_insulator, T_outside_hot, "heating")
    T_final_cold = calculate_temperature(t_insulator, k_insulator, T_outside_cold, "cooling")
    
    T_final_values_hot.append(T_final_hot)
    T_final_values_cold.append(T_final_cold)
    
    # Find the required insulation thicknesses
    if T_final_hot <= T_max and optimal_t_hot is None:
        optimal_t_hot = t_insulator
    if T_final_cold >= T_min and optimal_t_cold is None:
        optimal_t_cold = t_insulator

# Choose the thickest required insulation thickness
print(f"Optimal Insulation Thickness (Hot): {optimal_t_hot*1000:.3f} mm")
print(f"Optimal Insulation Thickness (Cold): {optimal_t_cold*1000:.3f} mm")
optimal_t = max(optimal_t_hot, optimal_t_cold)

# Plot results
plt.figure(figsize=(8, 6))
plt.plot(t_insulator_values, T_final_values_hot, label=f'Final Temperature vs. Insulation Thickness (Hot - {insulation_material})')
plt.plot(t_insulator_values, T_final_values_cold, label=f'Final Temperature vs. Insulation Thickness (Cold - {insulation_material})')
plt.axhline(y=T_max, color='r', linestyle='--', label='Maximum Functional Temperature')
plt.axhline(y=T_min, color='b', linestyle='--', label='Minimum Functional Temperature')
plt.xlabel('Insulation Thickness (m)')
plt.ylabel('Final Temperature (C)')
plt.legend()
plt.grid(True)
plt.show()

# Display result
if optimal_t:
    print(f"Minimum insulation thickness required ({insulation_material}): {optimal_t*1000:.3f} mm")
else:
    print("Even the maximum insulation thickness tested is insufficient.")
