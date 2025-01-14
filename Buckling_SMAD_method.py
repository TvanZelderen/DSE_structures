import numpy as np
import matplotlib.pyplot as plt

# Material Properties
# This dictionary contains mechanical properties for various materials.
# Each material includes values for:
# - Elastic modulus (E) in GPa
# - Ultimate tensile strength (UTS_MPa) in MPa
# - Yield strength (Yield_MPa) in MPa
# - Poisson's ratio (Poisson_ratio)
# - Density (Density_kg_m3) in kg/m³
materials = {
    "Aluminium 6061-T6": {
        "E": 68,
        "UTS_MPa": 290,
        "Yield_MPa": 240,
        "Poisson_ratio": 0.33,
        "Density_kg_m3": 2710
    },
    "Aluminium 2219-T62": {
        "E": 73.1,
        "UTS_MPa": 414,
        "Yield_MPa": 290,
        "Poisson_ratio": 0.33,
        "Density_kg_m3": 2840
    },
    "Aluminium 7075": {
        "E": 71,
        "UTS_MPa": 524,
        "Yield_MPa": 448,
        "Poisson_ratio": 0.33,
        "Density_kg_m3": 2800
    },
    "Steel 17-4PH": {
        "E": 196,
        "UTS_MPa": 660,
        "Yield_MPa": 970,
        "Poisson_ratio": 0.291,
        "Density_kg_m3": 7860
    },
    "Steel PH 15-7 Mo": {
        "E": 200,
        "UTS_MPa": 896,
        "Yield_MPa": 372,
        "Poisson_ratio": 0.28,
        "Density_kg_m3": 7804
    },
    "Ti6Al4V Grade": {
        "E": 114,
        "UTS_MPa": 1000,
        "Yield_MPa": 910,
        "Poisson_ratio": 0.342,
        "Density_kg_m3": 4420
    }
}

# Geometry and Mass
# Define the structural geometry and properties
material_type = "Aluminium 6061-T6"  # Select material from the dictionary
L = 0.79  # Length in meters
D = 0.29  # Diameter in meters
R = D / 2  # Radius in meters
m = 35  # Mass in kg
circumference = 2 * np.pi * R  # Circumference of the structure in m
No_stringers = 8  # Number of stringers
Stringer_spacing = circumference / No_stringers  # Spacing between stringers

# Load and Safety Factors
# Define loads, natural frequencies, and safety factors
Axial_acc = 13.8  # Axial acceleration in g
Lateral_acc = 3.1  # Lateral acceleration in g
fnat_ax = 25  # Axial natural frequency in Hz
fnat_lat = 10  # Lateral natural frequency in Hz
Yield_safety_factor = 1.5  # Safety factor for yield load

# Compute limit loads
Limit_load_lateral = Lateral_acc * 9.80665 * m  # Lateral load in N
Limit_load_axial = Axial_acc * 9.80665 * m  # Axial load in N

# Calculate bending moment and equivalent load
Bending_moment = Limit_load_lateral * (L / 2)  # Lateral load applied at CG
P_eq = Limit_load_axial + 2 * Bending_moment / R  # Equivalent load
P_safe = P_eq * Yield_safety_factor  # Safe ultimate load

# Material Constants
# Extract material properties for selected material
E = materials[material_type]["E"] * 10**9  # Elastic modulus in Pa
UTS = materials[material_type]["UTS_MPa"] * 10**6  # Ultimate tensile strength in Pa
Y = materials[material_type]["Yield_MPa"] * 10**6  # Yield strength in Pa
nu = materials[material_type]["Poisson_ratio"]  # Poisson's ratio

# Initial Skin Thickness
fixed_increment = 0.0000001  # Thickness increment in meters
MS_buckling = -1  # Initialize margin of safety below threshold

# Rigidity Requirements
A_req = (fnat_ax / 0.25)**2 * m * L / E  # Axial rigidity requirement
t_rigidity_ax = A_req / (np.pi * R * 2)  # Required thickness for axial rigidity
I_req_lat = (fnat_lat / 0.56)**2 * m * L**3 / E  # Lateral rigidity requirement
t_rigidity_lat = I_req_lat / (np.pi * R**3)  # Required thickness for lateral rigidity
t_rigidity = max(t_rigidity_ax, t_rigidity_lat)  # Maximum rigidity requirement

# Strength Requirements
t_req = P_safe / (2 * np.pi * R * Y)  # Required thickness for Yield

# Iteration for Buckling
t = 0.000001  # Initial thickness
thicknesses = []  # Store thicknesses during iterations
margins = []  # Store margin of safety values
boom_areas = []  # Store boom areas
cnt = 0  # Counter for iterations

# Iterative process to determine required skin thickness
while MS_buckling < 0:
    I_skin = np.pi * R**3 * t  # Skin moment of inertia
    I_stringers = max(I_req_lat - I_skin, 0)  # Stringer moment of inertia
    if I_stringers != 0:
        cnt += 1

    theta = np.linspace(0, 2 * np.pi, No_stringers, endpoint=False)
    y_coords = R * np.sin(theta)
    distance_squared_x = y_coords**2
    A_stringer = I_stringers / np.sum(distance_squared_x)  # Stringer area
    A_total = A_stringer * No_stringers + 2 * np.pi * R * t  # Total area

    # Buckling coefficient (K)
    Z = Stringer_spacing**2 / (R * t) * np.sqrt(1 - nu**2)
    if R / t < 500:
        K_value = 0.4292 * Z + 1.4337
    elif R / t < 700:
        K_value = 0.3174 * Z + 2.1836
    else:
        K_value = 0.1874 * Z + 4.1155

    Crippling_stress = K_value * np.pi**2 * E * (t / Stringer_spacing)**2 / (12 * (1 - nu**2))  # Crippling stress
    MS_buckling = Crippling_stress * A_total / P_safe - 1  # Margin of safety

    thicknesses.append(t * 1000)  # Convert thickness to mm
    margins.append(MS_buckling)
    boom_areas.append(A_stringer * 10**6)  # Convert area to mm²

    t += fixed_increment  # Increment thickness

# Check rigidity requirements
if t < t_rigidity:
    print("The thickness is not enough to meet the rigidity requirements:", t_rigidity)

# Output results
print("The area of the booms is 0 from iteration", cnt, "for a skin thickness of:", thicknesses[cnt])
print("Final Skin Thickness: {:.4f} mm".format(thicknesses[-1]))
print("Crippling Stress: {:.2f} Pa".format(Crippling_stress))

# Plot Results
plt.figure(figsize=(12, 6))

# Margin of Safety vs Thickness
plt.subplot(1, 2, 1)
plt.plot(thicknesses, margins, label="Margin of Safety")
plt.scatter(thicknesses, margins, color='blue', s=10)
plt.axhline(0, color='r', linestyle='--', label="Target MS = 0")
plt.xlabel("Skin Thickness (mm)")
plt.ylabel("Margin of Safety")
plt.title("Margin of Safety vs Skin Thickness")
plt.legend()
plt.grid()

# Boom Area vs Thickness
plt.subplot(1, 2, 2)
plt.plot(thicknesses, boom_areas, label="Boom Area", color='green')
plt.scatter(thicknesses, boom_areas, color='green', s=10)
plt.xlabel("Skin Thickness (mm)")
plt.ylabel("Boom Area (mm²)")
plt.title("Boom Area vs Skin Thickness")
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()