import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Material Properties Dictionary
# Contains properties of common aerospace materials such as:
# Young's Modulus (E), Ultimate Tensile Strength (UTS), Yield Strength, Poisson's ratio, and Density.
materials = {
    "Aluminium 6061-T6": {
        "E": 68,  # Young's Modulus (GPa)
        "UTS_MPa": 290,  # Ultimate Tensile Strength (MPa)
        "Yield_MPa": 240,  # Yield Strength (MPa)
        "Poisson_ratio": 0.33,  # Poisson's Ratio
        "Density_kg_m3": 2710  # Density (kg/m³)
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

# Geometry and Material Selection
material_type = "Aluminium 6061-T6"  # Selected material
L = 0.79  # Length of the structure (meters)
D = 0.29  # Diameter of the structure (meters)
R = D / 2  # Radius of the structure (meters)
t_skin = 0.0005  # Initial skin thickness (meters)
circumference = 2 * np.pi * R  # Circumference of the structure (meters)
No_stringers = 8  # Number of stringers
Stringer_spacing = circumference / No_stringers  # Spacing between stringers

# Length of each stage (center of gravity positions)
L_1_return_aft = 1.731 # Stage 1 length (meters)
L_2_payload = 1.382    # Stage 2 length (meters) 
L_3_return_fwd = 1.133 # Stage 3 length (meters) 

# Individual stage lengths
L_1 = 0.28   # Length of return aft stage (meters)
L_2 = 0.18   # Length of payload stage (meters)
L_3 = 0.13   # Length of return forward stage (meters)

# Mass of each stage
mass_1 = 32.5  # Stage 1 mass (kg)
mass_2 = 26.34  # Stage 2 mass (kg)
mass_3 = 18.42  # Stage 3 mass (kg)

# Center of gravity for each stage
cg1 = L_1_return_aft / 2
cg2 = L_2_payload / 2
cg3 = L_3_return_fwd / 2

# Load and Safety Parameters
safety_factor = 1.5  # Safety factor
Axial_acc = 13.8 * 9.81 * safety_factor  # Axial acceleration (m/s²)
Lateral_acc = 3.1 * 9.81 * safety_factor  # Lateral acceleration (m/s²)

# Generate stringer positions (circular arrangement)
theta = np.linspace(0, 2 * np.pi, No_stringers, endpoint=False)
x_coords = R * np.cos(theta)  # X-coordinates of stringers
y_coords = R * np.sin(theta)  # Y-coordinates of stringers

# Rotate to position one stringer at the bottom (-y axis)
offset = np.pi / 2
theta_rotated = (theta + offset) % (2 * np.pi)
x_coords_rotated = R * np.cos(theta_rotated)
y_coords_rotated = R * np.sin(theta_rotated)

# Calculate squared distance for parallel axis theorem
distance_squared_x = np.sum(y_coords_rotated**2)

# Plot stringer positions
plt.figure(figsize=(8, 8))
circle = plt.Circle((0, 0), R, color='blue', fill=False, linestyle='--', label='Circle Boundary')
plt.gca().add_artist(circle)
plt.scatter(x_coords_rotated, y_coords_rotated, color='red', label='Stringers')

for i, (x, y) in enumerate(zip(x_coords_rotated, y_coords_rotated)):
    plt.text(x, y, f'{i+1}', fontsize=12, ha='center', va='center', color='black')

plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlim(-R * 1.2, R * 1.2)
plt.ylim(-R * 1.2, R * 1.2)
plt.title("Stringer Positions Around the Circle")
plt.xlabel("X-coordinate (m)")
plt.ylabel("Y-coordinate (m)")
plt.legend()
plt.grid(True)
plt.show()

# Material constants for the selected material
E = materials[material_type]["E"] * 10**9  # Young's Modulus (Pa)
rho = materials[material_type]["Density_kg_m3"]  # Density (kg/m³)
Y = materials[material_type]["Yield_MPa"] * 10**6  # Yield Strength (Pa)
nu = materials[material_type]["Poisson_ratio"]  # Poisson's Ratio
Le_L_ratio = 2  # Effective length ratio (to be verified)

# Maximum moments for each stage

#Stage 1 
Max_Shear_1=Lateral_acc*mass_1
Max_Moment_1=Max_Shear_1*cg1
#Stage 2
Max_Shear_2=Lateral_acc*mass_2
Max_Moment_2=Max_Shear_2*cg2
#Stage 3
Max_Shear_3=Lateral_acc*mass_3 
Max_Moment_3=Max_Shear_3*cg3


# Initialize total stringers mass
total_stringers_mass = 0

# Function to calculate thickness based on required moment of inertia
def calculate_thickness(width, required_moment_of_inertia):
    """
    Calculate the thickness of a stringer based on the required moment of inertia.
    Ensures that thickness is non-negative and returns None if invalid.
    """
    try:
        t = (width - (width**4 - required_moment_of_inertia * 12)**0.25) / 2
        if t < 0:  # Invalid thickness
            return None
        return t
    except ValueError:  # Catch math domain errors
        return None

# STAGE 1 Calculations
a1 = 0.015  # Width of stringer (meters)
flag = False  # Flag for thickness adjustment

# Parameters for stage 1
B = np.pi**2 * E / (Le_L_ratio * L_1)**2  # Buckling coefficient
C = Axial_acc * mass_1 / No_stringers  # Axial force per stringer
D = Max_Moment_1 * R / distance_squared_x  # Moment-induced force
I = (C + D) / B  # Required moment of inertia for stability

t_0 = calculate_thickness(a1, I)  # Calculate thickness
if t_0 is None or t_0 < 0.0005:
    flag = True
    t_0 = 0.0005  # Minimum thickness constraint

Area_stringer = a1**2 - (a1 - 2 * t_0)**2  # Stringer cross-sectional area
Stringers_mass = Area_stringer * L_1 * rho * No_stringers  # Mass of stringers for stage 1
total_stringers_mass += Stringers_mass  # Add to total stringers mass

if flag:
    print("Thickness is less than 0.5mm. Stringer thickness set to 0.5mm")
print(f"Thickness for Stage 1 (t_0): {t_0 * 1000:.2f} mm")
print(f"Area of Stringer for Stage 1: {Area_stringer:.6f} m²")
print(f"Stringers Mass for Stage 1: {Stringers_mass:.2f} kg")

# STAGE 2 Calculations
a2 = 0.01  # Width of stringer (meters)
flag = False

# Parameters for stage 2
B = np.pi**2 * E / (Le_L_ratio * L_2)**2
C = Axial_acc * mass_2 / No_stringers
D = Max_Moment_2 * R / distance_squared_x
I = (C + D) / B

t_0 = calculate_thickness(a2, I)
if t_0 is None or t_0 < 0.0005:
    flag = True
    t_0 = 0.0005

Area_stringer = a2**2 - (a2 - 2 * t_0)**2
Stringers_mass = Area_stringer * L_2 * rho * No_stringers
total_stringers_mass += Stringers_mass

print("\n")
if flag:
    print("Thickness is less than 0.5mm. Stringer thickness set to 0.5mm")
print(f"Thickness for Stage 2 (t_0): {t_0 * 1000:.2f} mm")
print(f"Area of Stringer for Stage 2: {Area_stringer:.6f} m²")
print(f"Stringers Mass for Stage 2: {Stringers_mass:.2f} kg")

# STAGE 3 Calculations
a3 = 0.01  # Width of stringer (meters)
flag = False

# Parameters for stage 3
B = np.pi**2 * E / (Le_L_ratio * L_3)**2
C = Axial_acc * mass_3 / No_stringers
D = Max_Moment_3 * R / distance_squared_x
I = (C + D) / B

t_0 = calculate_thickness(a3, I)
if t_0 is None or t_0 < 0.0005:
    flag = True
    t_0 = 0.0005

Area_stringer = a3**2 - (a3 - 2 * t_0)**2
Stringers_mass = Area_stringer * L_3 * rho * No_stringers
total_stringers_mass += Stringers_mass

print("\n")
if flag:
    print("Thickness is less than 0.5mm. Stringer thickness set to 0.5mm")
print(f"Thickness for Stage 3 (t_0): {t_0 * 1000:.2f} mm")
print(f"Area of Stringer for Stage 3: {Area_stringer:.6f} m²")
print(f"Stringers Mass for Stage 3: {Stringers_mass:.2f} kg")

# Skin Mass Calculation
skin_mass = 2 * np.pi * R * t_skin * L * rho
total_structure_mass = total_stringers_mass + skin_mass

print("\n")
print(f"Total Stringers Mass: {total_stringers_mass:.2f} kg")
print(f"Skin Mass: {skin_mass:.2f} kg")
print(f"Total Structure Mass: {total_structure_mass:.2f} kg")
