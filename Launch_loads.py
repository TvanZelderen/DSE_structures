import math
import numpy as np
import matplotlib.pyplot as plt

pi = math.pi
g = 9.81

sigma_y = 240 * 10 ** 6  # alu 6061-T6
E = 68 * 10 ** 9  # alu 6061-T6
density = 2710

################## ADJUSTABLE:
safety_factor = 1.5
r = 0.29 / 2
t_cylinder = 0.5 / 1000

# Starting widths
w1 = 0.013
w2 = 0.01
w3 = 0.01

# Module lengths:
l_module1 = 312/1000
l_module2 = 167/1000
l_module3 = 102/1000

# Length of each stage (center of gravity positions)
L_1_return_aft = 1.831  # Stage 1 length (meters)
L_2_payload = 1.44994  # Stage 2 length (meters)
L_3_return_fwd = 1.21277  # Stage 3 length (meters)

# Mass of each stage
# Mass of each stage
mass_1 = 32.5  # Stage 1 mass (kg)
mass_2 = mass_1 - 7.384  # Stage 2 mass (kg)
mass_3 = mass_2 - 1.085 - 5.3152  # Stage 3 mass (kg)
print("Module masses: ", mass_1, mass_2, mass_3)

# Center of gravity for each stage
cg1 = L_1_return_aft / 2
cg2 = L_2_payload / 2
cg3 = L_3_return_fwd / 2

n_axial = 13.8
n_lateral = 3.1
I_cylinder = (pi * r ** 4) / 4 - (pi * (r - t_cylinder) ** 4) / 4  # skin

# moment for each module
M_max1 = mass_1 * g * n_lateral * cg1
M_max2 = mass_2 * g * n_lateral * cg2
M_max3 = mass_3 * g * n_lateral * cg3

mod_count = 0
error = False

d1_s = np.cos(np.deg2rad(0)) * r
d2_s = np.cos(np.deg2rad(45)) * r
d3_s = np.cos(np.deg2rad(90)) * r
d4_s = d2_s
print(d1_s, d2_s, d3_s, d4_s)


def stiffener_dimensions(w, M_x, l_module, m_supporting, mod_count=mod_count, error=error):
    mod_count += 1
    l_effective = 2 * l_module  # TODO: Check assumption!

    cnt = []
    I_square_list = []
    t_list = []
    P_list = []
    I_req_list = []

    cnt_nr = 0
    t_square = 0
    delta = 1
    while delta > 0.001 * 10 ** (-9) and error == False:
        cnt_nr += 1
        cnt.append(cnt_nr)

        t_square += 0.01 / 1000
        t_list.append(t_square)

        I_square = (w ** 4) / 12 - ((w - 2 * t_square) ** 4) / 12
        I_square_list.append(I_square)

        A_square = w * w - (w - 2 * t_square) ** 2
        I_stiffeners_combined = I_square * 8 + 2 * A_square * (d1_s ** 2 + d2_s ** 2 + d3_s ** 2 + d4_s ** 2)
        I_total = I_cylinder + I_stiffeners_combined
        sigma_b = (M_x * r) / (I_total)

        P_axial = (m_supporting * g * n_axial * safety_factor) / 8
        P_lateral = sigma_b * A_square
        P_eq = P_axial + P_lateral
        P_list.append(P_eq)

        I_req = (P_eq * l_effective ** 2) / (E * pi ** 2)
        I_req_list.append(I_req)

        delta = I_req - I_square
        if t_square * 2 > w:
            print("ERROR for module ", mod_count, ": SELECT LARGER WIDTH! I_square: ", I_square, "I_req: ", I_req)
            error = True

        test = sigma_b
    return cnt, I_square_list, t_list, P_list, I_req_list, test


def mass(w, t_square, l_module):
    A_square = w * w - (w - 2 * t_square) ** 2
    m_stiffeners = 8 * A_square * l_module * density
    return m_stiffeners, A_square


cnt_mod1, I_square_mod1, t_mod1, P_mod1, I_req_mod1, test1 = stiffener_dimensions(w=w1, M_x=M_max1, l_module=l_module1,
                                                                                  m_supporting=32.5)
if error == False:
    cnt_mod2, I_square_mod2, t_mod2, P_mod2, I_req_mod2, test2 = stiffener_dimensions(w=w2, M_x=M_max2,
                                                                                      l_module=l_module2,
                                                                                      m_supporting=26.34)
if error == False:
    cnt_mod3, I_square_mod3, t_mod3, P_mod3, I_req_mod3, test3 = stiffener_dimensions(w=w3, M_x=M_max3,
                                                                                      l_module=l_module3,
                                                                                      m_supporting=18.42)

# print(P_mod1[-1], "N ")
# print("Test ", test1)
print("w1: ", w1, "w2: ", w2, "w3: ", w3)
print("I_req 1 [m^4]: ", I_req_mod1[-1], "I_req 2 [m^4]: ", I_req_mod2[-1], "I_req 3 [m^4]: ", I_req_mod3[-1])
print("Required thicknesses: t1 [mm] = ", t_mod1[-1]*1000, "t2 [mm] = ", t_mod2[-1]*1000, "t3 [mm] = ", t_mod3[-1]*1000)

fig, axs = plt.subplots(3, 1, figsize=(8, 12))  # 3 rows, 1 column of subplots

# # Plot for the first dataset
axs[0].plot(cnt_mod1, I_req_mod1, label="I_req_mod1")
axs[0].plot(cnt_mod1, I_square_mod1, label="I_square_mod1")
axs[0].set_title("Module 1")
axs[0].legend()

# Plot for the second dataset
axs[1].plot(cnt_mod2, I_req_mod2, label="I_req_mod2")
axs[1].plot(cnt_mod2, I_square_mod2, label="I_square_mod2")
axs[1].set_title("Module 2")
axs[1].legend()

# Plot for the third dataset
axs[2].plot(cnt_mod3, I_req_mod3, label="I_req_mod3")
axs[2].plot(cnt_mod3, I_square_mod3, label="I_square_mod3")
axs[2].set_title("Module 3")
axs[2].legend()

plt.tight_layout()

fig, ax = plt.subplots(figsize=(8, 6))  # Single plot for Module 1

# Plot required and actual second moment of area
ax.plot(cnt_mod1, I_req_mod1, label="Required Moment of Inertia ($I$)", linestyle='--', color='r')
ax.plot(cnt_mod1, I_square_mod1, label="Calculated Moment of Inertia)", linestyle='-', color='b')

# Labels and legend
ax.set_xlabel("Iteration Step")  # X-axis label
ax.set_ylabel("Moment of Inertia ($m^4$)")  # Y-axis label
ax.legend()

plt.grid(True)  # Add grid for readability
plt.show()


thickness_module_1 = t_mod1[-1]
thickness_module_2 = t_mod2[-1]
thickness_module_3 = t_mod3[-1]

# width w1,w2,w3

# t_cylinder = 0.5/1000

m_min_stiffeners1, Area_min1 = mass(w=w1, t_square=thickness_module_1, l_module=l_module1)
m_min_stiffeners2, Area_min2 = mass(w=w2, t_square=thickness_module_2, l_module=l_module2)
m_min_stiffeners3, Area_min3 = mass(w=w3, t_square=thickness_module_3, l_module=l_module3)

# Generate stringer positions (circular arrangement)
theta = np.linspace(0, 2 * np.pi, 8, endpoint=False)
x_coords = r * np.cos(theta)  # X-coordinates of stringers
y_coords = r * np.sin(theta)  # Y-coordinates of stringers

# Rotate to position one stringer at the bottom (-y axis)
offset = np.pi / 2
theta_rotated = (theta + offset) % (2 * np.pi)
x_coords_rotated = r * np.cos(theta_rotated)
y_coords_rotated = r * np.sin(theta_rotated)

print(thickness_module_1, thickness_module_2, thickness_module_3)
print(m_min_stiffeners1, Area_min1)
print(m_min_stiffeners2, Area_min2)
print(m_min_stiffeners3, Area_min3)


def plot_section(radius, width, thickness, mass, area, color, title):
    scaling_factor = 1.5
    num_stringers = 8

    theta = np.linspace(0, 2 * np.pi, num_stringers, endpoint=False)
    x_coords = radius * np.cos(theta)
    y_coords = radius * np.sin(theta)


    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot the cylinder skin
    circle = plt.Circle((0, 0), radius, color='gray', fill=False, linestyle='--')
    ax.add_patch(circle)

    # Plot stringers as hollow squares
    for x, y in zip(x_coords, y_coords):
        side_length = width * scaling_factor
        inner_side_length = side_length - 2 * thickness * scaling_factor

        outer_rect = plt.Rectangle((x - side_length / 2, y - side_length / 2),
                                   side_length, side_length, color=color, fill=True, alpha=0.5)
        inner_rect = plt.Rectangle((x - inner_side_length / 2, y - inner_side_length / 2),
                                   inner_side_length, inner_side_length, color='white', fill=True)

        ax.add_patch(outer_rect)
        ax.add_patch(inner_rect)

    # Add text annotations
    ax.text(0, 0, f"Stringer Width: {width * 1000:.2f} mm\nStringer Thickness: {thickness * 1000:.2f} mm",
            fontsize=9, color='black', ha='center', va='center', fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))

    ax.set_aspect('equal')
    ax.set_xlim(-radius * 1.2, radius * 1.2)
    ax.set_ylim(-radius * 1.2, radius * 1.2)
    ax.set_title(f"{title}\nTotal Stiffener Mass: {mass:.3f} kg\nSingle Stiffener Area: {area * 1e6:.2f} mm²")
    ax.set_xlabel("X-coordinate (m)")
    ax.set_ylabel("Y-coordinate (m)")
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
    ax.grid(True)

    plt.show()


# Plot sections
plot_section(r, w1, thickness_module_1, m_min_stiffeners1, Area_min1, 'red', 'Section 1')
plot_section(r, w2, thickness_module_2, m_min_stiffeners2, Area_min2, 'green', 'Section 2')
plot_section(r, w3, thickness_module_3, m_min_stiffeners3, Area_min3, 'blue', 'Section 3')

print(sum([m_min_stiffeners1, m_min_stiffeners2, m_min_stiffeners3]))

t1 = float(input("Choose thicknes [mm] for module 1. Select value larger than required thickness: "))
t2 = float(input("Choose thicknes [mm] for module 2. Select value larger than required thickness: "))
t3 = float(input("Choose thicknes [mm] for module 3. Select value larger than required thickness: "))
m_stiffeners1, Area1 = mass(w = w1, t_square = t1/1000, l_module=l_module1)
m_stiffeners2, Area2 = mass(w = w2, t_square= t2/1000, l_module=l_module2)
m_stiffeners3, Area3 = mass(w = w3, t_square= t3/1000, l_module=l_module3)
print("Areas: ", Area1, Area2, Area3)
print("Mass 1 [kg]: ", m_stiffeners1, "Mass 2 [kg]: ", m_stiffeners2, "Mass 3 [kg]: ", m_stiffeners3, "Combined mass [kg]: ", m_stiffeners1+m_stiffeners2+m_stiffeners3)
