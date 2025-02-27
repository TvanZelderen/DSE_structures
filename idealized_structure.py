import numpy as np

diameter = 290  # mm
radius = diameter / 2
forward_return_module_length = 100  # mm
A = np.pi * radius**2

skin_thickness = 0.5  # mm
n_stringers = 8
unsupported_skin_length = (2 * np.pi * radius) / n_stringers
stringer_area = 0 # mm^2

offset = True
neutral_line_offset = np.sin(np.deg2rad(30)) * radius if offset else 0
d1 = np.sin(np.deg2rad(67.5)) * radius + neutral_line_offset
d2 = np.sin(np.deg2rad(22.5)) * radius + neutral_line_offset
d3 = np.sin(np.deg2rad(-22.5)) * radius + neutral_line_offset
d4 = np.sin(np.deg2rad(-67.5)) * radius + neutral_line_offset

# print(f"ds{d1, d2, d3, d4}")

d1_b = np.sin(np.deg2rad(67.5)) * radius + radius
d2_b = np.sin(np.deg2rad(22.5)) * radius + radius
d3_b = np.sin(np.deg2rad(-22.5)) * radius + radius
d4_b = np.sin(np.deg2rad(-67.5)) * radius + radius
# print(f"dbs{d1_b, d2_b, d3_b, d4_b}")

model_factor = 1.5
maximum_lift = 999  # N
maximum_drag = 153  # N
lift_moment = maximum_lift * model_factor * forward_return_module_length / 1000 / 4  # Nm
drag_moment = maximum_drag * model_factor * 0.3 # Nm, this assumes drag is a point force at the quarter span. Both wings have this force applied.


# Calculating the cross sectional area for the stringers
b1 = stringer_area + 1 * 115 / 6 * (2 + d1 / d1) + 1 * 115 / 6 * (2 + d2 / d1)  # mm^2
b2 = stringer_area + 1 * 115 / 6 * (2 + d1 / d2) + 1 * 115 / 6 * (2 + d3 / d2)  # mm^2
b3 = stringer_area + 1 * 115 / 6 * (2 + d2 / d3) + 1 * 115 / 6 * (2 + d4 / d3)  # mm^2
b4 = stringer_area + 1 * 115 / 6 * (2 + d3 / d4) + 1 * 115 / 6 * (2 + d4 / d4)  # mm^2

print(f"Areas: {b1, b2, b3, b4}")

# Centroid location from the bottom
centroid = (2 * b1 * d1_b + 2 * b2 * d2_b + 2 * b3 * d3_b + 2 * d4_b) / (
    2 * d1 + 2 * d2 + 2 * d3 + 2 * b4
)

# print(f"Centroid: {centroid}")

y1 = d1_b - centroid
y2 = d2_b - centroid
y3 = d3_b - centroid
y4 = d4_b - centroid
x1 = np.cos(np.deg2rad(67.5)) * radius
x2 = np.cos(np.deg2rad(22.5)) * radius

### Moment

ixx_1 = b1 * (y1) ** 2  # mm^4
ixx_2 = b2 * (y2) ** 2  # mm^4
ixx_3 = b3 * (y3) ** 2  # mm^4
ixx_4 = b4 * (y4) ** 2  # mm^4
ixx = (ixx_1 + ixx_2 + ixx_3 + ixx_4) * 2  # mm^4
print(f"Ixx: {ixx}")

iyy_1 = b1 * (x1) ** 2  # mm^4
iyy_2 = b2 * (x2) ** 2  # mm^4
iyy_3 = b3 * (x2) ** 2  # mm^4
iyy_4 = b4 * (x1) ** 2  # mm^4
iyy = (iyy_1 + iyy_2 + iyy_3 + iyy_4) * 2  # mm^4
print(f"Iyy: {iyy}")

σ1_moment = lift_moment * y1 / ixx + drag_moment * x1 / iyy # N/mm^2 aka MPa
σ2_moment = lift_moment * y2 / ixx + drag_moment * x2 / iyy # N/mm^2
σ3_moment = lift_moment * y3 / ixx + drag_moment * x2 / iyy # N/mm^2
σ4_moment = lift_moment * y4 / ixx + drag_moment * x1 / iyy # N/mm^2
σ = np.array([σ1_moment, σ2_moment, σ3_moment, σ4_moment])
σ_max = np.max(np.abs(σ))
print("y's, moments")
print(y1, y2, y3, y4)
print(σ)
print(σ_max)

### Shear

factor = -maximum_lift * model_factor / ixx

# Basic shear flow
q18 = 0
q12 = factor * b1 * y1
q23 = factor * b2 * y2 + q12
q34 = factor * b3 * y3 + q23
q45 = factor * b4 * y4 + q34
q56 = factor * b4 * y4 + q45
q67 = factor * b3 * y3 + q56
q78 = factor * b2 * y2 + q67
print("Shear flow")
print(q18, q12, q23, q34, q45, q56, q67, q78)

# Correction through moment
qs0 = (
    -q12 * (-x2 + x1) * y1
    + q12 * (y2 - y1) * x1
    - q23 * (-x2 + x2) * y2
    + q23 * (y3 - y2) * x2
    - q34 * (-x1 + x2) * y3
    + q34 * (y4 - y3) * x2
    - q45 * (x1 + x1) * y4
    + q45 * (y4 - y4) * x1
    - q56 * (x2 - x1) * y4
    + q56 * (y3 - y4) * x1
    - q67 * (x2 - x2) * y3
    + q67 * (y2 - y3) * x2
    - q78 * (x1 - x2) * y2
    + q78 * (y1 - y2) * x2
    - q18 * (x1 + x1) * y1
    + q18 * (y1 - y1) * x1
) / (2 * A)
print("Moment correction")
print(qs0)

qs = np.array([q18, q12, q23, q34, q45, q56, q67, q78])
qs = qs + qs0
print(qs)
max_q = np.max(np.abs(qs))

tau = max_q / skin_thickness  # N/mm^2

von_mises = np.sqrt(σ_max**2 + 3 * tau**2)  # N/mm^2
print(f"Skin stress: {von_mises:.2f} MPa")
print(f"Skin SM: {240 / von_mises:.2f}. Aluminium 6061 T6")
print(f"Skin thickness: {skin_thickness} mm\n")

print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

### Wing attachement stringers
print("Wing attachement stringers")

# stringer_t = 1 # mm
# stringer_r = 12 # mm
# stringer_A = np.pi * (stringer_r**2 - (stringer_r-stringer_t)**2) # mm^2
# stringer_J = np.pi * stringer_t * (stringer_r)**3 * 2 # mm^4
# stringer_I = stringer_J / 2 # mm^4

square_l = 10 # mm
square_t = 1 # mm
square_A = square_l**2 - (square_l-square_t*2)**2
square_I = (square_l**4 - (square_l-square_t*2)**4) / 12
square_J = (square_l**4 - (square_l-2*square_t)**4) / 6

# print(f"Area: {square_A/stringer_A}")
# print(f"Ixx: {square_I/stringer_I}")

# σ_normal = maximum_drag * model_factor / stringer_A
# σ_bending_lift = maximum_lift * model_factor * forward_return_module_length / 4 * stringer_r / stringer_I
# σ_bending_drag = maximum_drag * model_factor * forward_return_module_length / 4 * stringer_r / stringer_I
# σ_bending = np.sqrt(σ_bending_lift**2 + σ_bending_drag**2)
# τ = maximum_lift * model_factor * 0.25 * stringer_r / stringer_J

# von_mises_stringer = np.sqrt(σ_normal**2 - σ_normal*σ_bending + σ_bending**2 + 3*τ**2)

# print(f"Wing stringer stress: {von_mises_stringer:.2f} MPa")
# print(f"Stringer SF: {240 / von_mises_stringer:.2f}. Aluminium 6061 T6")
# print(f"Stringer thickness: {stringer_t} mm, stringer radius: {stringer_r} mm")

σ_normal = maximum_drag / 2 * model_factor / square_A
σ_bending_lift = maximum_lift / 2 * model_factor * forward_return_module_length / 4 * square_l / 2 / square_I
σ_bending_drag = (maximum_drag*0.3) * model_factor * square_l / 2 / square_I
σ_bending_wing_moment = (maximum_lift * 0.25 - 0.204) * model_factor * square_l / 2 / square_I                           # TODO: fill in torque for 1 wing, for fin comment out
σ_bending = np.sqrt((σ_bending_lift+σ_bending_wing_moment)**2 + σ_bending_drag**2)
τ = maximum_lift / 2 * model_factor * 0.25 * square_l / 2 / square_J

von_mises_stringer = np.sqrt(σ_normal**2 - σ_normal*σ_bending + σ_bending**2 + 3*τ**2)

# print(f"\nSquare:")
print(f"Wing stringer stress: {von_mises_stringer:.2f} MPa")
print(f"Stringer SM: {240 / von_mises_stringer:.2f}. Aluminium 6061 T6")
print(f"Stringer thickness: {square_t} mm, stringer side length: {square_l} mm")

