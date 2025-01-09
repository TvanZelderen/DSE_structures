import numpy as np

diameter = 290  # mm
radius = diameter / 2
forward_return_module_length = 100  # mm
A = np.pi * radius**2

skin_thickness = 1  # mm
n_stringers = 8
unsupported_skin_length = (2 * np.pi * radius) / n_stringers
stringer_area = 1  # mm^2

neutral_line_offset = np.sin(np.deg2rad(-30)) * radius
d1 = np.sin(np.deg2rad(67.5)) * radius - neutral_line_offset
d2 = np.sin(np.deg2rad(22.5)) * radius - neutral_line_offset
d3 = np.sin(np.deg2rad(-22.5)) * radius - neutral_line_offset
d4 = np.sin(np.deg2rad(-67.5)) * radius - neutral_line_offset

d1_b = np.sin(np.deg2rad(67.5)) * radius + radius
d2_b = np.sin(np.deg2rad(22.5)) * radius + radius
d3_b = np.sin(np.deg2rad(-22.5)) * radius + radius
d4_b = np.sin(np.deg2rad(-67.5)) * radius + radius

model_factor = 1.5
maximum_lift = 998.44  # N
maximum_drag = 152.01  # N
moment = maximum_lift * model_factor * forward_return_module_length / 1000 / 4  # Nm

# Calculating the cross sectional area for the stringers
b1 = stringer_area + 1 * 115 / 6 * (2 + d1 / d1) + 1 * 115 / 6 * (2 + d2 / d1)  # mm^2
b2 = stringer_area + 1 * 115 / 6 * (2 + d1 / d2) + 1 * 115 / 6 * (2 + d3 / d2)  # mm^2
b3 = stringer_area + 1 * 115 / 6 * (2 + d2 / d3) + 1 * 115 / 6 * (2 + d4 / d3)  # mm^2
b4 = stringer_area + 1 * 115 / 6 * (2 + d3 / d4) + 1 * 115 / 6 * (2 + d4 / d4)  # mm^2

# Centroid location from the bottom
centroid = (2 * b1 * d1_b + 2 * b2 * d2_b + 2 * b3 * d3_b + 2 * d4_b) / (
    2 * d1 + 2 * d2 + 2 * d3 + 2 * b4
)

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
ixx = ixx_1 + ixx_2 + ixx_3 + ixx_4  # mm^4

σ1_moment = moment * y1 / ixx  # N/mm^2 aka MPa
σ2_moment = moment * y2 / ixx  # N/mm^2
σ3_moment = moment * y3 / ixx  # N/mm^2
σ4_moment = moment * y4 / ixx  # N/mm^2
σ = np.array([σ1_moment, σ2_moment, σ3_moment, σ4_moment])
σ_max = np.max(np.abs(σ))

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

qs = np.array([q18, q12, q23, q34, q45, q56, q67, q78])
qs = qs + qs0
max_q = np.max(np.abs(qs))

tau = max_q / skin_thickness  # N/mm^2

von_mises = np.sqrt(σ_max**2 + 3 * tau**2)  # N/mm^2
print(von_mises)
