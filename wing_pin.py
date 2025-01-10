import numpy as np

forward_return_module_length = 100  # mm

model_factor = 1.5
maximum_lift = 998.44  # N
maximum_drag = 152.01  # N
moment = maximum_lift * model_factor * forward_return_module_length / 4  # Nmm

radius = 5
screw_size = 3
thickness = (radius * 2 - screw_size)/2

# Wing pin
theta = np.sqrt((moment / (np.pi * thickness * radius**2))**2 + (maximum_drag * model_factor / (np.pi*(radius**2-(radius-thickness)**2)))**2)
print(theta)
print(240/theta)



# Screw failure
screw_theta = maximum_lift * model_factor / ((screw_size * 0.8 / 2)**2*np.pi)
print(screw_theta)
print(640/screw_theta)
