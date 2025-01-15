import numpy as np

launch_ring_mass = 33  # kg
landing_g = 1.5
deceleration_g = 0.3
G = 9.81
SF = 1.5

rear_ski_length = 0.2  # m
rear_ski_width = 0.04  # m
front_ski_length = 0.11  # m
front_ski_width = 0.04  # m

rear_leg_length = 0.3265  # m
front_leg_length = 0.1  # m
drag_angle = np.deg2rad(120)  # rad
normal_angle = np.deg2rad(127.76)  # rad

yield_strength_aluminium = 240e6  # Pa
density_aluminium = 2710  # kg/m^3

maximum_moment = landing_g * G * SF * launch_ring_mass * rear_ski_length / 4  # Nm
print(f"Maximum moment: {maximum_moment:.2f} Nm")

required_ski_thickness = np.sqrt(
    6 / yield_strength_aluminium / rear_ski_width * maximum_moment
)  # m
print(f"Required ski thickness: {required_ski_thickness*1000:.2f} mm")

rear_ski_volume = rear_ski_length * rear_ski_width * required_ski_thickness
front_ski_volume = front_ski_length * front_ski_width * required_ski_thickness

rear_ski_mass = rear_ski_volume * density_aluminium
front_ski_mass = front_ski_volume * density_aluminium
print(f"Mass of ski's: {rear_ski_mass*2 + front_ski_mass*2:.3f} kg")

leg_normal_force = (
    np.cos(drag_angle) * deceleration_g * G + np.cos(normal_angle) * landing_g * G * SF
)  # N
deceleration_force = np.sin(drag_angle) * deceleration_g * G * SF
normal_force = np.sin(normal_angle) * landing_g * G * SF


t = 0.001  # m
d = 0.008  # m

Ixx = np.pi * t * d**3 / 8
A = np.pi * ((d / 2) ** 2 - (d / 2 - t) ** 2)

normal_stress = leg_normal_force / A
deceleration_stress = deceleration_force * rear_leg_length * d / 2 / Ixx
normal_stress = normal_force * rear_leg_length * d / 2 / Ixx

stress = normal_stress + np.sqrt(deceleration_stress**2 + normal_stress**2)  # Pa
print(f"Stress: {stress/10**6:.2f} MPa")
rear_leg_volume = A * rear_leg_length
rear_leg_mass = rear_leg_volume * density_aluminium
print(f"Leg mass: {rear_leg_mass:.2f} kg")
front_leg_volume = A * front_leg_length
front_leg_mass = front_leg_volume * density_aluminium

total_mass = (
    rear_ski_mass * 2 + front_ski_mass * 2 + rear_leg_mass * 2 + front_leg_mass * 2
)
print(f"Total landing system mass: {total_mass:.3f} kg")
