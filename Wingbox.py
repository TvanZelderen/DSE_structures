import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

# Assumptions: Thin walled, evenly distributed loads, buckling boundary conditions(!), constant shear (in wrong direction I believe)


# Inputs:
L_max = 998.44/2           # [N] TBC
D_max = 152.01/2           # [N] TBC
T_max = -8.68/2          # [Nm], pitch down. TBC

bh = 0.6                # [m]
w = 66.8/1000          # [m] TBC
h = 9/1000           # [m] TBC
d = 32.5/1000               # [m] TBC, this is de distance from wingbox center to Vy application point
# TODO: This is an assumption ^ (quarter chord)

# sigma_y = 787*10**6                 # composite
# tau_y = 128*10**6
# v = 0.3
# E = 125*10**9
# density = 1580

sigma_y = 1000*10**6                # steel 17-4ph https://protoxyz.com/materials/metal/Stainless_Steel_17-4_PH
tau_y = 827*10**6
v = 0.28
E = 193*10**9
density = 7800

# sigma_y = 276*10**6                # alu 6061-t6 https://asm.matweb.com/search/specificmaterial.asp?bassnum=ma6061t6
# tau_y = 207*10**6
# v = 0.33
# E = 68.9*10**9
# density = 2700

# TODO: Assume constant shear throughout z

# Safety factors
A = 1.5                 # Design factor, accounts for uncertainty in models, TBC
# TODO: Investigate gust loads for small UAVs


# Setup:
L = A*L_max
D = A*D_max
T = A*T_max + L*d     # Subtracts the torsion component of Vy, which acts at a certain distance from the wingbox center
# print("Torque_z: ", T)


# Reaction forces and moments:
R_y = L
X_m = L * (1/2) * bh
R_x = D
Y_m = R_x * (1/2) * bh

# Shear forces, bending moments, and torsion at root:
V_y = -R_y
V_x = R_x
M_x = X_m
M_y = -Y_m

# print("Shear and moment", V_y, V_x, M_x, M_y)

qT = -T/(2*w*h)         #negative, because shearflow is taken cw positive
# TODO: double check signs, especially qT!


# Functions:
def I(t):
    I_xx = (w * h ** 3) / (12) - ((w - 2*t) * (h - 2*t) ** 3) / 12
    I_yy = (h * w ** 3) / (12) - ((h - 2*t) * (w - 2*t) ** 3) / 12
    return I_xx, I_yy

dx = 0.01/1000
x = np.arange(-w/2, w/2 +dx, dx)
dy = 0.01/1000
y = np.arange(-h/2 ,h/2 +dy, dy)
def s_bending(I_xx, I_yy):
    s_bending_top = np.max((M_y * x) / I_yy + (M_x * h/2) / I_xx)
    s_bending_side = np.max((M_y * w/2) / I_yy + (M_x * y) / I_xx)
    s_bending_max = s_bending_side
    return max(s_bending_side, s_bending_top)

def s_shear(I_xx, I_yy, t):
    q12b = - (V_y / I_xx) * (0.5 * t * h * s1 - 0.5 * t * s1 ** 2) - (V_x * w * t * s1) / (2 * I_yy)
    q2b = q12b[-1]
    q23b = q2b + (V_y / I_xx) * (0.5 * t * h * s2) - (V_x / I_yy) * (0.5 * t * w * s2 - 0.5 * t * s2 ** 2)
    q3b = q23b[-1]
    q34b = q3b - (V_y / I_xx) * (-0.5 * t * h * s3 + 0.5 * t * s3 ** 2) + (V_x * w * t * s3) / (2 * I_yy)
    q4b = q34b[-1]
    q41b = q4b - (V_y / I_xx) * (0.5 * t * h * s4) - (V_x / I_yy) * (-0.5 * t * w * s4 + 0.5 * t * s4 ** 2)
    q1b = q41b[-1]

    Mb = np.sum(q12b * ds * 0.5 * w) + np.sum(q23b * ds * 0.5 * h) + np.sum(q34b * ds * 0.5 * w) + np.sum(
        q41b * ds * 0.5 * h)
    qs0 = (-Mb) / (2 * w * h)

    q12 = q12b + qs0 + qT
    q23 = q23b + qs0 + qT
    q34 = q34b + qs0 + qT
    q41 = q41b + qs0 + qT

    q_max = max(
        np.max(np.abs(q12)),
        np.max(np.abs(q23)),
        np.max(np.abs(q34)),
        np.max(np.abs(q41) )
    )

    s_shear = q_max / t
    return s_shear


def mass(t_selected1, t_selected2, t_selected3, t_selected4, bh = bh, w = w, h = h):
    A_cross_section = w * h - (w - t_selected1 - t_selected3) * (h - t_selected2 - t_selected4)
    m_stiffeners = A_cross_section * bh * density
    return m_stiffeners


# TODO: Assumption: Max stress is applied for both shear and compression buckling (conservative)
def cbuckling(t, b, C = 4):
    # Top wall: side = 0, Aft wall: side = 1
    I_xx, I_yy = I(t)
    s_b = s_bending(I_xx,I_yy)
    t_cb = math.sqrt((s_b * 12 * (1 - v ** 2) * b ** 2) / (C * math.pi ** 2 * E))
    if t_cb>t:
        while t_cb > t:
            t += 0.01/1000
            I_xx, I_yy = I(t)
            s_b = s_bending(I_xx, I_yy)
            t_cb = math.sqrt((s_b * 12 * (1 - v ** 2) * b ** 2) / (C * math.pi ** 2 * E))
    return t_cb


def sbuckling(t, b):
    k0 = 5.34 + 4 / (bh / b) ** 2
    I_xx, I_yy = I(t)
    s_s = s_shear(I_xx, I_yy, t)
    t_sb = math.sqrt( (s_s * 12 * (1 - v**2) * b**2) / (k0 * math.pi**2 * E))
    if t_sb>t:
        while t_sb > t:
            t += 0.01/1000
            I_xx, I_yy = I(t)
            s_s = s_shear(I_xx, I_yy, t)
            t_sb = math.sqrt( (s_s * 12 * (1 - v**2) * b**2) / (k0 * math.pi**2 * E))
    return t_sb


def combined_buckling(t, b, C=4):
    k0 = 5.34 + 4 / (bh / b) ** 2
    I_xx, I_yy = I(t)
    s_s = s_shear(I_xx, I_yy, t)
    s_b = s_bending(I_xx, I_yy)
    s_cr = (C * math.pi ** 2 * E * t ** 2) / (12 * (1 - v ** 2) * b ** 2)
    tau_cr = (k0 * math.pi**2 * E * t**2 )/( 12*(1-v**2) * b**2 )
    ratio = s_b/s_cr + (s_s/tau_cr)**2
    # print("Bending stress ", s_b, "Shear stress ", s_s, "s_cr ", s_cr, "tau_cr ", tau_cr)
    while ratio > 1:
        t += 0.01 / 1000
        I_xx, I_yy = I(t)
        s_s = s_shear(I_xx, I_yy, t)
        s_b = s_bending(I_xx, I_yy)
        s_cr = (C * math.pi ** 2 * E * t ** 2) / (12 * (1 - v ** 2) * b ** 2)
        tau_cr = (k0 * math.pi ** 2 * E * t ** 2) / (12 * (1 - v ** 2) * b ** 2)
        ratio = s_b / s_cr + (s_s / tau_cr) ** 2
    return t, ratio


# Bending calculations:
t_b = 0
s_b = sigma_y + 1

while sigma_y<s_b:
    t_b += (0.01/1000)
    I_xx, I_yy = I(t_b)
    s_b = s_bending(I_xx, I_yy)
print("Bending moment top skin thickness required = ", t_b*1000, "mm")





# Shear calculations:
ds = 1/1000000
s1 = np.arange(0,h+ds,ds)
s2 = np.arange(0,w+ds,ds)
s3 = np.arange(0,h+ds,ds)
s4 = np.arange(0,w+ds,ds)

t_s = 0
s_s = tau_y + 1

while tau_y<s_s:
    t_s+=(0.01/1000)
    I_xx, I_yy = I(t_s)
    s_s = s_shear(I_xx, I_yy, t_s)
print("Shear thickness required = ", t_s*1000, "mm")
# TODO: Double check signs of shear flow


# Buckling calculations:

#Compression buckling:

# Top wall: side = 0, Aft wall: side = 1
t_cb_top = cbuckling(t_b, w)
t_cb_side = cbuckling(t_b, h)
print("Top skin compression buckling t [mm]: ", t_cb_top*1000, ", Side skin compression buckling t [mm]: ", t_cb_side*1000)

t_sb_top = sbuckling(t_s , w)
t_sb_side = sbuckling(t_s , h)
print("Top skin shear buckling t [mm]:", t_sb_top*1000, ", Side skin shear buckling t [mm]: ", t_sb_side*1000)

# Top wall: side = 0, Aft wall: side = 1
t_combined_b_top, ratio_top = combined_buckling(t=max(t_sb_top, t_b, t_cb_top), b = w)
t_combined_b_side, ratio_side = combined_buckling(t=max(t_sb_side, t_b, t_cb_side), b = h)
print("Top skin combined buckling t [mm]: ", t_combined_b_top*1000, "Ratio top skin: ", ratio_top)
print("Side skin combined buckling t [mm]: ", t_combined_b_side*1000, "Ratio side skin: ", ratio_side)




# TODO: Add progression plot of thickness calculation

t_selected1 = float(input("Choose thicknes [mm] for top wall. Select value larger than required thickness: "))
t_selected2 = float(input("Choose thicknes [mm] for right wall. Select value larger than required thickness: "))
t_selected3 = float(input("Choose thicknes [mm] for bottom wall. Select value larger than required thickness: "))
t_selected4 = float(input("Choose thicknes [mm] for left wall. Select value larger than required thickness: "))

mass = mass(t_selected1/1000, t_selected2/1000, t_selected3/1000, t_selected4/1000)
print("Mass per wing [kg]: ", mass, "Total mass [kg]: ", 2*mass)





