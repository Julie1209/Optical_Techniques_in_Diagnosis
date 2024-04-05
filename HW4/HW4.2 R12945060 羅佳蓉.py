# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 22:12:26 2023

@author: cdpss
"""

import numpy as np
# import matplotlib.pyplot as plt

# n1 = 1.4
# n2 = 1
# critical_angle = np.arcsin(n2/n1)
# incident_angle = np.linspace(0, np.pi/2, num=90)
# def calculate_reflection_coefficient(n1, n2, incident_angle_array):
#     r_array = []
#     for incident_angle in incident_angle_array:
#         if 0 <= incident_angle and incident_angle < critical_angle: # 0 <= incident_angle < critical_angle
#             refractive_angle = np.arcsin(n1 / n2 * np.sin(incident_angle))
#             Rs = ((n1 * np.cos(incident_angle) - n2 * np.cos(refractive_angle)) / (n1 * np.cos(incident_angle) + n2 * np.cos(refractive_angle))) ** 2
#             Rp = ((n1 * np.cos(refractive_angle) - n2 * np.cos(incident_angle)) / (n1 * np.cos(refractive_angle) + n2 * np.cos(incident_angle))) ** 2
#             r = (Rs + Rp) / 2
#             r_array.append(r)
#         else:
#             r = 1 
#             r_array.append(r)
            
#     return r_array



# reflection_coefficient = calculate_reflection_coefficient(n1, n2, incident_angle)

# plt.plot(incident_angle, reflection_coefficient)
# plt.xlabel("incident angle (rad)")
# plt.ylabel("reflection coefficient")
# title = "Relationship between incident_angle and reflection_coefficient (n1={}, n2={})".format(n1, n2)
# plt.title(title)
# plt.grid()
# plt.show()


import scipy.integrate as integrate

def exclude_total_reflection(incident_angle, n1, n2):
    refractive_angle = np.arcsin(n1/n2 * np.sin(incident_angle))
    Rs = ((n1 * np.cos(incident_angle) - n2 * np.cos(refractive_angle)) / (n1 * np.cos(incident_angle) + n2 * np.cos(refractive_angle)))**2
    Rp = ((n1 * np.cos(refractive_angle) - n2 * np.cos(incident_angle)) / (n1 * np.cos(refractive_angle) + n2 * np.cos(incident_angle)))**2
    r = (Rs + Rp) / 2
    y = r * np.cos(incident_angle) * np.sin(incident_angle)
    return y

def include_total_reflection(incident_angle):
    y = np.cos(incident_angle) * np.sin(incident_angle)
    return y

def calculate_reflection_and_transmission(n1, n2):
    critical_angle = np.arcsin(n2 / n1)
    rd_half = integrate.quad(exclude_total_reflection, 0, critical_angle, args=(n1, n2))[0] + integrate.quad(include_total_reflection, critical_angle, np.pi/2)[0]
    rd = 2 * rd_half
    print("n1 =", n1, "n2 =", n2)
    print("Critical Angle:", critical_angle)
    print("rd:", rd)
    print("1 - rd:", 1 - rd)

n1_values = [1.4, 1.4]
n2_values = [1, 1.33]

for n1, n2 in zip(n1_values, n2_values):
    calculate_reflection_and_transmission(n1, n2)


 
