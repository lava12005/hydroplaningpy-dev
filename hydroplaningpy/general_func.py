import numpy as np
import math

__all__ =['deg_to_rad', 'rad_to_deg']

def _mat_rotation(theta_d):
    """Convert angle into rotation matrix. Input is in degree"""
    theta = theta_d*math.pi/180
    return np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def deg_to_rad(theta_d):
    return theta_d*math.pi/180

def rad_to_deg(theta):
    return theta*180/math.pi