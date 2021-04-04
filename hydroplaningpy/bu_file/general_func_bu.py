import numpy as np
import math

def mat_rotation(theta_d):
    """Convert angle into rotation matrix. Input is in degree"""
    theta = theta_d*math.pi/180
    return np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def mat_orientation(theta_d, pos_x, pos_z):
    """To make orientation matrix. Angle is in degree"""
    mat_orient = np.identity(3)
    mat_orient[0,2] = pos_x
    mat_orient[1,2] = pos_z
    mat_orient[:2,:2] = mat_rotation(theta_d)
    return mat_orient

def get_angle_from_orient(mat_orient):
    """To get the angle in degree from an orientation matrix"""
    return 180/math.pi*np.arccos(mat_orient[0,0])

def get_x_from_orient(mat_orient):
    """To get the angle in degree from an orientation matrix"""
    return mat_orient[0,2]

def get_z_from_orient(mat_orient):
    """To get the angle in degree from an orientation matrix"""
    return mat_orient[1,2]

def zcg_to_zmain(craft, z_cg, trim_d):
    """
    A function that transforms z_cg to the z of the main hull
    """
    cg_to_main = mat_orientation(0, -craft.lcg, -craft.vcg)
    cg_orient = mat_orientation(trim_d, 0, z_cg)
    main_orient = np.matmul(cg_orient, cg_to_main)
    return get_z_from_orient(main_orient)

def zmain_to_zcg(craft, z_main, trim_d):
    """
    A function that transforms z of the main hull to the z_cg
    """
    main_to_cg  = mat_orientation(0, craft.lcg, craft.vcg)
    main_orient = mat_orientation(trim_d, 0, z_main)
    cg_orient   = np.matmul(main_orient, main_to_cg)
    return get_z_from_orient(cg_orient)

def deg_to_rad(theta_d):
    return theta_d*math.pi/180

def rad_to_deg(theta):
    return theta*180/math.pi

