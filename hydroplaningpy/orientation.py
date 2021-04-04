"""
To handle all related to the orientation matrix.
This orientation matrix is basically a homogeneous transformation for 2D case
"""
import numpy as np
from .general_func import *

__all__ =['Orientation', 'orient_mul', 'zmain_to_zcg', 'zcg_to_zmain']

class Orientation(object):
    """
    Orientation here is a 2D homogeneous transformation matrix.
    It is useful to denote the attitude of a hull as well as transformation between them.
    """
    def __init__(self, theta_d, pos_x, pos_z):
        self.theta_d = theta_d
        self.theta   = deg_to_rad(theta_d)
        self.pos_x   = pos_x
        self.pos_z   = pos_z

    def __call__(self, array = None):
        """
        If there is no input then the state will be printed.
        If there is any input, then the state will be the output in an array format.
        """
        if array is None:
            print('The orientation is:')
            print('theta_d    = %6.2f deg' %self.theta_d)
            print('theta      = %6.2f rad' %self.theta)
            print('x position = %6.2f m'   %self.pos_x)
            print('z position = %6.2f m'   %self.pos_z)
        else:
            self._orient = np.identity(3)
            self._orient[0,2] = self.pos_x
            self._orient[1,2] = self.pos_z
            self._orient[:2,:2] = np.array([[np.cos(self.theta), -np.sin(self.theta)],
                                            [np.sin(self.theta), np.cos(self.theta)]])
            return self._orient

def orient_mul(orient_1, orient_2):
    """
    Basically a wrapper for matrix multiplication for orientation object
    orient_1: Orient2D, object
            The first orientation

    orient_2: Orient2D, object
            The second orientation
    :return:
    """
    theta_d = orient_1.theta_d + orient_2.theta_d
    orient_mat = np.matmul(orient_1(0), orient_2(0)) # Returns a 3x3 array
    pos_x = orient_mat[0,2]
    pos_z = orient_mat[1,2]
    return Orientation(theta_d, pos_x, pos_z)

def zcg_to_zmain(craft, z_cg, trim_d):
    """
    A function that transforms z_cg to the z of the main hull
    """
    cg_to_main = Orientation(0, -craft.lcg, -craft.vcg)
    cg_orient  = Orientation(trim_d, 0, z_cg)
    main_orient = orient_mul(cg_orient, cg_to_main)
    return main_orient.pos_z

def zmain_to_zcg(craft, z_main, trim_d):
    """
    A function that transforms z of the main hull to the z_cg
    """
    main_to_cg  = Orientation(0, craft.lcg, craft.vcg)
    main_orient = Orientation(trim_d, 0, z_main)
    cg_orient   = orient_mul(main_orient, main_to_cg)
    return cg_orient.pos_z