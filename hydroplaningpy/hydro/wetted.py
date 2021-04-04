import numpy as np
import math
from hydroplaningpy.general_func import *
from hydroplaningpy import constant

# To calculate wetted area with different method

## To extend into warped hull

def savitsky_wetted(hull):
    """
    This is a function based on Daniel Savitsky 1964 / 1976 paper to calculate wetted area
    :param hull: Hull Object
    :return: wetted area information
    """
    v    = hull.state.v
    beta_d = hull.beta_d
    beta   = deg_to_rad(beta_d)
    trim_d = hull.state.trim_d
    trim   = deg_to_rad(trim_d)

    cv = v/math.sqrt(constant.gravity*hull.beam)

    lambda_k     = -hull.state.z / np.sin(trim) / hull.beam

    if hull.method.wetted == 1: # Based on Savitsky 1964
        w        = np.tan(beta)/(math.pi*np.tan(trim))
        lambda_m = lambda_k - 0.5*w
        lambda_c = lambda_k - w

    elif hull.method.wetted == 2: # Based on Savitsky 1976
        w = (0.57 + beta_d/1000)*(math.tan(beta)/2/math.tan(trim) - beta_d/167)
        lambda_c = lambda_k - w

        if (lambda_c <1):
            lambda_c = (lambda_k - w) - 0.2*math.exp(-(lambda_k - w)/0.3)

        lambda_m = (lambda_k + lambda_c)/2 + 0.03
    else:
        print('Method of calculating wetted area does not exist.')
        return

    # Checking whether chines is dry or wet
    if lambda_c < 0:
        lambda_c = 0
        lambda_m = 0.5 * lambda_k
        beam_mul = lambda_k / w
    else:
        beam_mul = 1

    trim_k_d = hull.state.trim_d
    trim_e_d = trim_k_d  # FOR NOW
    beta_e_d = hull.beta_d

    wetted = {
        'cv' : cv,
        'trim_k_d' : trim_k_d,
        'trim_e_d' : trim_e_d,
        'beta_e_d' : beta_e_d,
        'lambda_k' : lambda_k,
        'lambda_m' : lambda_m,
        'lambda_c' : lambda_c,
        'beam_mul' : beam_mul
    }

    return wetted