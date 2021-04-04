import numpy as np
import math

from hydroplaningpy import constant
from hydroplaningpy.general_func import *

# Apart from the force and moment class, we still need to have many of the force related calculation here

class HydroForceMoment(object):
    """A class for hydro force and moment"""

    def __init__(self):
        self.lift_pres = 0
        self.lift_visc = 0
        self.lift_total = 0
        self.drag_pres = 0
        self.drag_visc = 0
        self.drag_total = 0
        self.moment = 0
        self.lcp = 0
        self.vcp = 0

    def __call__(self, array=None):
        """
        If there is no input, state will be printed.
        If there is any input, then the state will be the output in an array format.
        """
        if array is None:
            print('The hydro force is:')
            print('Lift by pressure  = %8.2f N' % self.lift_pres)
            print('Lift by viscosity = %8.2f N' % self.lift_visc)
            print('Total Lift        = %8.2f N' % self.lift_total)
            print('Drag by pressure  = %8.2f N' % self.drag_pres)
            print('Drag by viscosity = %8.2f N' % self.drag_visc)
            print('Total Drag        = %8.2f N' % self.drag_total)
            print('Total Moment      = %8.2f Nm'% self.moment)
            print('Longitudinal CP   = %8.2f m' % self.lcp)
            print('Vertical CP       = %8.2f m' % self.vcp)
        else:
            x = np.zeros(9)
            x[0] = self.lift_pres
            x[1] = self.lift_visc
            x[2] = self.lift_total
            x[3] = self.drag_pres
            x[4] = self.drag_visc
            x[5] = self.drag_total
            x[6] = self.moment
            x[7] = self.lcp
            x[8] = self.vcp
            return x

    def reset(self):
        self.lift_pres = 0
        self.lift_visc = 0
        self.lift_total = 0
        self.drag_pres = 0
        self.drag_visc = 0
        self.drag_total = 0
        self.moment = 0
        self.lcp = 0
        self.vcp = 0

class Wetted(object):
    """
    The non_dimensional state of a hull
    CV
    trim_k_d
    trim_e_d
    beta_e_d
    lambda_k
    lambda_m
    lambda_c
    beam_mul
    """

    def __init__(self, cv=0, trim_k_d=0, trim_e_d=0, beta_e_d=0,
                 lambda_k=0, lambda_m=0, lambda_c=0, beam_mul=0):
        self.cv = cv
        self.trim_k_d = trim_k_d  # deg
        self.trim_e_d = trim_e_d  # deg
        self.beta_e_d = beta_e_d  # deg
        self.lambda_k = lambda_k
        self.lambda_m = lambda_m
        self.lambda_c = lambda_c
        self.beam_mul = beam_mul

    def __call__(self, array=None):
        """
        If there is no input, state_nd will be printed.
        If there is any input, then the state will be the output in an array format.
        """
        if array is None:
            print('The state_nd is:')
            print('CV       = %8.2f' % self.cv)
            print('trim_k_d = %8.2f deg' % self.trim_k_d)
            print('trim_e_d = %8.2f deg' % self.trim_e_d)
            print('beta_e_d = %8.2f deg' % self.beta_e_d)
            print('lambda_k = %8.2f ' % self.lambda_k)
            print('lambda_m = %8.2f ' % self.lambda_m)
            print('lambda_c = %8.2f ' % self.lambda_c)
            print('beam_mul = %8.2f ' % self.beam_mul)
        else:
            x = np.zeros(8)
            x[0] = self.cv
            x[1] = self.trim_k_d
            x[2] = self.trim_e_d
            x[3] = self.beta_e_d
            x[4] = self.lambda_k
            x[5] = self.lambda_m
            x[6] = self.lambda_c
            x[7] = self.beam_mul
            return x

class HydroMethod(object):
    """
    List of hydrodynamic method

    HAS TO BE EXPLAINED HERE LATER
    """
    def __init__(self):
        # These are the default values
        self.wetted         = 1 # 1=Savistky 1964 ; 2=Savitsky 1976
        self.dry_chine      = 1 # 1=Savitsky 1964 correction ; 2=Milwitzky 1948 ; 3=Milwitzky CFD corrected
        self.warp_lambda    = 1 # 1 = Existing equivalent method ; 2=Wagner pi/2 wave rise factor
        self.symmetric_hull = 1 # 1=Half of symmetric hull; 2=Morabito 2011 formulation
        self.roughness = 0.0004 # 0.0004 is often used in full scale
        self.wake_profile   = 1 # 1=Flat from step; 2=Savitsky 2010; 3=Korvin ;4=CFD data; 5=TEST
        self.wake_method    = 1 # 1=Savitsky; 2=KHB; 3=Svahn Method; 4=KHB_2; 5=SVAHN_2 ;6=KHB_3
        self.hum_correction = 0 # 0=No correction  %1=Blount 1976 %2=Savitsky 1976(NOT YET IMPLEMENTED)
        self.whisker_spray  = 1 # 0=Not Calculated %1=Savitsky 2007
        self.spray_method   = 0 # 0=Ignored  %1=Savitsky Chine %2=CFD generated step
        self.spray_sf       = 1 # Safety factor for the spray height

    def __call__(self):
        print('This is to be expanded later!')

def get_hydroFM(hull):
    """
    To calculate the force and moment given hull with its wetted condition
    :param hull:
    :return: hydroFM dictionary
    """

    # Important info
    v    = hull.state.v
    beam = hull.beam

    cv = hull.wetted.cv
    trim_e_d = hull.wetted.trim_e_d
    trim_e   = deg_to_rad(trim_e_d)
    beta_e_d = hull.wetted.beta_e_d
    beta_e   = deg_to_rad(beta_e_d)
    lambda_m = hull.wetted.lambda_m

    roughness = hull.method.roughness

    # Force and Moment Calculation

    cl_0 = trim_e_d ** 1.1 * (0.012 * lambda_m ** 0.5 + 0.0055 * lambda_m ** 2.5 / cv ** 2)
    cl_beta = cl_0 - 0.0065 * beta_e_d * cl_0 ** 0.6
    lift_pres = 0.5 * constant.density_water * v ** 2 * beam ** 2 * cl_beta
    drag_pres = lift_pres * np.tan(trim_e)

    # Calculating viscous force
    cl_0_dyn = trim_e_d ** 1.1 * 0.012 * lambda_m ** 0.5
    cl_beta_dyn = cl_0_dyn - 0.0065 * beta_e_d * cl_0_dyn ** 0.6
    v_1 = v * math.sqrt(1 - cl_beta_dyn / lambda_m / np.cos(trim_e))
    reynolds_no = constant.density_water * v_1 * lambda_m * beam / 0.00122
    c_f = 0.075 / (math.log10(reynolds_no) - 2) ** 2 + roughness

    visc_force = 0.5 * 1025 * v_1 ** 2 * (lambda_m * beam ** 2) / np.cos(beta_e) * c_f
    lift_visc = -visc_force * np.sin(trim_e)
    drag_visc = visc_force * np.cos(trim_e)

    lift_total = lift_pres + lift_visc
    drag_total = drag_pres + drag_visc

    lcp = (0.75 - 1 / (5.21 * (cv / lambda_m) ** 2 + 2.39)) * lambda_m * beam
    vcp = 0.25 * np.tan(beta_e) * beam

    # Need to calculate the moment here USE TRIM_E AS REFERENCE
    # The moment is calculated w.r.t the transom of the hull
    moment = lift_total * (np.cos(trim_e) * lcp - np.sin(trim_e) * vcp) + drag_total * (np.sin(trim_e) * lcp - np.cos(trim_e) * vcp)

    hydroFM = {
        'lift_pres'  : lift_pres,
        'lift_visc'  : lift_visc,
        'lift_total' : lift_total,
        'drag_pres'  : drag_pres,
        'drag_visc'  : drag_visc,
        'drag_total' : drag_total,
        'moment'     : moment,
        'lcp'        : lcp,
        'vcp'        : vcp
    }

    return hydroFM