import numpy as np
import math

from .general_func import *
from .orientation import *

from . import hydro
from . import constant

class Craft(object):
    """
    Creates a craft object 
    
    Parameters
    -------------------------------------
    name -> name of the craft
    LCG -> Longitudinal CG of the craft
    VCG -> Vertical CG of the craft
    
    List of Variables that can be defined
    -------------------------------------
    mass
    i_yy
    """
    
    # The following value is to be initialized later on
    mass = 1000 # kg
    i_yy = 1000 * (0.25 * 10) ** 2 # kg*m^2
    
    def __init__(self, name='craft', lcg=0, vcg=0):
        self.name = name
        self.lcg  = lcg
        self.vcg  = vcg

        """
        The list of hull of the craft.
        hull[0] is reserved for calm water
        hull[1] will be the main hull, referred most of the time by 'main'
        """
        self.hull =[] # A list of hull
        self.hull.append(Hull('calm_water', 0, 0))
        
        # Initialize state and orientation
        self.state       = State()
        self.wetted      = hydro.Wetted()
        self.hydroFM     = hydro.HydroForceMoment()
        self.main_to_cg  = Orientation(0,  self.lcg,  self.vcg)
        self.cg_to_main  = Orientation(0, -self.lcg, -self.vcg)
        self.main_orient = Orientation(0, 0, 0)
        self.cg_orient   = Orientation(0, 0, 0)
        
    # Printing out the relevant information
    def __call__(self):
        print('The craft information:')
        print('Name = ',             self.name)
        print('LCG  = %8.2f m'      %self.lcg)
        print('VCG  = %8.2f m'      %self.vcg)
        print('Mass = %8.2f kg'     %self.mass)
        print('I_yy = %8.2f kg*m^2' % self.i_yy)
        
    # Adding hull
    def add_hull(self, name, beam, beta_d, warp_b = 0, pos_x = 0, pos_z = 0, relative_trim_d = 0, wake_of= 0):
        available_id = len(self.hull)
        
        # Make sure that each hull is a wake of other hull
        if wake_of > available_id:
            return 'The hull %d is not in the wake of a valid hull' %available_id
        else: 
            self.hull.append(Hull(name, beam, beta_d, warp_b, pos_x, pos_z, relative_trim_d, wake_of, hull_id = available_id))
    
    def init_orient(self):
        """Initialize the orientation of the craft and its hull"""
        
        # Need to get the main_orient first
        trim = math.pi/180*self.state.trim_d
        pos_x = self.lcg*np.cos(trim) - self.vcg*np.sin(trim)

        self.cg_orient = Orientation(self.state.trim_d, pos_x, self.state.z)
        self.main_orient = orient_mul(self.cg_orient, self.cg_to_main)

        no_of_hull = len(self.hull) - 1
        for i in range(1, no_of_hull+1): # Remember hull_id starts with 1
            self.hull[i].hull_orient = orient_mul(self.main_orient, self.hull[i].main_to_hull)
    
    def init_hull_state(self):
        """
        Initialize the state of the craft's hull. 
        Needs to be done after initializing the orientation
        """
        
        no_of_hull = len(self.hull) - 1
        for i in range(1, no_of_hull+1): # Remember hull_id starts with 1
            if self.hull[i].wake_of == 0:
                v      = self.state.v
                trim_d = self.hull[i].hull_orient.theta_d
                pos_z  = self.hull[i].hull_orient.pos_z
                self.hull[i].state.set_state([v, trim_d, pos_z, 0, 0])
    
    def get_hydro(self):
        """get hydro from all hull"""
        
        # First is to reset value to 0 first
        self.hydroFM.reset()
        self.init_orient()
        self.init_hull_state()
        
        no_of_hull = len(self.hull) - 1
        for i in range(1, no_of_hull+1): # Remember hull_id starts with 1
            self.hull[i].get_wetted()
            self.hull[i].get_hydro()
            
            self.hydroFM.lift_pres  += self.hull[i].hydroFM.lift_pres
            self.hydroFM.lift_visc  += self.hull[i].hydroFM.lift_visc
            self.hydroFM.lift_total += self.hull[i].hydroFM.lift_total
            self.hydroFM.drag_pres  += self.hull[i].hydroFM.drag_pres
            self.hydroFM.drag_visc  += self.hull[i].hydroFM.drag_visc
            self.hydroFM.drag_total += self.hull[i].hydroFM.drag_total
    
            x_hull = self.hull[i].hull_orient.pos_x
            z_hull = self.hull[i].hull_orient.pos_z
            z_main = self.main_orient.pos_z
            self.hydroFM.moment   += self.hull[i].hydroFM.moment + \
                                   self.hull[i].hydroFM.lift_total * x_hull + self.hull[i].hydroFM.drag_total * (z_hull - z_main)
        
        # This is not really important, but a representation is given
        self.hydroFM.vcp = self.hull[1].hydroFM.vcp
        self.hydroFM.lcp = (self.hydroFM.moment - self.hydroFM.vcp * self.hydroFM.drag_total) / self.hydroFM.lift_total
        
class State(object):
    """
    The hydrodynamic state is given by 
    v             (m/s)
    trim_d        (deg)
    z             (m)
    trim_d_dot    (deg/s)
    z_dot         (m/s)
    """
    
    def __init__(self):
        self.v          = 10   # m/s
        self.trim_d     = 4    # deg
        self.z          = 0    # m
        self.trim_d_dot = 0    # deg/s
        self.z_dot      = 0    # m/s
        
    def __call__(self, array = None):
        """
        If there is no input, state will be printed.
        If there is any input, then the state will be the output in an array format.
        """
        if array is None:
            print('The state is:')
            print('V         = %6.2f m/s'   %self.v)
            print('trim      = %6.2f deg'   %self.trim_d)
            print('z         = %6.2f m'     %self.z)
            print('trim_dot  = %6.2f deg/s' %self.trim_d_dot)
            print('z_dot     = %6.2f m/s'   %self.z_dot)
        else:
            x = np.zeros(5)
            x[0] = self.v
            x[1] = self.trim_d
            x[2] = self.z
            x[3] = self.trim_d_dot
            x[4] = self.z_dot
            return x
        
    def set_state(self, x):
        """
        Convert an array data to the state
        """
        self.v          = x[0]
        self.trim_d     = x[1]
        self.z          = x[2]
        self.trim_d_dot = x[3]
        self.z_dot      = x[4]
        
class Hull(object):
    """    
    Creates an hull object 
    
    Parameters
    -------------------------------------
    name -> name of the hull
    beam -> width of the beam (m)
    beta_d -> deadrise angle (deg) 
    warp_b -> deadrise warping (deg/beam)
    pos_x -> relative position in x-axis (horizontal) positive is in front
    pos_z -> relative position in z-axis (vertical) positive is upwards
    relative_trim_d -> the relative trim of this hull w.r.t the main hull (deg)
    hull_id -> id of the hull, 1 corresponds to the main hull
    wake_of_hull -> the hull right in front of the existing hull
    
    """
    
    def __init__(self, name, beam, beta_d, warp_b = 0, pos_x = 0, pos_z = 0, relative_trim_d = 0, wake_of = 0, hull_id = 0):
        self.name            = name
        self.beam            = beam
        self.beta_d          = beta_d
        self.warp_b          = warp_b
        self.pos_x           = pos_x
        self.pos_z           = pos_z
        self.relative_trim_d = relative_trim_d
        self.id              = hull_id
        self.wake_of         = wake_of
        
        # Initialize state and orientation
        self.state        = State()
        self.wetted       = hydro.Wetted()
        self.hydroFM      = hydro.HydroForceMoment()
        self.main_to_hull = Orientation(relative_trim_d, pos_x, pos_z)
        self.hull_orient  = Orientation(relative_trim_d, pos_x, pos_z)

        # Initialize Hydro Method
        self.method       = hydro.HydroMethod()

    def __call__(self):
        print('The hull information:')
        print('Name            = ',               self.name)
        print('id              = %8.2f '         %self.id)
        print('beam            = %8.2f m'        %self.beam)
        print('beta_d          = %8.2f deg'      %self.beta_d)
        print('warp_b          = %8.2f deg/beam' %self.warp_b)
        print('relative_trim_d = %8.2f deg'      %self.relative_trim_d)
        print('wake_of hull id = %8.0f'          %self.wake_of)
        
    def get_wetted(self):
        """
        To calculate the non-dimensional aspect of the wetted area.
        Input 
        
        NOW THERE IS NO WARP FOR NOW
        """

        _wetted = hydro.wetted.savitsky_wetted(self)

        self.wetted.cv       = _wetted['cv']
        self.wetted.trim_k_d = _wetted['trim_k_d']
        self.wetted.trim_e_d = _wetted['trim_e_d']
        self.wetted.beta_e_d = _wetted['beta_e_d']
        self.wetted.lambda_k = _wetted['lambda_k']
        self.wetted.lambda_m = _wetted['lambda_m']
        self.wetted.lambda_c = _wetted['lambda_c']
        self.wetted.beam_mul = _wetted['beam_mul']
        
    def get_hydro(self):
        """
        To calculate the force and moment on the hull
        """

        _hydroFM = hydro.get_hydroFM(self)

        self.hydroFM.lift_pres  = _hydroFM['lift_pres']
        self.hydroFM.lift_visc  = _hydroFM['lift_visc']
        self.hydroFM.lift_total = _hydroFM['lift_total']
        self.hydroFM.drag_pres  = _hydroFM['drag_pres']
        self.hydroFM.drag_visc  = _hydroFM['drag_visc']
        self.hydroFM.drag_total = _hydroFM['drag_total']
        self.hydroFM.moment     = _hydroFM['moment']
        self.hydroFM.lcp        = _hydroFM['lcp']
        self.hydroFM.vcp        = _hydroFM['vcp']

def craft_eom(state, craft):
    """
    The EOM of a craft with multi hull

    state here is an array, not an object (MAYBE IN THE FUTURE CAN BE BOTH?)
    """

    # ===========
    # The state
    # ===========
    v          = state[0]
    trim_d     = state[1]
    z          = state[2]
    trim_d_dot = state[3]
    z_dot      = state[4]

    # ==============
    # Final values
    # ==============
    # NOT NOW, LATER MAYBE

    # TO CALCULATE THE HYDRO FORCE
    craft.state.set_state(state)
    craft.get_hydro()

    # ===================
    # Equation of Motion
    # ===================
    cg_x = craft.cg_orient.pos_x
    cg_z = craft.cg_orient.pos_z
    force_z   = craft.hydroFM.lift_total - craft.mass * 9.81
    moment_cg = craft.hydroFM.moment - craft.hydroFM.lift_total * cg_x - craft.hydroFM.drag_total * cg_z

    # Artificial Damping to stabilize
    zeta_f = 20 * max(1, abs(craft.hydroFM.lift_total / craft.mass))
    zeta_m = 8 * max(1, abs(moment_cg / craft.i_yy))

    dx0 = 0  # Speed is constant
    dx1 = trim_d_dot
    dx2 = z_dot
    dx3 = -zeta_m * (trim_d_dot * math.pi / 180) + moment_cg / craft.i_yy
    dx4 = -zeta_f * z_dot + force_z / craft.mass

    return np.array((dx0, dx1, dx2, dx3*180/math.pi, dx4))