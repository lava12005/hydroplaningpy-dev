import numpy as np
import math

from .general_func_bu import *

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
        hull[1] will be the main hull, referred most of the time by 'main
        """
        self.hull =[] # A list of hull
        self.hull.append(Hull('calm_water', 0, 0))
        
        # Initialize state and orientation
        self.state       = State()
        self.state_nd    = StateNd()
        self.hydro       = HydroForceMoment()
        self.main_to_cg  = mat_orientation(0,  self.lcg,  self.vcg)
        self.cg_to_main  = mat_orientation(0, -self.lcg, -self.vcg)
        self.main_orient = mat_orientation(0, 0, 0)
        self.cg_orient   = mat_orientation(0, 0, 0)
        
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

        self.cg_orient = mat_orientation(self.state.trim_d, pos_x, self.state.z)
        self.main_orient = np.matmul(self.cg_orient, self.cg_to_main)

        no_of_hull = len(self.hull) - 1
        for i in range(1, no_of_hull+1): # Remember hull_id starts with 1
            self.hull[i].hull_orient = np.matmul(self.main_orient, self.hull[i].main_to_hull)
    
    def init_hull_state(self):
        """
        Initialize the state of the craft's hull. 
        Needs to be done after initializing the orientation
        """
        
        no_of_hull = len(self.hull) - 1
        for i in range(1, no_of_hull+1): # Remember hull_id starts with 1
            if self.hull[i].wake_of == 0:
                v      = self.state.v
                trim_d = get_angle_from_orient(self.hull[i].hull_orient)
                pos_z  = get_z_from_orient(self.hull[i].hull_orient)
                self.hull[i].state.set_state([v, trim_d, pos_z, 0, 0])
    
    def get_hydro(self):
        """get hydro from all hull"""
        
        # First is to reset value to 0 first
        self.hydro.reset()
        self.init_orient()
        self.init_hull_state()
        
        no_of_hull = len(self.hull) - 1
        for i in range(1, no_of_hull+1): # Remember hull_id starts with 1
            self.hull[i].get_wetted()
            self.hull[i].get_hydro()
            
            self.hydro.lift_pres  += self.hull[i].hydro.lift_pres
            self.hydro.lift_visc  += self.hull[i].hydro.lift_visc
            self.hydro.lift_total += self.hull[i].hydro.lift_total
            self.hydro.drag_pres  += self.hull[i].hydro.drag_pres
            self.hydro.drag_visc  += self.hull[i].hydro.drag_visc
            self.hydro.drag_total += self.hull[i].hydro.drag_total
    
            x_hull = get_x_from_orient(self.hull[i].hull_orient)
            z_hull = get_z_from_orient(self.hull[i].hull_orient)
            z_main = get_z_from_orient(self.main_orient)
            self.hydro.moment     += self.hull[i].hydro.moment + \
                                     self.hull[i].hydro.lift_total * x_hull + self.hull[i].hydro.drag_total * (z_hull - z_main)
        
        # This is not really important, but a representation is given
        self.hydro.vcp = self.hull[1].hydro.vcp
        self.hydro.lcp = (self.hydro.moment - self.hydro.vcp * self.hydro.drag_total) / self.hydro.lift_total
        
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
    
class StateNd(object):
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
    
    def __init__(self, cv = 0, trim_k_d = 0, trim_e_d = 0, beta_e_d = 0,
                 lambda_k = 0, lambda_m = 0, lambda_c = 0, beam_mul = 0):
        self.cv       = cv
        self.trim_k_d = trim_k_d # deg
        self.trim_e_d = trim_e_d # deg
        self.beta_e_d = beta_e_d # deg
        self.lambda_k = lambda_k
        self.lambda_m = lambda_m
        self.lambda_c = lambda_c
        self.beam_mul = beam_mul
        
    def __call__(self, array = None):
        """
        If there is no input, state_nd will be printed.
        If there is any input, then the state will be the output in an array format.
        """
        if array is None:
            print('The state_nd is:')
            print('CV       = %8.2f'     %self.cv)
            print('trim_k_d = %8.2f deg' %self.trim_k_d)
            print('trim_e_d = %8.2f deg' %self.trim_e_d)
            print('beta_e_d = %8.2f deg' %self.beta_e_d)
            print('lambda_k = %8.2f '    %self.lambda_k)
            print('lambda_m = %8.2f '    %self.lambda_m)
            print('lambda_c = %8.2f '    %self.lambda_c)
            print('beam_mul = %8.2f '    %self.beam_mul)
        else:
            x = np.zeros(7)
            x[0] = self.cv
            x[1] = self.trim_k_d
            x[2] = self.trim_e_d
            x[3] = self.beta_e_d
            x[4] = self.lambda_k
            x[5] = self.lambda_m
            x[6] = self.lambda_c
            x[7] = self.beam_mul
            return x
        
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
        self.name = name
        self.beam = beam
        self.beta_d = beta_d
        self.warp_b = warp_b
        self.pos_x  = pos_x
        self.pos_z  = pos_z
        self.relative_trim_d = relative_trim_d
        self.id = hull_id
        self.wake_of = wake_of
        
        # Initialize state and orientation
        self.state        = State()
        self.state_nd     = StateNd()
        self.hydro        = HydroForceMoment()
        self.main_to_hull = mat_orientation(relative_trim_d, pos_x, pos_z)
        self.hull_orient  = mat_orientation(relative_trim_d, pos_x, pos_z)
        
    def __call__(self):
        print('The hull information:')
        print('Name            = ',               self.name)
        print('id              = %8.2f '         %self.id)
        print('beam            = %8.2f m'        %self.beam)
        print('beta_d          = %8.2f deg'      %self.beta_d)
        print('warp_b          = %8.2f deg/beam' %self.warp_b)
        print('relative_trim_d = %8.2f deg'      %self.relative_trim_d)
        print('wake_of hull id = %8.0f'          %self.wake_of)
        
    def get_state_nd(self):
        """
        To calculate the non-dimensional aspect of the wetted area.
        Input 
        
        NOW THERE IS NO WARP
        
        This function will populate the state_nd object of the hull
        """

        v    = self.state.v
        beta = math.pi/180*self.beta_d
        trim = math.pi/180*self.state.trim_d
        
        cv = v/math.sqrt(9.81*self.beam)
        
        lambda_k     = -self.state.z / np.sin(trim) / self.beam
        delta_lambda = np.tan(beta)/(math.pi*np.tan(trim))
        lambda_m     = lambda_k - 0.5*delta_lambda
        lambda_c     = lambda_k - delta_lambda
        
        # Checking whether chines is dry or wet
        if lambda_c < 0:
            lambda_c = 0
            lambda_m = 0.5*lambda_k
            beam_mul = lambda_k/delta_lambda
        else:
            beam_mul = 1
        
        trim_k_d = self.state.trim_d
        trim_e_d = trim_k_d # FOR NOW
        beta_e_d = self.beta_d
        
        self.state_nd.cv       = cv
        self.state_nd.trim_k_d = trim_k_d
        self.state_nd.trim_e_d = trim_e_d
        self.state_nd.beta_e_d = beta_e_d
        self.state_nd.lambda_k = lambda_k
        self.state_nd.lambda_m = lambda_m
        self.state_nd.lambda_c = lambda_c
        self.state_nd.beam_mul = beam_mul
        
    def get_hydro(self):
        """
        To calculate the force and moment on the hull
        
        This function will populate the hydro_forceMoment object of the hull
        """
        v        = self.state.v
        beam     = self.beam
        
        cv       = self.state_nd.cv
        trim_e_d = self.state_nd.trim_e_d
        trim_e   = trim_e_d*math.pi/180   # Convert to radian
        beta_e_d = self.state_nd.beta_e_d
        beta_e   = beta_e_d*math.pi/180   # Convert to radian
        lambda_m = self.state_nd.lambda_m
        
        cl_0      = trim_e_d**1.1 * (0.012*lambda_m**0.5 + 0.0055*lambda_m**2.5/cv**2)
        cl_beta   = cl_0 - 0.0065*beta_e_d*cl_0**0.6
        lift_pres = 0.5*1025*v**2*beam**2*cl_beta
        drag_pres = lift_pres*np.tan(trim_e)
        
        #Calculating viscous force
        cl_0_dyn    = trim_e_d**1.1*0.012*lambda_m**0.5
        cl_beta_dyn = cl_0_dyn - 0.0065*beta_e_d*cl_0_dyn**0.6
        v_1         = v*math.sqrt(1-cl_beta_dyn/lambda_m/np.cos(trim_e))
        reynolds_no = 1025*v_1*lambda_m*beam/0.00122
        c_f         = 0.075/(math.log10(reynolds_no) - 2)**2 + 0.0004
        visc_force  = 0.5*1025*v_1**2*(lambda_m*beam**2)/np.cos(beta_e)*c_f
        lift_visc   = -visc_force*np.sin(trim_e)
        drag_visc   =  visc_force*np.cos(trim_e)
        
        lift_total  = lift_pres + lift_visc
        drag_total  = drag_pres + drag_visc
        
        lcp = (0.75 - 1/(5.21*(cv/lambda_m)**2 + 2.39))*lambda_m*beam
        vcp = 0.25*np.tan(beta_e)*beam
        
        # Need to calculate the moment here USE TRIM_E AS REFERENCE
        # The moment is calculated w.r.t the transom of the hull
        moment = lift_total*(np.cos(trim_e)*lcp - np.sin(trim_e)*vcp) + drag_total*(np.sin(trim_e)*lcp - np.cos(trim_e)*vcp)
        
        self.hydro.lift_pres  = lift_pres
        self.hydro.lift_visc  = lift_visc
        self.hydro.lift_total = lift_pres + lift_visc
        self.hydro.drag_pres  = drag_pres
        self.hydro.drag_visc  = drag_visc
        self.hydro.drag_total = drag_pres + drag_visc
        self.hydro.moment     = moment
        self.hydro.lcp        = lcp
        self.hydro.vcp        = vcp

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
    cg_x = get_x_from_orient(craft.cg_orient)
    cg_z = get_z_from_orient(craft.cg_orient)
    force_z   = craft.hydro.lift_total - craft.mass * 9.81
    moment_cg = craft.hydro.moment - craft.hydro.lift_total * cg_x - craft.hydro.drag_total * cg_z

    # Damping to stabilize
    zeta_f = 20 * max(1, abs(craft.hydro.lift_total / craft.mass))
    zeta_m = 8 * max(1, abs(moment_cg / craft.i_yy))

    dx0 = 0  # Speed is constant
    dx1 = trim_d_dot
    dx2 = z_dot
    dx3 = -zeta_m * (trim_d_dot * math.pi / 180) + moment_cg / craft.i_yy
    dx4 = -zeta_f * z_dot + force_z / craft.mass

    return np.array((dx0, dx1, dx2, dx3*180/math.pi, dx4))