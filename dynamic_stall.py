import numpy as np
import time

from BEM_code import Blade, BladeElement, loads, c_thrust
from dynamic_inflow import pitt_peters, larsen_madsen, oye


'''
Put the dynamic stall modules here
'''

## *** Nomenclature ***
#   AoA     Angle of attack
#   TE      Trailing edge
#   CSF     Carlos Simoa Ferreira


## *** Dynamic Stall Model - Module 1 ***
def Module1(t,dt,U_0,c,X_lag_old,Y_lag_old,b1,b2,A1,A2,dalpha_qs,alpha_qs,\
            dCl_dalpha,alpha_0,D_nc_old,a_s0,Kalpha,dalpha_qs_dt_new,dalpha_qs_dt_old):
    # === Unsteady attached flow ===
    
    # Inputs
    #   t                   [s]     Time
    #   dt                  [s]     Step in time
    #   U_0                 [m/s]   Freestream velocity
    #   c                   [m]     Airfoil chord
    #   X_lag_old           [rad?]  First lag state in AoA; comes from first term in Wagner expression (previous time)
    #   Y_lag_old           [rad?]  Second lag state in AoA; comes from second term in Wagner expression (previous time)
    #   b1                  [-]     Wagner expression coefficient; Given to be b1 = 0.045
    #   b2                  [-]     Wagner expression coefficient; Given to be b2 = 0.3
    #   A1                  [-]     Wagner expression coefficient; Given to be A1 = 0.165
    #   A2                  [-]     Wagner expression coefficient; Given to be A2 = 0.335
    #   dalpha_qs           [-(?)]  Increment in quasi-steady Aoa
    #   dCl_dalpha          [-]     Lift slope
    #   alpha_0             [rad]   Zero-liftangle of attack
    #   D_nc_old            [-]     Non-circulatory deficiency function (previous time)
    #   a_s0                [m/s]   Speed of sound
    #   Kalpha              [-]     Constant (?) ; Given to be 0.75
    #   dalpha_qs_dt_new    [-]     Time change in quasi-steady AoA (current time)
    #   dalpha_qs_dt_old    [-]     Time change in quasi-steady AoA (previous time)
    
    # Output
    #   X_lag_new   [rad?]  First lag state in AoA (current time)
    #   Y_lag_new   [rad?]  First lag state in AoA (current time)
    #   alpha_eq    [rad]   Equivalent AoA
    #   D_nc_new    [-]     Non-circulatory deficiency function (current time)
    #   C_n_c       [-]     Circulatory normal force coefficient
    #   C_n_nc      [-]     Non-circulatory normal force coefficient
    
    # Evaluate the reduced time and step in reduced time
    s = 2*t*U_0/c
    ds = 2*dt*U_0/c
    
    # If first time, then initialise the lag states
    if s == 0:
        X_lag_new = 0
        Y_lag_new = 0
        
    # Otherwise evaluate the current lag states
    else:
        X_lag_new = X_lag_old*np.exp(-b1*ds) + dalpha_qs*A1*np.exp(-b1*ds/2)
        Y_lag_new = Y_lag_old*np.exp(-b2*ds) + dalpha_qs*A2*np.exp(-b2*ds/2)
        
    # Evaluate the equivalent angle of attack
    alpha_eq = alpha_qs - X_lag_new - Y_lag_new
    
    # Evaluate the circulatory normal force coefficient
    C_n_c = dCl_dalpha * (alpha_eq - alpha_0)
    
    # Evaluate the deficiency function for the non-circulatory normal force coefficient
    D_nc_new = D_nc_old * np.exp(-a_s0*dt/(Kalpha*c)) + (dalpha_qs_dt_new - dalpha_qs_dt_old) * np.exp(-a_s0*dt/(2*Kalpha*c))
    
    # Evaluate the non-circulatory normal force coefficient
    C_n_nc = 4*Kalpha*c/U_0 * (dalpha_qs_dt_new - D_nc_new)
    
    return X_lag_new, Y_lag_new, alpha_eq, D_nc_new, C_n_c, C_n_nc

## *** Dynamic Stall Model - Module 2 ***
def Module2(D_pf_old,ds,C_n_p_old,C_n_p_new,dCl_dalpha,alpha_0,D_bl_old,f_p_old,C_n_nc):
    # === Non-linear trailing edge flow separation ===
    
    # Inputs
    #   D_pf_old        [-?]    Pressure lag deficiency function (previous time)
    #   ds              [-]     Reduced time increment
    #   C_n_p_old       [-]     Total normal force coefficent (previous time)
    #   C_n_p_new       [-]     Total normal force coefficent (current time)
    #   dCl_dalpha      [-]     Lift slope
    #   alpha_0         [rad]   Zero-liftangle of attack
    #   D_bl_old        [-?]    Boundary-layer deficit function (previous time)
    #   f_p_old         [-]     Separation location (previous time)
    #   C_n_nc          [-]     Non-circulatory normal force coefficient
    
    # Output
    #   f_p_new     [-]     Separation location (current time)
    #   D_bl_new    [-]     Boundary-layer deficit function (current time)
    #   C_n_f       [-]     Unsteady non-linear normal force coefficient including TE separation
    
    
    # Define user's parameters
    Tp = 1.7                # [-]       (?) "Mostly related on Mach number and in small amount on the airfoil shape"
    Tf = 3.0                # [-]       (?) "Depends on the airfoil shape"
    AoA_f_1 = 7*np.pi/180   # [rad]     Control AoA to change the evaluation of f_sep
    AoA_f_2 = 15*np.pi/180  # [rad]     Control AoA to change the evaluation of f_sep
    AoA_f_3 = 21*np.pi/180  # [rad]     Control AoA to change the evaluation of f_sep
    
    # Evaluate the pressure lag deficiency function
    D_pf_new = D_pf_old * np.exp(-ds/Tp) + (C_n_p_old - C_n_p_new) * np.exp(-ds/(2*Tp))
    
    # Evaluate the lagged potential flow load
    C_n_p_prime = C_n_p_new - D_pf_new
    
    # Evaluate the equivalent angle of attack
    alpha_f = C_n_p_prime/dCl_dalpha + alpha_0
    
    # Evaluate the pseudo-location of the separation point over the airfoil section (Taken from CSF)    # NOT SURE IF NEEDED; NOT IN THE ALGORITHMS DERIVED ON THE BOARD IN CLASS BUT SEEM REQUIRED
    if alpha_f <= AoA_f_1:
        f_sep = 1   # f = 1 : Separation at TE
    elif alpha_f <= AoA_f_2:
        f_set = 1 - 0.8 * ((alpha_f-AoA_f_1)/(AoA_f_2 - AoA_f_1))           
    elif alpha_f <= AoA_f_3:
        f_set = 0.2 * (1 - ((alpha_f-AoA_f_2)/(AoA_f_3 - AoA_f_2))**0.3)
    else:
        f_sep = 0   # f = 0 : Separation at LE
    
    # Evaluate the non-linear normal force coefficient in steady flow including TE separation
    C_n_st = dCl_dalpha * ((1 + np.sqrt(f_sep))/2)**2 * (alpha_f - alpha_0)             # NOT SURE IF USING CORRECT EQUIVALENT AOA (alpha_f OR alpha_eq) ???
    
    # Evaluate the separation location for the lagged potential flow
    f_p_new = (2 * (C_n_st/(dCl_dalpha * (alpha_f - alpha_0)))**(1/2) - 1)**2           # NOT SURE IF USING CORRECT EQUIVALENT AOA (alpha_f OR alpha_eq) ???
    
    # Evaluate the deficiency function for the bondary-layer development
    D_bl_new = D_bl_old * np.exp(-ds/Tf) + (f_p_new - f_p_old) * np.exp(-ds/(2*Tf))
    
    # Evaluate the separation location for the boundary-layer
    f_bl = f_p_new - D_bl_new
    
    # Evaluate the unsteady non-linear normal force coefficient including TE separation
    C_n_f = dCl_dalpha * ((1 + np.sqrt(f_bl))/2)**2 * (alpha_eq - alpha_0) + C_n_nc     # NOT SURE IF USING CORRECT EQUIVALENT AOA (alpha_f OR alpha_eq) ???
    
    return f_p_new, D_bl_new, C_n_f


class DSAirfoil:
    """
    Replacement for the DU95W150 airfoil class for the dynamic stall model
    """
    def __init__(self):
        self.t = None
        self.delta_t = None

    def set_time_data(self, t, delta_t):
        self.t = t
        self.delta_t = delta_t

    def cl(self, alpha):
        """
        I will program the dynamic stall loop here
        :param alpha: the quasi-steady angle of attack
        :return: The dynamic stall lift coefficient
        """
        _ = self.t + alpha
        return None

    def cd(self, alpha):
        _ = self.t + alpha
        return 0


class DSBladeElement(BladeElement):
    """
    Class to represent blade elements in the dynamic stall model. This one has all the functions from BladeElement with
    the addition that it keeps track of time for the dynamic stall model.
    """
    def __init__(self, pos_r: float, chord: float, twist: float, airfoil: DSAirfoil):
        super().__init__(pos_r, chord, twist, airfoil)

        self.t = None
        self.delta_t = None

    def __repr__(self): # print blade element
        return f"<DS Blade Element at r={self.r}, c={self.c}, beta={self.twist}>"

    def set_time_data(self, t, delta_t):
        self.t = t
        self.delta_t = delta_t
        self.airfoil.set_time_data(t, delta_t)


class DSBlade(Blade):
    """
    Class to represent the blade in the dynamic stall model. This one has all the functions from Blade with
    the addition that it keeps track of time for the dynamic stall model.
    """
    def __init__(self, n_blades, airfoil, r_start, r_end, blade_pitch, n_elements):
        super().__init__(n_blades, airfoil, r_start, r_end, blade_pitch, n_elements, element_class=DSBladeElement)

    def set_time_data(self, t, delta_t):
        [be.set_time_data(t, delta_t) for be in self.blade_elements]


class DSTurbine:
    """
    New turbine class for the dynamic stall model. Basically the same as the Turbine class for dynamic inflow, but
    restricted to sinusoidal velocity signals instead of the whole array of signals we needed in the previous assignment
    """
    def __init__(self, n_annuli):
        self.blade = DSBlade(3, DSAirfoil, .2 * 50, 50, -2, n_annuli)

    def u_inf_func(self, u_inf_0, delta_u_inf, reduced_freq, v0, tsr, model='pp'):
        """
        Determine and plot the time evolution of the turbine properties given a step in thrust coefficient
        :param u_inf_0: Mean inflow velocity
        :param delta_u_inf: Amplitude of the inflow velocity variation
        :param reduced_freq: Reduced frequency of the dynamic inflow
        :param v0: The incoming velocity
        :param tsr: The turbine tip-speed ratio
        :param model: Selection of the dynamic inflow model (pp: Pitt-Peters, lm: Larsen-Madsen, oye: Oye)
        :return: None
        """
        if model not in ('pp', 'lm', 'oye'):
            raise ValueError("Unknown model, please enter one of the following: 'pp', 'lm', 'oye'.")

        # Initialise a timer list to check compute time
        timer = [time.time(), ]

        # Initialise the time parameters: time step, start and final time
        delta_t = .04 * self.blade.r / v0
        t_0 = -.2 * self.blade.r / v0
        t_final = 4 * np.pi / reduced_freq * self.blade.r / v0
        t_list = np.round(np.arange(t_0, t_final + delta_t, delta_t), 9)

        # Extract the radial positions of the blade elements and the radial length of each
        r_list = self.blade.r_list[1:-1]
        dr = r_list[1] - r_list[0]

        # Generate the sinusoid and determine the corresponding pitch time series
        u_inf = u_inf_0 + delta_u_inf / v0 * np.cos(reduced_freq * v0 / self.blade.r * t_list)
        u_inf *= v0

        # Initialise the output value arrays: induction, AoA, thrust coefficient.
        # The shape is (time series x spanwise distribution).
        a = np.empty((t_list.size, r_list.size))
        alpha = np.empty((t_list.size, r_list.size))
        phi = np.empty((t_list.size, r_list.size))
        ctr = np.empty((t_list.size, r_list.size))
        cqr = np.empty((t_list.size, r_list.size))

        # Initialise the intermediate induced velocity array.
        # The shape is (time series x spanwise distribution).
        v_int = np.empty((t_list.size, r_list.size))

        # Initialise the quasi-steady value arrays: induction, AoA, thrust coefficient.
        # The shape is (time series x spanwise distribution).
        a_qs = np.empty((t_list.size, r_list.size))
        alpha_qs = np.empty((t_list.size, r_list.size))
        phi_qs = np.empty((t_list.size, r_list.size))
        ctr_qs = np.empty((t_list.size, r_list.size))
        cqr_qs = np.empty((t_list.size, r_list.size))

        # Loop over time, with index 'n' and time 't'
        for n, t in enumerate(t_list):
            # Just some stuff to print the status every once in a while and to monitor compute time.
            if not n:
                print(f't = {t}s\t\t(t_final = {round(t_final, 1)}s)\t(Preparation computed in {round(time.time() - timer[-1], 3)} s)')
                timer.append(time.time())
            elif t % 5 == 0:
                print(f't = {t}s\t\t(t_final = {round(t_final, 1)}s)\t(Last 5 seconds computed in {round(time.time() - timer[-1], 3)} s)')
                timer.append(time.time())

            # Some stuff for efficiency
            # In case the pitch does not change (I used this check because there sometimes is a machine error here)
            # Also ensure the first time step gets correct values by ignoring it in this check with the 2nd condition
            if abs(u_inf[n] - u_inf[n - 1]) < 1e-15 and n != 0:
                # Just reuse the thrust coefficient distribution from the previous time step
                ctr_qs[n, :] = ctr_qs[n - 1, :]
                cqr_qs[n, :] = cqr_qs[n - 1, :]

            # In case the pitch has changed since last time step
            else:
                # Run the BEM code for this pitch angle
                self.blade.reset()
                self.blade.determine_cp_ct(u_inf[n], tsr * v0 / u_inf[n], 0)
                # Get the new qs thrust coefficient distribution
                ctr_qs[n, :] = c_thrust(self.blade.p_n_list[1:-1], u_inf[n], r_list, self.blade.b, dr)
                cqr_qs[n, :] = c_thrust(self.blade.p_t_list[1:-1], u_inf[n], r_list, self.blade.b, dr) * r_list / self.blade.r

            # Set the time parameters inside the blade elements at each time step for the dynamics stall model
            # Do this here because the self.blade.reset() resets the timekeeping and we don't want that
            self.blade.set_time_data(t, delta_t)

            # Loop over the blade elements
            for i, be in enumerate(self.blade.blade_elements[1:-1]):
                # Set a tuple with parameters that the loads() function will need inside the different models
                params = (be.r, be.twist, be.c, self.blade.r, 0, be.airfoil, u_inf[n], tsr * v0 / self.blade.r, 0, 0)

                a_qs[n, i] = be.a
                alpha_qs[n, i] = be.alpha
                phi_qs[n, i] = be.phi

                # At the first time step, just initialise the output and intermediate value arrays
                if n == 0:
                    a[0, i] = be.a
                    pn, pt, alpha[0, i], phi[0, i] = loads(a[0, i], *params)
                    ctr[0, i] = c_thrust(pn, u_inf[n], be.r, self.blade.b, dr)
                    cqr[0, i] = c_thrust(pt, u_inf[n], be.r, self.blade.b, dr) * be.r / self.blade.r
                    v_int[0, i] = -a[0, i] * v0

                # If the model is Pitt-Peters
                elif model == 'pp':
                    # Propagate the AoA and induction factor of this blade element with pitt_peters()
                    a[n, i] = pitt_peters(ctr_qs[n, i], a[n - 1, i], delta_t, params, dr, self.blade.b)

                # In case of Larsen-Madsen
                elif model == 'lm':
                    # Propagate the AoA and induction factor of this blade element with larsen_madsen().
                    # be.a is the quasi-steady induction factor that L-M requires
                    a[n, i] = larsen_madsen(be.a, a[n - 1, i], delta_t, params)

                elif model == 'oye':
                    # Propagate the AoA and induction factor of this blade element with oye().
                    a[n, i], v_int[n, i] = oye(a_qs[n, i], a_qs[n - 1, i], a[n - 1, i], delta_t, params,
                                               v_int[n - 1, i])

                pn, pt, alpha[n, i], phi[n, i] = loads(a[n, i], *params)
                ctr[n, i] = c_thrust(pn, u_inf[n], be.r, self.blade.b, dr)
                cqr[n, i] = c_thrust(pt, u_inf[n], be.r, self.blade.b, dr) * be.r / self.blade.r

        # Just a happy print statement bacause the code is done running :D
        print(f'Done! (Entire time series computed in {round(timer[-1] - timer[0], 3)} s)')

        # Return the outputs for later plotting
        return r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs