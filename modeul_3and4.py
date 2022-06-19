"""
Created on Thu Jun 16 20:26:24 2022

@author: Sowmya
"""

import numpy as np


## *** Dynamic Stall Model - Module 1 ***
def Module1(t, dt, U_0, c, X_lag_old, Y_lag_old, b1, b2, A1, A2, dalpha_qs, alpha_qs,
            dCl_dalpha, alpha_0, D_nc_old, a_s0, Kalpha, dalpha_qs_dt_new, dalpha_qs_dt_old):
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
    s = 2 * t * U_0 / c
    ds = 2 * dt * U_0 / c

    # If first time, then initialise the lag states
    if s == 0:
        X_lag_new = 0
        Y_lag_new = 0

    # Otherwise evaluate the current lag states
    else:
        X_lag_new = X_lag_old * np.exp(-b1 * ds) + dalpha_qs * A1 * np.exp(-b1 * ds / 2)
        Y_lag_new = Y_lag_old * np.exp(-b2 * ds) + dalpha_qs * A2 * np.exp(-b2 * ds / 2)

    # Evaluate the equivalent angle of attack
    alpha_eq = alpha_qs - X_lag_new - Y_lag_new

    # Evaluate the circulatory normal force coefficient
    C_n_c = dCl_dalpha * (alpha_eq - alpha_0)

    # Evaluate the deficiency function for the non-circulatory normal force coefficient
    D_nc_new = D_nc_old * np.exp(-a_s0 * dt / (Kalpha * c)) + (dalpha_qs_dt_new - dalpha_qs_dt_old) * np.exp(
        -a_s0 * dt / (2 * Kalpha * c))

    # Evaluate the non-circulatory normal force coefficient
    C_n_nc = 4 * Kalpha * c / U_0 * (dalpha_qs_dt_new - D_nc_new)

    return X_lag_new, Y_lag_new, alpha_eq, D_nc_new, C_n_c, C_n_nc


## *** Dynamic Stall Model - Module 2 ***
def Module2(D_pf_old, ds, C_n_p_old, C_n_p_new, dCl_dalpha, alpha_qs, alpha_0, D_bl_old, f_p_old, C_n_nc):
    # === Non-linear trailing edge flow separation ===

    # Inputs
    #   D_pf_old        [-?]    Pressure lag deficiency function (previous time)
    #   ds              [-]     Reduced time increment
    #   C_n_p_old       [-]     Total normal force coefficent (previous time)
    #   C_n_p_new       [-]     Total normal force coefficent (current time)
    #   dCl_dalpha      [-]     Lift slope
    #   alpha_qs        [rad]   Quasi-steady AoA
    #   alpha_0         [rad]   Zero-liftangle of attack
    #   D_bl_old        [-?]    Boundary-layer deficit function (previous time)
    #   f_p_old         [-]     Separation location (previous time)
    #   C_n_nc          [-]     Non-circulatory normal force coefficient

    # Output
    #   f_p_new     [-]     Separation location (current time)
    #   D_bl_new    [-]     Boundary-layer deficit function (current time)
    #   C_n_f       [-]     Unsteady non-linear normal force coefficient including TE separation

    # Define user's parameters
    Tp = 1.7  # [-]       (?) "Mostly related on Mach number and in small amount on the airfoil shape"
    Tf = 3.0  # [-]       (?) "Depends on the airfoil shape"
    AoA_f_1 = 7 * np.pi / 180  # [rad]     Control AoA to change the evaluation of f_sep
    AoA_f_2 = 15 * np.pi / 180  # [rad]     Control AoA to change the evaluation of f_sep
    AoA_f_3 = 21 * np.pi / 180  # [rad]     Control AoA to change the evaluation of f_sep

    # Evaluate the pressure lag deficiency function
    D_pf_new = D_pf_old * np.exp(-ds / Tp) + (C_n_p_old - C_n_p_new) * np.exp(-ds / (2 * Tp))

    # Evaluate the lagged potential flow load
    C_n_p_prime = C_n_p_new - D_pf_new

    # Evaluate the equivalent angle of attack
    alpha_f = C_n_p_prime / dCl_dalpha + alpha_0

    # Evaluate the pseudo-location of the separation point over the airfoil section (Taken from CSF)    # NOT SURE IF NEEDED; NOT IN THE ALGORITHMS DERIVED ON THE BOARD IN CLASS BUT SEEM REQUIRED
    if alpha_f <= AoA_f_1:
        f_sep = 1  # f = 1 : Separation at TE
    elif alpha_f <= AoA_f_2:
        f_set = 1 - 0.8 * ((alpha_f - AoA_f_1) / (AoA_f_2 - AoA_f_1))
    elif alpha_f <= AoA_f_3:
        f_set = 0.2 * (1 - ((alpha_f - AoA_f_2) / (AoA_f_3 - AoA_f_2)) ** 0.3)
    else:
        f_sep = 0  # f = 0 : Separation at LE

    # Evaluate the non-linear normal force coefficient in steady flow including TE separation
    C_n_st = dCl_dalpha * ((1 + np.sqrt(f_sep)) / 2) ** 2 * (alpha_qs - alpha_0)

    # Evaluate the separation location for the lagged potential flow
    f_p_new = (2 * (C_n_st / (dCl_dalpha * (alpha_f - alpha_0))) ** (1 / 2) - 1) ** 2

    # Evaluate the deficiency function for the bondary-layer development
    D_bl_new = D_bl_old * np.exp(-ds / Tf) + (f_p_new - f_p_old) * np.exp(-ds / (2 * Tf))

    # Evaluate the separation location for the boundary-layer
    f_bl = f_p_new - D_bl_new

    # Evaluate the unsteady non-linear normal force coefficient including TE separation
    C_n_f = dCl_dalpha * ((1 + np.sqrt(f_bl)) / 2) ** 2 * (
                alpha_eq - alpha_0) + C_n_nc  # NOT SURE IF USING CORRECT EQUIVALENT AOA (alpha_f OR alpha_eq) ???

    return f_p_new, D_bl_new, C_n_f


def leading_edge_flow_separation(CNf, Uinf, delta_alpha_E, s, delta_s, tau_v_previous,CN1=1.0093):
    """
    Calculate the location of the separation on the airfoil as a reduced time
    :param CNf: Unsteady normal foDrcee coefficient
    :param Uinf: Freestream velocity
    :param delta_alpha_E: Step of equivalent quasi-steady angle of attack
    :param s: reduced time
    :param delta_s: reduced time step
    :param tau_v_previous: previous location as reduced time (s-delta s)
    :param CN1: Some critical value --> i do not know what this is
    :return: tau_v_current: current location as reduced time (s)
    """
    if CNf>CN1:
        #initialise tau_v_previous for first iteration
        tau_v_current = tau_v_previous + 0.45*delta_s
        #0.45Uinf is the vortex convection speed
        #Uinf is not used here so I am not sure if we need it or if it is normalised somehow
    else:
        if (delta_alpha_E<0 and tau_v_previous>0):
            tau_v_current = tau_v_previous + 0.45*delta_s
        else:
            tau_v_current = 0 
            
    tau_v_previous = tau_v_current
    #this is not in his slides but i think its the logic? 
    #or is this a matlab thing and python can take it as an input
    
    return tau_v_current


def vortex_shedding_module(CNp,fbl,s,delta_s,tau_v_current, Cv_previous, CNv_previous,Tv = 6.0,Tvl = 5.0):
    """
    Calculate contribution of vortex to lift
    :param CNp: Circulatory and Non-circulatory force sum 
    :param fbl: separation position function
    :param s: reduced time
    :param delta_s: reduced time step
    :param Cv_previous: Forcing term
    :param CNv_previous: Coefficient of normal force from vortex at reduced time (s-delta s)
    :param Tv: Reference value
    :param Tvl: Reference value
    :return: CNv_previous: Coefficient of normal force from vortex at reduced time (s)
    """

    #initialise CNv_previous, Cv_previous for first iteration
    Cv_current = CNp*(1 - ((1+np.sqrt(fbl))/2)**2)
    
    if (tau_v_current<Tvl and tau_v_current>0):
        CNv_current = CNv_previous*np.exp(-delta_s/Tv) + (Cv_current - Cv_previous)*np.exp(-delta_s/2/Tv)
    else:
        CNv_current = CNv_previous*np.exp(-delta_s/Tv)
        
        
    CNv_previous = CNv_current    
    Cv_previous = Cv_current 
    #this is not in his slides but i think its the logic? 
    #or is this a matlab thing and python can take it as an input

    return CNv_current