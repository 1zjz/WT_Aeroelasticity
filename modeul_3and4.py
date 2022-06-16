"""
Created on Thu Jun 16 20:26:24 2022

@author: Sowmya
"""

import numpy as np

def leading_edge_flow_separation(CNf, CN1, Uinf, delta_alpha_E, s, delta_s, tau_v_previous):
    """
    Calculate the location of the separation on the airfoil as a reduced time
    :param CNf: Unsteady normal foDrcee coefficient
    :param CN1: Some critical value --> i do not know what this is
    :param Uinf: Freestream velocity
    :param delta_alpha_E: Step of equivalent quasi-steady angle of attack
    :param s: reduced time
    :param delta_s: reduced time step
    :param tau_v_previous: previous location as reduced time (s-delta s)
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


def vortex_shedding_module(CNp,fbl,s,delta_s,tau_v_current,Tv = 6.0,Tvl = 5.0, Cv_previous, CNv_previous):
    """
    Calculate contribution of vortex to lift
    :param CNp: Circulatory and Non-circulatory force sum 
    :param fbl: separation position function
    :param s: reduced time
    :param delta_s: reduced time step
    :param Tv: Reference value
    :param Tvl: Reference value
    :param Cv_previous: Forcing term
    :param CNv_previous: Coefficient of normal force from vortex at reduced time (s-delta s)
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