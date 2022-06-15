import numpy as np
import time

from BEM_code import Blade, BladeElement, loads, c_thrust
from dynamic_inflow import pitt_peters, larsen_madsen, oye


'''
Put the dynamic stall modules here
'''


class DSAirfoil:
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
    def __init__(self, pos_r: float, chord: float, twist: float, airfoil: DSAirfoil):
        super().__init__(pos_r, chord, twist, airfoil)

        self.t = None
        self.delta_t = None

    def set_time_data(self, t, delta_t):
        self.t = t
        self.delta_t = delta_t
        self.airfoil.set_time_data(t, delta_t)


class DSBlade(Blade):
    def __init__(self, n_blades, airfoil, r_start, r_end, blade_pitch, n_elements):
        super().__init__(n_blades, airfoil, r_start, r_end, blade_pitch, n_elements, element_class=DSBladeElement)

    def set_time_data(self, t, delta_t):
        [be.set_time_data(t, delta_t) for be in self.blade_elements]


class DSTurbine:
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
        u_inf = u_inf_0 + delta_u_inf / v0 * np.sin(reduced_freq * v0 / self.blade.r * t_list)
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
            self.blade.set_time_data(t, delta_t)
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