import numpy as np
import time

from BEM_code import DU95W150, Blade, interpolate, loads, c_thrust
from read_write import read_from_file, write_to_file


ct_steps = ((.5, .4), (.9, -.4), (.2, .9), (1.1, -.7),)
ct_sins = ((.5, .5), (.9, .3), (.2, .7),)

u_inf_steps = ((1., .5), (1., -.3), (1., .2), (1., -.1))
u_inf_sins = ((1., .5), (.7, .3), (1.2, .5))


class Turbine:
    def __init__(self, n_annuli):
        self.blade = Blade(3, DU95W150, .2 * 50, 50, -2, n_annuli)

    def ct_pitch(self, v0, tsr):
        """
        Determine the thrust coefficient vs. thrust curve
        :param v0: The incoming velocity
        :param tsr: The turbine tip-speed ratio
        :return: None
        """
        pitch = np.round(np.arange(-10, 15 + .01, .01), 2)
        ct = np.empty(pitch.shape)
        for i, theta in enumerate(pitch):
            print(theta)
            self.blade.determine_cp_ct(v0, tsr, theta)
            ct[i] = self.blade.c_thrust

        out_array = np.array([pitch, ct])
        write_to_file(out_array, f'ct_pitch_{v0}_{tsr}.csv')

    @staticmethod
    def pitch(ct_in, v0, tsr):
        """
        Determine the pitch required to achieve a given thrust coefficient
        MAKE SURE TO HAVE THE THRUST-PITCH CURVE DETERMINED USING ct_pitch()
        :param ct_in: Input thrust coefficient
        :param v0: The incoming velocity
        :param tsr: The turbine tip-speed ratio
        :return: The pitch angle in degrees
        """
        pitch, ct = read_from_file(f'ct_pitch_{v0}_{tsr}.csv')

        ct1, ct2 = ct[ct > ct_in][-1], ct[ct <= ct_in][0]
        pitch1, pitch2 = pitch[ct > ct_in][-1], pitch[ct <= ct_in][0]

        return interpolate(pitch1, pitch2, ct1, ct2, ct_in)

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
        t_final = 10 * self.blade.r / v0 if reduced_freq is None else 4 * np.pi / reduced_freq * self.blade.r / v0
        t_list = np.round(np.arange(t_0, t_final + delta_t, delta_t), 9)

        # Extract the radial positions of the blade elements and the radial length of each
        r_list = self.blade.r_list[1:-1]
        dr = r_list[1] - r_list[0]

        # Set the ct and pitch time series depending on whether the case is a step function or a sinusoidal function
        if reduced_freq is None:
            # In case of a step function, start with an empty array for both ct and pitch
            u_inf = np.empty((t_list.size,))

            # Fill all the values before t=0 with the initial ct and the corresponding pitch
            u_inf[t_list <= 0] = u_inf_0
            # Fill all the values after t=0 with the final ct and the corresponding pitch
            u_inf[t_list > 0] = u_inf_0 + delta_u_inf / v0
        else:
            # In case of a sinusoidal function, generate the sinusoid and determine the corresponding pitch time series
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

    def ct_func(self, ct0, delta_ct, reduced_freq, v0, tsr, model='pp'):
        """
        Determine and plot the time evolution of the turbine properties given a step in thrust coefficient
        :param ct0: Mean thrust coefficient
        :param delta_ct: Amplitude of the thrust coefficient variation
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
        t_final = 10 * self.blade.r / v0 if reduced_freq is None else 4 * np.pi / reduced_freq * self.blade.r / v0
        t_list = np.round(np.arange(t_0, t_final + delta_t, delta_t), 9)

        # Extract the radial positions of the blade elements and the radial length of each
        r_list = self.blade.r_list[1:-1]
        dr = r_list[1] - r_list[0]

        # Set the ct and pitch time series depending on whether the case is a step function or a sinusoidal function
        if reduced_freq is None:
            # In case of a step function, start with an empty array for both ct and pitch
            ct, pitch = np.empty((2, t_list.size))

            # Fill all the values before t=0 with the initial ct and the corresponding pitch
            ct[t_list <= 0] = ct0
            pitch[t_list <= 0] = self.pitch(ct0, v0, tsr)
            # Fill all the values after t=0 with the final ct and the corresponding pitch
            ct[t_list > 0] = ct0 + delta_ct
            pitch[t_list > 0] = self.pitch(ct0 + delta_ct, v0, tsr)
        else:
            # In case of a sinusoidal function, generate the sinusoid and determine the corresponding pitch time series
            ct = ct0 + delta_ct * np.sin(reduced_freq * v0 / self.blade.r * t_list)
            pitch = np.array([self.pitch(ctn, v0, tsr) for ctn in ct])

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
            if abs(pitch[n] - pitch[n-1]) < 1e-15 and n != 0:
                # Just reuse the thrust coefficient distribution from the previous time step
                ctr_qs[n, :] = ctr_qs[n - 1, :]
                cqr_qs[n, :] = cqr_qs[n - 1, :]

            # In case the pitch has changed since last time step
            else:
                # Run the BEM code for this pitch angle to set qs values
                self.blade.reset()
                self.blade.determine_cp_ct(v0, tsr, pitch[n])
                # Get the new qs thrust coefficient distribution
                ctr_qs[n, :] = c_thrust(self.blade.p_n_list[1:-1], v0, r_list, self.blade.b, dr)
                cqr_qs[n, :] = c_thrust(self.blade.p_t_list[1:-1], v0, r_list, self.blade.b, dr) * r_list / self.blade.r

            # Loop over the blade elements
            for i, be in enumerate(self.blade.blade_elements[1:-1]):
                # Set a tuple with parameters that the loads() function will need inside the different models
                params = (be.r, be.twist, be.c, self.blade.r, pitch[n], be.airfoil, v0, tsr * v0 / self.blade.r, 0, 0)

                a_qs[n, i] = be.a
                alpha_qs[n, i] = be.alpha
                phi_qs[n, i] = np.degrees(be.phi)

                # At the first time step, just initialise the output and intermediate value arrays
                if n == 0:
                    a[0, i] = be.a
                    pn, pt, alpha[0, i], phi[0, i] = loads(a[0, i], *params)
                    ctr[0, i] = c_thrust(pn, v0, be.r, self.blade.b, dr)
                    cqr[0, i] = c_thrust(pt, v0, be.r, self.blade.b, dr) * be.r / self.blade.r
                    v_int[0, i] = -a[0, i] * v0

                # If the model is Pitt-Peters
                elif model == 'pp':
                    # Propagate the AoA and induction factor of this blade element with pitt_peters()
                    a[n, i] = pitt_peters(ctr_qs[n, i], a[n-1, i], delta_t, params, dr, self.blade.b)

                # In case of Larsen-Madsen
                elif model == 'lm':
                    # Propagate the AoA and induction factor of this blade element with larsen_madsen().
                    # be.a is the quasi-steady induction factor that L-M requires
                    a[n, i] = larsen_madsen(be.a, a[n-1, i], delta_t, params)

                elif model == 'oye':
                    # Propagate the AoA and induction factor of this blade element with oye().
                    a[n, i], v_int[n, i] = oye(a_qs[n, i], a_qs[n-1, i], a[n-1, i], delta_t, params, v_int[n-1, i])

                pn, pt, alpha[n, i], phi[n, i] = loads(a[n, i], *params)
                ctr[n, i] = c_thrust(pn, v0, be.r, self.blade.b, dr)
                cqr[n, i] = c_thrust(pt, v0, be.r, self.blade.b, dr) * be.r / self.blade.r

        # Just a happy print statement bacause the code is done running :D
        print(f'Done! (Entire time series computed in {round(timer[-1] - timer[0], 3)} s)')

        # Return the outputs for later plotting
        return r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs


def pitt_peters(c_thrust_current, a_previous, dt, be_params, dr, b):
    """
    Calculate the new induction factor using Pitts-Peters
    :param c_thrust_current: Thrust coefficient at this time step
    :param a_previous: Induction factor at previous time step
    :param dt: Time step
    :param be_params: Parameters of the blade element required for the loads() function
        (r, twist, c, r_blade, pitch, airfoil, v_0, omega, yaw, azimuth)
    :param dr: Radial length of the blade element
    :param b: Number of turbine blades
    :return: The current time step: angle of attack, induction factor
    """
    # Determine the thrust loading on the blade element based on the previous time step induction factor
    p_n, _, _, _ = loads(a_previous, *be_params)
    # Use the thrust loading to determine the local thrust coefficient
    c_thrust_ind = c_thrust(p_n, be_params[6], be_params[0], b, dr)
    # print(c_thrust_current - c_thrust_ind)

    # Calculate the time derivative of the induction factor
    da_dt = (c_thrust_current - c_thrust_ind) / (16 / (3 * np.pi)) * (
             be_params[6] ** 2 / be_params[3]) / be_params[6]
    # Calculate the new induction factor with time propagation
    a_current = a_previous - da_dt * dt
    return a_current


def larsen_madsen(a_qs_current, a_previous, dt, be_params):
    """
    Calculate the new induction factor using the Larsen-Madsen model
    :param a_qs_current: Steady-state induction factor at this time step
    :param a_previous: Induction factor at previous time step
    :param dt: Time step
    :param be_params: Parameters of the blade element required for the loads() function
        (r, twist, c, r_blade, pitch, airfoil, v_0, omega, yaw, azimuth)
    :return: The current time step: angle of attack, induction factor
    """
    # # Determine the thrust loading on the blade element based on the previous time step induction factor
    # p_n, _, _ = loads(a_previous, *be_params)

    # Evaluate the wake velocity
    v_wake = be_params[6] * (1 - a_previous)

    # Evaluate the time scale time scale of the model
    tau = 0.5 * be_params[3] / v_wake

    # Evaluate the transient and quasi steady induction factors
    a_transient = a_previous * np.exp(-dt/tau)
    a_quasteady = a_qs_current * (1 - np.exp(-dt/tau))

    # Evaluate the new induction factor
    a_current = a_transient + a_quasteady

    # Evaluate the time rate of change of the induction factor
    _ = (a_previous - a_current)/dt

    return a_current


def oye(a_qs_current, a_qs_previous, a_previous, dt, be_params, v_int_previous):
    """
    Calculate the new induction factor using the Oye model
    :param a_qs_current: Steady-state induction factor at this time step
    :param a_qs_previous: Steady-state induction factor at previous time step
    :param a_previous: Induction factor at previous time step
    :param dt: Time step
    :param be_params: Parameters of the blade element required for the loads() function
        (r, twist, c, r_blade, pitch, airfoil, v_0, omega, yaw, azimuth)
    :param v_int_previous: intermediate induced velocity at the previous time step
    :return: The current time step: angle of attack, induction factor and intermediate induced velocity
    """
    # # Determine the thrust loading on the blade element based on the previous time step induction factor
    # _, _, alpha = loads(a_previous, *be_params)

    # calculate quasi-steady induction velocity
    v_qs_previous = -a_qs_previous * be_params[6]
    # calculate induction velocity of the previous time step
    v_ind_previous = a_previous * be_params[6]

    # calculate time scales of the model
    t1 = 1.1 / (1 - 1.3 * a_previous) * be_params[3] / be_params[6]
    t2 = (0.39 - 0.26 * (be_params[0] / be_params[3])**2) * t1

    # calculate next-time-step quasi-steady induction velocity
    v_qs_current = -a_qs_current * be_params[6]

    # calculate time derivative of intermediate velocity
    dvint_dt = (v_qs_previous + (v_qs_current - v_qs_previous) / dt * 0.6 * t1 - v_int_previous) / t1

    # calculate new intermediate velocity
    v_int_current = v_int_previous + dvint_dt * dt

    # calculate time derivaive of the induced velocity
    dvz_dt = ((v_int_previous + v_int_current) / 2 + v_ind_previous) / t2

    # calculate new induced velocity
    a_current = (v_ind_previous - dvz_dt * dt) / be_params[6]
    return a_current, v_int_current


def generate_data():
    # Create the turbine with 25 blade elements
    turbine = Turbine(25)

    for i, model in enumerate(('pp', 'lm', 'oye')):
        print(f'============== MODEL = {model.upper()} ==============')
        for case in ct_steps:
            print(f'----- ct0 = {case[0]}, d_ct = {case[1]}, rf = {0} -----')
            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = turbine.ct_func(*case, None, 10, 10, model=model)
            write_to_file([r_list, ], f'./{model}/ct_step/{case[0]}_{case[1]}_r_list.csv')
            write_to_file([t_list, ], f'./{model}/ct_step/{case[0]}_{case[1]}_t_list.csv')
            write_to_file(ctr, f'./{model}/ct_step/{case[0]}_{case[1]}_ctr.csv')
            write_to_file(cqr, f'./{model}/ct_step/{case[0]}_{case[1]}_cqr.csv')
            write_to_file(a, f'./{model}/ct_step/{case[0]}_{case[1]}_a.csv')
            write_to_file(alpha, f'./{model}/ct_step/{case[0]}_{case[1]}_alpha.csv')
            write_to_file(phi, f'./{model}/ct_step/{case[0]}_{case[1]}_phi.csv')
            write_to_file(ctr_qs, f'./{model}/ct_step/{case[0]}_{case[1]}_ctr_qs.csv')
            write_to_file(cqr_qs, f'./{model}/ct_step/{case[0]}_{case[1]}_cqr_qs.csv')
            write_to_file(a_qs, f'./{model}/ct_step/{case[0]}_{case[1]}_a_qs.csv')
            write_to_file(alpha_qs, f'./{model}/ct_step/{case[0]}_{case[1]}_alpha_qs.csv')
            write_to_file(phi_qs, f'./{model}/ct_step/{case[0]}_{case[1]}_phi_qs.csv')
            print()

        for case in ct_sins:
            for rf in np.round(np.arange(0.05, 0.35, 0.05), 2):
                print(f'----- ct0 = {case[0]}, d_ct = {case[1]}, rf = {rf} -----')
                r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = turbine.ct_func(*case, rf, 10, 10, model=model)
                write_to_file([r_list, ], f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_r_list.csv')
                write_to_file([t_list, ], f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_t_list.csv')
                write_to_file(ctr, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_ctr.csv')
                write_to_file(cqr, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_cqr.csv')
                write_to_file(a, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_a.csv')
                write_to_file(alpha, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_alpha.csv')
                write_to_file(phi, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_phi.csv')
                write_to_file(ctr_qs, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_ctr_qs.csv')
                write_to_file(cqr_qs, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_cqr_qs.csv')
                write_to_file(a_qs, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_a_qs.csv')
                write_to_file(alpha_qs, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_alpha_qs.csv')
                write_to_file(phi_qs, f'./{model}/ct_sin/{case[0]}_{case[1]}_{rf}_phi_qs.csv')
                print()

        for case in u_inf_steps:
            print(f'----- u0 = {case[0]}, d_u = {case[1]}, rf = {0} -----')
            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = turbine.u_inf_func(*case, None, 10, 10, model=model)
            write_to_file([r_list, ], f'./{model}/u_inf_step/{case[0]}_{case[1]}_r_list.csv')
            write_to_file([t_list, ], f'./{model}/u_inf_step/{case[0]}_{case[1]}_t_list.csv')
            write_to_file(ctr, f'./{model}/u_inf_step/{case[0]}_{case[1]}_ctr.csv')
            write_to_file(cqr, f'./{model}/u_inf_step/{case[0]}_{case[1]}_cqr.csv')
            write_to_file(a, f'./{model}/u_inf_step/{case[0]}_{case[1]}_a.csv')
            write_to_file(alpha, f'./{model}/u_inf_step/{case[0]}_{case[1]}_alpha.csv')
            write_to_file(phi, f'./{model}/u_inf_step/{case[0]}_{case[1]}_phi.csv')
            write_to_file(ctr_qs, f'./{model}/u_inf_step/{case[0]}_{case[1]}_ctr_qs.csv')
            write_to_file(cqr_qs, f'./{model}/u_inf_step/{case[0]}_{case[1]}_cqr_qs.csv')
            write_to_file(a_qs, f'./{model}/u_inf_step/{case[0]}_{case[1]}_a_qs.csv')
            write_to_file(alpha_qs, f'./{model}/u_inf_step/{case[0]}_{case[1]}_alpha_qs.csv')
            write_to_file(phi_qs, f'./{model}/u_inf_step/{case[0]}_{case[1]}_phi_qs.csv')
            print()

        for case in u_inf_sins:
            for rf in np.round(np.arange(0.05, 0.35, 0.05), 2):
                print(f'----- u0 = {case[0]}, d_u = {case[1]}, rf = {rf} -----')
                try:
                    r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = turbine.u_inf_func(*case, rf, 10, 10, model=model)
                except ValueError:
                    global relaxation
                    relaxation = .25
                    r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = turbine.u_inf_func(*case, rf, 10, 10, model=model)
                    relaxation = 1

                write_to_file([r_list, ], f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_r_list.csv')
                write_to_file([t_list, ], f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_t_list.csv')
                write_to_file(ctr, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_ctr.csv')
                write_to_file(cqr, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_cqr.csv')
                write_to_file(a, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_a.csv')
                write_to_file(alpha, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_alpha.csv')
                write_to_file(phi, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_phi.csv')
                write_to_file(ctr_qs, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_ctr_qs.csv')
                write_to_file(cqr_qs, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_cqr_qs.csv')
                write_to_file(a_qs, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_a_qs.csv')
                write_to_file(alpha_qs, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_alpha_qs.csv')
                write_to_file(phi_qs, f'./{model}/u_inf_sin/{case[0]}_{case[1]}_{rf}_phi_qs.csv')
                print()


if __name__ == '__main__':
    raise RuntimeError("Do not run this file, it has no use.")
