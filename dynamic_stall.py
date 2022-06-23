import matplotlib.pyplot as plt
import numpy as np
import time

from BEM_code import Blade, BladeElement, loads, c_thrust
from dynamic_inflow import Turbine, pitt_peters, larsen_madsen, oye
from read_write import read_from_file, write_to_file


a_s0 = 343                  # [m/s]
k_alpha = .75               # [-]
Tp = 1.7                    # [-]
Tf = 3.0                    # [-]
alpha_f1 = 7 * np.pi / 180  # [rad]
alpha_f2 = 15 * np.pi / 180 # [rad]
alpha_f3 = 21 * np.pi / 180 # [rad]
cn_1 = 1.0093               # [-]
tv = 6.0                    # [-]
tvl = 5.0                   # [-]

# cases = {'Dyn1': (1, .5, .05, 10, 8), 'Dyn2': (1, .5, .3, 10, 8), 'no_les': (1, .5, .05, 10, 8, True)}
cases = {'no_les': (1, .5, .05, 10, 8, True)}

## *** Nomenclature ***
#   AoA     Angle of attack
#   TE      Trailing edge
#   CSF     Carlos Simoa Ferreira


class DSAirfoil:
    """
    Replacement for the DU95W150 airfoil class for the dynamic stall model
    """
    def __init__(self, chord):
        data = read_from_file('DU95W150.csv')
        self.alpha_lst = data[:, 0]
        self.cd_lst = data[:, 2]

        self.t = None
        self.delta_t = None
        self.s = 0
        self.delta_s = None

        self.velocity_scale = None

        self.alpha_0 = -2 * np.pi / 180
        self.dcl_dalpha = 2 * np.pi
        self.c = chord
        self.no_les = False

        # Initialise "old" values for module 1
        self.x_lag = 0
        self.y_lag = 0
        self.d_nc = 0
        self.alpha = None
        self.dalpha_qs_dt = 0
        self.cn_p = None
        self.alpha_eq = None

        # Initialise "old" values for module 2
        self.d_pf = 0
        self.d_bl = 0
        self.f_p = 0

        # Initialise "old" values for module 3
        self.tau_v = 0

        # Initialise "old" values for module 4
        self.c_v = None
        self.cn_v = 0

        # Wagner function coefficients
        self.a1 = .165
        self.a2 = .335
        self.b1 = .045
        self.b2 = .3

    def set_time_data(self, t, delta_t, delta_s, velocity_scale, no_les=False):
        self.t = t
        self.delta_t = delta_t

        self.delta_s = delta_s
        self.velocity_scale = velocity_scale

        self.no_les = no_les

    def cl(self, alpha, propagate: bool = False, **kwargs):
        """
        I will program the dynamic stall loop here
        :param alpha: the quasi-steady angle of attack in degrees
        :param propagate: Indicate whether or not to propagate the airfoil in time
        :return: The dynamic stall lift coefficient
        """
        # Convert angle of attack to radians, as that is needed by all further equations
        alpha = np.radians(alpha)
        # Initialise the angle of attack at the first time step
        if self.s == 0:
            self.alpha = alpha

        # Run through module 1 and determine the total cn of this module
        cn_c, cn_nc, alpha_eq = self._module1(alpha, propagate)
        cn_p = cn_c + cn_nc
        # Initialise the cn_p at the first time step
        if self.s == 0:
            self.cn_p = cn_p

        # Run through module 2
        cn_f, f_bl = self._module2(alpha, cn_p, cn_nc, alpha_eq, propagate)
        # Initialise the equivalent angle of attack at the first time step
        if self.s == 0:
            self.alpha_eq = alpha_eq

        if not self.no_les:
            # Run through modules 3 and 4
            tau_v = self._module3(cn_f, alpha_eq)
            cn_v = self._module4(cn_p, f_bl, tau_v, propagate)
        else:
            # Ignore leading edge separation effects if so desired
            cn_v = 0
            tau_v = 0

        # If time propogation is required, set the current values into the "old" values
        if propagate:
            self.alpha = alpha
            self.alpha_eq = alpha_eq
            self.cn_p = cn_p
            self.cn_v = cn_v
            self.tau_v = tau_v
            self.s = self.s + self.delta_s

        # Return the total force coefficient, which in this case is equivalent to the lift coefficient
        return cn_f + cn_v

    def _module1(self, alpha, propagate):
        """
        Unsteady attached flow module
        :param alpha: quasi-steady angle of attack
        :param propagate: indication whether to propagate
        :return: circulatory, non-circulatory normal force and equivalent angle of attack
        """
        # Evaluate the quasi-steady angle of attack difference and time derivative
        dalpha_qs = self.alpha - alpha
        dalpha_qs_dt_new = dalpha_qs / self.delta_t

        # Evaluate the current lag states
        x_lag_new = self.x_lag * np.exp(-self.b1 * self.delta_s) +\
            dalpha_qs * self.a1 * np.exp(-self.b1 * self.delta_s / 2)
        y_lag_new = self.y_lag * np.exp(-self.b2 * self.delta_s) +\
            dalpha_qs * self.a2 * np.exp(-self.b2 * self.delta_s / 2)

        # Evaluate the equivalent angle of attack
        alpha_eq = alpha - x_lag_new - y_lag_new

        # Evaluate the circulatory normal force coefficient
        cn_c = self.dcl_dalpha * (alpha_eq - self.alpha_0)

        # Evaluate the deficiency function for the non-circulatory normal force coefficient
        d_nc_new = self.d_nc * np.exp(-a_s0 * self.delta_t / (k_alpha * self.c)) +\
            (dalpha_qs_dt_new - self.dalpha_qs_dt) * np.exp(-a_s0 * self.delta_t / (2 * k_alpha * self.c))

        # Evaluate the non-circulatory normal force coefficient
        cn_nc = 4 * k_alpha * self.c / self.velocity_scale * (dalpha_qs_dt_new - d_nc_new)

        # If this is the final run at this timestep, propagate the values required for the next time step
        if propagate:
            self.x_lag = x_lag_new
            self.y_lag = y_lag_new
            self.d_nc = d_nc_new
            self.dalpha_qs_dt = dalpha_qs_dt_new

        # Return the circulatory and non-circulatory normal force coefficients
        return cn_c, cn_nc, alpha_eq

    def _module2(self, alpha, cn_p, cn_nc, alpha_eq, propagate):
        """
        Non-linear trailing edge separation module
        :param alpha: quasi-steady angle of attack
        :param cn_p: sum of cn_c and cn_nc from module 1
        :param cn_nc: non-circulatory normal force from module 1
        :param alpha_eq: equivalent angle of attack from module 1
        :param propagate: indication whether to propagate
        :return: normal force and boundary layer function value
        """
        # Evaluate the pressure lag deficiency function
        d_pf_new = self.d_pf * np.exp(-self.delta_s / Tp) + (self.cn_p - cn_p) * np.exp(-self.delta_s / (2 * Tp))
        # Evaluate the lagged potential flow load
        cn_p_prime = cn_p - d_pf_new
        # Evaluate the equivalent angle of attack
        alpha_f = cn_p_prime / self.dcl_dalpha + self.alpha_0

        # Evaluate the pseudo-location of the separation point over the airfoil section (Taken from CSF)
        if alpha_f <= alpha_f1:
            f_sep = 1  # f = 1 : Separation at TE
        elif alpha_f <= alpha_f2:
            f_sep = 1 - 0.8 * ((alpha_f - alpha_f1) / (alpha_f2 - alpha_f1))
        elif alpha_f <= alpha_f3:
            f_sep = 0.2 * (1 - ((alpha_f - alpha_f2) / (alpha_f3 - alpha_f2)) ** 0.3)
        else:
            f_sep = 0  # f = 0 : Separation at LE

        # Evaluate the non-linear normal force coefficient in steady flow including TE separation
        cn_st = self.dcl_dalpha * ((1 + np.sqrt(f_sep)) / 2) ** 2 * (alpha - self.alpha_0)
        # Evaluate the separation location for the lagged potential flow
        f_p_new = (2 * (cn_st/(self.dcl_dalpha * (alpha_f - self.alpha_0)))**(1/2) - 1)**2
        # Initialise the value for f_p at the first time step since it is non-zero initialy
        if self.s == 0:
            self.f_p = f_p_new

        # Evaluate the deficiency function for the bondary-layer development
        d_bl_new = self.d_bl * np.exp(-self.delta_s / Tf) + (f_p_new - self.f_p) * np.exp(-self.delta_s / (2 * Tf))
        # Evaluate the separation location for the boundary-layer
        f_bl = f_p_new - d_bl_new
        # Evaluate the unsteady non-linear normal force coefficient including TE separation
        cn_f = self.dcl_dalpha * ((1 + np.sqrt(f_bl)) / 2) ** 2 * (alpha_eq - self.alpha_0) + cn_nc

        # If this is the final run at this timestep, propagate the values required for the next time step
        if propagate:
            self.d_pf = d_pf_new
            self.d_bl = d_bl_new
            self.f_p = f_p_new

        return cn_f, f_bl

    def _module3(self, cn_f, alpha_eq):
        """
        Leading edge flow separation module
        :param cn_f: normal force from module 2
        :param alpha_eq: equivalent angle of attack from module 1
        :return: the time scale of the leading edge flow separation
        """
        if cn_f > cn_1 or self.alpha_eq - alpha_eq <= 0:
            tau_v = self.tau_v + .45 * self.delta_s
        else:
            tau_v = 0

        return tau_v

    def _module4(self, cn_p, f_bl, tau_v, propagate):
        """
        Vortex shedding module
        :param cn_p: sum of cn_c and cn_nc from module 1
        :param f_bl: boundary layer function value from module 2
        :param tau_v: the time scale from module 3
        :param propagate: indication whether to propagate
        :return: vortex normal force
        """
        c_v = cn_p * (1 - ((1 + np.sqrt(f_bl)) / 2) ** 2)
        if self.s == 0:
            self.c_v = c_v

        if 0 < tau_v < tvl:
            cn_v = self.cn_v * np.exp(-self.delta_s / tv) + (c_v - self.c_v) * np.exp(-self.delta_s / 2 / tv)
        else:
            cn_v = self.cn_v * np.exp(-self.delta_s / tv)

        if propagate:
            self.c_v = c_v

        return cn_v

    def cd(self, alpha): return np.interp(alpha, self.alpha_lst, self.cd_lst)


class DSBladeElement(BladeElement):
    """
    Class to represent blade elements in the dynamic stall model. This one has all the functions from BladeElement with
    the addition that it keeps track of time for the dynamic stall model.
    """
    def __init__(self, pos_r: float, chord: float, twist: float, airfoil: DSAirfoil, r_blade):
        super().__init__(pos_r, chord, twist, airfoil, r_blade)

        self.t = None
        self.delta_t = None
        self.delta_s = None

    def __repr__(self): # print blade element
        return f"<DS Blade Element at r={self.r}, c={self.c}, beta={self.twist}>"

    def set_time_data(self, t, delta_t, v0, tsr, no_les=False):
        self.t = t
        self.delta_t = delta_t

        velocity_scale = v0 * np.sqrt(1 + (tsr * tsr * self.r * self.r) / (self.r_blade * self.r_blade))
        self.delta_s = 2 * delta_t * velocity_scale / self.c

        self.airfoil.set_time_data(t, delta_t, self.delta_s, velocity_scale, no_les=no_les)


class DSBlade(Blade):
    """
    Class to represent the blade in the dynamic stall model. This one has all the functions from Blade with
    the addition that it keeps track of time for the dynamic stall model.
    """
    def __init__(self, n_blades, airfoil, r_start, r_end, blade_pitch, n_elements):
        super().__init__(n_blades, airfoil, r_start, r_end, blade_pitch, n_elements, create_elements=False)

        self.r_list = []
        self.blade_elements = list()

        # Divide the blade up in n_elements pieces;
        for i in range(n_elements + 1):
            r = r_start + (r_end - r_start) / n_elements * i
            self.r_list.append(r)
            # Sorry for hardcoding the equations below- taken from the assignment description :)
            twist = 14 * (1 - r / r_end)
            chord = (3 * (1 - r / r_end) + 1)

            relative_pitch = blade_pitch + twist

            self.blade_elements.append(DSBladeElement(r, chord, relative_pitch, airfoil(chord), self.r))

        self.r_list = np.array(self.r_list)

    def set_time_data(self, t, delta_t, v0, tsr, no_les=False):
        [be.set_time_data(t, delta_t, v0, tsr, no_les=no_les) for be in self.blade_elements]


class DSTurbine:
    """
    New turbine class for the dynamic stall model. Basically the same as the Turbine class for dynamic inflow, but
    restricted to sinusoidal velocity signals instead of the whole array of signals we needed in the previous assignment
    """
    def __init__(self, n_annuli):
        self.blade = DSBlade(3, DSAirfoil, .2 * 50, 50, -2, n_annuli)

    def u_inf_func(self, u_inf_0, delta_u_inf, reduced_freq, v0, tsr, no_les=False, model='pp'):
        """
        Determine and plot the time evolution of the turbine properties given a step in thrust coefficient
        :param u_inf_0: Mean inflow velocity
        :param delta_u_inf: Amplitude of the inflow velocity variation
        :param reduced_freq: Reduced frequency of the dynamic inflow
        :param v0: The incoming velocity
        :param tsr: The turbine tip-speed ratio
        :param no_les: Trigger to turn of the leading edge separation (modules 3-4)
        :param model: Selection of the dynamic inflow model (pp: Pitt-Peters, lm: Larsen-Madsen, oye: Oye)
        :return: None
        """
        if model not in ('pp', 'lm', 'oye'):
            raise ValueError("Unknown model, please enter one of the following: 'pp', 'lm', 'oye'.")

        # Initialise a timer list to check compute time
        timer = [time.time(), ]

        # Initialise the time parameters: time step, start and final time,
        # based on the reduced frequency and reduced time
        delta_t = .04 * self.blade.r / v0
        t_0 = -.2 * self.blade.r / v0
        t_final = 4 * np.pi / reduced_freq * self.blade.r / v0
        t_list = np.round(np.arange(t_0, t_final + delta_t, delta_t), 9)

        # Extract the radial positions of the blade elements and the radial length
        r_list = self.blade.r_list[1:-1]
        dr = r_list[1] - r_list[0]

        # Generate the sinusoid
        u_inf = u_inf_0 + delta_u_inf / v0 * np.cos(reduced_freq * v0 / self.blade.r * t_list)
        u_inf *= v0

        # Initialise the output value arrays: induction, AoA, inflow angle, thrust coefficient and torque coefficient.
        # The shape is (time series x spanwise distribution).
        a = np.empty((t_list.size, r_list.size))
        alpha = np.empty((t_list.size, r_list.size))
        phi = np.empty((t_list.size, r_list.size))
        ctr = np.empty((t_list.size, r_list.size))
        cqr = np.empty((t_list.size, r_list.size))
        # Initialise an extra set of airfoils for determining the unsteady response with the dynamic inflow models
        # this is done because the unsteady time evolution will be different from
        # the quasi-steady dynamic stall response
        airfoils_us = [DSAirfoil(be.c) for be in self.blade.blade_elements]

        # Initialise the intermediate induced velocity array.
        # The shape is (time series x spanwise distribution).
        v_int = np.empty((t_list.size, r_list.size))

        # Initialise the quasi-steady value arrays: induction, AoA, inflow angle, thrust coefficient and
        # torque coefficient. I call it quasi-steady, but it really is the dynamic stall without dynamic inflow.
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

            # Run the BEM code for the current wind speed. Time variables of the airfoils are set here as well.
            self.blade.reset()
            self.blade.set_time_data(t, delta_t, u_inf[n], tsr * v0 / u_inf[n], no_les=no_les)
            self.blade.determine_cp_ct(u_inf[n], tsr * v0 / u_inf[n], 0)
            # Get the new qs thrust and torque coefficient distribution
            ctr_qs[n, :] = c_thrust(self.blade.p_n_list[1:-1], u_inf[n], r_list, self.blade.b, dr)
            cqr_qs[n, :] = c_thrust(self.blade.p_t_list[1:-1], u_inf[n], r_list, self.blade.b, dr) * r_list / self.blade.r

            # Loop over the blade elements
            for i, be in enumerate(self.blade.blade_elements[1:-1]):
                # Set the unsteady airfoil's time data as well, since that hasn't happened above.
                airfoils_us[i].set_time_data(t, delta_t, u_inf[n], tsr * v0 / u_inf[n], no_les=no_les)
                # Set a tuple with parameters that the loads() function will need inside the different models
                params_us = (be.r, be.twist, be.c, self.blade.r, 0, airfoils_us[i],
                             u_inf[n], tsr * v0 / self.blade.r, 0, 0)

                # Set the quasi-steady induction, AoA and inflow angle of this blade element
                a_qs[n, i] = be.a
                alpha_qs[n, i] = be.alpha
                phi_qs[n, i] = be.phi

                # At the first time step, initialise the output and intermediate value arrays with
                # the quasi-steady values
                if n == 0:
                    a[0, i] = be.a
                    pn, pt, alpha[0, i], phi[0, i] = loads(a[0, i], *params_us)
                    ctr[0, i] = c_thrust(pn, u_inf[n], be.r, self.blade.b, dr)
                    cqr[0, i] = c_thrust(pt, u_inf[n], be.r, self.blade.b, dr) * be.r / self.blade.r
                    v_int[0, i] = -a[0, i] * v0

                # If the model is Pitt-Peters
                elif model == 'pp':
                    # Propagate the induction factor of this blade element with pitt_peters()
                    a[n, i] = pitt_peters(ctr_qs[n, i], a[n - 1, i], delta_t, params_us, dr, self.blade.b)

                # In case of Larsen-Madsen
                elif model == 'lm':
                    # Propagate the induction factor of this blade element with larsen_madsen().
                    # be.a is the quasi-steady induction factor that L-M requires
                    a[n, i] = larsen_madsen(be.a, a[n - 1, i], delta_t, params_us)

                elif model == 'oye':
                    # Propagate the induction factor of this blade element with oye().
                    a[n, i], v_int[n, i] = oye(a_qs[n, i], a_qs[n - 1, i], a[n - 1, i], delta_t, params_us,
                                               v_int[n - 1, i])

                # Determine the loads based on the unsteady induction factor above.
                pn, pt, alpha[n, i], phi[n, i] = loads(a[n, i], *params_us, propagate=True)
                ctr[n, i] = c_thrust(pn, u_inf[n], be.r, self.blade.b, dr)
                cqr[n, i] = c_thrust(pt, u_inf[n], be.r, self.blade.b, dr) * be.r / self.blade.r

        # Just a happy print statement bacause the code is done running :D
        print(f'Done! (Entire time series computed in {round(timer[-1] - timer[0], 3)} s)')

        # Return the outputs for later plotting
        return r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs


def generate_data():
    turbine_dynamic_stall = DSTurbine(25)
    turbine_dynamic_inflow = Turbine(25)

    for name, case in cases.items():
        print(f'============== CASE:\t{name.upper()} ==============')

        for i, model in enumerate(('pp', 'lm', 'oye')):
            print(f' ------------- Model:\t{model}\t ------------- ')
            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_ds, cqr_ds, a_ds, alpha_ds, phi_ds = \
                turbine_dynamic_stall.u_inf_func(*case, model=model)
            *_, ctr_di, cqr_di, a_di, alpha_di, phi_di, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = \
                turbine_dynamic_inflow.u_inf_func(*case[:5], model=model, ds=True)

            # Save the discrete time and blade nodes
            write_to_file([r_list, ], f'./{model}/{name}/r_list.csv')
            write_to_file([t_list, ], f'./{model}/{name}/t_list.csv')

            # Save fully unsteady
            write_to_file(ctr, f'./{model}/{name}/ctr.csv')
            write_to_file(cqr, f'./{model}/{name}/cqr.csv')
            write_to_file(a, f'./{model}/{name}/a.csv')
            write_to_file(alpha, f'./{model}/{name}/alpha.csv')
            write_to_file(phi, f'./{model}/{name}/phi.csv')

            # Save dynamic stall
            write_to_file(ctr_ds, f'./{model}/{name}/ctr_ds.csv')
            write_to_file(cqr_ds, f'./{model}/{name}/cqr_ds.csv')
            write_to_file(a_ds, f'./{model}/{name}/a_ds.csv')
            write_to_file(alpha_ds, f'./{model}/{name}/alpha_ds.csv')
            write_to_file(phi_ds, f'./{model}/{name}/phi_ds.csv')

            # Save dynamic inflow
            write_to_file(ctr_di, f'./{model}/{name}/ctr_di.csv')
            write_to_file(cqr_di, f'./{model}/{name}/cqr_di.csv')
            write_to_file(a_di, f'./{model}/{name}/a_di.csv')
            write_to_file(alpha_di, f'./{model}/{name}/alpha_di.csv')
            write_to_file(phi_di, f'./{model}/{name}/phi_di.csv')

            # Save quasi-steady
            write_to_file(ctr_qs, f'./{model}/{name}/ctr_qs.csv')
            write_to_file(cqr_qs, f'./{model}/{name}/cqr_qs.csv')
            write_to_file(a_qs, f'./{model}/{name}/a_qs.csv')
            write_to_file(alpha_qs, f'./{model}/{name}/alpha_qs.csv')
            write_to_file(phi_qs, f'./{model}/{name}/phi_qs.csv')


if __name__ == '__main__':
    raise RuntimeError("Do not run this file, it has no use.")
