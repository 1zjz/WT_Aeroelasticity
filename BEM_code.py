import numpy as np
import scipy.integrate as spig
import matplotlib.pyplot as plt
import time


relaxation = 0.1
rho = 1.225
p_atm = 101325


class DU95W150:
    def __init__(self):
        data = read_from_file('DU95W150.csv')
        self.alpha_lst = data[:, 0]
        self.cl_lst = data[:, 1]
        self.cd_lst = data[:, 2]
        self.cm_lst = data[:, 3]

    def cl(self, alpha): return np.interp(alpha, self.alpha_lst, self.cl_lst)

    def cd(self, alpha): return np.interp(alpha, self.alpha_lst, self.cd_lst)

    def cm(self, alpha): return np.interp(alpha, self.alpha_lst, self.cm_lst)

    def plot_polars(self, axes):
        axes[0].plot(self.alpha_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)], 'k')
        axes[1].plot(self.cd_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)], 'k')

        optimal = np.argmax(self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)] /
                            self.cd_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)])

        axes[0].plot(self.alpha_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal], 'ro')
        axes[1].plot(self.cd_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal], 'ro')


class BladeElement: #one of blade elements, the self.r stores radial position in 
#the blade element function 
    def __init__(self, pos_r: float, chord: float, relative_pitch: float, airfoil):
        self.r = pos_r
        self.c = chord
        self.beta = relative_pitch

        self.a = None #at some point we will also get this
        # self.a_prime = None
        self.axial_induction = None
        self.azimuthal_induction = None
        self.phi = None
        self.alpha = None
        self.p_n = None
        self.p_t = None
        self.u_tangential = None
        self.u_normal = None

        self.af = airfoil
        self.airfoil = airfoil()

    def __repr__(self): #print blade element
        return f"<Blade Element at r={self.r}, c={self.c}, beta={self.beta}>"

    def determine_loads(self, v_0, omega, theta_p, b, r_blade, r_root, yaw=0, azimuth=0, loss=True):
        # bem code for a single blade element, loss is for tip/root loss
        yaw = np.radians(yaw)
        azimuth = np.radians(azimuth)
        # Set initial loop values
        self.a = 0
        error_a = 1
        i = 0
        # Iterative solver for a and a_prime until the difference between the iterations becomes very small
        while True:
            self.u_tangential = omega * self.r # * (1 + self.a_prime)
            self.u_normal = v_0 * (1 - self.a)

            # For the previous a and a_prime, find the flow angle and angle of attack
            self.phi = np.arctan2(self.u_normal, self.u_tangential)
            self.alpha = np.degrees(self.phi) - self.beta - theta_p

            # With the angle of attack, determine the lift and drag coefficient from airfoil data interpolation
            cl = self.airfoil.cl(self.alpha)
            cd = self.airfoil.cd(self.alpha)

            # Use these to find the normal and tangential force coefficients
            cn = cl * np.cos(self.phi) + cd * np.sin(self.phi)
            ct = cl * np.sin(self.phi) - cd * np.cos(self.phi)

            # Break conditions for the a-loop
            if error_a <= 1e-9: # and error_a_dash <= 1e-9:
                break
            elif i > 1e3:
                raise ValueError(f"r={self.r}: Solution for a not converging. a={self.a}.")

            # Determine the solidity and Prandtl’s tip loss correction
            solidity = self.c * b / (2 * np.pi * self.r)
            f_tip = (2/np.pi) * np.arccos(np.exp(-(b * (r_blade - self.r) / (2 * self.r * np.sin(abs(self.phi)))))) if loss else 1
            f_root = (2 / np.pi) * np.arccos(np.exp(-(b * (self.r - r_root) / (2 * self.r * np.sin(abs(self.phi))))))
            f = f_root * f_tip

            # Determine the new a and a_prime
            # If it's higher than 0.33, use a glauert correction
            if self.a >= 0.33:
                c_thrust = ((1 - self.a) ** 2 * cn * solidity) / (np.sin(self.phi) ** 2)

                a_star = c_thrust / (4 * f * (1 - 0.25*(5 - 3 * self.a) * self.a))
                a_new = relaxation * a_star + (1-relaxation) * self.a

            else:
                a_new = 1 / ((4 * f * np.sin(self.phi)**2) / (solidity * cn) + 1)

            # Determine the difference between this and the previous iteration
            error_a = abs(a_new - self.a)

            # Get ready for the next iteration
            self.a = a_new
            i += 1

        self.p_n, self.p_t, _ = loads(self.a, self.r, self.beta, self.c, r_blade, theta_p, self.airfoil,
                                   v_0, omega, yaw, azimuth)

    def get_loads(self):
        if self.p_t is None or self.p_n is None:
            raise ValueError(f"Loads have not been determined. Run .determine_loads() first.")
        else:
            return self.p_n, self.p_t

    def reset(self):
        self.__init__(self.r, self.c, self.beta, self.af)


class Blade:
    def __init__(self, no_blades, airfoil, r_start, r_end, blade_pitch, n_elements):
        self.b = no_blades

        self.power = None
        self.thrust = None
        self.c_power = None
        self.c_thrust = None

        self.r_list = []
        self.p_n_list = None
        self.p_t_list = None

        self.blade_elements = list()
        # Divide the blade up in n_elements pieces;
        for i in range(n_elements + 1):
            r = r_start + (r_end - r_start)/n_elements * i
            self.r_list.append(r)
            # Sorry for hardcoding the equations below- taken from the assignment description :)
            twist = 14*(1-r/r_end)
            chord = (3*(1-r/r_end)+1)

            # BladeElement takes in argument relative_pitch, I assume that this means total? So offset with the blade pitch
            relative_pitch = blade_pitch + twist

            self.blade_elements.append(BladeElement(r, chord, relative_pitch, airfoil))

        self.r_list = np.array(self.r_list)
        self.r = r_end

    def find_pn_pt(self, v_0, theta_p, omega, yaw=0, azimuth=0, loss=True):
        # Initialise the lists for p_n and p_t
        p_n_list, p_t_list = list(), list()
        for blade_element in self.blade_elements:
            if self.r_list[0] < blade_element.r < self.r:
                blade_element.determine_loads(v_0, omega, theta_p, self.b, self.r, self.r_list[0], yaw, azimuth, loss)
                p_n, p_t = blade_element.get_loads()

                p_n_list.append(p_n)
                p_t_list.append(p_t)

            else:
                # Add zero load at the blade tip and root
                p_n_list.append(0)
                p_t_list.append(0)

        return np.array(p_n_list), np.array(p_t_list), self.r_list

    def determine_cp_ct(self, v_0, lamda, theta_p, yaw=0, azimuth=0, loss=True):
        """
        Let the BEM code Determine the power and thrust coefficient of the turbine
        :param v_0: Incoming velocity
        :param lamda: Tip speed ratio
        :param theta_p: pitch angle
        :param yaw: yaw angle
        :param azimuth: azimuthal position in the turbine disk
        :param loss: turn on or off the tip loss factor
        :return: None
        """
        # Determine the rotational speed of the turbine
        omega = lamda * v_0 / self.r
        # Get the loads on the blade elements
        self.p_n_list, self.p_t_list, r_list = self.find_pn_pt(v_0, theta_p, omega, yaw, azimuth, loss)

        # Determine the thrust and power of the turbine
        self.thrust = self.b * spig.trapz(self.p_n_list, self.r_list)
        self.power = omega * self.b * spig.trapz(self.p_t_list * self.r_list, self.r_list)

        # Determine the thrust and power coefficient
        self.c_thrust = self.thrust / (0.5 * rho * np.pi * self.r**2 * v_0**2)
        self.c_power = self.power / (0.5 * rho * np.pi * self.r**2 * v_0**3)

    def reset(self):
        for be in self.blade_elements:
            be.reset()


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

    def ct_step(self, ct1, ct2, v0, tsr):
        """
        Determine and plot the time evolution of the turbine properties given a step in thrust coefficient
        :param ct1: Initial thrust coefficient
        :param ct2: Final thrust coefficient
        :param v0: The incoming velocity
        :param tsr: The turbine tip-speed ratio
        :return: None
        """
        # Initialise the time parameters: time step, start and final time
        delta_t = 0.01
        t_0, t_final = -5 * v0 / self.blade.r, 250 * v0 / self.blade.r
        # Extract the radial positions of the blade elements and the radial length of each
        r_list = self.blade.r_list
        dr = r_list[1] - r_list[0]

        # Determine the pitch angles required for both thrust coefficients
        pitch_1 = self.pitch(ct1, v0, tsr)
        pitch_2 = self.pitch(ct2, v0, tsr)

        # Use the BEM code to determine the local thrust coefficients based on the thrust loading distribution
        # of the blade at the initial pitch angle.
        # Also determine the induction factor distribution for the initial thrust coefficient
        self.blade.determine_cp_ct(v0, tsr, pitch_1)
        c_thrust_r_1 = c_thrust(self.blade.p_n_list, v0, r_list, self.blade.b, dr)

        a_0 = np.empty(r_list.shape)
        alpha_0 = np.empty(r_list.shape)
        blade_params = []
        for i, be in enumerate(self.blade.blade_elements):
            a_0[i] = be.a
            alpha_0[i] = be.alpha
            # Also extract the blade element properties required for the loads() function now that we are
            # looping over the blade elements
            blade_params.append([be.r, be.beta, be.c, self.blade.r, pitch_1, be.airfoil,
                                 v0, tsr * v0 / self.blade.r, 0, 0])

        blade_params = np.array(blade_params)

        # Same thing but for the final thrust coefficient (and thus its corresponding pitch angle
        self.blade.determine_cp_ct(v0, tsr, pitch_2)
        c_thrust_r_2 = c_thrust(self.blade.p_n_list, v0, r_list, self.blade.b, dr)

        a_f = np.empty(r_list.shape)
        alpha_f = np.empty(r_list.shape)
        for i, be in enumerate(self.blade.blade_elements):
            a_f[i] = be.a
            alpha_f[i] = be.alpha

        # Initialise the time propagation of the induction factor distribution by setting the thrust, a counter,
        # the initial time and initial induction factor distribution.
        c_thrust_r = c_thrust_r_1
        i = 0
        t = [t_0, ]
        a = [a_0, ]
        alpha = [alpha_0, ]
        # Loop until the final time step is reached
        while t[-1] <= t_final:
            # At t=0s, step from the initial to the final thrust coefficient distribution
            if t[-1] == 0:
                c_thrust_r = c_thrust_r_2
                blade_params[:, 4] = pitch_2

            # At this time step, apply pitts-peters to each blade element to determine the distribution of
            # induction factor for this time step
            a_current = np.empty(a_0.shape)
            alpha_current = np.empty(alpha_0.shape)
            for j, be in enumerate(self.blade.blade_elements):
                alpha_current[j], a_current[j], da_dt = pitt_peters(c_thrust_r[j], a[-1][j], delta_t, blade_params[j], dr, self.blade.b)

            # At this to the huge matrix of results
            a.append(a_current)
            alpha.append(alpha_current)
            # Update the counter and time
            i += 1
            t.append(round(t[-1] + delta_t, 3))

        # Plot a set of distributions with equal time spacing to show change in the whole blade
        a = np.array(a)
        alpha = np.array(alpha)
        for j, tme in enumerate(t):
            if round(tme, 0) == round(tme, 3):
                plt.plot(r_list, alpha[j], 'k')

        # Also plot the initial and final induction factor distributions
        plt.plot(r_list, alpha_0, 'r')
        plt.plot(r_list, alpha_f, 'r')
        plt.show()

        # Plot the time evolution of the induction factor for each blade element
        for j, be in enumerate(self.blade.blade_elements):
            plt.plot(t, alpha[:, j], 'k')
        plt.show()


def c_thrust(p_n, v0, r, b, dr):
    """
    Determine the local thrust coefficient based on the local loading
    :param p_n: Local thrust loading in [N/m]
    :param v0: Turbine incoming velocity
    :param r: Radial position
    :param b: Number of turbine blades
    :return: the local thrust coefficient [-]
    """
    return b * p_n / (.5 * rho * v0 ** 2 * np.pi * r * dr)


def pitt_peters(c_thrust_current, a_previous, dt, be_params, dr, b):
    """
    Calculate the new induction factor using Pitts-Peters
    :param c_thrust_current: Thrust coefficient at this time step
    :param a_previous: Induction factor at previous time step
    :param dt: Time step
    :param be_params: Parameters of the blade element required for the loads() function
    :param dr: Radial length of the blade element
    :param b: Number of turbine blades
    :return: The current time step induction factor and its derivative
    """
    # Determine the thrust loading on the blade element based on the previous time step induction factor
    p_n, _, alpha = loads(a_previous, *be_params)
    # Use the thrust loading to determine the local thrust coefficient
    c_thrust_ind = c_thrust(p_n, be_params[6], be_params[0], b, dr)

    # Calculate the time derivative of the induction factor
    da_dt = (c_thrust_current - c_thrust_ind) / (16 / (3 * np.pi)) * (
             be_params[6] ** 2 / be_params[3]) / be_params[6]
    # Calculate the new induction factor with time propagation
    a_current = a_previous - da_dt * dt
    return alpha, a_current, da_dt


def xi(a, yaw):
    # Using the approximation given in slides 2.2.2:12.
    return (0.6 * a + 1) * yaw


def loads(a, r, twist, c, r_blade, pitch, airfoil, v_0, omega, yaw, azimuth):
    """
    Determine the local loading based on geometry and induction
    :param a: local induction factor
    :param r: radial position
    :param twist: local blade twist angle in degrees
    :param c: local blade chord
    :param r_blade: blade radius
    :param pitch: global blade pitch
    :param airfoil: the used airfoil
    :param v_0: incoming velocity
    :param omega: turbine rotational speed
    :param yaw: yaw angle in degrees
    :param azimuth: azimuthal position in degrees
    :return: the thrust and power loading
    """
    # Determining skew angle of outgoing flow
    x = xi(a, yaw)

    # Using Coleman's model for vortex cylinder in yaw
    K_xi = 2 * np.tan(x / 2)

    # Using Glauert theory for yawed motion, determine separate induction factors. (slides 2.2.2:9)
    axial_induction = a * (1 + K_xi * r * np.sin(azimuth - np.pi / 2) / r_blade)
    # self.azimuthal_induction = self.a_prime

    u_tangential = (omega * r - v_0 * np.sin(yaw) * np.sin(azimuth)) # * (1 + self.a_prime)
    u_normal = v_0 * (np.cos(yaw) - axial_induction)

    # For the previous a and a_prime, find the flow angle and angle of attack
    phi = np.arctan2(u_normal, u_tangential)
    alpha = np.degrees(phi) - twist - pitch

    # With the angle of attack, determine the lift and drag coefficient from airfoil data interpolation
    cl = airfoil.cl(alpha)
    cd = airfoil.cd(alpha)

    # Use these to find the normal and tangential force coefficients
    cn = cl * np.cos(phi) + cd * np.sin(phi)
    ct = cl * np.sin(phi) - cd * np.cos(phi)

    # Determine the relative velocity with the velocity triangle
    v_rel = np.sqrt(u_normal**2 + u_tangential**2)

    # Using the previous calculations, find the forces on the blade element
    p_n = 0.5 * rho * v_rel ** 2 * c * cn
    p_t = 0.5 * rho * v_rel ** 2 * c * ct

    return p_n, p_t, alpha


def solve_a(cp):
    a_lst = np.arange(0, .33, .001)
    cp_lst = 4 * a_lst * (1 - a_lst)**2

    cp_lower = cp_lst[cp_lst < cp][-1]
    cp_upper = cp_lst[cp_lst >= cp][0]

    a_lower = a_lst[cp_lst < cp][-1]
    a_upper = a_lst[cp_lst >= cp][0]

    return interpolate(a_lower, a_upper, cp_lower, cp_upper, cp)


def interpolate(value1, value2, co1, co2, co_interpolation):
    dy_dx = (value2 - value1) / (co2 - co1)
    return dy_dx * (co_interpolation - co1) + value1


def write_to_file(array, path):
    lines = []
    for row in array:
        line = ''
        for num in row:
            line = line + f'{num},'

        lines.append(line[:-1] + '\n')

    f = open(path, 'w')
    f.writelines(lines)
    f.close()


def read_from_file(path):
    f = open(path)
    lines = f.readlines()
    out_list = [[float(num) for num in line.strip('\n').split(',')] for line in lines]
    return np.array(out_list)


if __name__ == '__main__':
    # Create the turbine with 25 blade elements
    turbine = Turbine(25)

    turbine.ct_step(.5, .9, 10, 10)
