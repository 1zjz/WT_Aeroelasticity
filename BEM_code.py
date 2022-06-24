import numpy as np
import scipy.integrate as spig
import matplotlib.pyplot as plt

from read_write import read_from_file


relaxation = 1
rho = 1.225
p_atm = 101325


class DU95W150:
    def __init__(self):
        data = read_from_file('DU95W150.csv')
        self.alpha_lst = data[:, 0]
        self.cl_lst = data[:, 1]
        self.cd_lst = data[:, 2]
        self.cm_lst = data[:, 3]

    def cl(self, alpha, **kwargs): return np.interp(alpha, self.alpha_lst, self.cl_lst)

    def cd(self, alpha): return np.interp(alpha, self.alpha_lst, self.cd_lst)

    def plot_polars(self):
        fig, axes = plt.subplots(1, 3, figsize=(9, 3.5))
        axes[0].plot(self.alpha_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)], 'k')
        axes[1].plot(self.cd_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)], 'k')
        axes[2].plot(self.alpha_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)] /
                     self.cd_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)], 'k')

        optimal = np.argmax(self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)] /
                            self.cd_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)])

        axes[0].plot(self.alpha_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal], 'ro')
        axes[1].plot(self.cd_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal], 'ro')
        axes[2].plot(self.alpha_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal],
                     self.cl_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal] /
                     self.cd_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal], 'ro')

        axes[0].set_xlabel('$\\alpha$ ($^{\\circ}$)')
        axes[2].set_xlabel('$\\alpha$ ($^{\\circ}$)')
        axes[1].set_xlabel('$C_d$ (-)')

        axes[0].set_ylabel('$C_l$ ($^{\\circ}$)')
        # axes[1].set_ylabel('$C_l$ ($^{\\circ}$)')
        axes[2].set_ylabel('$C_l/C_d$ (-)')
        fig.set_tight_layout(True)
        plt.savefig('airfoil_polars.pdf')
        plt.show()
        print(self.alpha_lst[np.logical_and(self.alpha_lst >= -6, self.alpha_lst <= 12)][optimal])


class BladeElement:
    def __init__(self, pos_r: float, chord: float, twist: float, airfoil, r_blade):
        # One of blade elements
        # Fixed values
        self.r = pos_r
        self.c = chord
        self.twist = twist
        self.r_blade = r_blade
        # Values to be determined with other functions
        self.a = None
        self.axial_induction = None
        self.azimuthal_induction = None
        self.phi = None
        self.alpha = None
        self.p_n = None
        self.p_t = None
        self.u_tangential = None
        self.u_normal = None

        self.airfoil = airfoil

    def __repr__(self): # print blade element
        return f"<Blade Element at r={self.r}, c={self.c}, beta={self.twist}>"

    def determine_loads(self, v_0, omega, pitch, b, r_blade, r_root, yaw=0, azimuth=0, loss=True):
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
            self.alpha = np.degrees(self.phi) - self.twist - pitch

            # With the angle of attack, determine the lift and drag coefficient from airfoil data interpolation
            cl = self.airfoil.cl(self.alpha)
            cd = self.airfoil.cd(self.alpha)

            # Use these to find the normal and tangential force coefficients
            cn = cl * np.cos(self.phi) + cd * np.sin(self.phi)
            ct = cl * np.sin(self.phi) - cd * np.cos(self.phi)

            # Break conditions for the a-loop
            if error_a <= 1e-3: # 1e-9:
                break
            elif i > 1e2:
                raise ValueError(f"r={self.r}: Solution for a not converging. a={self.a}. Last delta: {error_a}")

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

        self.p_n, self.p_t, _, _ = loads(self.a, self.r, self.twist, self.c, r_blade, pitch, self.airfoil,
                                         v_0, omega, yaw, azimuth, propagate=True)

    def get_loads(self):
        if self.p_t is None or self.p_n is None:
            raise ValueError(f"Loads have not been determined. Run .determine_loads() first.")
        else:
            return self.p_n, self.p_t

    def reset(self):
        self.__init__(self.r, self.c, self.twist, self.airfoil, self.r_blade)


class Blade:
    def __init__(self, n_blades, airfoil, r_start, r_end, blade_pitch, n_elements, create_elements=True):
        self.b = n_blades
        self.r = r_end

        self.power = None
        self.thrust = None
        self.c_power = None
        self.c_thrust = None
        self.p_n_list = None
        self.p_t_list = None

        if create_elements:
            self.r_list = []
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

                self.blade_elements.append(BladeElement(r, chord, relative_pitch, airfoil(), self.r))

            self.r_list = np.array(self.r_list)

    def find_pn_pt(self, v_0, pitch, omega, yaw=0, azimuth=0, loss=True):
        # Initialise the lists for p_n and p_t
        p_n_list, p_t_list = list(), list()
        for blade_element in self.blade_elements:
            if self.r_list[0] < blade_element.r < self.r:
                blade_element.determine_loads(v_0, omega, pitch, self.b, self.r, self.r_list[0], yaw, azimuth, loss)
                p_n, p_t = blade_element.get_loads()

                p_n_list.append(p_n)
                p_t_list.append(p_t)

            else:
                # Add zero load at the blade tip and root
                p_n_list.append(0)
                p_t_list.append(0)

        return np.array(p_n_list), np.array(p_t_list), self.r_list

    def determine_cp_ct(self, v_0, tsr, pitch, yaw=0, azimuth=0, loss=True):
        """
        Let the BEM code Determine the power and thrust coefficient of the turbine
        :param v_0: Incoming velocity
        :param tsr: Tip speed ratio
        :param pitch: pitch angle
        :param yaw: yaw angle
        :param azimuth: azimuthal position in the turbine disk
        :param loss: turn on or off the tip loss factor
        :return: None
        """
        # Determine the rotational speed of the turbine
        omega = tsr * v_0 / self.r
        # Get the loads on the blade elements
        self.p_n_list, self.p_t_list, r_list = self.find_pn_pt(v_0, pitch, omega, yaw, azimuth, loss)

        # Determine the thrust and power of the turbine
        self.thrust = self.b * spig.trapz(self.p_n_list, self.r_list)
        self.power = omega * self.b * spig.trapz(self.p_t_list * self.r_list, self.r_list)

        # Determine the thrust and power coefficient
        self.c_thrust = self.thrust / (0.5 * rho * np.pi * self.r**2 * v_0**2)
        self.c_power = self.power / (0.5 * rho * np.pi * self.r**2 * v_0**3)

    def reset(self):
        for be in self.blade_elements:
            be.reset()


def xi(a, yaw):
    # Using the approximation given in slides 2.2.2:12.
    return (0.6 * a + 1) * yaw


def loads(a, r, twist, c, r_blade, pitch, airfoil, v_0, omega, yaw, azimuth, propagate=False):
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
    :param propagate: indicate whether to do the time propagation of the dynamic stall model
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
    cl = airfoil.cl(alpha, propagate=propagate)
    cd = airfoil.cd(alpha)

    # Use these to find the normal and tangential force coefficients
    cn = cl * np.cos(phi) + cd * np.sin(phi)
    ct = cl * np.sin(phi) - cd * np.cos(phi)

    # Determine the relative velocity with the velocity triangle
    v_rel = np.sqrt(u_normal**2 + u_tangential**2)

    # Using the previous calculations, find the forces on the blade element
    p_n = 0.5 * rho * v_rel ** 2 * c * cn
    p_t = 0.5 * rho * v_rel ** 2 * c * ct

    return p_n, p_t, alpha, np.degrees(phi)


def c_thrust(p_n, v0, r, b, dr):
    """
    Determine the local thrust coefficient based on the local loading
    :param p_n: Local thrust loading in [N/m]
    :param v0: Turbine incoming velocity
    :param r: Radial position
    :param b: Number of turbine blades
    :param dr: Blade element length
    :return: the local thrust coefficient [-]
    """
    return b * p_n / (.5 * rho * v0 ** 2 * np.pi * r * dr)


def interpolate(value1, value2, co1, co2, co_interpolation):
    """
    Interpolate linearly between two points
    :param value1: f(x1)
    :param value2: f(x2)
    :param co1: x1
    :param co2: x2
    :param co_interpolation: x
    :return: f(x)
    """
    df_dx = (value2 - value1) / (co2 - co1)
    return df_dx * (co_interpolation - co1) + value1


if __name__ == '__main__':
    raise RuntimeError("Do not run this file, it has no use.")
