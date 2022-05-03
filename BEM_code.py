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

    def determine_loads(self, v_0, omega, theta_p, b, r_blade, r_root, yaw, azimuth, loss=True):
        #bem code for a single blade element, loss is for tip/root loss
        yaw = np.radians(yaw)
        azimuth = np.radians(azimuth)
        # Set initial loop values
        self.a = 0
        # self.a_prime = 0
        error_a = 1
        error_a_dash = 1
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
                # raise ValueError(f"r={self.r}: Solution for a and a' not converging. a={self.a}, a' = {self.a_prime}.")
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

            # a_prime_new = 1 / ((4 * f * np.sin(self.phi) * np.cos(self.phi)) / (solidity * ct) - 1)

            # Determine the difference between this and the previous iteration
            error_a = abs(a_new - self.a)
            # error_a_dash = abs(a_prime_new - self.a_prime)

            # Get ready for the next iteration
            self.a = a_new
            # self.a_prime = a_prime_new
            i += 1

        # Determining skew angle of outgoing flow
        x = xi(self.a, yaw)

        # Using Coleman's model for vortex cylinder in yaw
        K_xi = 2 * np.tan(x / 2)

        # Using Glauert theory for yawed motion, determine separate induction factors. (slides 2.2.2:9)
        self.axial_induction = self.a * (1 + K_xi * self.r * np.sin(azimuth - np.pi / 2) / r_blade)
        # self.azimuthal_induction = self.a_prime

        self.u_tangential = (omega * self.r - v_0 * np.sin(yaw) * np.sin(azimuth)) # * (1 + self.a_prime)
        self.u_normal = v_0 * (np.cos(yaw) - self.axial_induction)

        # For the previous a and a_prime, find the flow angle and angle of attack
        self.phi = np.arctan2(self.u_normal, self.u_tangential)
        self.alpha = np.degrees(self.phi) - self.beta - theta_p

        # With the angle of attack, determine the lift and drag coefficient from airfoil data interpolation
        cl = self.airfoil.cl(self.alpha)
        cd = self.airfoil.cd(self.alpha)

        # Use these to find the normal and tangential force coefficients
        cn = cl * np.cos(self.phi) + cd * np.sin(self.phi)
        ct = cl * np.sin(self.phi) - cd * np.cos(self.phi)

        # Determine the relative velocity with the velocity triangle
        v_rel = np.sqrt(self.u_normal**2 + self.u_tangential**2)

        # Using the previous calculations, find the forces on the blade element
        self.p_n = 0.5 * rho * v_rel ** 2 * self.c * cn
        self.p_t = 0.5 * rho * v_rel ** 2 * self.c * ct

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

    def find_pn_pt(self, v_0, theta_p, omega, yaw, azimuth, loss=True):
        # Initialise the lists for p_n and p_t
        p_n_list, p_t_list = list(), list()
        for blade_element in self.blade_elements:
            if self.r_list[0] < blade_element.r < self.r:
                blade_element.determine_loads(v_0, omega, theta_p, self.b, self.r, self.r_list[0], yaw, azimuth, loss=loss)
                p_n, p_t = blade_element.get_loads()

                p_n_list.append(p_n)
                p_t_list.append(p_t)

            else:
                # Add zero load at the blade tip and root
                p_n_list.append(0)
                p_t_list.append(0)

        return np.array(p_n_list), np.array(p_t_list), self.r_list

    def determine_cp_ct(self, v_0, lamda, theta_p, yaw, azimuth, loss=True):
        # Determine the rotational speed of the turbine
        omega = lamda * v_0 / self.r
        # Get the loads on the blade elements
        self.p_n_list, self.p_t_list, r_list = self.find_pn_pt(v_0, theta_p, omega, yaw, azimuth, loss=loss)

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


def xi(a, yaw):
    # Using the approximation given in slides 2.2.2:12.
    return (0.6 * a + 1) * yaw


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


def read_from_file(path):
    f = open(path)
    lines = f.readlines()
    out_list = [[float(num) for num in line.strip('\n').split(',')] for line in lines]
    return np.array(out_list)


if __name__ == '__main__':
    turbine = Turbine(50)
#turbine for 50 blade elements

