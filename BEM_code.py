import numpy as np
import scipy.integrate as spig
import matplotlib.pyplot as plt
import time


relaxation = 0.1
rho = 1.225
p_atm = 101325 #pa


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


class BladeElement:
    def __init__(self, pos_r: float, chord: float, relative_pitch: float, airfoil):
        self.r = pos_r
        self.c = chord
        self.beta = relative_pitch

        self.a = None
        self.a_prime = None
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

    def __repr__(self):
        return f"<Blade Element at r={self.r}, c={self.c}, beta={self.beta}>"

    def determine_loads(self, v_0, omega, theta_p, b, r_blade, r_root, yaw, azimuth, loss=True):
        yaw = np.radians(yaw)
        azimuth = np.radians(azimuth)
        # Set initial loop values
        self.a = 0
        self.a_prime = 0
        error_a = 1
        error_a_dash = 1
        i = 0
        # Iterative solver for a and a_prime until the difference between the iterations becomes very small
        while True:
            self.u_tangential = omega * self.r * (1 + self.a_prime)
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
            if error_a <= 1e-9 and error_a_dash <= 1e-9:
                break
            elif i > 1e3:
                raise ValueError(f"r={self.r}: Solution for a and a' not converging. a={self.a}, a' = {self.a_prime}.")

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

            a_prime_new = 1 / ((4 * f * np.sin(self.phi) * np.cos(self.phi)) / (solidity * ct) - 1)

            # Determine the difference between this and the previous iteration
            error_a = abs(a_new - self.a)
            error_a_dash = abs(a_prime_new - self.a_prime)

            # Get ready for the next iteration
            self.a = a_new
            self.a_prime = a_prime_new
            i += 1

        # Determining skew angle of outgoing flow
        x = xi(self.a, yaw)

        # Using Coleman's model for vortex cylinder in yaw
        K_xi = 2 * np.tan(x / 2)

        # Using Glauert theory for yawed motion, determine separate induction factors. (slides 2.2.2:9)
        self.axial_induction = self.a * (1 + K_xi * self.r * np.sin(azimuth - np.pi / 2) / r_blade)
        self.azimuthal_induction = self.a_prime

        self.u_tangential = (omega * self.r - v_0 * np.sin(yaw) * np.sin(azimuth)) * (1 + self.a_prime)
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

    def get_static_pressure_before_rotor(self):
        # Using bernoulli
        return p_atm - 0.5 * rho * (self.u_normal**2 * (1-(1-self.a)**2))

    def get_static_pressure_after_rotor(self):
        u02 = self.u_normal**2
        return p_atm - 0.5 * rho * ( (u02*(1-self.a)**2 + u02*(1-2*self.a)**2) )

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

    def cp_lamda(self):
        tsr = np.round(np.arange(4, 12.1, 0.1), 1)
        cp = np.zeros(tsr.shape)
        thrust = np.zeros(tsr.shape)
        torque = np.zeros(tsr.shape)

        for i, lamda in enumerate(tsr):
            self.blade.determine_cp_ct(10, lamda, 0, 0, 0)
            cp[i] = self.blade.c_power
            thrust[i] = self.blade.thrust
            torque[i] = self.blade.power / (lamda * 10 / self.blade.r)

            if lamda in (6, 8, 10):
                plt.plot(lamda, self.blade.c_power, 'k^')

            self.blade.reset()

        plt.xlabel("$\\lambda\\ [-]$")
        plt.ylabel("$C_P\\ [-]$")
        plt.tight_layout()
        plt.plot(tsr, cp, 'k')
        plt.savefig('cp-lambda.pdf')
        plt.grid()
        plt.show()

        fig, ax1 = plt.subplots()

        ax1.plot(tsr, thrust / 1e3, color='tab:blue')
        ax1.set_xlabel('$\\lambda$ [-]')
        ax1.set_ylabel('$T$ [kN]', color='tab:blue')
        ax1.tick_params(axis='y', colors='tab:blue')
        ax1.set_yticks(np.linspace(0, 450, 10))

        ax2 = ax1.twinx()
        ax2.plot(tsr, torque / 1e3, color='tab:orange')
        ax2.set_ylabel('$Q$ [kNm]', color='tab:orange')
        ax2.tick_params(axis='y', colors='tab:orange')
        ax2.set_yticks(np.linspace(600, 1500, 10))

        plt.tight_layout()
        ax1.grid()
        ax2.grid()
        plt.savefig('cp-thrust-torque.pdf')

        plt.show()

    def spanwise_distributions(self):
        linestyles = ('dashed', 'solid', 'dotted')
        for j, tsr in enumerate((6, 8, 10)):
            self.blade.determine_cp_ct(10, tsr, 0, 0, 0)
            pn, pt = self.blade.p_n_list, self.blade.p_t_list
            alpha, phi, a, a_prime, twist = np.zeros((5, len(self.blade.blade_elements)))
            for i, be in enumerate(self.blade.blade_elements):
                alpha[i] = be.alpha
                phi[i] = be.phi
                a[i] = be.a
                a_prime[i] = be.a_prime
                twist[i] = be.beta

            plt.figure(1, figsize=(5, 5))
            plt.plot(self.blade.r_list, alpha, linestyle=linestyles[j], color='tab:blue', label=f'Angle of Attack ($\\alpha$) ($\\lambda={tsr}$)')
            plt.plot(self.blade.r_list, np.degrees(phi), linestyle=linestyles[j], color='tab:orange', label=f'Inflow Angle ($\\phi$) ($\\lambda={tsr}$)')
            # plt.plot(self.blade.r_list, twist, label=f'Twist Angle ($\\beta$)')
            plt.xlabel('$r$ [m]')
            plt.ylabel('$Angle$ [$^{\\circ}$]')

            plt.figure(2, figsize=(5, 5))
            plt.plot(self.blade.r_list, a, linestyle=linestyles[j], color='tab:blue', label=f'Axial Induction ($a$) ($\\lambda={tsr}$)')
            plt.plot(self.blade.r_list, a_prime, linestyle=linestyles[j], color='tab:orange', label=f"Azimuthal Induction ($a'$) ($\\lambda={tsr}$)")
            plt.xlabel('$r$ [m]')
            plt.ylabel('$Induction\\ factor$ [-]')
            plt.yticks(np.arange(0, 0.6, 0.1))

            plt.figure(3, figsize=(5, 5))
            plt.plot(self.blade.r_list, pn, linestyle=linestyles[j], color='tab:blue', label=f'Thrust Loading ($p_n$) ($\\lambda={tsr}$)')
            plt.plot(self.blade.r_list, pt, linestyle=linestyles[j], color='tab:orange', label=f'Azimuthal Loading ($p_t$) ($\\lambda={tsr}$)')
            plt.xlabel('$r$ [m]')
            plt.ylabel('$p$ [N/m]')

            plt.figure(4, figsize=(5, 5))
            # plt.plot(self.blade.r_list, pt * self.blade.r_list, linestyle=linestyles[j], color='tab:green', label=f'Torque Loading ($p_q$) ($\\lambda={tsr}$)')
            plt.plot(self.blade.r_list, 4 * a * (1-a), linestyle=linestyles[j], color='tab:blue', label=f'Thrust Loading ($c_n$) ($\\lambda={tsr}$)')
            plt.plot(self.blade.r_list, 4 * a_prime * (1-a) * tsr, linestyle=linestyles[j], color='tab:orange', label=f'Azimuthal Loading ($c_t$) ($\\lambda={tsr}$)')
            plt.plot(self.blade.r_list, 4 * a_prime * (1-a) * tsr * self.blade.r_list / self.blade.r, linestyle=linestyles[j], color='tab:green', label=f'Torque Loading ($c_q$) ($\\lambda={tsr}$)')
            plt.xlabel('$r$ [m]')
            plt.ylabel('Coefficient [-]')

        plt.figure(1, figsize=(5, 5))
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig('Angles_no_yaw.pdf')

        plt.figure(2, figsize=(5, 5))
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig('Inductions_no_yaw.pdf')

        plt.figure(3, figsize=(5, 5))
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig('Forces_no_yaw.pdf')

        plt.figure(4, figsize=(5, 5))
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig('Coefficients_no_yaw.pdf')

        plt.show()

    def yaw_polar_plots(self):
        yaw = (0, 15, 30,)
        azimuth = np.arange(0, 360 + 1, 1)
        r_grid, az_grid = np.meshgrid(self.blade.r_list[1:-1], np.radians(azimuth))
        fig1, (ax1, ax2) = create_axes(1)
        fig2, (ax3, ax4) = create_axes(2)
        fig3, (ax5, ax6) = create_axes(3)
        fig4, (ax7, ax8) = create_axes(4)

        for i, y in enumerate(yaw):
            alpha, phi, a, a_prime, pn, pt, un, ut = np.zeros((8, azimuth.size, len(self.blade.blade_elements[1:-1])))
            for j, az in enumerate(azimuth):
                # print(az)
                self.blade.determine_cp_ct(10, 8, 0, y, az)

                pn[j, :], pt[j, :] = self.blade.p_n_list[1:-1], self.blade.p_t_list[1:-1]
                for k, be in enumerate(self.blade.blade_elements[1:-1]):
                    alpha[j, k] = be.alpha
                    phi[j, k] = be.phi
                    a[j, k] = be.axial_induction
                    a_prime[j, k] = be.azimuthal_induction
                    un[j, k] = be.u_normal
                    ut[j, k] = be.u_tangential

                self.blade.reset()

            color = 'k'
            cmap = None
            linewidth = .5

            contour_plot(az_grid, r_grid, alpha, ax1, fig1,
                         (i, y, cmap, color, linewidth, 0, 'Angle of Attack $\\alpha$ [$^{\\circ}$]'))
            contour_plot(az_grid, r_grid, np.degrees(phi), ax2, fig1,
                         (i, y, cmap, color, linewidth, 1, 'Inflow Angle $\\phi$ [$^{\\circ}$]'))

            contour_plot(az_grid, r_grid, a, ax3, fig2,
                         (i, y, cmap, color, linewidth, 2, 'Axial Induction $a_{total}$ [-]'))
            contour_plot(az_grid, r_grid, a_prime, ax4, fig2,
                         (i, y, cmap, color, linewidth, 3, "Azimuthal Induction $a'$ [-]"))

            contour_plot(az_grid, r_grid, pn, ax5, fig3,
                         (i, y, cmap, color, linewidth, 4, 'Normal Force $p_n$ [N/m]'))
            contour_plot(az_grid, r_grid, pt, ax6, fig3,
                         (i, y, cmap, color, linewidth, 5, 'Tangential Force $p_t$ [N/m]'))

            contour_plot(az_grid, r_grid, un, ax7, fig4,
                         (i, y, cmap, color, linewidth, 6, 'Normal Velocity $u_n$ [m/s]'))
            contour_plot(az_grid, r_grid, ut, ax8, fig4,
                         (i, y, cmap, color, linewidth, 7, 'Tangential Velocity $u_t$ [m/s]'))

            # fig1.tight_layout()
            # fig2.tight_layout()
            # fig3.tight_layout()
            # fig4.tight_layout()

        fig1.savefig('Angles.png')
        fig2.savefig('Inductions.png')
        fig3.savefig('Forces.png')
        fig4.savefig('Velocities.png')
        plt.show()

    def loss_comparison(self):
        linestyles = ('solid', 'dashed')
        for j in range(2):
            self.blade.determine_cp_ct(10, 8, 0, 0, 0, loss=bool(j))

            alpha, phi, a, a_prime, twist = np.zeros((5, len(self.blade.blade_elements[1:-1])))
            for i, be in enumerate(self.blade.blade_elements[1:-1]):
                a[i] = be.a
                a_prime[i] = be.a_prime

            if j:
                txt = '(with tip correction)'
            else:
                txt = '(without tip correction)'

            plt.figure(1)
            plt.plot(self.blade.r_list[1:-1], a, linestyle=linestyles[j], color='tab:blue', label=f'Axial Induction ($a$) {txt}')
            plt.plot(self.blade.r_list[1:-1], a_prime, linestyle=linestyles[j], color='tab:orange', label=f"Azimuthal Induction ($a'$) {txt}")
            plt.figure(2)
            plt.plot(self.blade.r_list[1:-1], self.blade.p_n_list[1:-1], linestyle=linestyles[j], color='tab:blue', label=f'Thrust Loading ($p_n$) {txt}')
            plt.plot(self.blade.r_list[1:-1], self.blade.p_t_list[1:-1], linestyle=linestyles[j], color='tab:orange', label=f'Azimuthal Loading ($p_t$) {txt}')
            self.blade.reset()

        plt.figure(1)
        plt.ylim(0, .5)
        plt.grid()
        plt.xlabel('$r$ [m]')
        plt.ylabel('$Induction\\ factor$ [-]')
        plt.legend()
        plt.tight_layout()
        plt.savefig('Tip_loss_induction.pdf')

        plt.figure(2)
        plt.ylim(0, 4500)
        plt.grid()
        plt.xlabel('$r$ [m]')
        plt.ylabel('$p$ [N/m]')
        plt.legend()
        plt.tight_layout()
        plt.savefig('Tip_loss_forces.pdf')

        plt.show()

    def enthalpy_distributions(self, v_0):
        self.blade.determine_cp_ct(10, 8, 0, 0, 0)
        # row = station, column is azimuthal pos
        enthalpies = np.zeros((4, len(self.blade.blade_elements)))

        for i, be in enumerate(self.blade.blade_elements):
            # Ignore first and last one, not set because tip & root loss factors
            if (i == 0 or i == len(self.blade.blade_elements)-1):
                continue;

            # At the end and the start, the static pressure is just the atmospheric pressure.
            enthalpies[0,i] = p_atm/rho + 0.5*v_0**2
            pressure_before_rotor = be.get_static_pressure_before_rotor()
            enthalpies[1,i] = pressure_before_rotor/rho + 0.5*v_0**2 * (1-be.a)**2
            pressure_after_rotor = be.get_static_pressure_after_rotor()
            enthalpies[2,i] = pressure_after_rotor/rho + 0.5*v_0**2 * (-(1-be.a)**2 + (1-2*be.a)**2)
            enthalpies[3,i] = p_atm/rho + 0.5 * v_0**2*(1-2*be.a)**2

        # Skipping the first & last index, so remove those columns
        enthalpies = np.delete(enthalpies, [0, -1], 1)
        r_sliced = np.delete(self.blade.r_list, [0, -1])

        # Do some nice plotssss
        plt.figure(1)
        plt.plot(r_sliced, enthalpies[0], label='upwind', color='b')
        plt.plot(r_sliced, enthalpies[1], label='upwind rotor', color='tab:orange')
        plt.plot(r_sliced, enthalpies[2], label='downwind rotor', color='tab:brown')
        plt.plot(r_sliced, enthalpies[3], label='downwind', color='c')
        plt.xlabel('$r$ [m]')
        plt.ylabel('relative Specific enthalpy')
        plt.grid()
        plt.legend()
        plt.show()

    def enthalpy_trial(self):
        self.blade.determine_cp_ct(10, 8, 0, 0, 0)

        a_lst = np.zeros(len(self.blade.blade_elements[1:-1]))
        for i, be in enumerate(self.blade.blade_elements[1:-1]):
            a_lst[i] = be.a

        a = spig.trapz(a_lst, self.blade.r_list[1:-1]) / (self.blade.r_list[1:-1][-1] - self.blade.r_list[1:-1][0])

        p01 = p_atm + .5 * 1.225 * 10 ** 2
        p2 = p01 - .5 * 1.225 * (10 * (1 - a)) ** 2
        p04 = p_atm + .5 * 1.225 * (10 * (1 - 2 * a)) ** 2
        p3 = p04 - .5 * 1.225 * (10 * (1 - a)) ** 2

        hs1 = (p_atm / 1.225 + .5 * 10 ** 2) * np.ones(len(self.blade.blade_elements[1:-1]))
        hs2 = p2 / 1.225 + .5 * (10 * (1 - a_lst)) ** 2
        hs3 = p3 / 1.225 + .5 * (10 * (1 - a_lst)) ** 2
        hs4 = p_atm / 1.225 + .5 * (10 * (1 - 2 * a_lst)) ** 2

        plt.plot(self.blade.r_list[1:-1], hs1, label='Infinity Upwind')
        plt.plot(self.blade.r_list[1:-1], hs2, label='Upwind Rotor')
        plt.plot(self.blade.r_list[1:-1], hs3, label='Downwind Rotor')
        plt.plot(self.blade.r_list[1:-1], hs4, label='Infinity Downwind')
        plt.legend()
        plt.grid()
        plt.show()


def create_axes(num):
    fig, axes = plt.subplots(3, 2, num=num, subplot_kw=dict(projection='polar'), figsize=(9, 12), sharey='all')
    ax1, ax2 = axes.T
    fig.subplots_adjust(0.07, 0.01, .98, .97, .12, .1)

    return fig, (ax1, ax2)


def contour_plot(az_grid, r_grid, values, axes, figure, options: tuple):
    def contour_levels():
        base_levels = np.linspace(0, 1, 1000)

        minima = (2, 2, 0, 0, 500, 50, 3, 10)
        maxima = (12, 20, .6, .06, 4500, 450, 8, 90)
        nticks = (11, 10, 11, 11, 9, 9, 11, 9)

        return (base_levels * (maxima[qtt] - minima[qtt]) + minima[qtt],
                np.linspace(minima[qtt], maxima[qtt], nticks[qtt]))

    idx, theta, cmap, color, linewidth, qtt, label, *_ = options

    axis = axes[idx]

    levels, ticks = contour_levels()

    ctr = axis.contourf(az_grid, r_grid, values, levels, cmap=cmap)
    ctr1 = axis.contour(az_grid, r_grid, values, ticks, colors=color, linewidths=linewidth)
    axis.set_yticks(np.arange(0, 50 + 10, 10))
    cbar = figure.colorbar(ctr, ax=axis, ticks=ticks)

    if not idx:
        axis.set_title(label)
    if not qtt % 2:
        axis.set_ylabel(f'$\\theta = {theta}' + '^{\\circ}$')
        axis.yaxis.set_label_coords(-0.13, 0.5)


def xi(a, yaw):
    # Using the approximation given in slides 2.2.2:12.
    return (0.6 * a + 1) * yaw
#     val = yaw.copy()
#     diff = 1
#     c = 0
#     relax = 0.1
#     while diff > 1e-3 and c < 1e2:
#         # new = np.arctan2(2 * np.tan(val / 2), 1 - np.tan(val / 2)**2)
#         new = relax * np.arctan2(np.sin(yaw) - a * np.tan(val/2), np.cos(yaw) - a) + (1 - relax) * val
#         diff = abs(new - val)
#         val = new
#         c += 1
#
#     # print(c)
#     if c < 1e3:
#         return val
#     else:
#         print(f'Not converged for yaw={yaw}, a={a}.')
#         return yaw


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


def convergence():
    thrust = []
    annuli = np.arange(5, 1000 + 5, 5)
    times = []

    for n in annuli:
        blade = Blade(3, DU95W150, .2 * 50, 50, -2, n)
        t0 = time.time()
        blade.determine_cp_ct(10, 8, 0, 0, 0)
        times.append(time.time() - t0)
        thrust.append(blade.c_thrust)

    plt.figure(1)
    plt.hlines(max(thrust), 0, annuli[-1], linestyles='dotted')
    plt.plot(annuli, thrust)
    plt.xlabel('$N$ [-]')
    plt.ylabel('$C_T$ [-]')
    plt.grid()
    plt.savefig('Convergence_History.pdf')

    plt.figure(2)
    plt.plot(annuli, times)

    plt.figure(3)
    plt.plot(annuli, (max(thrust) - np.array(thrust)) / np.array(times))
    plt.show()


def airfoil_polars():
    airfoil = DU95W150()
    fig, axes = plt.subplots(1, 2, sharey='all')
    airfoil.plot_polars(axes)
    blade = Blade(3, DU95W150, .2 * 50, 50, -2, 50)
    blade.determine_cp_ct(10, 8, 0, 0, 0)

    alpha = []
    chord = []
    radius = []
    for be in blade.blade_elements[1:-1]:
        alpha.append(be.alpha)
        chord.append(be.c)
        radius.append(be.r)

    alpha_mm = [min(alpha), max(alpha)]
    axes[0].plot(alpha_mm, [airfoil.cl(a) for a in alpha_mm], 'k^')
    axes[0].set_xlabel('$\\alpha$ [$^{\\circ}$]')
    axes[0].set_ylabel('$C_l$ [-]')
    axes[1].plot([airfoil.cd(a) for a in alpha_mm], [airfoil.cl(a) for a in alpha_mm], 'k^')
    axes[1].set_xlabel('$C_d$ [-]')

    plt.tight_layout()
    plt.savefig('polars.pdf')
    axes[0].grid()
    axes[1].grid()
    plt.show()

    fig, axes = plt.subplots(2, 1, sharex='all')
    axes[0].plot(radius, [airfoil.cl(a) for a in alpha])
    axes[1].plot(radius, chord)
    plt.show()


if __name__ == '__main__':
    # convergence()
    # airfoil_polars()

    turbine = Turbine(50)
    # turbine.cp_lamda()
    turbine.spanwise_distributions()
    # turbine.yaw_polar_plots()
    # turbine.loss_comparison()
    # turbine.enthalpy_trial()

    # a = .82
    # yaw = np.radians(30)
    # print(np.degrees((.6 * a + 1) * yaw))
    # print(np.degrees(xi(a, yaw)))
