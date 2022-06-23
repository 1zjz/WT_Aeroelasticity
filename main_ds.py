import matplotlib.pyplot as plt
import numpy as np

from dynamic_stall import generate_data, cases
from read_write import read_from_file


def read_data(case, model):
    """
    Read the generated data
    :param case: Selection of the assignment case ('Dyn1' or 'Dyn2')
    :param model: Selection of the dynamic inflow model ('pp': Pitt-Peters, 'lm': Larsen-Madsen, 'oye': Oye)
    :return: (r, t) followed by (ct, cq, a, alpha, phi) for fully unsteady, dynamic stall only, dynamic inflow only and fully
        quasi-steady respectively
    """
    return ((read_from_file(f'./{model}/{case}/r_list.csv')[0],
             read_from_file(f'./{model}/{case}/t_list.csv')[0]),

            (read_from_file(f'./{model}/{case}/ctr.csv'),
             read_from_file(f'./{model}/{case}/cqr.csv'),
             read_from_file(f'./{model}/{case}/a.csv'),
             read_from_file(f'./{model}/{case}/alpha.csv'),
             read_from_file(f'./{model}/{case}/phi.csv')),

            (read_from_file(f'./{model}/{case}/ctr_ds.csv'),
             read_from_file(f'./{model}/{case}/cqr_ds.csv'),
             read_from_file(f'./{model}/{case}/a_ds.csv'),
             read_from_file(f'./{model}/{case}/alpha_ds.csv'),
             read_from_file(f'./{model}/{case}/phi_ds.csv')),

            (read_from_file(f'./{model}/{case}/ctr_di.csv'),
             read_from_file(f'./{model}/{case}/cqr_di.csv'),
             read_from_file(f'./{model}/{case}/a_di.csv'),
             read_from_file(f'./{model}/{case}/alpha_di.csv'),
             read_from_file(f'./{model}/{case}/phi_di.csv')),

            (read_from_file(f'./{model}/{case}/ctr_qs.csv'),
             read_from_file(f'./{model}/{case}/cqr_qs.csv'),
             read_from_file(f'./{model}/{case}/a_qs.csv'),
             read_from_file(f'./{model}/{case}/alpha_qs.csv'),
             read_from_file(f'./{model}/{case}/phi_qs.csv'))
            )


def example_read():
    colors = ('r', 'g', 'b')
    for name in cases.keys():
        for i, model in enumerate(('pp', 'lm', 'oye')):
            ((r_list, t_list), (ctr, cqr, a, alpha, phi), (ctr_ds, cqr_ds, a_ds, alpha_ds, phi_ds),
             (ctr_di, cqr_di, a_di, alpha_di, phi_di), (ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs)) = read_data(name, model)

            for j in (0, 8, -2, -1):
                plt.figure(1)
                plt.plot(t_list, a[:, j], colors[i])
                plt.plot(t_list, a_ds[:, j], colors[i], linestyle='dashdot')
                plt.plot(t_list, a_di[:, j], colors[i], linestyle='dashed')
                plt.plot(t_list, a_qs[:, j], colors[i], linestyle='dotted')

                plt.figure(2)
                plt.plot(t_list, ctr[:, j], colors[i])
                plt.plot(t_list, ctr_ds[:, j], colors[i], linestyle='dashdot')
                plt.plot(t_list, ctr_di[:, j], colors[i], linestyle='dashed')
                plt.plot(t_list, ctr_qs[:, j], colors[i], linestyle='dotted')

                plt.figure(3)
                plt.plot(t_list, alpha[:, j], colors[i])
                plt.plot(t_list, alpha_ds[:, j], colors[i], linestyle='dashdot')
                plt.plot(t_list, alpha_di[:, j], colors[i], linestyle='dashed')
                plt.plot(t_list, alpha_qs[:, j], colors[i], linestyle='dotted')

                plt.figure(4)
                plt.plot(t_list, cqr[:, j], colors[i])
                plt.plot(t_list, cqr_ds[:, j], colors[i], linestyle='dashdot')
                plt.plot(t_list, cqr_di[:, j], colors[i], linestyle='dashed')
                plt.plot(t_list, cqr_qs[:, j], colors[i], linestyle='dotted')

                plt.figure(5)
                plt.plot(t_list, phi[:, j], colors[i])
                plt.plot(t_list, np.degrees(phi_ds[:, j]), colors[i], linestyle='dashdot')
                plt.plot(t_list, phi_di[:, j], colors[i], linestyle='dashed')
                plt.plot(t_list, np.degrees(phi_qs[:, j]), colors[i], linestyle='dotted')

        plt.figure(1)
        plt.xlabel('$t$ (s)')
        plt.ylabel('$a$ (-)')
        plt.tight_layout()

        plt.figure(2)
        plt.xlabel('$t$ (s)')
        plt.ylabel('$C_T$ (-)')
        plt.tight_layout()

        plt.figure(3)
        plt.xlabel('$t$ (s)')
        plt.ylabel('$\\alpha$ (deg)')
        plt.tight_layout()

        plt.figure(4)
        plt.xlabel('$t$ (s)')
        plt.ylabel('$C_Q$ (-)')
        plt.tight_layout()

        plt.figure(5)
        plt.xlabel('$t$ (s)')
        plt.ylabel('$\\phi$ (deg)')
        plt.tight_layout()

        plt.show()


if __name__ == '__main__':
    generate_data()
    example_read()
