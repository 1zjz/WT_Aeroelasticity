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

def plot_combined_subplot(y_label, set_1_tag, set_2_tag, ax1, ax2, ax3, ax4, x_lst, y_1_mat, y_2_mat, blade_loc_id, blade_loc_tag, color, line, i):
    if set_2_tag == 'Quasi-steady':
        ax1.plot(x_lst, y_1_mat[:, blade_loc_id[0]], color=color, linestyle=line, label=set_1_tag)
        ax1.set_title(blade_loc_tag[0])
        ax1.set_ylabel(y_label)
        ax1.grid()
        ax2.plot(x_lst, y_1_mat[:, blade_loc_id[1]], color=color, linestyle=line)
        ax2.set_title(blade_loc_tag[1])
        ax2.set_ylabel(y_label)
        ax2.grid()
        ax3.plot(x_lst, y_1_mat[:, blade_loc_id[2]], color=color, linestyle=line)
        ax3.set_title(blade_loc_tag[2])
        ax3.set_ylabel(y_label)
        ax3.grid()
        ax4.plot(x_lst, y_1_mat[:, blade_loc_id[3]], color=color, linestyle=line)
        ax4.set_title(blade_loc_tag[3])
        ax4.set_ylabel(y_label)
        ax4.set_xlabel('Time [s]')
        ax4.grid()
        if not i:
            ax1.plot(x_lst, y_2_mat[:, blade_loc_id[0]], color=qs_color, linestyle='solid', label=set_2_tag)
            ax2.plot(x_lst, y_2_mat[:, blade_loc_id[1]], color=qs_color, linestyle='solid')
            ax3.plot(x_lst, y_2_mat[:, blade_loc_id[2]], color=qs_color, linestyle='solid')
            ax4.plot(x_lst, y_2_mat[:, blade_loc_id[3]], color=qs_color, linestyle='solid')
        
    else:
        ax1.plot(x_lst, y_1_mat[:, blade_loc_id[0]], color=color, linestyle=line, label=set_1_tag)
        ax1.plot(x_lst, y_2_mat[:, blade_loc_id[0]], color=color, linestyle='solid', label=set_2_tag)
        ax1.set_title(blade_loc_tag[0])
        ax1.set_ylabel(y_label)
        ax1.grid()
        ax2.plot(x_lst, y_1_mat[:, blade_loc_id[1]], color=color, linestyle=line)
        ax2.plot(x_lst, y_2_mat[:, blade_loc_id[1]], color=color, linestyle='solid')
        ax2.set_title(blade_loc_tag[1])
        ax2.set_ylabel(y_label)
        ax2.grid()
        ax3.plot(x_lst, y_1_mat[:, blade_loc_id[2]], color=color, linestyle=line)
        ax3.plot(x_lst, y_2_mat[:, blade_loc_id[2]], color=color, linestyle='solid')
        ax3.set_title(blade_loc_tag[2])
        ax3.set_ylabel(y_label)
        ax3.grid()
        ax4.plot(x_lst, y_1_mat[:, blade_loc_id[3]], color=color, linestyle=line)
        ax4.plot(x_lst, y_2_mat[:, blade_loc_id[3]], color=color, linestyle='solid')
        ax4.set_title(blade_loc_tag[3])
        ax4.set_ylabel(y_label)
        ax4.set_xlabel('Time [s]')
        ax4.grid()
    return
        

def plot_combined_subplot_elem(y_label,ax1,ax2,ax3,x_lst,y_mat,y_qs_mat,row_interest,model_tags,color,line,line_label,qs_color,model_i,time_step_counter):
    if model_i == 0 and time_step_counter == 0:
        # Plot the initial state of quasi-steady solution
        ax1.plot(x_lst, y_qs_mat[1, :], color='#069AF3', linestyle='solid', label='Quasi-steady initial')
        ax2.plot(x_lst, y_qs_mat[1, :], color='#069AF3', linestyle='solid')
        ax3.plot(x_lst, y_qs_mat[1, :], color='#069AF3', linestyle='solid')
        # Plot the final state of quasi-steady solution
        ax1.plot(x_lst, y_qs_mat[-1, :], color='#F97306', linestyle='solid', label='Quasi-steady final')
        ax2.plot(x_lst, y_qs_mat[-1, :], color='#F97306', linestyle='solid')
        ax3.plot(x_lst, y_qs_mat[-1, :], color='#F97306', linestyle='solid')
        ax1.grid()
        ax2.grid()
        ax3.grid()
    if model_i == 0:
        ax1.plot(x_lst, y_mat[row_interest, :], color=color, linestyle=line, label=line_label)
        ax1.set_title(model_tags[model_i])
        ax1.set_ylabel(y_label)
        ax1.set_xlabel('Blade radial position [m]')
    elif model_i == 1:
        ax2.plot(x_lst, y_mat[row_interest, :], color=color, linestyle=line)
        ax2.set_title(model_tags[model_i])
        ax2.set_xlabel('Blade radial position [m]')
    else:
        ax3.plot(x_lst, y_mat[row_interest, :], color=color, linestyle=line)
        ax3.set_title(model_tags[model_i])
        ax3.set_xlabel('Blade radial position [m]')
    return

def plot_save_figure(fig_tag, case_tag, response_tag, comparison_tag, folder_name):
    fig_tag.tight_layout()
    fig_tag.subplots_adjust(bottom=0.35)
    fig_tag.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.21), ncol=2)
    fig_name = case_tag + '_' + response_tag + '_' + comparison_tag + '_time.pdf'
    fig_tag.savefig(folder_name + '\\' + fig_name)
    return

def plot_save_figure_elem(fig_tag, case_tag, response_tag, folder_name):
    fig_tag.tight_layout()
    fig_tag.subplots_adjust(bottom=0.2)
    fig_tag.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.1), ncol=4)
    fig_name = case_tag + '_' + response_tag + '_blade_elem.pdf'
    fig_tag.savefig(folder_name + '\\' + fig_name)


if __name__ == '__main__':
    # To generate the data, uncomment the following two lines
    # generate_data()
    # example_read()
    
    # Define the range of case tags considered
    case_tag_range = ('Dyn1', 'Dyn2')
    
    # Define blade locations of interest and plotting styles
    blade_loc_id = (0, 8, -2, -1)                               # Index used correspond to radial distribution indexing
    blade_loc_tag = ('0.2 R', '0.5 R', '0.9 R', '1.0 R')        # Labelling follows the ordering defined in blade_loc_id
    blade_loc_line = ('solid', 'dotted', 'dashed', 'dashdot')   # Line style follows the ordering defined in blade_loc_id
    
    # Define the color range (one per model)
    model_colors = ('#15B01A', '#EF4026', '#F97306')
    model_marker = ('o', 's', 'p')
    model_line   = ('dashed', 'dotted', 'dashdot')
    model_tag    = ('Pitt-Peters', 'Larsen-Madsen', 'Oye')
    qs_color     = '#069AF3'
    
    # Define what model to be compared; 
    #       USDI = Unsteady-Dynamic Inflow
    comparison_range = ('USDS',)
    
    # Define the folder name to store figures
    figure_folder_tag = 'Figures_Assignment_2' 
    
    # Close and clear all plots
    plt.close('all')
    
    for case_tag_i, case_tag in enumerate(case_tag_range):
        print('=== Case {} ==='.format(case_tag))

        # Plotting the five responses over time for the three models (and the 6 reduced frequencies for case A2 and B2)
        print('== Plotting responses over time ==')
    
        # Define the set of the two flow models to be plotted
        for comparison_i, comparison_tag in enumerate(comparison_range):
            print(comparison_tag)
             # Initialise the plots
            fig_a, (ax_a1, ax_a2, ax_a3, ax_a4) = plt.subplots(4, 1, sharex='all',figsize=(9, 5))            # a: Induction factor
            fig_ct, (ax_ct1, ax_ct2, ax_ct3, ax_ct4) = plt.subplots(4, 1, sharex='all',figsize=(9, 5))       # ct: Thrust coefficient
            fig_cq, (ax_cq1, ax_cq2, ax_cq3, ax_cq4) = plt.subplots(4, 1, sharex='all',figsize=(9, 5))       # cq: Torque coefficient
            fig_aoa, (ax_aoa1, ax_aoa2, ax_aoa3, ax_aoa4) = plt.subplots(4, 1, sharex='all', figsize=(9, 5)) # aoa: Angle of attack (alpha)
            fig_phi, (ax_phi1, ax_phi2, ax_phi3, ax_phi4) = plt.subplots(4, 1, sharex='all',figsize=(9, 5))  # phi: Inflow angle

            # Loop over each model
            for i, model in enumerate(('pp', 'lm', 'oye')):
                print('= Model {} ='.format(model))
                ((r_list, t_list), (ctr, cqr, a, alpha, phi), (ctr_ds, cqr_ds, a_ds, alpha_ds, phi_ds),
                 (ctr_di, cqr_di, a_di, alpha_di, phi_di), (ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs)) = read_data(case_tag, model)

                if comparison_tag == 'USDI':   # Set 1: Fully unsteady; Set 2: Dynamic inflow         
                    set_1_tag = model_tag[i] + ' | Fully unsteady'
                    set_2_tag = model_tag[i] + ' | Dynamic inflow'
                    a_1 = a
                    a_2 = a_di
                    ctr_1 = ctr
                    ctr_2 = ctr_di
                    cqr_1 = cqr
                    cqr_2 = cqr_di
                    alpha_1 = alpha
                    alpha_2 = alpha_di
                    phi_1 = phi
                    phi_2 = phi_di
                elif comparison_tag == 'USDS':  # Set 1: Unsteady; Set 2: Dynamic-stall
                    set_1_tag = model_tag[i] + ' | Fully unsteady'
                    set_2_tag = model_tag[i] + ' | Dynamic stall'
                    a_1 = a
                    a_2 = a_ds
                    ctr_1 = ctr
                    ctr_2 = ctr_ds
                    cqr_1 = cqr
                    cqr_2 = cqr_ds
                    alpha_1 = alpha
                    alpha_2 = alpha_ds
                    phi_1 = phi
                    phi_2 = phi_ds
                elif  comparison_tag == 'USQS': # Set 1: Unsteady; Set 2: Quasi-steady
                    set_1_tag = model_tag[i] + ' | Fully unsteady'
                    set_2_tag = 'Quasi-steady'
                    a_1 = a
                    a_2 = a_qs
                    ctr_1 = ctr
                    ctr_2 = ctr_qs
                    cqr_1 = cqr
                    cqr_2 = cqr_qs
                    alpha_1 = alpha
                    alpha_2 = alpha_qs
                    phi_1 = phi
                    phi_2 = phi_qs
                else:
                    set_1_tag = model_tag[i] + ' | Fully unsteady'
                    set_2_tag = 'Quasi-steady'
                    a_1 = a
                    a_2 = a_qs
                    ctr_1 = ctr
                    ctr_2 = ctr_qs
                    cqr_1 = cqr
                    cqr_2 = cqr_qs
                    alpha_1 = alpha
                    alpha_2 = alpha_qs
                    phi_1 = phi
                    phi_2 = phi_qs
                
                    
                # Assemble the plots
                plot_combined_subplot('a [-]', set_1_tag, set_2_tag, ax_a1, ax_a2, ax_a3, ax_a4, t_list, a_1, a_2, blade_loc_id,
                                      blade_loc_tag, model_colors[i], model_line[i], i)
                plot_combined_subplot('$C_t$ [-]', set_1_tag, set_2_tag, ax_ct1, ax_ct2, ax_ct3, ax_ct4, t_list, ctr_1, ctr_2,
                                      blade_loc_id, blade_loc_tag, model_colors[i], model_line[i], i)
                plot_combined_subplot('$C_q$ [-]', set_1_tag, set_2_tag, ax_cq1, ax_cq2, ax_cq3, ax_cq4, t_list, cqr_1, cqr_2,
                                      blade_loc_id, blade_loc_tag, model_colors[i], model_line[i], i)
                plot_combined_subplot('$\\alpha$ [deg]', set_1_tag, set_2_tag, ax_aoa1, ax_aoa2, ax_aoa3, ax_aoa4, t_list, alpha_1,
                                      alpha_2, blade_loc_id, blade_loc_tag, model_colors[i], model_line[i], i)
                plot_combined_subplot('$\\phi$ [deg]', set_1_tag, set_2_tag, ax_phi1, ax_phi2, ax_phi3, ax_phi4, t_list, phi_1, phi_2,
                                      blade_loc_id, blade_loc_tag, model_colors[i], model_line[i], i)
        
            # Save the plots to .pdf
            plot_save_figure(fig_a, case_tag, 'a', comparison_tag, figure_folder_tag)
            plot_save_figure(fig_ct, case_tag, 'ct', comparison_tag, figure_folder_tag)
            plot_save_figure(fig_cq, case_tag, 'cq', comparison_tag, figure_folder_tag)
            plot_save_figure(fig_aoa, case_tag, 'aoa', comparison_tag, figure_folder_tag)
            plot_save_figure(fig_phi, case_tag, 'phi', comparison_tag, figure_folder_tag)
    
        
        # Plotting responses over blade radial position
        print('== Plotting responses over radial positions ==')
        # Initialise the plots
        fig_a_elem, (ax_a1_elem, ax_a2_elem, ax_a3_elem) = plt.subplots(1, 3, sharey='all',figsize=(9, 5))          # a: Induction factor
        fig_ct_elem, (ax_ct1_elem, ax_ct2_elem, ax_ct3_elem) = plt.subplots(1, 3, sharey='all', figsize=(9, 5))     # ct: Thrust coefficient
        fig_cq_elem, (ax_cq1_elem, ax_cq2_elem, ax_cq3_elem) = plt.subplots(1, 3, sharey='all', figsize=(9, 5))     # cq: Torque coefficient
        fig_aoa_elem, (ax_aoa1_elem, ax_aoa2_elem, ax_aoa3_elem) = plt.subplots(1, 3, sharey='all', figsize=(9, 5)) # aoa: Angle of attack (alpha)
        fig_phi_elem, (ax_phi1_elem, ax_phi2_elem, ax_phi3_elem) = plt.subplots(1, 3, sharey='all', figsize=(9, 5)) # phi: Inflow angle
    
        # Loop over each model
        for model_i, model in enumerate(('pp', 'lm', 'oye')):
            print('= Model {} ='.format(model))
            ((r_list, t_list), (ctr, cqr, a, alpha, phi), (ctr_ds, cqr_ds, a_ds, alpha_ds, phi_ds),
             (ctr_di, cqr_di, a_di, alpha_di, phi_di), (ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs)) = read_data(case_tag, model)
    
            # Initialise the counter of time steps performed
            time_step_counter = 0
    
            # Loop over all time steps of interest
            for row_of_a in range(a.shape[0]):
                time_sampling = 10  # Define the time step for which to retrieve the responses
                if t_list[row_of_a] % time_sampling == 0:
                    # Assemble the plots
                    time_step_grayscale = str((
                                                  time_step_counter) / 6)  # Divide by 6 to have a grayscale with depth, need to update this value if more time steps are considered
                    plot_combined_subplot_elem('a [-]', ax_a1_elem, ax_a2_elem, ax_a3_elem, r_list, a, a_qs,
                                               row_of_a, model_tag, time_step_grayscale, '--',
                                               't [s] = ' + str(t_list[row_of_a]), qs_color, model_i,
                                               time_step_counter)
                    plot_combined_subplot_elem('$C_t$ [-]', ax_ct1_elem, ax_ct2_elem, ax_ct3_elem, r_list, ctr,
                                               ctr_qs, row_of_a, model_tag, time_step_grayscale, '--',
                                               't [s] = ' + str(t_list[row_of_a]), qs_color, model_i,
                                               time_step_counter)
                    plot_combined_subplot_elem('$C_q$ [-]', ax_cq1_elem, ax_cq2_elem, ax_cq3_elem, r_list, cqr,
                                               cqr_qs, row_of_a, model_tag, time_step_grayscale, '--',
                                               't [s] = ' + str(t_list[row_of_a]), qs_color, model_i,
                                               time_step_counter)
                    plot_combined_subplot_elem('$\\alpha$ [deg]', ax_aoa1_elem, ax_aoa2_elem, ax_aoa3_elem,
                                               r_list, alpha, alpha_qs, row_of_a, model_tag,
                                               time_step_grayscale, '--', 't [s] = ' + str(t_list[row_of_a]),
                                               qs_color, model_i, time_step_counter)
                    plot_combined_subplot_elem('$\\phi$ [deg]', ax_phi1_elem, ax_phi2_elem, ax_phi3_elem,
                                               r_list, phi, phi_qs, row_of_a, model_tag, time_step_grayscale,
                                               '--', 't [s] = ' + str(t_list[row_of_a]), qs_color, model_i,
                                               time_step_counter)
    
                    # Increment the counter of time steps performed
                    time_step_counter += 1
    
        # Save the plots to .pdf
        plot_save_figure_elem(fig_a_elem, case_tag, 'a', figure_folder_tag)
        plot_save_figure_elem(fig_ct_elem, case_tag, 'ct', figure_folder_tag)
        plot_save_figure_elem(fig_cq_elem, case_tag, 'cq', figure_folder_tag)
        plot_save_figure_elem(fig_aoa_elem, case_tag, 'aoa', figure_folder_tag)
        plot_save_figure_elem(fig_phi_elem, case_tag, 'phi', figure_folder_tag)
    
    plt.show()
