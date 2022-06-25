import matplotlib.pyplot as plt
import numpy as np

from dynamic_inflow import Turbine
from dynamic_stall import DSTurbine, generate_data, cases
from read_write import read_from_file


def read_data(case, model):
    """
    Read the generated data
    :param case: Selection of the assignment case ('Dyn1' or 'Dyn2' or 'no_les')
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
             np.degrees(read_from_file(f'./{model}/{case}/phi_ds.csv'))),

            (read_from_file(f'./{model}/{case}/ctr_di.csv'),
             read_from_file(f'./{model}/{case}/cqr_di.csv'),
             read_from_file(f'./{model}/{case}/a_di.csv'),
             read_from_file(f'./{model}/{case}/alpha_di.csv'),
             read_from_file(f'./{model}/{case}/phi_di.csv')),

            (read_from_file(f'./{model}/{case}/ctr_qs.csv'),
             read_from_file(f'./{model}/{case}/cqr_qs.csv'),
             read_from_file(f'./{model}/{case}/a_qs.csv'),
             read_from_file(f'./{model}/{case}/alpha_qs.csv'),
             np.degrees(read_from_file(f'./{model}/{case}/phi_qs.csv')))
            )


def airfoil_hist_plots():
    plots = (((0, 'Dyn1', 'lm'), (2, 'Dyn1', 'lm'), '1.1.1'),
             ((0, 'Dyn1', 'lm'), (0, 'no_les', 'lm'), '1.1.2'),
             )

    models = {'pp': 'Pitt-Peters', 'lm': 'Larsen-Madsen', 'oye': 'Øye'}
    case_names = {'Dyn1': 'Dyn1', 'Dyn2': 'Dyn2', 'no_les': 'Dyn1 w/o Leading Edge Separation'}
    data_names = ('Fully unsteady', 'Dynamic stall', 'Dynamic inflow', 'Quasi-steady')
    be_loc = ('0.2 R', '0.5 R', '0.9 R', '1.0 R')

    for ip0, (*plot, fig_name) in enumerate(plots):
        fig0, axes0 = plt.subplots(2, 2, num=0)

        for ip1, (data_idx, name, model) in enumerate(plot):
            (r_list, t_list), *data = read_data(name, model)

            ctr, cqr, a, alpha, phi = data[data_idx]

            u_inf_0, delta_u_inf, freq, v0, tsr, *_ = cases[name]
            t_select = t_list >= .5 * t_list[-1]

            turbine = DSTurbine(25) if data_idx in (0, 1) else Turbine(25)
            no_les = True if name == 'no_les' else False
            delta_t = t_list[1] - t_list[0]
            u_inf = u_inf_0 + delta_u_inf / v0 * np.cos(freq * v0 / turbine.blade.r * t_list)
            u_inf *= v0

            cl = np.empty(alpha.shape)
            for n, t in enumerate(t_list):
                if data_idx in (0, 1):
                    turbine.blade.set_time_data(t, delta_t, u_inf[n], tsr * v0 / u_inf[n], no_les=no_les)
                for bi, be in enumerate(turbine.blade.blade_elements[1:-1]):
                    cl[n, bi] = be.airfoil.cl(alpha[n, bi], propagate=True)

            for ii, be_idx in enumerate((0, 8, -1, -2)):
                plt_idx0, plt_idx1 = ii // 2, ii % 2
                if not ip1:
                    axes0[plt_idx0, plt_idx1].set_title(be_loc[ii])
                    axes0[plt_idx0, plt_idx1].set_ylabel('$c_l$ (-)')
                    axes0[plt_idx0, plt_idx1].set_xlabel('$\\alpha$ (-)')

                if not ii:
                    axes0[plt_idx0, plt_idx1].plot(alpha[t_select, be_idx], cl[t_select, be_idx],
                                                   label=f'{case_names[name]}, {models[model]} | {data_names[data_idx]}.')
                else:
                    axes0[plt_idx0, plt_idx1].plot(alpha[t_select, be_idx], cl[t_select, be_idx])

        plt.tight_layout()
        fig0.subplots_adjust(bottom=0.25)
        fig0.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.15), ncol=1)
        plt.savefig(f'./Figures_Assignment_2/loops/{fig_name}_cl_alpha.pdf')

        plt.show()


def ct_a_plots():
    plots = (((0, 'Dyn1', 'lm'), (2, 'Dyn1', 'lm'), '1.1.1'),
             ((0, 'Dyn1', 'lm'), (0, 'no_les', 'lm'), '1.1.2'),
             ((0, 'Dyn2', 'lm'), (2, 'Dyn2', 'lm'), '1.2.2'),
             )

    qtt = ('ct', 'cq')
    models = {'pp': 'Pitt-Peters', 'lm': 'Larsen-Madsen', 'oye': 'Øye'}
    case_names = {'Dyn1': 'Dyn1', 'Dyn2': 'Dyn2', 'no_les': 'Dyn1 w/o Leading Edge Separation'}
    data_names = ('Fully unsteady', 'Dynamic stall', 'Dynamic inflow', 'Quasi-steady')
    be_loc = ('0.2 R', '0.5 R', '0.9 R', '1.0 R')

    for ip0, (*plot, fig_name) in enumerate(plots):
        fig0, axes0 = plt.subplots(2, 2, num=0)
        fig1, axes1 = plt.subplots(2, 2, num=1)

        for ip1, (data_idx, name, model) in enumerate(plot):
            (r_list, t_list), *data = read_data(name, model)

            ctr, cqr, a, alpha, phi = data[data_idx]

            _, _, freq, v0, _, *_ = cases[name]
            t_select = t_list >= .5 * t_list[-1]

            for ii, be_idx in enumerate((0, 8, -1, -2)):
                plt_idx0, plt_idx1 = ii // 2, ii % 2
                if not ip1:
                    axes0[plt_idx0, plt_idx1].set_title(be_loc[ii])
                    axes0[plt_idx0, plt_idx1].set_xlabel('$C_T$ (-)')
                    axes0[plt_idx0, plt_idx1].set_ylabel('$a$ (-)')

                    axes1[plt_idx0, plt_idx1].set_title(be_loc[ii])
                    axes1[plt_idx0, plt_idx1].set_xlabel('$C_Q$ x 1e3 (-)')
                    axes1[plt_idx0, plt_idx1].set_ylabel('$a$ (-)')

                if not ii:
                    axes0[plt_idx0, plt_idx1].plot(ctr[t_select, be_idx], a[t_select, be_idx],
                                                   label=f'{case_names[name]}, {models[model]} | {data_names[data_idx]}.')
                    axes1[plt_idx0, plt_idx1].plot(cqr[t_select, be_idx] * 1e3, a[t_select, be_idx],
                                                   label=f'{case_names[name]}, {models[model]} | {data_names[data_idx]}.')
                else:
                    axes0[plt_idx0, plt_idx1].plot(ctr[t_select, be_idx], a[t_select, be_idx])
                    axes1[plt_idx0, plt_idx1].plot(cqr[t_select, be_idx] * 1e3, a[t_select, be_idx])

        for fi in (0, 1):
            plt.figure(fi)
            plt.tight_layout()

        fig0.subplots_adjust(bottom=0.25)
        fig0.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.15), ncol=1)
        fig1.subplots_adjust(bottom=0.25)
        fig1.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.15), ncol=1)

        for fi in (0, 1):
            plt.figure(fi)
            plt.savefig(f'./Figures_Assignment_2/loops/{fig_name}_{qtt[fi]}.pdf')

        plt.show()



def example_read():
    colors = ('r', 'g', 'b')
    model = 'lm'
    for name in cases.keys():
        ((r_list, t_list), (ctr, cqr, a, alpha, phi), (ctr_ds, cqr_ds, a_ds, alpha_ds, phi_ds),
         (ctr_di, cqr_di, a_di, alpha_di, phi_di), (ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs)) = read_data(name, model)

        for j in (0, 8, -2, -1)[-1:]:
            plt.figure(1)
            plt.plot(t_list, a[:, j], 'k')
            # plt.plot(t_list, a_ds[:, j], 'k', linestyle='dashdot')
            plt.plot(t_list, a_di[:, j], 'k', linestyle='dashed')
            # plt.plot(t_list, a_qs[:, j], 'k', linestyle='dotted')

            plt.figure(2)
            plt.plot(t_list, ctr[:, j], 'k')
            # plt.plot(t_list, ctr_ds[:, j], 'k', linestyle='dashdot')
            plt.plot(t_list, ctr_di[:, j], 'k', linestyle='dashed')
            # plt.plot(t_list, ctr_qs[:, j], 'k', linestyle='dotted')

            plt.figure(3)
            plt.plot(t_list, alpha[:, j], 'k')
            # plt.plot(t_list, alpha_ds[:, j], 'k', linestyle='dashdot')
            plt.plot(t_list, alpha_di[:, j], 'k', linestyle='dashed')
            # plt.plot(t_list, alpha_qs[:, j], 'k', linestyle='dotted')

            plt.figure(4)
            plt.plot(t_list, cqr[:, j], 'k')
            # plt.plot(t_list, cqr_ds[:, j], 'k', linestyle='dashdot')
            plt.plot(t_list, cqr_di[:, j], 'k', linestyle='dashed')
            # plt.plot(t_list, cqr_qs[:, j], 'k', linestyle='dotted')

            plt.figure(5)
            plt.plot(t_list, phi[:, j], 'k')
            # plt.plot(t_list, phi_ds[:, j], 'k', linestyle='dashdot')
            plt.plot(t_list, phi_di[:, j], 'k', linestyle='dashed')
            # plt.plot(t_list, phi_qs[:, j], 'k', linestyle='dotted')

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

def plot_combined_subplot(y_label, set_1_tag, set_2_tag, ax1, ax2, ax3, ax4, x_lst, y_1_mat, y_2_mat, blade_loc_id, blade_loc_tag, model_set1_color, model_set2_color, line):
    if set_2_tag == 'Quasi-steady':
        ax1.plot(x_lst, y_1_mat[:, blade_loc_id[0]], color=model_set1_color, linestyle=line, label=set_1_tag)
        ax1.set_title(blade_loc_tag[0])
        ax1.set_ylabel(y_label)
        ax1.grid()
        ax2.plot(x_lst, y_1_mat[:, blade_loc_id[1]], color=model_set1_color, linestyle=line)
        ax2.set_title(blade_loc_tag[1])
        ax2.set_ylabel(y_label)
        ax2.grid()
        ax3.plot(x_lst, y_1_mat[:, blade_loc_id[2]], color=model_set1_color, linestyle=line)
        ax3.set_title(blade_loc_tag[2])
        ax3.set_ylabel(y_label)
        ax3.grid()
        ax4.plot(x_lst, y_1_mat[:, blade_loc_id[3]], color=model_set1_color, linestyle=line)
        ax4.set_title(blade_loc_tag[3])
        ax4.set_ylabel(y_label)
        ax4.set_xlabel('Time [s]')
        ax4.grid()
        if not i:
            ax1.plot(x_lst, y_2_mat[:, blade_loc_id[0]], color=model_set2_color, linestyle='dashdot', label=set_2_tag)
            ax2.plot(x_lst, y_2_mat[:, blade_loc_id[1]], color=model_set2_color, linestyle='dashdot')
            ax3.plot(x_lst, y_2_mat[:, blade_loc_id[2]], color=model_set2_color, linestyle='dashdot')
            ax4.plot(x_lst, y_2_mat[:, blade_loc_id[3]], color=model_set2_color, linestyle='dashdot')
        
    else:
        ax1.plot(x_lst, y_1_mat[:, blade_loc_id[0]], color=model_set1_color, linestyle=line, label=set_1_tag)
        ax1.plot(x_lst, y_2_mat[:, blade_loc_id[0]], color=model_set1_color, linestyle='dashdot', label=set_2_tag)
        ax1.set_title(blade_loc_tag[0])
        ax1.set_ylabel(y_label)
        ax1.grid()
        ax2.plot(x_lst, y_1_mat[:, blade_loc_id[1]], color=model_set1_color, linestyle=line)
        ax2.plot(x_lst, y_2_mat[:, blade_loc_id[1]], color=model_set1_color, linestyle='dashdot')
        ax2.set_title(blade_loc_tag[1])
        ax2.set_ylabel(y_label)
        ax2.grid()
        ax3.plot(x_lst, y_1_mat[:, blade_loc_id[2]], color=model_set1_color, linestyle=line)
        ax3.plot(x_lst, y_2_mat[:, blade_loc_id[2]], color=model_set1_color, linestyle='dashdot')
        ax3.set_title(blade_loc_tag[2])
        ax3.set_ylabel(y_label)
        ax3.grid()
        ax4.plot(x_lst, y_1_mat[:, blade_loc_id[3]], color=model_set1_color, linestyle=line)
        ax4.plot(x_lst, y_2_mat[:, blade_loc_id[3]], color=model_set1_color, linestyle='dashdot')
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

def plot_combined_subplot_elem_one_model(y_label,ax1,x_lst,y_set1_mat,y_set_qs_mat,row_interest,model_tags,color,line,line_label,qs_color,time_step_counter):
    if time_step_counter == 0:
        # Plot the initial state of quasi-steady solution
        ax1.plot(x_lst, y_set_qs_mat[1, :], color='#069AF3', linestyle='solid', label='Quasi-steady initial')
        # Plot the final state of quasi-steady solution
        ax1.plot(x_lst, y_set_qs_mat[-1, :], color='#F97306', linestyle='solid', label='Quasi-steady final')
        ax1.grid()
    else:
        ax1.plot(x_lst, y_set1_mat[row_interest, :], color=color, linestyle=line, label=line_label)
        ax1.set_title(model_tags)
        ax1.set_ylabel(y_label)
        ax1.set_xlabel('Blade radial position [m]')
    return

def plot_combined_subplot_elem_one_model_meanmaxmin(y_label,ax1,x_lst,y_set1_mean,y_set1_max,y_set1_min,y_set_qs_mean,y_set_qs_max,y_set_qs_min,model_tags):
    # Plot the mean,max,min states of quasi-steady solution
    ax1.plot(x_lst, y_set_qs_mean, color='#069AF3', linestyle='solid', label='Quasi-steady | Mean')
    ax1.plot(x_lst, y_set_qs_max, color='#069AF3', linestyle='dashed', label='Quasi-steady | Max')
    ax1.plot(x_lst, y_set_qs_min, color='#069AF3', linestyle='dotted', label='Quasi-steady | Min')
    # Plot the mean,max,min states of solution set 1
    ax1.plot(x_lst, y_set1_mean, color='#F97306', linestyle='solid', label=model_tags+ ' | Mean')
    ax1.plot(x_lst, y_set1_max, color='#F97306', linestyle='dashed', label=model_tags+' | Max')
    ax1.plot(x_lst, y_set1_min, color='#F97306', linestyle='dotted', label=model_tags+' | Min')
    ax1.set_ylabel(y_label)
    ax1.set_xlabel('Blade radial position [m]')
    ax1.grid()
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
    
def plot_save_figure_elem_meanmaxmin(fig_tag, case_tag, response_tag, folder_name):
    fig_tag.tight_layout()
    fig_tag.subplots_adjust(bottom=0.35)
    fig_tag.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.21), ncol=2)
    fig_name = case_tag + '_' + response_tag + '_MeanMaxMin_blade_elem.pdf'
    fig_tag.savefig(folder_name + '\\' + fig_name)

def getMeanMaxMinAllColumnsDarray(mat, AxisInterest):
    mat_mean_array = np.mean(mat,axis=AxisInterest)
    mat_max_array  = np.max(mat,axis=AxisInterest)
    mat_min_array  = np.min(mat,axis=AxisInterest)
    return mat_mean_array, mat_max_array, mat_min_array

        
        
if __name__ == '__main__':
    # Close and clear all plots
    plt.close('all')
    
    # To generate the data, uncomment the following two lines
    # generate_data()
    # example_read()
    
    # Plot the ct-alpha hysteresis
    #ct_a_plots()
    #airfoil_hist_plots()
    
    # Define the range of case tags considered
    case_tag_range = ('Dyn1','Dyn2','no_les')
    
    # Define blade locations of interest and plotting styles
    blade_loc_id = (0, 8, -2, -1)                               # Index used correspond to radial distribution indexing
    blade_loc_tag = ('0.2 R', '0.5 R', '0.9 R', '1.0 R')        # Labelling follows the ordering defined in blade_loc_id
    blade_loc_line = ('solid', 'dotted', 'dashed', 'dashdot')   # Line style follows the ordering defined in blade_loc_id
    
    # Define model for which to plot responses over radial position
    model = 'lm'
    
    # Define the plotting format for the model
    model_line   = 'dashed'
    model_tag    = 'Larsen-Madsen'
    model_set1_color  = '#EF4026'
    model_set2_color  = '#069AF3'
    
    # Define what model to be compared; 
    #       USDI = Unsteady-Dynamic Inflow
    comparison_range = ('USDI',)
    
    # Define the folder name to store figures
    figure_folder_tag = 'Figures_Assignment_2' 
    
    # Define the number of lines that shall be hown in the blade element plots
    n_time_lines = 6
    
    for case_tag_i, case_tag in enumerate(case_tag_range):
        print('=== Case {} ==='.format(case_tag))
#    
#        # Define the set of the two flow models to be plotted
#        for comparison_i, comparison_tag in enumerate(comparison_range):
#            print('-- Comparison {} --'.format(comparison_tag))
#            
#            # Plotting the five responses over time for the three models (and the 6 reduced frequencies for case A2 and B2)
#            print('== Plotting responses over time ==')
#             # Initialise the plots
#            fig_a, (ax_a1, ax_a2, ax_a3, ax_a4) = plt.subplots(4, 1, sharex='all',figsize=(9, 5))            # a: Induction factor
#            fig_ct, (ax_ct1, ax_ct2, ax_ct3, ax_ct4) = plt.subplots(4, 1, sharex='all',figsize=(9, 5))       # ct: Thrust coefficient
#            fig_cq, (ax_cq1, ax_cq2, ax_cq3, ax_cq4) = plt.subplots(4, 1, sharex='all',figsize=(9, 5))       # cq: Torque coefficient
#            fig_aoa, (ax_aoa1, ax_aoa2, ax_aoa3, ax_aoa4) = plt.subplots(4, 1, sharex='all', figsize=(9, 5)) # aoa: Angle of attack (alpha)
#            fig_phi, (ax_phi1, ax_phi2, ax_phi3, ax_phi4) = plt.subplots(4, 1, sharex='all',figsize=(9, 5))  # phi: Inflow angle
#
#            # Retreive the model of interest
#            print('Model {} '.format(model))
#            
#            # Retrieve the responses of the the model of interest
#            ((r_list, t_list), (ctr, cqr, a, alpha, phi), (ctr_ds, cqr_ds, a_ds, alpha_ds, phi_ds),
#             (ctr_di, cqr_di, a_di, alpha_di, phi_di), (ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs)) = read_data(case_tag, model)
#
#            if comparison_tag == 'USDI':   # Set 1: Fully unsteady; Set 2: Dynamic inflow         
#                set_1_tag = model_tag + ' | Fully unsteady'
#                set_2_tag = model_tag + ' | Dynamic inflow'
#                a_1 = a
#                a_2 = a_di
#                ctr_1 = ctr
#                ctr_2 = ctr_di
#                cqr_1 = cqr
#                cqr_2 = cqr_di
#                alpha_1 = alpha
#                alpha_2 = alpha_di
#                phi_1 = phi
#                phi_2 = phi_di
#            elif comparison_tag == 'USDS':  # Set 1: Unsteady; Set 2: Dynamic-stall
#                set_1_tag = model_tag + ' | Fully unsteady'
#                set_2_tag = model_tag + ' | Dynamic stall'
#                a_1 = a
#                a_2 = a_ds
#                ctr_1 = ctr
#                ctr_2 = ctr_ds
#                cqr_1 = cqr
#                cqr_2 = cqr_ds
#                alpha_1 = alpha
#                alpha_2 = alpha_ds
#                phi_1 = phi
#                phi_2 = phi_ds
#            elif  comparison_tag == 'USQS': # Set 1: Unsteady; Set 2: Quasi-steady
#                set_1_tag = model_tag + ' | Fully unsteady'
#                set_2_tag = 'Quasi-steady'
#                a_1 = a
#                a_2 = a_qs
#                ctr_1 = ctr
#                ctr_2 = ctr_qs
#                cqr_1 = cqr
#                cqr_2 = cqr_qs
#                alpha_1 = alpha
#                alpha_2 = alpha_qs
#                phi_1 = phi
#                phi_2 = phi_qs
#            else:
#                set_1_tag = model_tag + ' | Fully unsteady'
#                set_2_tag = 'Quasi-steady'
#                a_1 = a
#                a_2 = a_qs
#                ctr_1 = ctr
#                ctr_2 = ctr_qs
#                cqr_1 = cqr
#                cqr_2 = cqr_qs
#                alpha_1 = alpha
#                alpha_2 = alpha_qs
#                phi_1 = phi
#                phi_2 = phi_qs
#            
#                
#            # Assemble the plots
#            plot_combined_subplot('a [-]', set_1_tag, set_2_tag, ax_a1, ax_a2, ax_a3, ax_a4, t_list, a_1, a_2, blade_loc_id,
#                                  blade_loc_tag, model_set1_color, model_set2_color, model_line)
#            plot_combined_subplot('$C_t$ [-]', set_1_tag, set_2_tag, ax_ct1, ax_ct2, ax_ct3, ax_ct4, t_list, ctr_1, ctr_2,
#                                  blade_loc_id, blade_loc_tag, model_set1_color, model_set2_color, model_line)
#            plot_combined_subplot('$C_q$ [-]', set_1_tag, set_2_tag, ax_cq1, ax_cq2, ax_cq3, ax_cq4, t_list, cqr_1, cqr_2,
#                                  blade_loc_id, blade_loc_tag, model_set1_color, model_set2_color, model_line)
#            plot_combined_subplot('$\\alpha$ [deg]', set_1_tag, set_2_tag, ax_aoa1, ax_aoa2, ax_aoa3, ax_aoa4, t_list, alpha_1,
#                                  alpha_2, blade_loc_id, blade_loc_tag, model_set1_color, model_set2_color, model_line)
#            plot_combined_subplot('$\\phi$ [deg]', set_1_tag, set_2_tag, ax_phi1, ax_phi2, ax_phi3, ax_phi4, t_list, phi_1, phi_2,
#                                  blade_loc_id, blade_loc_tag, model_set1_color, model_set2_color, model_line)
#    
#        # Save the plots to .pdf
#        plot_save_figure(fig_a, case_tag, 'a', comparison_tag, figure_folder_tag)
#        plot_save_figure(fig_ct, case_tag, 'ct', comparison_tag, figure_folder_tag)
#        plot_save_figure(fig_cq, case_tag, 'cq', comparison_tag, figure_folder_tag)
#        plot_save_figure(fig_aoa, case_tag, 'aoa', comparison_tag, figure_folder_tag)
#        plot_save_figure(fig_phi, case_tag, 'phi', comparison_tag, figure_folder_tag)
#    
#        
#        # Plotting responses over blade radial position
#        print('== Plotting responses over radial positions ==')
#        
#        # Initialise the plots
#        fig_a_elem, ax_a1_elem = plt.subplots(1, 1, sharey='all',figsize=(9, 5))          # a: Induction factor
#        fig_ct_elem, ax_ct1_elem = plt.subplots(1, 1, sharey='all', figsize=(9, 5))     # ct: Thrust coefficient
#        fig_cq_elem, ax_cq1_elem = plt.subplots(1, 1, sharey='all', figsize=(9, 5))     # cq: Torque coefficient
#        fig_aoa_elem, ax_aoa1_elem = plt.subplots(1, 1, sharey='all', figsize=(9, 5)) # aoa: Angle of attack (alpha)
#        fig_phi_elem, ax_phi1_elem = plt.subplots(1, 1, sharey='all', figsize=(9, 5)) # phi: Inflow angle
#
#        # Retreive the model of interest
#        print('Model {}'.format(model))
#        
#        # Retrieve the responses of the the model of interest
#        ((r_list, t_list), (ctr, cqr, a, alpha, phi), (ctr_ds, cqr_ds, a_ds, alpha_ds, phi_ds),
#         (ctr_di, cqr_di, a_di, alpha_di, phi_di), (ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs)) = read_data(case_tag, model)
#
#        # Initialise the counter of time steps performed
#        time_step_counter = 0
#        
#        # Loop over all time steps of interest
#        for row_of_a in range(a.shape[0]):
#            # Evaluate the sampling time required to plot only 'n_time_line' number of time lines
#            time_sampling = round(t_list[-1]/100)*100 * 1/(n_time_lines-1)
#            
#            # If time range too short, then sampling time will be zero, thus force it to be every 1 second
#            if time_sampling == 0:
#                time_sampling = 1
#                
#            # If the time equals to a multiple of the sampling time, then plot it
#            if t_list[row_of_a] % time_sampling == 0:
#                # Define the gray-scale color
#                time_step_grayscale = str((time_step_counter+2)/(n_time_lines*2))
#
#                # Assemble the plots
#                plot_combined_subplot_elem_one_model('a [-]', ax_a1_elem,r_list, a, a_qs,row_of_a, 'Larsen-Madsen | Fully unsteady', time_step_grayscale, '--','t [s] = ' + str(t_list[row_of_a]), model_set1_color,time_step_counter)
#                plot_combined_subplot_elem_one_model('$C_t$ [-]', ax_ct1_elem,r_list, ctr, ctr_qs,row_of_a, 'Larsen-Madsen | Fully unsteady', time_step_grayscale, '--','t [s] = ' + str(t_list[row_of_a]), model_set1_color,time_step_counter)
#                plot_combined_subplot_elem_one_model('$C_q$ [-]', ax_cq1_elem,r_list, cqr, cqr_qs,row_of_a, 'Larsen-Madsen | Fully unsteady', time_step_grayscale, '--','t [s] = ' + str(t_list[row_of_a]), model_set1_color,time_step_counter)
#                plot_combined_subplot_elem_one_model('$\\alpha$ [deg]', ax_aoa1_elem,r_list, alpha, alpha_qs,row_of_a, 'Larsen-Madsen | Fully unsteady', time_step_grayscale, '--','t [s] = ' + str(t_list[row_of_a]), model_set1_color,time_step_counter)
#                plot_combined_subplot_elem_one_model('$\\phi$ [deg]', ax_phi1_elem,r_list, phi, phi_qs,row_of_a, 'Larsen-Madsen | Fully unsteady', time_step_grayscale, '--','t [s] = ' + str(t_list[row_of_a]), model_set1_color,time_step_counter)
#            
#                # Increment the counter of time steps performed
#                time_step_counter += 1
#    
#        # Save the plots to .pdf
#        plot_save_figure_elem(fig_a_elem, case_tag, 'a', figure_folder_tag)
#        plot_save_figure_elem(fig_ct_elem, case_tag, 'ct', figure_folder_tag)
#        plot_save_figure_elem(fig_cq_elem, case_tag, 'cq', figure_folder_tag)
#        plot_save_figure_elem(fig_aoa_elem, case_tag, 'aoa', figure_folder_tag)
#        plot_save_figure_elem(fig_phi_elem, case_tag, 'phi', figure_folder_tag)
        
        # Plotting responses over blade radial position
        print('== Plotting responses over radial positions with mean, max, min ==')
        
        # Initialise the plots
        fig_a_elem, ax_a1_elem = plt.subplots(1, 1, sharey='all',figsize=(9, 5))          # a: Induction factor
        fig_ct_elem, ax_ct1_elem = plt.subplots(1, 1, sharey='all', figsize=(9, 5))     # ct: Thrust coefficient
        fig_cq_elem, ax_cq1_elem = plt.subplots(1, 1, sharey='all', figsize=(9, 5))     # cq: Torque coefficient
        fig_aoa_elem, ax_aoa1_elem = plt.subplots(1, 1, sharey='all', figsize=(9, 5)) # aoa: Angle of attack (alpha)
        fig_phi_elem, ax_phi1_elem = plt.subplots(1, 1, sharey='all', figsize=(9, 5)) # phi: Inflow angle

        # Retreive the model of interest
        print('Model {}'.format(model))
        
        # Retrieve the responses of the the model of interest
        ((r_list, t_list), (ctr, cqr, a, alpha, phi), (ctr_ds, cqr_ds, a_ds, alpha_ds, phi_ds),
         (ctr_di, cqr_di, a_di, alpha_di, phi_di), (ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs)) = read_data(case_tag, model)
        
        a_mean,a_max,a_min              = getMeanMaxMinAllColumnsDarray(a, 0)
        ctr_mean,ctr_max,ctr_min        = getMeanMaxMinAllColumnsDarray(ctr, 0)
        cqr_mean,cqr_max,cqr_min        = getMeanMaxMinAllColumnsDarray(cqr, 0)
        alpha_mean,alpha_max,alpha_min  = getMeanMaxMinAllColumnsDarray(alpha, 0)
        phi_mean,phi_max,phi_min        = getMeanMaxMinAllColumnsDarray(phi, 0)
        
        a_qs_mean,a_qs_max,a_qs_min              = getMeanMaxMinAllColumnsDarray(a_qs, 0)
        ctr_qs_mean,ctr_qs_max,ctr_qs_min        = getMeanMaxMinAllColumnsDarray(ctr_qs, 0)
        cqr_qs_mean,cqr_qs_max,cqr_qs_min        = getMeanMaxMinAllColumnsDarray(cqr_qs, 0)
        alpha_qs_mean,alpha_qs_max,alpha_qs_min  = getMeanMaxMinAllColumnsDarray(alpha_qs, 0)
        phi_qs_mean,phi_qs_max,phi_qs_min        = getMeanMaxMinAllColumnsDarray(phi_qs, 0)

        # Assemble the plots
        plot_combined_subplot_elem_one_model_meanmaxmin('a [-]', ax_a1_elem,r_list,a_mean,a_max,a_min,a_qs_mean,a_qs_max,a_qs_min,'Fully unsteady | Larsen-Madsen')
        plot_combined_subplot_elem_one_model_meanmaxmin('$C_t$ [-]', ax_ct1_elem,r_list,ctr_mean,ctr_max,ctr_min,ctr_qs_mean,ctr_qs_max,ctr_qs_min,'Fully unsteady | Larsen-Madsen')
        plot_combined_subplot_elem_one_model_meanmaxmin('$C_q$ [-]', ax_cq1_elem,r_list,cqr_mean,cqr_max,cqr_min,cqr_qs_mean,cqr_qs_max,cqr_qs_min,'Fully unsteady | Larsen-Madsen')
        plot_combined_subplot_elem_one_model_meanmaxmin('$\\alpha$ [deg]', ax_aoa1_elem,r_list,alpha_mean,alpha_max,alpha_min,alpha_qs_mean,alpha_qs_max,alpha_qs_min,'Fully unsteady | Larsen-Madsen')
        plot_combined_subplot_elem_one_model_meanmaxmin('$\\phi$ [deg]', ax_phi1_elem,r_list,phi_mean,phi_max,phi_min,phi_qs_mean,phi_qs_max,phi_qs_min,'Fully unsteady | Larsen-Madsen')

        # Save the plots to .pdf
        plot_save_figure_elem_meanmaxmin(fig_a_elem, case_tag, 'a', figure_folder_tag)
        plot_save_figure_elem_meanmaxmin(fig_ct_elem, case_tag, 'ct', figure_folder_tag)
        plot_save_figure_elem_meanmaxmin(fig_cq_elem, case_tag, 'cq', figure_folder_tag)
        plot_save_figure_elem_meanmaxmin(fig_aoa_elem, case_tag, 'aoa', figure_folder_tag)
        plot_save_figure_elem_meanmaxmin(fig_phi_elem, case_tag, 'phi', figure_folder_tag)
    
    plt.show()
