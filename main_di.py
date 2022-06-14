import matplotlib.pyplot as plt

from read_write import read_from_file
from dynamic_inflow import Turbine, generate_data, ct_steps, ct_sins, u_inf_steps, u_inf_sins


def read_data(select, initial, delta, reduced_freq, model):
    if reduced_freq is None:
        return (read_from_file(f'./{model}/{select}_step/{initial}_{delta}_r_list.csv')[0],
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_t_list.csv')[0],
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_ctr.csv'),
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_cqr.csv'),
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_a.csv'),
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_alpha.csv'),
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_phi.csv'),
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_ctr_qs.csv'),
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_cqr_qs.csv'),
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_a_qs.csv'),
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_alpha_qs.csv'),
                read_from_file(f'./{model}/{select}_step/{initial}_{delta}_phi_qs.csv'))
    else:
        return (read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_r_list.csv')[0],
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_t_list.csv')[0],
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_ctr.csv'),
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_cqr.csv'),
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_a.csv'),
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_alpha.csv'),
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_phi.csv'),
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_ctr_qs.csv'),
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_cqr_qs.csv'),
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_a_qs.csv'),
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_alpha_qs.csv'),
                read_from_file(f'./{model}/{select}_sin/{initial}_{delta}_{reduced_freq}_phi_qs.csv'))

def plot_combined_subplot(y_label,ax1,ax2,ax3,ax4,x_lst,y_mat,y_qs_mat,blade_loc_id,blade_loc_tag,color,line,tag,i):
    ax1.plot(x_lst, y_mat[:, blade_loc_id[0]], color=color, linestyle=line, label=tag)
    ax1.set_title(blade_loc_tag[0])
    ax1.set_ylabel(y_label)
    ax1.grid()
    ax2.plot(x_lst, y_mat[:, blade_loc_id[1]], color=color, linestyle=line)
    ax2.set_title(blade_loc_tag[1])
    ax2.set_ylabel(y_label)
    ax2.grid()
    ax3.plot(x_lst, y_mat[:, blade_loc_id[2]], color=color, linestyle=line)
    ax3.set_title(blade_loc_tag[2])
    ax3.set_ylabel(y_label)
    ax3.grid()
    ax4.plot(x_lst, y_mat[:, blade_loc_id[3]], color=color, linestyle=line)
    ax4.set_title(blade_loc_tag[3])
    ax4.set_ylabel(y_label)
    ax4.set_xlabel('Time [s]')
    ax4.grid()
    if not i:
        ax1.plot(x_lst, y_qs_mat[:, blade_loc_id[0]], color=qs_color, linestyle='solid', label='Quasi-steady')
        ax2.plot(x_lst, y_qs_mat[:, blade_loc_id[1]], color=qs_color, linestyle='solid')
        ax3.plot(x_lst, y_qs_mat[:, blade_loc_id[2]], color=qs_color, linestyle='solid')
        ax4.plot(x_lst, y_qs_mat[:, blade_loc_id[3]], color=qs_color, linestyle='solid')
    return

def plot_combined_subplot_red_freq(y_label,ax1,ax2,ax3,ax4,x_lst,y_mat,y_qs_mat,blade_loc_id,blade_loc_tag,color,line,tag,i):
    # Plot once the quasi steady solution with its label
    if not i:
        ax1.plot(x_lst, y_qs_mat[:, blade_loc_id[0]], color=qs_color, linestyle='solid', label='Quasi-steady')
        ax2.plot(x_lst, y_qs_mat[:, blade_loc_id[1]], color=qs_color, linestyle='solid')
        ax3.plot(x_lst, y_qs_mat[:, blade_loc_id[2]], color=qs_color, linestyle='solid')
        ax4.plot(x_lst, y_qs_mat[:, blade_loc_id[3]], color=qs_color, linestyle='solid')
        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax4.grid()
    # Otherwise just plot the quasi steady solution, without label
    else:
        ax1.plot(x_lst, y_qs_mat[:, blade_loc_id[0]], color=qs_color, linestyle='solid')
        ax2.plot(x_lst, y_qs_mat[:, blade_loc_id[1]], color=qs_color, linestyle='solid')
        ax3.plot(x_lst, y_qs_mat[:, blade_loc_id[2]], color=qs_color, linestyle='solid')
        ax4.plot(x_lst, y_qs_mat[:, blade_loc_id[3]], color=qs_color, linestyle='solid')

    ax1.plot(x_lst, y_mat[:, blade_loc_id[0]], color=color, linestyle=line, label=tag)
    ax1.set_title(blade_loc_tag[0])
    ax1.set_ylabel(y_label)
    ax1.set_xlim(0,628.4)
    ax2.plot(x_lst, y_mat[:, blade_loc_id[1]], color=color, linestyle=line)
    ax2.set_title(blade_loc_tag[1])
    ax2.set_ylabel(y_label)
    ax2.set_xlim(0,628.4)
    ax3.plot(x_lst, y_mat[:, blade_loc_id[2]], color=color, linestyle=line)
    ax3.set_title(blade_loc_tag[2])
    ax3.set_ylabel(y_label)
    ax3.set_xlim(0,628.4)
    ax4.plot(x_lst, y_mat[:, blade_loc_id[3]], color=color, linestyle=line)
    ax4.set_title(blade_loc_tag[3])
    ax4.set_ylabel(y_label)
    ax4.set_xlabel('Time [s]')
    ax4.set_xlim(0,628.4)
    return

def plot_combined_subplot_red_freq_norm(y_label,ax1,ax2,ax3,ax4,x_lst,y_mat,y_qs_mat,blade_loc_id,blade_loc_tag,color,line,tag,i):
    # Plot once the quasi steady solution with its label
    if not i:
        ax1.plot(x_lst, y_qs_mat[:, blade_loc_id[0]], color=qs_color, linestyle='solid', label='Quasi-steady')
        ax2.plot(x_lst, y_qs_mat[:, blade_loc_id[1]], color=qs_color, linestyle='solid')
        ax3.plot(x_lst, y_qs_mat[:, blade_loc_id[2]], color=qs_color, linestyle='solid')
        ax4.plot(x_lst, y_qs_mat[:, blade_loc_id[3]], color=qs_color, linestyle='solid')
        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax4.grid()
    # Otherwise just plot the quasi steady solution, without label
    else:
        ax1.plot(x_lst, y_qs_mat[:, blade_loc_id[0]], color=qs_color, linestyle='solid')
        ax2.plot(x_lst, y_qs_mat[:, blade_loc_id[1]], color=qs_color, linestyle='solid')
        ax3.plot(x_lst, y_qs_mat[:, blade_loc_id[2]], color=qs_color, linestyle='solid')
        ax4.plot(x_lst, y_qs_mat[:, blade_loc_id[3]], color=qs_color, linestyle='solid')

    ax1.plot(x_lst, y_mat[:, blade_loc_id[0]], color=color, linestyle=line, label=tag)
    ax1.set_title(blade_loc_tag[0])
    ax1.set_ylabel(y_label)
    ax2.plot(x_lst, y_mat[:, blade_loc_id[1]], color=color, linestyle=line)
    ax2.set_title(blade_loc_tag[1])
    ax2.set_ylabel(y_label)
    ax3.plot(x_lst, y_mat[:, blade_loc_id[2]], color=color, linestyle=line)
    ax3.set_title(blade_loc_tag[2])
    ax3.set_ylabel(y_label)
    ax4.plot(x_lst, y_mat[:, blade_loc_id[3]], color=color, linestyle=line)
    ax4.set_title(blade_loc_tag[3])
    ax4.set_ylabel(y_label)
    ax4.set_xlabel('Non-dimensional time: $k \\cdot s$  $\\left(\\omega t \\right)$ [-]')
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


def plot_save_figure(fig_tag,case_tag,case_ID,response_tag,freq_red_tag,folder_name):
    fig_tag.tight_layout()
    if freq_red_tag == 0.0:
        fig_tag.subplots_adjust(bottom=0.15)
        fig_tag.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.07), ncol=4)
        fig_name = case_tag + '_' + str(case_ID) + '_' + response_tag + '_time.pdf'
    else:
        fig_tag.subplots_adjust(bottom=0.28)
        fig_tag.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.14), ncol=4)
        fig_name = case_tag + '_' + str(case_ID) + '_' + response_tag + '_k' + str(freq_red_tag) + '_time.pdf'
    fig_tag.savefig(folder_name + '\\' + fig_name)
    return

def plot_save_figure_elem(fig_tag,case_tag,case_ID,response_tag,folder_name):
    fig_tag.tight_layout()
    fig_tag.subplots_adjust(bottom=0.2)
    fig_tag.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.1), ncol=4)
    fig_name = case_tag + '_' + str(case_ID) + '_' + response_tag + '_blade_elem.pdf'
    fig_tag.savefig(folder_name + '\\' + fig_name)
    return

if __name__ == '__main__':
    # Create the turbine with 25 blade elements
    turbine = Turbine(25)
    generate_data()

#    # User inputs: Uncomment the following line to manually select the case and condition to be plotted
#    # Select the case to be plotted; Either: A1, A2, B1, B2
#    case_tag = 'B2'
#
#    # Select the condition number to be plotted (i.e. row number of interest in the table of the assignment);
#    #   For A1 : 1-4
#    #   For A2 : 1-3
#    #   For B1 : 1-4
#    #   For B2 : 1-3
#    case_ID = 3

    # Define the range of case tags considered
    case_tag_range = ('A1','A2','B1','B2')

    # Define the range of conditions considered under each case (4 conditions for A1 and B1, and 3 conditions under A2 and B2)
    case_cond_range_A1B1 = (1,2,3,4)
    case_cond_range_A2B2 = (1,2,3)

    # Define blade locations of interest and plotting styles
    blade_loc_id = (0, 8, -2, -1)
    blade_loc_tag = ('0.2 R','0.5 R','0.9 R','1.0 R')
    blade_loc_line = ('solid', 'dotted', 'dashed', 'dashdot')

    # Define the color range (one per model)
    model_colors = ('#15B01A', '#EF4026', '#F97306')
    model_marker = ('o','s','p')
    model_line = ('dashed', 'dotted', 'dashdot')
    model_tag = ('Pitt-Peters', 'Larsen-Madsen', 'Oye')
    qs_color = '#069AF3'

    # Close and clear all plots
    plt.close('all')

    for case_tag_i, case_tag in enumerate(case_tag_range):
        # Define the range of redced frequency for sinusodal perturbations
        if case_tag == 'A2' or case_tag == 'B2':
            freq_red_range = (0.05,0.10,0.15,0.20,0.25,0.30)    # Reduced frequency range required from assignment
            case_ID_range = case_cond_range_A2B2
        else:
            freq_red_range = [0.00]      # Dummy variable, as long as it is one value, step change responses will only be plotted once
            case_ID_range = case_cond_range_A1B1

        for case_ID_i, case_ID in enumerate(case_ID_range):
            # Plotting the five responses over time for the three models (and the 6 reduced frequencies for case A2 and B2)
            print('Plotting responses over time.')
            for freq_red_index, freq_red in enumerate(freq_red_range):
                print('k = ',freq_red)
                # Initialise the plots
                fig_a, (ax_a1,ax_a2,ax_a3,ax_a4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5))           # a: Induction factor
                fig_ct, (ax_ct1,ax_ct2,ax_ct3,ax_ct4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5))      # ct: Thrust coefficient
                fig_cq, (ax_cq1,ax_cq2,ax_cq3,ax_cq4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5))      # cq: Torque coefficient
                fig_aoa, (ax_aoa1,ax_aoa2,ax_aoa3,ax_aoa4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5)) # aoa: Angle of attack (alpha)
                fig_phi, (ax_phi1,ax_phi2,ax_phi3,ax_phi4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5)) # phi: Inflow angle

                # Loop over each model
                for i, model in enumerate(('pp', 'lm', 'oye')):
                    print(model)
                    if case_tag == 'A1':
                        r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('ct', *ct_steps[case_ID-1], None, model=model) # NB: Use case_ID-1 to comply with Python indexing convention
                    elif case_tag == 'A2':
                        r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('ct', *ct_sins[case_ID-1], freq_red, model=model)
                    elif case_tag == 'B1':
                        r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('u_inf', *u_inf_steps[case_ID-1], None, model=model)
                    elif case_tag == 'B2':
                        r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('u_inf', *u_inf_sins[case_ID-1], freq_red, model=model)
                    else:
                        print('Warning: Invalid case tag enterred.')

                    # OLD CODE; Used to run the simulations
                    # r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = turbine.ct_func(.5, .4, None, 10, 10, model=model)
                    # r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = turbine.ct_func(.5, .5, .3, 10, 10, model=model)
                    # r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = turbine.u_inf_func(1., .5, None, 10, 10, model=model)
                    # r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = turbine.u_inf_func(1., .5, .3, 10, 10, model=model)

                    # Assemble the plots
                    plot_combined_subplot('a [-]',ax_a1,ax_a2,ax_a3,ax_a4,t_list,a,a_qs,blade_loc_id,blade_loc_tag,model_colors[i],model_line[i],model_tag[i],i)
                    plot_combined_subplot('$C_t$ [-]',ax_ct1,ax_ct2,ax_ct3,ax_ct4,t_list,ctr,ctr_qs,blade_loc_id,blade_loc_tag,model_colors[i],model_line[i],model_tag[i],i)
                    plot_combined_subplot('$C_q$ [-]',ax_cq1,ax_cq2,ax_cq3,ax_cq4,t_list,cqr,cqr_qs,blade_loc_id,blade_loc_tag,model_colors[i],model_line[i],model_tag[i],i)
                    plot_combined_subplot('$\\alpha$ [deg]',ax_aoa1,ax_aoa2,ax_aoa3,ax_aoa4,t_list,alpha,alpha_qs,blade_loc_id,blade_loc_tag,model_colors[i],model_line[i],model_tag[i],i)
                    plot_combined_subplot('$\\phi$ [deg]',ax_phi1,ax_phi2,ax_phi3,ax_phi4,t_list,phi,phi_qs,blade_loc_id,blade_loc_tag,model_colors[i],model_line[i],model_tag[i],i)


                # Save the plots to .pdf
                plot_save_figure(fig_a,case_tag,case_ID,'a',freq_red,'Figures')
                plot_save_figure(fig_ct,case_tag,case_ID,'ct',freq_red,'Figures')
                plot_save_figure(fig_cq,case_tag,case_ID,'cq',freq_red,'Figures')
                plot_save_figure(fig_aoa,case_tag,case_ID,'aoa',freq_red,'Figures')
                plot_save_figure(fig_phi,case_tag,case_ID,'phi',freq_red,'Figures')

            # Plot each model on separate plot; plot range of reduced frequencies
            if case_tag == 'A2' or case_tag == 'B2':    # Loop over each model
                print('Plotting responses over time over range of reduced frequencies.')
                for i, model in enumerate(('pp', 'lm', 'oye')):
                    print(model)
                    # Initialise the plots
                    fig_a_k, (ax_ak1,ax_ak2,ax_ak3,ax_ak4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5))           # a: Induction factor
                    fig_ct_k, (ax_ctk1,ax_ctk2,ax_ctk3,ax_ctk4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5))      # ct: Thrust coefficient
                    fig_cq_k, (ax_cqk1,ax_cqk2,ax_cqk3,ax_cqk4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5))      # cq: Torque coefficient
                    fig_aoa_k, (ax_aoak1,ax_aoak2,ax_aoak3,ax_aoak4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5)) # aoa: Angle of attack (alpha)
                    fig_phi_k, (ax_phik1,ax_phik2,ax_phik3,ax_phik4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5)) # phi: Inflow angle

                    for freq_red_index, freq_red in enumerate(freq_red_range):  # Loop over each frequency
                        print(freq_red)
                        if case_tag == 'A1':
                            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('ct', *ct_steps[case_ID-1], None, model=model) # NB: Use case_ID-1 to comply with Python indexing convention
                        elif case_tag == 'A2':
                            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('ct', *ct_sins[case_ID-1], freq_red, model=model)
                        elif case_tag == 'B1':
                            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('u_inf', *u_inf_steps[case_ID-1], None, model=model)
                        elif case_tag == 'B2':
                            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('u_inf', *u_inf_sins[case_ID-1], freq_red, model=model)
                        else:
                            print('Warning: Invalid case tag enterred.')

                        # Assemble the plots
                        freq_red_grayscale = str(1 - (freq_red_index+1)/len(freq_red_range))    # Plot the different reduced frequency lines with a shade of grey [0 = Black; 1 = White]
                        plot_combined_subplot_red_freq('a [-]',ax_ak1,ax_ak2,ax_ak3,ax_ak4,t_list,a,a_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)
                        plot_combined_subplot_red_freq('$C_t$ [-]',ax_ctk1,ax_ctk2,ax_ctk3,ax_ctk4,t_list,ctr,ctr_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)
                        plot_combined_subplot_red_freq('$C_q$ [-]',ax_cqk1,ax_cqk2,ax_cqk3,ax_cqk4,t_list,cqr,cqr_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)
                        plot_combined_subplot_red_freq('$\\alpha$ [deg]',ax_aoak1,ax_aoak2,ax_aoak3,ax_aoak4,t_list,alpha,alpha_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)
                        plot_combined_subplot_red_freq('$\\phi$ [deg]',ax_phik1,ax_phik2,ax_phik3,ax_phik4,t_list,phi,phi_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)


                    # Save the plots to .pdf
                    plot_save_figure(fig_a_k,case_tag,case_ID,'a','_'+str(model),'Figures')
                    plot_save_figure(fig_ct_k,case_tag,case_ID,'ct','_'+str(model),'Figures')
                    plot_save_figure(fig_cq_k,case_tag,case_ID,'cq','_'+str(model),'Figures')
                    plot_save_figure(fig_aoa_k,case_tag,case_ID,'aoa','_'+str(model),'Figures')
                    plot_save_figure(fig_phi_k,case_tag,case_ID,'phi','_'+str(model),'Figures')

            # Plot each model on separate plot; plot range of reduced frequencies normalised
            if case_tag == 'A2' or case_tag == 'B2':    # Loop over each model
                print('Plotting responses over time over range of reduced frequencies normalised.')
                for i, model in enumerate(('pp', 'lm', 'oye')):
                    print(model)
                    # Initialise the plots
                    fig_a_k, (ax_ak1,ax_ak2,ax_ak3,ax_ak4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5))           # a: Induction factor
                    fig_ct_k, (ax_ctk1,ax_ctk2,ax_ctk3,ax_ctk4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5))      # ct: Thrust coefficient
                    fig_cq_k, (ax_cqk1,ax_cqk2,ax_cqk3,ax_cqk4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5))      # cq: Torque coefficient
                    fig_aoa_k, (ax_aoak1,ax_aoak2,ax_aoak3,ax_aoak4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5)) # aoa: Angle of attack (alpha)
                    fig_phi_k, (ax_phik1,ax_phik2,ax_phik3,ax_phik4) = plt.subplots(4, 1,sharex=True, figsize=(9, 5)) # phi: Inflow angle

                    for freq_red_index, freq_red in enumerate(freq_red_range):  # Loop over each frequency
                        print(freq_red)
                        if case_tag == 'A1':
                            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('ct', *ct_steps[case_ID-1], None, model=model) # NB: Use case_ID-1 to comply with Python indexing convention
                        elif case_tag == 'A2':
                            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('ct', *ct_sins[case_ID-1], freq_red, model=model)
                        elif case_tag == 'B1':
                            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('u_inf', *u_inf_steps[case_ID-1], None, model=model)
                        elif case_tag == 'B2':
                            r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('u_inf', *u_inf_sins[case_ID-1], freq_red, model=model)
                        else:
                            print('Warning: Invalid case tag enterred.')

                        # Assemble the plots
                        V0 = 10                             # [m/s] : Wind speed
                        R = 50                              # [m] : Blade length
                        norm_time = t_list*V0/R*freq_red    # [-] : Normalised time scale
                        freq_red_grayscale = str(1 - (freq_red_index+1)/len(freq_red_range))    # Plot the different reduced frequency lines with a shade of grey [0 = Black; 1 = White]
                        plot_combined_subplot_red_freq_norm('a [-]',ax_ak1,ax_ak2,ax_ak3,ax_ak4,norm_time,a,a_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)
                        plot_combined_subplot_red_freq_norm('$C_t$ [-]',ax_ctk1,ax_ctk2,ax_ctk3,ax_ctk4,norm_time,ctr,ctr_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)
                        plot_combined_subplot_red_freq_norm('$C_q$ [-]',ax_cqk1,ax_cqk2,ax_cqk3,ax_cqk4,norm_time,cqr,cqr_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)
                        plot_combined_subplot_red_freq_norm('$\\alpha$ [deg]',ax_aoak1,ax_aoak2,ax_aoak3,ax_aoak4,norm_time,alpha,alpha_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)
                        plot_combined_subplot_red_freq_norm('$\\phi$ [deg]',ax_phik1,ax_phik2,ax_phik3,ax_phik4,norm_time,phi,phi_qs,blade_loc_id,blade_loc_tag,freq_red_grayscale,'--','k = '+str(freq_red),freq_red_index)


                    # Save the plots to .pdf
                    plot_save_figure(fig_a_k,case_tag,case_ID,'a','_norm_'+str(model),'Figures')
                    plot_save_figure(fig_ct_k,case_tag,case_ID,'ct','_norm_'+str(model),'Figures')
                    plot_save_figure(fig_cq_k,case_tag,case_ID,'cq','_norm_'+str(model),'Figures')
                    plot_save_figure(fig_aoa_k,case_tag,case_ID,'aoa','_norm_'+str(model),'Figures')
                    plot_save_figure(fig_phi_k,case_tag,case_ID,'phi','_norm_'+str(model),'Figures')

            # Plotting responses over blade radial position
            if case_tag == 'A1' or case_tag == 'B1':
                print('Plotting responses over radial positions.')
                # Initialise the plots
                fig_a_elem, (ax_a1_elem,ax_a2_elem,ax_a3_elem) = plt.subplots(1, 3, sharey=True, figsize=(9, 5))            # a: Induction factor
                fig_ct_elem, (ax_ct1_elem,ax_ct2_elem,ax_ct3_elem) = plt.subplots(1, 3, sharey=True, figsize=(9, 5))        # ct: Thrust coefficient
                fig_cq_elem, (ax_cq1_elem,ax_cq2_elem,ax_cq3_elem) = plt.subplots(1, 3,sharey=True, figsize=(9, 5))         # cq: Torque coefficient
                fig_aoa_elem, (ax_aoa1_elem,ax_aoa2_elem,ax_aoa3_elem) = plt.subplots(1, 3,sharey=True, figsize=(9, 5))     # aoa: Angle of attack (alpha)
                fig_phi_elem, (ax_phi1_elem,ax_phi2_elem,ax_phi3_elem) = plt.subplots(1, 3,sharey=True, figsize=(9, 5))     # phi: Inflow angle


                # Loop over each model
                for model_i, model in enumerate(('pp', 'lm', 'oye')):
                    print(model)
                    if case_tag == 'A1':
                        r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('ct', *ct_steps[case_ID-1], None, model=model) # NB: Use case_ID-1 to comply with Python indexing convention
                    elif case_tag == 'B1':
                        r_list, t_list, ctr, cqr, a, alpha, phi, ctr_qs, cqr_qs, a_qs, alpha_qs, phi_qs = read_data('u_inf', *u_inf_steps[case_ID-1], None, model=model)
                    else:
                        print('Warning: Invalid case tag enterred.')

                    # Initialise the counter of time steps performed
                    time_step_counter = 0

                    # Loop over all time steps of interest
                    for row_of_a in range(a.shape[0]):
                        time_sampling = 10  # Define the time step for which to retrieve the responses
                        if t_list[row_of_a]%time_sampling == 0:
                            # Assemble the plots
                            time_step_grayscale = str((time_step_counter)/6)    # Divide by 6 to have a grayscale with depth, need to update this value if more time steps are considered
                            plot_combined_subplot_elem('a [-]',ax_a1_elem,ax_a2_elem,ax_a3_elem,r_list,a,a_qs,row_of_a,model_tag,time_step_grayscale,'--','t [s] = '+str(t_list[row_of_a]),qs_color,model_i,time_step_counter)
                            plot_combined_subplot_elem('$C_t$ [-]',ax_ct1_elem,ax_ct2_elem,ax_ct3_elem,r_list,ctr,ctr_qs,row_of_a,model_tag,time_step_grayscale,'--','t [s] = '+str(t_list[row_of_a]),qs_color,model_i,time_step_counter)
                            plot_combined_subplot_elem('$C_q$ [-]',ax_cq1_elem,ax_cq2_elem,ax_cq3_elem,r_list,cqr,cqr_qs,row_of_a,model_tag,time_step_grayscale,'--','t [s] = '+str(t_list[row_of_a]),qs_color,model_i,time_step_counter)
                            plot_combined_subplot_elem('$\\alpha$ [deg]',ax_aoa1_elem,ax_aoa2_elem,ax_aoa3_elem,r_list,alpha,alpha_qs,row_of_a,model_tag,time_step_grayscale,'--','t [s] = '+str(t_list[row_of_a]),qs_color,model_i,time_step_counter)
                            plot_combined_subplot_elem('$\\phi$ [deg]',ax_phi1_elem,ax_phi2_elem,ax_phi3_elem,r_list,phi,phi_qs,row_of_a,model_tag,time_step_grayscale,'--','t [s] = '+str(t_list[row_of_a]),qs_color,model_i,time_step_counter)

                            # Increment the counter of time steps performed
                            time_step_counter += 1

                # Save the plots to .pdf
                plot_save_figure_elem(fig_a_elem,case_tag,case_ID,'a','Figures')
                plot_save_figure_elem(fig_ct_elem,case_tag,case_ID,'ct','Figures')
                plot_save_figure_elem(fig_cq_elem,case_tag,case_ID,'cq','Figures')
                plot_save_figure_elem(fig_aoa_elem,case_tag,case_ID,'aoa','Figures')
                plot_save_figure_elem(fig_phi_elem,case_tag,case_ID,'phi','Figures')

    plt.show()