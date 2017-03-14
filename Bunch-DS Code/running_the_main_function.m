clc
clear all
clf

ti_size = [64 64]; %   TI_A size=[64 64]; TI_A size=[75 88];
sim_dim = [50 50 1];
search_radius_x = 15;
search_radius_y = 15;
fract_of_ti_to_scan = 0.1;
dist_threshold = 0.0;
max_no_nodes_in_data_event = 40;
window_type = 2; % circular
bunch_size = 1; % 25 nodes in the bunch
Condition_mode = 1;
ti_path_mode = 3; % unilateral with starting point
sim_path_mode = 2; % unilateral
No_Realizations = 1;
post_proc_mode = 0;
output_parfile = 1;
simulation_show_mode = 1; % show while doing the simulation

ti_name = 'TI_A'; % or it can be TI_D for the second example
cond_data_name = 'TI_A_Cond'; % or it can be TI_D_Cond for the second example
                             
[Realizations] = Bunch_DS (ti_name,ti_size,cond_data_name,sim_dim,search_radius_x,search_radius_y,fract_of_ti_to_scan,...
                           dist_threshold,max_no_nodes_in_data_event,window_type,bunch_size ,...
                           Condition_mode,ti_path_mode,sim_path_mode,No_Realizations ,...
                           post_proc_mode,output_parfile,simulation_show_mode);
                       