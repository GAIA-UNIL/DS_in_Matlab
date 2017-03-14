Paper: "Multiple-Point Geostatistical Simulation Using Bunch-Pasting Direct Sampling Method"



Authors: Hassan Rezaee, M.Sc. student of Mining Eng., at University of Tehran, Iran, email: h.rezaee@ut.ac.ir

         Gregoire Mariethoz, School of Civil and Environmental Engineering, University of New South Wales, 
         Sydney, Australia, e-mail: gregoire.mariethoz@minds.ch



Thank you for your interest in Bunch-DS algorithm:

This is a Matlab (2012) code written for the purpose of simulation using bunch-pasting multiple-point geostatistical
simulation. It is for both categorical and continuous training images. The inputs are as follow:



%   ti_size                     : the size of ti along x, y and z directions
%                                   for example [10 10 1] introduces a 2D ti
%   sim_dim                     : number of simulation grid nodes along x, y and z directions
%                                   for example [10 10 1] is 2D
%   search_radius_x &
%   search_radius_y             : in the case of rectangular window type, these are the dimension of window
%                                   if disc shape window is selected then: 
%                                   outer_search_radius = search_radius_x
%                                   inner_search_radius = search_radius_y
%                                   and if circular window is selected then:
%                                   radius = search_radius_x
%   fract_of_ti_to_scan         : the fraction of ti to be scanned; a number is to be chosen between [0-1]
%                                   for example 0.1 or 0.5
%   dist_threshold              : distance threshold between dev_ti and dev_cond, default = 0.05
%   max_no_nodes_in_data_event  : maximum number of nodes in data event
%                                   which is formed by both conditioning and previously-simulated data
%                                   a number of 30-40 would deliver realistic results
%                                   higher numbers exponentially decreases the computational efficiency of the algorithm
%                                 
%   window_type                 : consider the following codes to specify the search window type
%                                   1. circular
%                                   2. squared
%                                   3. disc Shape
%   bunch_size                  : if bunch_size > 0 the simulation mode is through bunch pasting manner
%                                   otherwise (bunch_size=0), pixel-based (traditional) DS is done
%                                   for example if "bunch_size = 1", then a number of (2*bunch_size+1)^2 = 9 nodes 
%                                   in the bunch are to be pasted
%   Condition_mode              : consider the following codes to do:
%                                   0. non-Conditional Simulation
%                                   1. conditional Simulation.
%   cond_data                   : If the conditional simulation is selected then the 
%                                   algorithm will ask you the name of conditioning data file which should be a
%                                   text file located in the directory of the matlab code
%                                   this data file should be in the form [x_coord y_coord z_coord variable]
%                                   otherwise the algorithm will consider cond_data = [-1 -1 -1 0] to do the non-conditional simulation
%                                   which carries no primary conditioning data while doing simulation
%   ti_path_mode                : use the following codes to determine the way through which the ti is scanned:
%                                 For this purpose the sub-function sim_path_mode_2d is employed.
%                                   1. random
%                                   3. unilateral with a starting point
%   No_Realizations             : number of realizations to be produced
%   post_proc_mode              : use the following codes to carry out post-processing steps:
%                                   0. no post-processing step is taken after the simulation
%                                   1. E-type, Inter Quartile Range (IQR) and Quantile maps are calculated after the simulation
%   output_parfile              : 1. the input parameter file is saved under the name of "Parfile.txt" in the directory
%   sim_path_mode_2d               : use the following code to specify the simulation path
%                                   1. random
%                                   2. unilateral
%                                   3. unilateral with a starting point
%                                   4. bisector 45
%                                   5. bisector 135
%                                   6. circle starting from center
%                                   7. circle starting from outside
%   simulation_show_mode        : to illustrate the realization while simulation consider 1 otherwise any other code


All parameters are explained enough to run the algorithm. This algorithm can be easily changed into the original DS
algorithm by considering bunch_size=0 which pastes only one node at a time to the simulation grid.

The output of this function is one sole matrice which is the final realization produced by the Bunch_DS code. Any extra 
output can be easily got through changing the code in the relevant way. 

All parts of the code are tried to be self-explained enough which is done by the adopting the meaningful names for the
variables employed, however extra explanations are associated in different steps of the code.