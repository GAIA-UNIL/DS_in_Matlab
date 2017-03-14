%---------------------------------------------------------------------------------------------
% Authors: Hassan Rezaee, Gregoire Mariethoz
% Paper: "Multiple-Point Geostatistical Simulation Using Bunch-Pasting Direct Sampling Method"
%---------------------------------------------------------------------------------------------

function [Realizations] = Bunch_DS (ti_name,ti_size,cond_data_name,sim_dim,search_radius_x,search_radius_y,fract_of_ti_to_scan,...
                                 dist_threshold,max_no_nodes_in_data_event,window_type,bunch_size ,...
                                 Condition_mode,ti_path_mode,sim_path_mode,No_Realizations ,...
                                 post_proc_mode,output_parfile,simulation_show_mode)
                             
%--------------------------------------------------------------------------------------------------------------------
% Conditional/Non-conditional Simulation of Categorical/Continuous Variables Using Bunch-Pasting/Pixel-Based DS Method
%--------------------------------------------------------------------------------------------------------------------
%
% USE: [Realizations] = Bunch_DS (ti_size,sim_dim,search_radius_x,search_radius_y,fract_of_ti_to_scan,...
%                                  dist_threshold,max_no_nodes_in_data_event,window_type,bunch_size ,...
%                                  Condition_mode,ti_path_mode,sim_path_mode,No_Realizations ,...
%                                  post_proc_mode,output_parfile,simulation_show_mode)
%
% INPUT:
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
%-----------------------------------------------------------------------------------------------------------
% OUPUT: 
%   Realizations                : final realizations produced by this method in the form [x_coord y_coord variable]

%-----------------------------------------------------------------------------------------------------------
% This program uses the following subroutines:
%   find_node_in_grid_2d.m       : find a node's id in a grid
%   dist_cal_2d.m                : calculate the distance between data events from simulation grid and ti 
%   dlmcell.m                    : write cell array to text file (Roland Pfister, Version: 01.06.2010)
%   cond_dev_extract_2d.m        : extracts the data event from the neighborhood of the node under simulation
%   sim_post_proc.m              : calculates the post-processing maps of E-Type, Interquartile and Quantile maps

%-----------------------------------------------------------------------------------------------------------
% An Example of Input Parameters for Simulation:

% sim_dim = [50 50 1];
% search_radius_x = 8;
% search_radius_y = 8;
% fract_of_ti_to_scan = 0.3;
% dist_threshold = 0.05;
% max_no_nodes_in_data_event = 35;
% window_type = 2;
% bunch_size = 2;
% Condition_mode = 1;
% post_proc_mode = 0;
% No_Realizations = 1;
% output_parfile = 1;
% simulation_show_mode = 0;
% sim_path_mode = 2;

% that forms the following input for the program:

% [Realizations] = DSSIM_2D ([64 64 1],[25 25 1],8,6,0.5,0.05,35,2,2,1,3,2,1,0,0,1)
% please enter the name of ti data file that is to be a text file (put the name between Single Quotes) : 'TI_A'
% please enter the name of conditioning data file that is to be a text file (put the name between Single Quotes) : 'TI_A_Cond'

%-----------------------------------------------------------------------------------------------------------

% the grid of the training image [x_coord y_coord z_coord variable] should be introduced to the algorithm
%T = input('please enter the name of ti data file that is to be a text file (put the name between Single Quotes) : ');
ti_grid = textread([ti_name '.txt']);
ti=zeros(ti_size(1),ti_size(2));
counter=0;
for i=1:ti_size(1)
    for j=1:ti_size(2)
        counter=counter+1;
        ti(i,j) = ti_grid(counter,3);
    end
end

% check for the type of variable; it is supposed that the number of facies in categorical TI is less than 10!
[B]=sort(ti_grid(:,3),1);
No_facies=0;
for i=2:length(B)
    if (B(i-1)~=B(i))
        No_facies = No_facies+1;
    end
end
if No_facies > 10
    disp('Continuous Variable is selected')
    pause(2)
    Distance_type=2;
else
    disp('Categorical Variable is selected')
    pause(2)
    Distance_type=1;
end

%% creating empty simulation grid
simulation = nan(sim_dim(1)+2*floor(search_radius_x/2),...
                 sim_dim(2)+2*floor(search_radius_y/2),sim_dim(3));

% Creating an outer shell around Simulation Grid
counter = 0;
sim_grid=zeros((sim_dim(1)+2*floor(search_radius_x/2))*(sim_dim(2)+2*floor(search_radius_y/2)),4);
for i=1:sim_dim(1)+2*floor(search_radius_x/2)
    for j=1:sim_dim(2)+2*floor(search_radius_y/2)
        for k=1:sim_dim(3)
            counter=counter+1;
            sim_grid(counter,1:3) = [i j k];
            sim_grid(counter,4) = simulation(i,j,k);
        end
    end
end
x_active = and(sim_grid(:,1)>floor(search_radius_x/2),sim_grid(:,1)<=max(sim_grid(:,1))-floor(search_radius_x/2));
y_active = and(sim_grid(:,2)>floor(search_radius_y/2),sim_grid(:,2)<=max(sim_grid(:,2))-floor(search_radius_y/2));
%f=and(x_active==1,y_active==1);
%active = sim_grid(f==1,:);

counter=0;
for i=1:sim_dim(1)+2*floor(search_radius_x/2)
    for j=1:sim_dim(2)+2*floor(search_radius_y/2)
        for k=1:sim_dim(3)
            counter=counter+1;
            sim_grid(counter,1:3) = [i j k];
            sim_grid(counter,4) = simulation(i,j,k);
        end
    end
end
sim_grid(or((x_active==0),y_active==0),4)=nan;
SIM_GRID = sim_grid;

%% Realizations over sim_grid are to be produced

% extra variables to show the contribution of different states that a node can be simulated
c_empty_dv=0;
c_count_thre=0;
c_distance=0;

Realization=zeros(size(sim_grid,1),No_Realizations);
for real=1:No_Realizations
    
    active = sim_grid(and((x_active==1),y_active==1),:);
    active(:,5)=0;
    diff=size(sim_grid,1)-sim_dim(1)*sim_dim(2);

% CONDITIONAL SIMULATION
if Condition_mode==1
    
    % 1. First to read the conditioning data file
    if real==1
    %cond_data_name = input('please enter the name of conditioning data file that is to be a text file (put the name between Single Quotes) : ');
    cond_data = textread([cond_data_name '.txt']);
    Hard=cond_data;
    elseif real>1
        cond_data=Hard;
    end
    % 2. Due to the outer shell around sim_grid, hard data are re-located.
    cond_data(:,1) = cond_data(:,1)+floor(search_radius_x/2);
    cond_data(:,2) = cond_data(:,2)+floor(search_radius_y/2);
    
    % 3. Just the "cond data" within the simulation grid should be identified
    a = and(cond_data(:,1)<=sim_dim(1)+floor(search_radius_x/2),...
            cond_data(:,2)<=sim_dim(2)+floor(search_radius_y/2));
    cond_data = cond_data(a,:);
    
    % 4. Inserting cond data in simulation grid
    for i=1:size(cond_data,1)
        a= and(sim_grid(:,1)==cond_data(i,1),sim_grid(:,2)==cond_data(i,2));
        sim_grid(a,4)=cond_data(i,4);
    end
    %clear sim_path
    % 5. Nodes are visited along a pre-determined path: random or unilateral

    sim_path(:,1)=sim_path_type_2d(sim_path_mode,sim_dim(1),sim_dim(2),0);

    
    % 6. Eliminating hard data point from the simulation path
    for i=1:size(sim_path,1)
        a=find(and(active(i,1)==cond_data(:,1),active(i,2)==cond_data(:,2)), 1);
       
        if ~isempty(a)
        sim_path(i,1)=-1;
        end
        
    end
    sim_path(sim_path == -1) = [];
    hard_data = cond_data;
    
% NON-CONDITIONAL SIMULATION    
elseif Condition_mode==0
    
    cond_data = [-1 -1 -1 0];
    
    % Nodes are visited along a pre-determined path mentioned first
    
    sim_path(:,1)=sim_path_type_2d(sim_path_mode,sim_dim(1),sim_dim(2),0);
    
end

tic

while 1==1
    
% for every simulation node do
for i = 1:size(sim_path,1)
    %Active_size=size(active,1);
    
    % The coordinates of the nodes under simulation
    sim_x=active(sim_path(i,1),1);
    sim_y=active(sim_path(i,1),2);
    sim_node_id = and(sim_grid(:,1)==sim_x , sim_grid(:,2)==sim_y);
    
if isnan(sim_grid(sim_node_id==1,4)) == 1
    
    % To show the progress of simulation
    disp(['Realization No.  ' num2str(real) '  is being produced   '...
           num2str(100-((sum(isnan(sim_grid(:,4)))-diff)/(sim_dim(1)*sim_dim(2)))*100),'  % completed']);

    % To extract the bunch-shape space to which the bunch from TI is to be pasted  

    sim_bunch_counter = 0;
    points_sim=zeros((2*bunch_size+1)^2,2);
    bunch_ids_sim_coord=zeros((2*bunch_size+1)^2,2);
    
        for x_id_sim = -bunch_size:bunch_size
            for y_id_sim = -bunch_size:bunch_size
                p_x=sim_x + x_id_sim; p_y=sim_y + y_id_sim;
                
                if and(isnan(sim_grid(and(sim_grid(:,1) == p_x , sim_grid(:,2) == p_y),4)) == 1,...
                       and(and(p_x <= max(active(:,1)) , p_x >= (min(active(:,1)))),...
                           and(p_y <= max(active(:,2)) , p_y >= (min(active(:,2))))))
                       
                    sim_bunch_counter = sim_bunch_counter+1;
                    points_sim(sim_bunch_counter,1:2)=[p_x   p_y];
                    bunch_ids_sim_coord(sim_bunch_counter,1:2) = [x_id_sim   y_id_sim];
                    
                end
            end
        end
        bunch_coord_sim = points_sim;
 
        s=max(sim_grid(:,2))-min(sim_grid(:,2))+1; 
        
        bunch_ids_sim = (s).*(points_sim(1:sim_bunch_counter,1)-1)+...
                         (points_sim(1:sim_bunch_counter,2));
                    
        bunch_ids_sim = bunch_ids_sim(bunch_ids_sim<=size(sim_grid,1));

% conditioing data from both primary cond data and previously simulated data are extracted
[data_event_cond informed_nodes_id] = cond_dev_extract_2d(window_type,[sim_x sim_y],cond_data,search_radius_x,search_radius_y);

[c] = size(data_event_cond,1);

% maximum number of nodes in search window. In circular case the n closest nodes are considered for further calculations
if c > max_no_nodes_in_data_event
    dist_to_central_node = sqrt((data_event_cond(:,1)-sim_x).^2+(data_event_cond(:,2)-sim_y).^2);
    [~,IX] = sort(dist_to_central_node);
    data_event_cond = cond_data(informed_nodes_id(IX(1:max_no_nodes_in_data_event)),1:4);
end

[c] = size(data_event_cond,1);

    % data-event is defined based on the lag vectors h
    h_x = data_event_cond(:,1)-sim_x;  
    h_y = data_event_cond(:,2)-sim_y;   
    
%% Define a search window in the TI grid by using the dimensions of the data event.
    h_x_pos=h_x(h_x>0);h_x_neg=h_x(h_x<0);
    h_y_pos=h_y(h_y>0);h_y_neg=h_y(h_y<0);
    
    if isempty(h_x_pos)==1
        h_x_pos=0;
    end
    if isempty(h_x_neg)==1
        h_x_neg=0;
    end
    if isempty(h_y_pos)==1
        h_y_pos=0;
    end
    if isempty(h_y_neg)==1
        h_y_neg=0;
    end
    
    x_max = max(h_x_pos); x_min = abs(min(h_x_neg));
    y_max = max(h_y_pos); y_min = abs(min(h_y_neg));
    
    if isempty(x_max)==1
        x_max=0;
    end
    if isempty(x_min)==1
        x_min=0;
    end
    if isempty(y_max)==1
        y_max=0;
    end
    if isempty(y_min)==1
        y_min=0;
    end

    ti_active_id_x = and(ti_grid(:,1) > (min(ti_grid(:,1))+x_min) , ti_grid(:,1) < (max(ti_grid(:,1))-x_max));
    ti_active_id_y = and(ti_grid(:,2) > (min(ti_grid(:,2))+y_min) , ti_grid(:,2) < (max(ti_grid(:,2))-y_max));
    
    ti_active_id = and(ti_active_id_x==1 , ti_active_id_y==1);
    
    ti_active = ti_grid(ti_active_id==1,:);
    
    x_edge=min(ti_active(:,1));y_edge=min(ti_active(:,2));
    
    s1=reshape(ti_active(:,3),1+max(ti_active(:,2))-min(ti_active(:,2)),1+max(ti_active(:,1))-min(ti_active(:,1)));
    
   % ti_nan=s1;

%% 2. scan the TI for matching pattern

    [ti_path] = sim_path_type_2d(ti_path_mode,ti_size(1),ti_size(2),0);

    % matching criteria is reached when Distance<=distance threshold

    count_thre = floor(size(ti_path,2)*fract_of_ti_to_scan);
    counter = 0;
    mindist = inf;

        while 1==1
            
        counter = counter+1;
        
        if ti_path(counter)>size(ti_active,1)
           ti_path=1:size(ti_active,1);
        end
        
        ti_node_id = ti_path(counter);
        point_ti(1,1:2) = ti_active(ti_node_id,1:2);

        % Extracting a bunch of nodes around the central node in the ti data event
        
        clear Points
        
        Points(:,1)=point_ti(1,1) + bunch_ids_sim_coord(1:sim_bunch_counter,1);
        Points(:,2)=point_ti(1,2) + bunch_ids_sim_coord(1:sim_bunch_counter,2);
         
        nodes_id_ti = (size(s1,1)).*(Points(:,1)-x_edge)+(Points(:,2)-y_edge)+1;        

        bunch_ids_ti = nodes_id_ti(and(nodes_id_ti <= size(ti_active,1) , nodes_id_ti>0),1);
        ti_bunch_counter=size(bunch_ids_ti,1);
        
        if ti_bunch_counter==0
           break
        end
    
        min_size=min(ti_bunch_counter , sim_bunch_counter);
    
        ti_bunch = ti_active(bunch_ids_ti(:,1),3); 
        
        % the case when data event from simulation grid is empty: a random pattern is pasted
            if or(c < 1 , sim_bunch_counter == 0)
               c_empty_dv = c_empty_dv+1;
               sim_grid(bunch_ids_sim(1:min_size,1),4) = ti_bunch(1:min_size,1);
               s=max(active(:,2))-min(active(:,2))+1;
                          
               nodes_id_sim = (s).*(sim_grid(bunch_ids_sim,1)-(floor(search_radius_x/2)+1))+...
                              (sim_grid(bunch_ids_sim,2)-(floor(search_radius_y/2)));
               
               nodes_id_sim = nodes_id_sim((nodes_id_sim)>0,1);
               active(nodes_id_sim,5) = -1;
               break
            end

            data_event_ti_coord(1:c,1) = ti_active(ti_node_id,1)+h_x;
            data_event_ti_coord(1:c,2) = ti_active(ti_node_id,2)+h_y;
                
            % based on the lag vectors defined above the data_event_ti is made then
            data_event_ti = find_node_in_grid_2d(data_event_ti_coord(1:c,:) , ti , c);
           
            % evaluating the distance between the data event in the simulation and that of ti
            Distance=dist_cal(data_event_cond , data_event_ti , c , Distance_type);  
            
            % The minimun distance calculated so fat is updated
            if Distance < mindist
                mindist = Distance;
                best_bunch = ti_bunch;
                %best_point = ti_grid(ti_path(counter),3);
            end
            
            % the case when the pre-specified fraction of TI is scanned and the matching pattern is not found yet
            if counter > count_thre
               c_count_thre = c_count_thre+1;
               s=max(active(:,2))-min(active(:,2))+1;
               id = (s).*(bunch_coord_sim(1:size(best_bunch,1),1)-(floor(search_radius_x/2)+1))+...
                        (bunch_coord_sim(1:size(best_bunch,1),2)-(floor(search_radius_y/2)));
               
               id=id(id<=size(active,1));
               active(id,5) = -1;
               sim_grid(bunch_ids_sim(1:size(best_bunch,1),1),4) = best_bunch(:,1);
            break
            
            % the case when the matching pattern is found that satisfies the distance threshold
            elseif Distance <= dist_threshold
               c_distance = c_distance+1;
               s=max(active(:,2))-min(active(:,2))+1;
               id2 = (s).*(bunch_coord_sim(1:size(ti_bunch,1),1)-(floor(search_radius_x/2)+1))+...
                          (bunch_coord_sim(1:size(ti_bunch,1),2)-(floor(search_radius_y/2)));
               
               id2=id2(and(id2<=size(active,1),id2>0));
               
               active(id2,5) = -1;
               sim_grid(bunch_ids_sim(1:size(ti_bunch,1),1),4) = ti_bunch(:,1);
            break
            end

        end
        
    % with the addittion of simulated nodes the conditioning data is increased each time
    cond_data = sim_grid(isnan(sim_grid(:,4))==0,:);
    
    % The illustration of realization while doing the simulation
    if simulation_show_mode == 1
    sim = reshape(sim_grid(:,4),sim_dim(2)+2*floor(search_radius_y/2),...
                  sim_dim(1)+2*floor(search_radius_x/2));
    sim(isnan(sim)==1) = -2;
    shg
    imagesc(sim');axis xy;axis equal;xlim([1 sim_dim(1)+2*floor(search_radius_x/2)]);
    ylim([1 sim_dim(2)+2*floor(search_radius_y/2)]);
    title('Bunch DS simulation')
    end
end
end

    active=active(active(:,5)~=-1,:);
    sim = reshape(sim_grid(:,4),sim_dim(1)+2*floor(search_radius_x/2),...
                  sim_dim(2)+2*floor(search_radius_y/2));
    sim(isnan(sim)==1) = -2;
    
% checks for the criteria that determines the end of the simulation 
if or( sum(active(:,5)~=-1)==0 , sum((isnan(sim_grid(:,4)))==0) == (sim_dim(1)*sim_dim(2)) )
    toc
    disp(['--------------------------','END OF SIMULATION USING BUNCH-DS-Rezaee et al., 2012'])
    break
end

clear sim_path

end

Realization(:,real)=sim_grid(:,4);
clear sim_grid
sim_grid=SIM_GRID;

end

%% Illustration of outputs
a=reshape(sim,(sim_dim(1)+2*floor(search_radius_x/2))*(sim_dim(2)+2*floor(search_radius_y/2)),1);
b=a~=-2;
s=a(b==1);
SIM=reshape(s,sim_dim(2),sim_dim(1));
imagesc(SIM);axis xy;axis equal;xlim([1 sim_dim(1)]);ylim([1 sim_dim(2)]);
xlabel('X Coordinate');ylabel('Y Coordinate')
title('Bunch DS simulation');

if Distance_type == 1
    colormap('gray')
end

if Condition_mode == 1
    hold on 
    hard_data = Hard;
    scatter(hard_data(:,1),hard_data(:,2) , 'o' , 'MarkerEdgeColor','k');
    legend('Conditioning Data', 'Location','SouthOutside');
end

nan_id = isnan(Realization(:,1))==0;
Realizations(:,1:2)=sim_grid(nan_id==1,[1 2]);

for r=1:No_Realizations
    Realizations(:,r+2) = Realization(nan_id==1,r);
end
Realizations(:,1)=Realizations(:,1)-(min(Realizations(:,1))-1);
Realizations(:,2)=Realizations(:,2)-(min(Realizations(:,2))-1);

disp(' ')
disp([num2str(c_empty_dv/(c_empty_dv+c_count_thre+c_distance)*100) ' % of nodes are simulated by pasting a random pattern from TI'])
disp([num2str(c_count_thre/(c_empty_dv+c_count_thre+c_distance)*100) ' % of nodes are simulated by drawing the pattern of minimum distance found so far'])
disp([num2str(c_distance/(c_empty_dv+c_count_thre+c_distance)*100) ' % of nodes are simulated by pasting the matching pattern satisfying the distance threshold'])


%% Carrying out the post processing steps on the realizations produced
if and(post_proc_mode == 1 , No_Realizations==1)
    errordlg('The number of realization should be more than 1 to produce post-processing maps');
end
if and(post_proc_mode == 1 , No_Realizations>1)

[post] = sim_post_proc(Realization , real , 3);
a1=isnan(post(:,1))==0;
post=post(a1==1,1:3);

etype = reshape(post(:,1), sim_dim(1) , sim_dim(2));
imagesc(etype);axis xy;axis equal;xlim([1 sim_dim(2)]);ylim([1 sim_dim(1)])
hold on
scatter(hard_data(:,1),hard_data(:,2),'o');title('E-type')

IQR = reshape(post (:,2), sim_dim(1) , sim_dim(2));
figure;imagesc(IQR);axis xy;axis equal;xlim([1 sim_dim(2)]);ylim([1 sim_dim(1)])
hold on
scatter(hard_data(:,1),hard_data(:,2),'o');title('IQR')

Quantile = reshape(post (:,3), sim_dim(1) , sim_dim(2));
figure;imagesc(Quantile);axis xy;axis equal;xlim([1 sim_dim(2)]);ylim([1 sim_dim(1)])
hold on
scatter(hard_data(:,1),hard_data(:,2),'o');title('Quantile')

end

%% To save the input parameter file
parfile=cell(11,4);
if output_parfile == 1
                             
    parfile(1,1) = {'Simulation Grid Dimension (x,y,z)'};
    parfile(1,2:4) = num2cell(sim_dim);  %size of the simulation: y x z
    parfile(2,1) = {'Search Radii (X , Y)'};
    parfile(2,2) = num2cell(search_radius_x); % In the case of rectangular window 
    parfile(2,3) = num2cell(search_radius_y); % In the case of rectangular window 
    parfile(3,1) = {'Fraction of ti which is scanned'};
    parfile(3,2) = num2cell(fract_of_ti_to_scan);  %fraction of ti to be scanned
    parfile(4,1) = {'Distance type'};
    if Distance_type==1
        parfile(4,2) = {'Categorical variable'};
    elseif Distance_type==2
        parfile(4,2) = {'Continuous variable'};
    end
    
    parfile(5,1) = {'Distance threshold'};
    parfile(5,2) = num2cell(dist_threshold);
    parfile(6,1) = {'Maximum number of nodes in dataevent'};
    parfile(6,2) = num2cell(max_no_nodes_in_data_event);
    
    parfile(7,1) = {'Window type'};
    if window_type==1
    parfile(7,2) = {'Circular'};
    elseif window_type==2 
        parfile(7,2) = {'Squared'};
    elseif window_type == 3
        parfile(7,2) = {'Disc Shape'};
    end
    parfile(8,1) = {'Bunch size'};
    parfile(8,2) = num2cell(bunch_size); 
    
     parfile(9,1) = {'Conditioning mode'};
    if Condition_mode == 1
    parfile(9,2) = {'Conditional Simulation'};
    elseif Condition_mode == 0
        parfile(9,2) = {'Non-Conditional Simulation'};
    end
    
    parfile(10,1) = {'Simulation path mode'};
    if sim_path_mode == 1
        parfile(10,2) = {'Nodes are visited "Randomly"'};
            elseif sim_path_mode == 2
        parfile(10,2) = {'The Uni-Lateral or Raster path is adopted'};
            elseif sim_path_mode == 3
        parfile(10,2) = {'The Uni-Lateral with a starting point is adopted'};
            elseif sim_path_mode == 4
        parfile(10,2) = {'The bisector 45 path is adopted'};
            elseif sim_path_mode == 5
        parfile(10,2) = {'The bisector 135 path is adopted'};
            elseif sim_path_mode == 6
        parfile(10,2) = {'The circle starting from center path is adopted'};
            elseif sim_path_mode == 7
        parfile(10,2) = {'The circle starting from outside path is adopted'};
    end
    
    parfile(11,1) = {'Number of realizations to be produced'};
    parfile(11,2) = num2cell(No_Realizations); 
    
end
dlmcell('Parfile.txt',parfile)
