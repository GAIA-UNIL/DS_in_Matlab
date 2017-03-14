function [data_event_cond informed_nodes_id] = cond_dev_extract(window_type , point , hard_data1 , R_x , R_y)

sim_x1 = point(1,1);
sim_y1 = point(1,2);

search_radius = R_x; 

search_radius_x = R_x;
search_radius_y = R_y;

outer_search_radius = R_x;
inner_search_radius = R_y;


    if window_type == 1
    % 1.  Extracting conditioning data within a CIRCLE
    dist=sqrt((sim_x1-hard_data1(:,1)).^2+(sim_y1-hard_data1(:,2)).^2);   
    informed_nodes_id = find(and(dist <= search_radius , dist > 0));
    data_event_cond = hard_data1(informed_nodes_id,1:4);

    elseif window_type == 2
    % 2.  Extracting conditioning data within a SQUARE
        informed_nodes_id = find(and(and(hard_data1(:,1)>(sim_x1-floor(search_radius_x/2)),hard_data1(:,1)<(sim_x1+floor(search_radius_x/2))),...
        and(hard_data1(:,2)>(sim_y1-floor(search_radius_y/2)),hard_data1(:,2)<(sim_y1+floor(search_radius_y/2)))));
        data_event_cond = hard_data1(informed_nodes_id,1:4);
        
    elseif window_type == 3
    % 3.  Extracting from a DISC shape area
        dist_disc = sqrt((sim_x1-hard_data1(:,1)).^2+(sim_y1-hard_data1(:,2)).^2);
        informed_nodes_id = find(and(dist_disc>=inner_search_radius,dist_disc<=outer_search_radius));
        data_event_cond = hard_data1(informed_nodes_id,:);
    end
    