function [path]=sim_path_type(path_type,sim_size_x,sim_size_y,show_mode)

% Available Simulation Pathes:

% random: code = 1
% unilateral: code = 2
% unilateral with a starting point: code = 3
% bisector 45: code = 4
% bisector 135: code = 5
% circle starting from center: code = 6
% circle starting from outside: code = 7
if sim_size_x==sim_size_y
X=meshgrid(1:sim_size_x);
Y=meshgrid(1:sim_size_y);
Xg=reshape(X,[],1);
Yg=reshape(Y',[],1);
sim_grid(:,1)=Xg;
sim_grid(:,2)=Yg;
sim_grid(:,3)=1;
all_points=sim_grid;
else
counter=0;
for i=1:sim_size_x
for j=1:sim_size_y
counter=counter+1;
sim_grid(counter,1:3)=[i,j,1];
end
end 
all_points=sim_grid;
end

if path_type == 1
    %disp('random simulation path is selected')
    path_size = sim_size_x*sim_size_y;
    path = randperm(path_size);
    
elseif path_type == 2
    %disp('unilateral simulation path is selected')
    path_size = sim_size_x*sim_size_y;
    path = 1:path_size;
elseif path_type == 3
    %disp('unilateral with a starting point simulation path is selected')
    path_size = sim_size_x*sim_size_y;
    path = 1:path_size;
    no = randperm(path_size);
    seed_no = no(1);
    path = path+seed_no;
    id = path <= length(sim_grid);
    path(1:sum(id)) = path(id==1);
    path(sum(id)+1:end)=1:(length(path)-sum(id));
elseif path_type == 4
    %disp('bisector 45 simulation path is selected')
    counter=0;
    slope = 1;
    intercept = -1*(sim_size_y-1):(sim_size_y-1);
    intercept=intercept(intercept~=0==1)-0.1;
    for i=1:length(intercept)-1
        intercept(length(intercept)-i+1);
        intercept(length(intercept)-i);
        points_sign_above = slope.*(sim_grid(:,1)) - (sim_grid(:,2)-intercept(length(intercept)-i+1));
        points_sign_below = slope.*(sim_grid(:,1)) - (sim_grid(:,2)-intercept(length(intercept)-i));
        pts = find(and(points_sign_above>0 , points_sign_below<0));
        if ~isempty(pts)
        for i=1:length(pts)
            counter=counter+1;
            path(counter) = pts(i);
        end
    end
    end
    
elseif path_type == 5
    disp('bisector 135 simulation path is selected')
    counter=0;
    slope = -1;
    intercept = 1:(2*(sim_size_y)+1);
    intercept=intercept(intercept~=0==1)-0.1;
    for i=1:length(intercept)-1
        intercept(i+1);
        intercept(i);
        points_sign_above = slope.*(sim_grid(:,1)) - (sim_grid(:,2)-intercept(i+1));
        points_sign_below = slope.*(sim_grid(:,1)) - (sim_grid(:,2)-intercept(i));
        pts = find(and(points_sign_above>0 , points_sign_below<0));
        if ~isempty(pts)
        for i=1:length(pts)
            counter=counter+1;
            path(counter) = pts(i); %#ok<SAGROW>
        end
    end
    end
    
elseif path_type == 6
    disp('circular simulation path starting from center is selected')
    if sim_size_x/2-ceil(sim_size_x/2)~=0
       x0=sim_size_x/2;
       y0=sim_size_y/2;
       r=1:sim_size_x/2+3;
    else
       x0=sim_size_x/2+0.5;
       y0=sim_size_y/2+0.5;
       r=1:sim_size_x/2+1.5;
    end
                
       all_pts=ones(length(sim_grid),1);  
       counter=0;
    for i=1:length(r)
        radius = r(i);
        dist = (sim_grid(:,1)-x0).^2+(sim_grid(:,2)-y0).^2-radius^2;
        pts = find(and(dist<0 , all_pts~=-1));
        if ~isempty(pts)
           all_pts(pts)=-1;
           for j=1:length(pts)
               counter=counter+1;
               pts(j)
               path(counter) = pts(j);
            end
        end
     end

%--------------------------------------------------------------------------
elseif path_type == 7
    disp('circular simulation path starting from outside is selected')
    if sim_size_x/2-ceil(sim_size_x/2)~=0
       x0=sim_size_x/2;
       y0=sim_size_y/2;
       r=0:sim_size_x/2;
    else
       x0=sim_size_x/2+0.5;
       y0=sim_size_y/2+0.5;
       r=0:sim_size_x/2;
    end
        
       all_pts=ones(length(sim_grid),1);  
       counter=0;
    for i=1:length(r)-1
        radius = r(length(r)-i);
        dist = (sim_grid(:,1)-x0).^2+(sim_grid(:,2)-y0).^2-radius^2;
        pts = find(and(dist>0 , all_pts~=-1));
        if ~isempty(pts)
           all_pts(pts)=-1;
           for j=1:length(pts)
               counter=counter+1;
               pts(j)
               path(counter) = pts(j);
            end
        end
    end
        
end

if show_mode == 1;
scatter(sim_grid(:,1),sim_grid(:,2),1,'.');
hold on

for i=1:length(path)
    pt=sim_grid(path(i),:);
    shg;scatter(pt(1,1),pt(1,2),'o');
    xlim([1 sim_size_x]);ylim([1 sim_size_y]);
end    
    hold off
end
    
    
    