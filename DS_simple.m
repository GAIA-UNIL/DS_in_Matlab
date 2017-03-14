%DO NOT USE FR BENCHMARKING

%% Initialization
clear; home;
%Performs a multiple-points simulation by Direct Sampling.
%The training image is generated using objects.
%Parameters are:
simul_size = [80 80];       %size of the simulation: y x
ti_size = [100 100];        %size of the ti: y x
template = [9 9];           %size of the template: y x
fract_of_ti_to_scan = 0.1; %maximum fraction of the ti to scan
thr = 0.0;                  %threshold position (between 0 and 1)

%% Generation of a training image and display
ti=zeros(ti_size(1),ti_size(2)); [xx,yy]=meshgrid(1:ti_size(1),1:ti_size(2));
for i=1:35
    obj_param=[random('unif',1,100) random('unif',1,100),random('unif',4,8)];
    d=sqrt((xx-obj_param(1)).^2 + (yy-obj_param(2)).^2);
    ti(d <= obj_param(3))=1;
end

figure(1); clf; subplot(1,2,1); colormap gray;
imagesc(ti); title('training image'); axis equal tight xy; drawnow;
H=gca; set(H,'position',[0, 0.25, 0.5 0.5]); tic;

%% DS Simulation
%defining shifts related to the template size
yshift = floor(template(1)/2); xshift = floor(template(2)/2);
%reducing the size of the ti to avoid scanning outside of it
ti_size(1) = size(ti,1)-(2*yshift); ti_size(2) = size(ti,2)-(2*xshift);
%creating empty simulation with a wrapping of NaNs of the size "shift"
simul = nan(simul_size(1)+(2*yshift),simul_size(2)+(2*xshift));
%defining path in ti and in simulation
path_ti = randperm(ti_size(1)*ti_size(2));
path_sim = randperm(simul_size(1)*simul_size(2));
sizesim = size(path_sim,2);
progress = 0; tinod = 0;

%looping simulation nodes
for simnod = 1:sizesim

    progress_current=ceil((simnod*100)/sizesim);
    if progress_current>progress
        progress=progress_current; disp([num2str(progress),' % completed'])
    end
   
    %find node in the simulation grid
    xsim = ceil(path_sim(simnod)/simul_size(1));
    ysim = path_sim(simnod)-((xsim-1)*simul_size(1));
 
    %define the point and shifting it to avoid scanning NaNs
    point_sim = [ysim+yshift xsim+xshift];
 
    %define data event at simulated point
    data_event_sim = simul(point_sim(1)-yshift:point_sim(1)+yshift ,...
        point_sim(2)-xshift:point_sim(2)+xshift);
 
    %scan the ti
    mindist = inf;  %initial best distance
    tries = 0;      %counter of attempts
    max_scan = size(path_ti,2)*fract_of_ti_to_scan; %max number of scans
 
    %reducing the data event to its informed nodes
    no_data_indicator = isfinite(data_event_sim);
    data_event_sim = data_event_sim(no_data_indicator);
 
    while 1==1 %scan the ti
 
        tinod = tinod+1; tries = tries+1;
        %if arriving at the end of the path, restart from the beginning
        if (tinod > size(path_ti,2))
            tinod = 1;
        end
 
        %find the point in the ti
        xti = ceil(path_ti(tinod)/ti_size(1));
        yti = path_ti(tinod)-((xti-1)*ti_size(1));
 
        %find scanned point and data event
        point_ti = [yti+yshift xti+xshift];
        data_event_ti = ti(point_ti(1)-yshift:point_ti(1)+yshift ,...
            point_ti(2)-xshift:point_ti(2)+xshift);
 
        %if template is totally unknown, take the first value
        if sum(no_data_indicator(:)) == 0;
            simul(point_sim(1),point_sim(2)) = ti(point_ti(1),point_ti(2));
            break
        end
 
        %find the data event at this point in the ti
        data_event_ti = data_event_ti(no_data_indicator);
 
        %evaluate the distance between both data events
        %using another distance here allows using continuous variable
        distance=mean(data_event_sim~=data_event_ti);
 
        %if distance under threshold, the point is accepted
        if distance <= thr
            simul(point_sim(1),point_sim(2)) = ti(point_ti(1),point_ti(2));
            break
        else
            %check if the distance is under the minimum found so far
            if distance < mindist
                mindist = distance;
                bestpoint = point_ti;
            end
            %if max_scan nodes have been scanned, take the best point found
            if tries > max_scan
                simul(point_sim(1),point_sim(2))=ti(bestpoint(1),bestpoint(2));
                break
            end
        end
    end
end
 
% remove the "wrapping" of NaNs
toc; simul = simul(yshift+1:end-yshift,xshift+1:end-xshift);
 
%show simul
subplot(1,2,2); imagesc(simul); title('simulation'); axis equal tight xy;
H=gca; set(H,'position',[0.5, 0.25, 0.5*simul_size(1)/ti_size(1) ...
    0.5*simul_size(2)/ti_size(2)])