%DO NOT USE FR BENCHMARKING

%Performs a multiple-point realizations using Direct Sampling (DS).
%The DS method is described here:
% Mariethoz G, Renard P, Straubhaar J (2010). The Direct Sampling method to perform 
% multiple-point geostatistical simulations. Water resources Research 46(W11536), 
% doi:10.1029/2008WR007621. 
%Parameterization of the method is described here:
% Meerschman E, Pirot G, Mariethoz G, Straubhaar J, Van Merivenne M, Renard P (2013). 
% A Practical Guide to Performing Multiple-Point Geostatistical Simulations 
% with the Direct Sampling Algorithm, Computers & Geosciences 52, 307-324, 
% doi: 10.1016/j.cageo.2012.09.019. 

clear;home

% Simulation parameters
simul_size = [50 50 1];%size of the simulation: y x z
search_radius = 20;    %maximum extension of the data events
n = 20;                %maximum number of points in the data event
f = 0.5;               %maximum fraction of scanned training image
t = 0.05;               %distance threshold (between 0 and 1)
distance_type=1;       %type of variable (0=categorical, 1= continuous with 
                       %Manhattan distance, 2= continuous with euclidean distance)

%DO NOT USE FR BENCHMARKING

%loading training image
tifile='ti.txt';
%loading conditioning data
datafile = 'cond.txt';

ti=load(tifile);
data = load(datafile);
ti_size = [size(ti,1) size(ti,2) size(ti,3)]; %to make sure size_ti is a vector of size 3
sizeyxti = ti_size(1)*ti_size(2);
sizeyxsim = simul_size(1)*simul_size(2);
tic

%defining path in simulation
path_sim = randperm(simul_size(1)*simul_size(2)*simul_size(3));

%creating empty simulation
simul = nan(simul_size(1)*simul_size(2)*simul_size(3),1);

%inserting conditioning data in realization
for i = 1:size(data,1)
    %finding ids of data points
    dataid = (simul_size(1).*simul_size(2)).*(data(i,3)-1)+simul_size(1).*(data(i,2)-1)+data(i,1);
    simul(dataid) = data(i,4);
    %eliminating data points in the path
    ind = find(path_sim==dataid);
    path_sim(ind) = [];
end

bestmin = zeros(size(path_sim,2),1);
nbtries = zeros(size(path_sim,2),1);
thr = inf;

%looping simulation nodes
for simnod = 1:size(path_sim,2)

    %finding where we are in the simulation
    [ysim,xsim,zsim] = findcoord(path_sim(simnod),simul_size(1),simul_size(2));

    %finding distance with all previously simulated points
    all_sim_pts = find(isfinite(simul)==1);
    [yt,xt,zt] = findcoord(all_sim_pts,simul_size(1),simul_size(2));
    d = sqrt((ysim-yt).^2+(xsim-xt).^2+(zsim-zt).^2);

    %sorting
    [d,s] = sort(d);
    all_sim_pts = all_sim_pts(s);

    %checking how many points are within outer search radius
    nb_in_search_radius = sum(d<=search_radius);

    %taking the closest n or less
    if nb_in_search_radius < n
        informed_nodes_counter = nb_in_search_radius;
    else
        informed_nodes_counter = n;
    end

    if informed_nodes_counter == 0
        vect_to_informed_nodesy = 0;
        vect_to_informed_nodesx = 0;
        vect_to_informed_nodesz = 0;
        data_event_sim = 0;
    else
        %finding vectors to each point and data event
        vect_to_informed_nodesx=zeros(1,informed_nodes_counter);
        vect_to_informed_nodesy=zeros(1,informed_nodes_counter);
        vect_to_informed_nodesz=zeros(1,informed_nodes_counter);
        data_event_sim=zeros(1,informed_nodes_counter);
        for i = 1:informed_nodes_counter
            [yt,xt,zt] = findcoord(all_sim_pts(i),simul_size(1),simul_size(2));
            vect_to_informed_nodesy(i) = yt-ysim;
            vect_to_informed_nodesx(i) = xt-xsim;
            vect_to_informed_nodesz(i) = zt-zsim;
            
            %indice of the data point in the simul
            id = (sizeyxsim).*(zt-1)+simul_size(1).*(xt-1)+yt;
            data_event_sim(i) = simul(id);
        end
    end

    %defining dark zone related to the extension of the data event
    yshiftmax = ti_size(1)-max(vect_to_informed_nodesy);
    yshiftmin = -min(vect_to_informed_nodesy);
    xshiftmax = ti_size(2)-max(vect_to_informed_nodesx);
    xshiftmin = -min(vect_to_informed_nodesx);
    zshiftmax = ti_size(3)-max(vect_to_informed_nodesz);
    zshiftmin = -min(vect_to_informed_nodesz);

    %taking care not to go out of the training image
    if yshiftmax > ti_size(1)
        yshiftmax = ti_size(1);
    end
    if yshiftmin < 0
        yshiftmin = 0;
    end
    if xshiftmax > ti_size(2)
        xshiftmax = ti_size(2);
    end
    if xshiftmin < 0
        xshiftmin = 0;
    end
    if zshiftmax > ti_size(3)
        zshiftmax = ti_size(3);
    end
    if zshiftmin < 0
        zshiftmin = 0;
    end

    %defining a random path
    window = [(yshiftmax-yshiftmin) (xshiftmax-xshiftmin) (zshiftmax-zshiftmin)];
    window_size = window(1)*window(2)*window(3);
    window_sizeyx = window(1)*window(2);
    path_ti = randperm(window_size);

    %scanning ti
    mindist = inf;  %initial best distance is set to inf. Updated with every best distance encountered
    nb_of_tries = ceil(window_size*f); %number of tries in the ti

    for i = 1:nb_of_tries

        %finding coordinates of the point in the random path and shifting
        zti = ceil(path_ti(i)./window_sizeyx);
        id2d = (path_ti(i)-(window_sizeyx.*(zti-1)));
        xti = ceil(id2d./window(1));
        yti = id2d-((xti-1).*window(1));

        yti = yti+yshiftmin;
        xti = xti+xshiftmin;
        zti = zti+zshiftmin;

        %vector of indices of informed points in the ti
        id = 1+(sizeyxti.*((zti+vect_to_informed_nodesz)-1)+ti_size(1).*((xti+vect_to_informed_nodesx)-1)+(yti+vect_to_informed_nodesy-1));
        data_event_ti = ti(id);

        %evaluating the distance between the data event in the simulation and the one in the ti.
        if distance_type == 0;              %for categorial variable
            distance = sum(data_event_sim~=data_event_ti);  %using sum instead of mean: takes less time
        elseif distance_type == 1;          %for continuous variable with euclidean distance
            distance = mean((data_event_sim-data_event_ti).^2);
        elseif distance_type == 2;          %continuous variable, weighting with inverse square distance
            distance = sum(((data_event_sim-data_event_ti).^2).*weights);
        end

        %checking if the distance is under the minimum distance found for thes node in the simulation
        if distance < mindist
            mindist = distance;
            bestpoint = sizeyxti.*(zti-1)+ti_size(1).*(xti-1)+yti;
        end
        %if distance under t, the best point is accepted.
        if mindist <= thr*(i/nb_of_tries)
            break
        end
    end

    nbtries(simnod) = i;
    bestmin(simnod) = mindist;
    thr = mean(bestmin(1:simnod))*t;
    %if enough of the ti has been scanned, take the best point so far (minimum distance)
    simul(path_sim(simnod)) = ti(bestpoint);

end

evalquality = mean(bestmin);
evaltries = mean(nbtries);
simul=reshape(simul,simul_size(1),simul_size(2),simul_size(3));
toc

%% visualizing
figure(1); clf;
subplot(1,2,1);hold on
imagesc(simul)
axis equal tight
colormap gray
scatter(data(:,2),data(:,1),20,data(:,4),'o','filled')
scatter(data(:,2),data(:,1),'ok')
title('simulation')

subplot(1,2,2);hold on
imagesc(ti)
axis equal tight
colormap gray
title('training image')


