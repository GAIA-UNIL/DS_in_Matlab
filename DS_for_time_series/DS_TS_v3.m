%DO NOT USE FR BENCHMARKING

%Performs time series reconstruction based on a multiple-point realizations using Direct Sampling (DS).
%The general DS method is described here:
% Mariethoz G, Renard P, Straubhaar J (2010). The Direct Sampling method to perform
% multiple-point geostatistical simulations. Water resources Research 46(W11536),
% doi:10.1029/2008WR007621.
%Its application to times series is described here:
% Mariethoz G, Linde N, Jougnot D, Rezaee H (2015). Feature-preserving interpolation
% and filtering of environmental time series. Environmental Modelling & Software, 72, 71-76,
% doi: 10.1016/j.envsoft.2015.07.001.

clear;home

%input file should be a dime series as an ascii vector of equally spaced
%values, with -999999999 or nan for unknown values
datafile='incomplete_TS.txt';
ti=load(datafile);

ti(ti==-999999999)=nan;
ti=ti(5000:10000);

% Simulation parameters
simul_size = size(ti,1);%size of the simulation: y x z
search_radius = 50;    %maximum radius of the data events
n = 25;                %maximum number of points in the data event
f = 1;               %maximum fraction of scanned training image
t_DS = 0.001;               %distance threshold (between 0 and 1)
distance_type=2;       %type of variable (0=categorical, 1= continuous with
%Manhattan distance, 2= continuous with euclidean distance)
postprocessing=0;

t=t_DS*( max(ti(isfinite(ti))) - min(ti(isfinite(ti))) );

ti_size=simul_size;
tic

%defining path in simulation
path_sim = randperm(simul_size);

%creating empty simulation
simul=ti;

bestmin = zeros(size(path_sim,2),1);
nbtries = zeros(size(path_sim,2),1);
thr = inf;

for postpro=1:postprocessing+1
    %looping simulation nodes
    for simnod = 1:size(path_sim,2)
        
        xsim=path_sim(simnod);
        
        if ( isfinite(ti(xsim)))
            continue
        end
        
        if ( isfinite(ti(xsim)) & (postpro==1))
            continue
        end        
        
        %finding distance with all previously simulated points
        all_sim_pts = find(isfinite(simul)==1);
        d = abs(xsim-all_sim_pts);
        
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
            vect_to_informed_nodesx = 0;
            data_event_sim = 0;
        else
            %finding vectors to each point and data event
            vect_to_informed_nodesx=zeros(informed_nodes_counter,1);
            data_event_sim=zeros(informed_nodes_counter,1);
            
            xt=all_sim_pts(1:informed_nodes_counter);
            vect_to_informed_nodesx(1:informed_nodes_counter) = xt-xsim;
            data_event_sim(1:informed_nodes_counter) = simul(all_sim_pts(1:informed_nodes_counter));
            
            %         for i = 1:informed_nodes_counter
            %             xt=all_sim_pts(i);
            %             vect_to_informed_nodesx(i) = xt-xsim;
            %
            %             %indice of the data point in the simul
            %             data_event_sim(i) = simul(all_sim_pts(i));
            %         end
        end
        
        %defining dark zone related to the extension of the data event
        xshiftmax = ti_size-max(vect_to_informed_nodesx);
        xshiftmin = -min(vect_to_informed_nodesx);
        
        %taking care not to go out of the training image
        if xshiftmax > ti_size
            xshiftmax = ti_size;
        end
        if xshiftmin < 0
            xshiftmin = 0;
        end
        
        %defining a random path
        window_size = (xshiftmax-xshiftmin);
        scanned_node_id=round(unifrnd(1,window_size));
        %path_ti = randperm(window_size);   %%%%%
        
        
        
        %scanning ti
        mindist = inf;  %initial best distance is set to inf. Updated with every best distance encountered
        nb_of_tries = ceil(window_size*f); %number of tries in the ti
        
        if distance_type == 0;              %for categorial variable
            t_now=t;
        elseif distance_type == 1;          %for continuous variable with euclidean distance
            t_now=t;
        elseif distance_type == 2;          %continuous variable, weighting with inverse square distance
            t_now=t^2;
        end
        
        for i = 1:nb_of_tries
            
            %finding point in the random path and shifting
            scanned_node_id=scanned_node_id+1;
            if (scanned_node_id >window_size)
                scanned_node_id=1;
            end
            %xti = path_ti(i);
            xti = scanned_node_id+xshiftmin;
            
            if isnan(ti(xti))
                continue
            end
            
            %vector of indices of informed points in the ti
            id = xti+vect_to_informed_nodesx;
            data_event_ti = ti(id);
            ind=isfinite(data_event_ti);
            nelm=sum(ind);
            
            %evaluating the distance between the data event in the simulation and the one in the ti.
            if distance_type == 0;              %for categorial variable
                distance = sum(data_event_sim(ind)~=data_event_ti(ind))/nelm;
            elseif distance_type == 1;          %for continuous variable with euclidean distance
                distance = (sum(abs(data_event_sim(ind)-data_event_ti(ind))))/nelm;
            elseif distance_type == 2;          %continuous variable, weighting with inverse square distance
                distance = (sum(data_event_sim(ind)-data_event_ti(ind)).^2)/nelm;
            end
            
            %checking if the distance is under the minimum distance found for thes node in the simulation
            if distance < mindist
                mindist = distance;
                bestpoint = xti;
            end
            %if distance under t, the best point is accepted.
            if mindist <= t_now;
                break
            end
        end
        
        nbtries(simnod) = i;
        bestmin(simnod) = mindist;
        %if enough of the ti has been scanned, take the best point so far (minimum distance)
        simul(path_sim(simnod)) = ti(bestpoint);
        
    end
    
end

evalquality = mean(bestmin);
evaltries = mean(nbtries);
toc

%% visualizing
figure(11);clf;hold on

plot(simul,'r')
plot(ti,'b')

