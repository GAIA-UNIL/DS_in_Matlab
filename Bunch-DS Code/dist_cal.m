function Distance=dist_cal(data_event_cond , data_event_ti , data_event_size , Distance_type)    
c = data_event_size;

         if Distance_type == 1;              %for categorial variable
            Distance = mean(data_event_cond(1:c,4)~=data_event_ti(1:c,1));  %using sum instead of mean: takes less time
         
         elseif Distance_type == 2;          %for continuous variable with euclidean Distance
         % First Option
%             d_max = max(abs(data_event_ti(1:c,1)-data_event_cond(1:c,4)));%abs(max(data_event_ti)-min(data_event_ti));
%          if d_max ~= 0
%             Distance = (sum(abs(data_event_ti(1:c,1)-data_event_cond(1:c,4))./d_max))/(c);
%          elseif d_max == 0
%             Distance = 0;
%          end
            % Second Option
            Distance = mean((data_event_cond(1:c,4)-data_event_ti(1:c,1)).^2);
        
         elseif Distance_type == 3;          %continuous variable, weighting with inverse square Distance
                clear d
                d(:,1) = sqrt((sim_x-data_event_ti_coord(1:c,1)).^2+((sim_x-data_event_ti_coord(1:c,1)).^2+...
                ((sim_x-data_event_ti_coord(1:c,1)).^2)));
                weights = 1./(d(:,1).^2);
                Distance = sum(((data_event_cond(:,4)-data_event_ti).^2).*weights);
         end 