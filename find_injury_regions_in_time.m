% Function: find_injury_regions
%
% Purpose: to calculate injury regions from flow surfaces
%
% Input parameters: 
%       all_rats: struct
%       dist: double vector
%       plot_bool: logical
%
% Output parameters:
%       pre_values: struct
%       post_values: struct
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function [pre_values, post_values] = find_injury_regions_in_time(all_rats, dist, plot_bool)

    if nargin <3
        plot_bool = false;
    end
        
    % initialize pre- and post-injury structure
    post_values = struct('Time',[], ...
                         'InjuryHeight',[], ...
                         'CranPenHeight',[], ...
                         'CaudPenHeight',[], ...
                         'CranNormHeight',[], ...
                         'CaudNormHeight',[], ...
                         'CranPenWidth',[], ...
                         'CaudPenWidth',[], ...
                         'UmbraLoc',[], ...
                         'CranLoc',[], ...
                         'CaudLoc',[], ...
                         'Transition',[]);

    pre_values = struct('InjuryHeight',[], ...
                        'CranPenHeight',[], ...
                        'CaudPenHeight',[], ...
                        'CranNormHeight',[], ...
                        'CaudNormHeight',[]);
    
    % perform processing for each rat
    for rat = 1:length(all_rats)

        % load relevant variable from structure
        time = [all_rats(rat).TimePost];
        time = time-all_rats(rat).InjuryTime;
        avpost = [all_rats(rat).FlowPost];
        avpre = [all_rats(rat).FlowPre];        

        % individual rat troubleshooting (cuts out artifactual recording)
        if rat == 2
            avpost = avpost(1:end-200,:);
            time = time(1:end-200);
        end

        if rat == 8
            avpost = avpost(1:266,:);
            time = time(1:266);
        end

        % skip if no data post-injury
        if isempty(avpost) continue
        end
        
        % find injury regions 
        [pre_rat, post_rat] = find_injury_regions(avpre, avpost, time, dist);

        % store values in structure
        post_values(rat) = post_rat; 
        
        % pre values are averaged over time
        pre_values(rat) = pre_rat;
        
        % plot to evaluate code/results
        if plot_bool
%             subplot(ceil(length(all_rats)/3),3,rat)
            [X,Y] = meshgrid(dist, (time)/60);
            surf(X,Y,avpost)
            shading(gca,'interp')
            hold on
            plot3(dist(post_rat.CranLoc), (post_rat.Time), post_values(rat).CranPenHeight(1:length(post_rat.CranLoc)), 'LineWidth',4, 'Color','g');
            plot3(dist(post_rat.CaudLoc), (post_rat.Time), post_values(rat).CaudPenHeight(1:length(post_rat.CaudLoc)), 'LineWidth',4, 'Color','g');
            plot3(dist(post_rat.UmbraLoc), (post_rat.Time), post_values(rat).InjuryHeight(1:length(post_rat.UmbraLoc)), 'LineWidth',4, 'Color','r');
            plot3(dist(post_rat.Transition), (post_rat.Time), ones(size(post_rat.InjuryHeight)), 'LineWidth',4, 'Color','y');
        
            title(all_rats(rat).Date)
            view(2)
            hold off
        end
        
    end

end