clear; close all;
% Purpose: To find changes in the spatial distribution after probe movement
% (Supplementary Fig. S2)
% This script was written for use with blood flow data stored as .mkvs
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% 



%% Step 1: Initialize data to analyze

root_path = 'D:\Data\';

% generate list of all rats with probe movement
a = dir(root_path);
rootexpsets = {a(contains({a.name}, '230505') | ...
                 contains({a.name}, '230512') | ...
                 contains({a.name}, '230530') | ...
                 contains({a.name}, '230717') | ...
                 contains({a.name}, '230829')).name};

%% Step 2: Load spatial distribution data from selected rats

figure(1)
pre = struct();

for rat = 1:length(rootexpsets)
    
    warning('') % Clear last warning message   
    
    % most data stored as USpre
    load(fullfile(root_path, rootexpsets{rat}, '\US\Spat. Dist. -12.0 to 12.0 by 0.1 mm\spatial_analysis_v1.mat'), 'USpre', 't_USpre', 'inj_time')
    
    pre(rat).t = t_USpre;
    pre(rat).US = USpre;
    
    % detect if data stored as Ipre
    [warnMsg, warnId] = lastwarn;
    if ~isempty(warnMsg)
        load(fullfile(root_path, rootexpsets{rat}, '\US\Spat. Dist. -12.0 to 12.0 by 0.1 mm\spatial_analysis_v1.mat'), 't_USpre', 'Ipre', 'inj_time')        
        pre(rat).US = Ipre;
    end

end


%% Step 3: Detect image hops (to determine bounds in nest step)

% plot each rat to determine probe change time
for rat = 1:length(pre)
    
    figure(2), subplot(2,3,rat)
    % flow throughout all distance
%     [X,Y] = meshgrid(dist, 1:length(pre(rat).t));
    [X,Y] = meshgrid(dist, pre(rat).t);
    
    Spre = surf(X(1:100:end,:), Y(1:100:end,:), pre(rat).US(1:100:end,:), 'EdgeColor','flat');
    view(2)

    title(rootexpsets{rat})
end


%% Step 4: Calculate spatial distribution on the images

% bounds determined from previous step
bounds = [9000 11000 50000 52000;...
          15000 17000 79000 81000;...
          9000 11000 67000 69000;...
          9000 11000 70000 72000;...
          9000 11000 40000 42000];

% Initialize structures for storing information
Sval_post  = struct('Time',[], ...
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

Sval_pre = struct('InjuryHeight',[], ...
                    'CranPenHeight',[], ...
                    'CaudPenHeight',[], ...
                    'CranNormHeight',[], ...
                    'CaudNormHeight',[]);

err_dist = zeros(length(pre), length(dist));


for rat = 1:length(pre)
    
    % set pre and post movement data
    avpre = mean(pre(rat).US(bounds(rat,1):bounds(rat,2),:), 1);
    avpost = mean(pre(rat).US(bounds(rat,3):bounds(rat,4),:), 1);

    err_dist(rat, :) = 100*(avpost-avpre)./avpre;

    % values for pre and post separately
    [pre_values, post_values] = find_injury_regions(avpre, avpost, 0, dist,0.03);

    % plot to verify
    figure(3), subplot(2,3,rat)
    plot(dist, avpre);
    hold on
    plot(dist, avpost);
    xline(dist(post_values.UmbraLoc), 'Color','r')
    xline(dist(post_values.CranLoc), 'Color','g')
    xline(dist(post_values.CaudLoc), 'Color','g')
    hold off
    xlim([-7.5 7.5])
    ylim([0 0.8])
    title(rootexpsets{rat})

    % store values
    Sval_pre(rat) = pre_values;
    Sval_post(rat) = post_values;
end

% store injury width
injwid = num2cell(diff(dist(cell2mat({Sval_post.Transition}')),1, 2));
[Sval_post.InjWidth] = injwid{:};

%% Step 5: Save values for Prism

save("C:\Users\Denis\OneDrive - Johns Hopkins\Lab Personal\230824 Spatial Distribution Paper\Figures\6. Spatial bar graphs\231016 Control Rats\controls.mat",...
     "Sval_pre", "Sval_post", 'pre')


