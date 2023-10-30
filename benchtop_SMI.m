clear; close all;
% Purpose: To test SMI VI at different flow and recording parameters
% This script was written for use with flow data stored as .mkvs
% Please email me for questions about data format.
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% 


%% Step 1: Read flowrate over time and convert csv times to datetime

A = readmatrix("D:\Data\231013 Benchtop SMI Test\231013_flowrate_times.csv");

flow_rate = A(:,1);
t_fr = datetime('2023-10-13','InputFormat','yyyy-MM-dd')+hours(A(:,2))+minutes(A(:,3))+seconds(A(:,4));

load('D:\Data\231013 Benchtop SMI Test\US\Phantom\spatial_analysis_v1.mat')

% % plot flow rate over time
% figure(1)
% yyaxis left
% plot(t_USpost, Ipost);
% yyaxis right
% plot(t_fr, flow_rate)

%% Step 2: Initialize data to analyze

% set flow variables and sampling intervals
US = Ipost;
t_US = t_USpost;
ind_step = ceil(seconds(15)/(t_US(2)-t_US(1)));
ind_buff = ceil(seconds(3)/(t_US(2)-t_US(1)));

% this specific time did not have good readings, clear from data
delete_list = [17];

% initialize storage variables
[US_fr, US_std] = deal(nan(size(flow_rate)));

% for each time/flowrate, calculate average flow
for fr = 1:length(t_fr)
    % removes data points with unsuccessful recording
    if ismember(fr, delete_list)
        continue;
    end
    
    % needed to refill Doppler fluid syringe, cut out truncated recordings,
    % see lab notebook pictures for details.
    [~, min_ind] = min(abs(t_US-t_fr(fr)));
    if fr == 37
        [~, min_ind] = min(abs(t_US-datetime("13-Oct-2023 12:37:23")));
    end

    if fr == 76
        [~, min_ind] = min(abs(t_US-datetime("13-Oct-2023 13:32:40")));
    end
    
    % calculate mean and SD of flow
    US_fr(fr) = mean(US(min_ind+ind_buff:min_ind+ind_step),'omitnan');
    US_std(fr) = std(US(min_ind+ind_buff:min_ind+ind_step),'omitnan'); 

end


% % plot average flow rate over time
% figure
% yyaxis left
% plot(t_fr, US_fr);
% yyaxis right
% plot(t_fr, flow_rate)

%% Step 3: Reorganize data matrix to flowrate vs tube diameter/SMI scale

US_fr = reshape(US_fr,[], 6);
US_std = reshape(US_std,[], 6);
flow_rate = reshape(flow_rate, [],6);
t_fr = reshape(t_fr, [],6);

% calculate velocity from flow rate
vel = flow_rate./(3600*pi*([0.0506, 0.0506,0.0506, 0.0408, 0.0408, 0.0408]/2).^2);

%% Step 4: Figure 4 (new): Plot flow velocity against US values

% choose diamter
trials = [1,4];
% trials = [2,5];
% trials = [3,6];

figure(2)
set(figure(2), 'Position', [367,403,298,268])

% set appropriate markers for consistent plots
markers = {'square'; 'diamond';'^'; 'square'; 'diamond';'^'};
scale = {'Maximum Scale', 'Scale: 0.9', 'Scale: 0.4'};
mfc = {[0 0.85 0]; [0 0.85 0]; [0 0.85 0]; 'k'; 'k'; 'k'};
legend_titles = {'506 \mum, Scale: 1.8',...
                 '506 \mum',...
                 '506 \mum',...
                 '408 \mum, Scale: 1.2',...
                 '408 \mum',...
                 '408 \mum'};


pl = errorbar(vel(:,trials), US_fr(:,trials), US_std(:,trials), 'LineStyle', 'none', 'MarkerSize', 10, 'MarkerEdgeColor', 'none', 'LineWidth', 2);
[pl.Marker] = markers{trials};
[pl.MarkerFaceColor] = mfc{trials};
[pl.Color] = mfc{trials};

legend(legend_titles(trials), 'Location', 'nw')
title(scale{trials(1)})

ylabel({'Area-Adjusted', 'Velocity Index (AU)'})
xlabel('Flow Velocity (cm/s)')
xlim([0 10])

set(gca,'FontWeight','bold')
set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)

exportgraphics(figure(2),...
               fullfile("C:\Users\Denis\OneDrive - Johns Hopkins\Lab Personal\230824 Spatial Distribution Paper\Figures\Supp. Benchtop\",...
               sprintf('231016 subplots %s.pdf', strrep(scale{trials(1)}, ':', ' -'))),...
               'ContentType','vector');

%% Step 5: Statistics (Spearman's correlation)

% Calculate correlations and fit linear models
[rho, pval] = corr(vel,US_fr,'rows','complete','type','Spearman');

for rr = 1:size(rho,2)
    fprintf('%s; size %s rho: %.02f; p: %.02f\n', scale{mod(rr-1,3)+1}, legend_titles{rr}, rho(1,rr), pval(1,rr))
end
