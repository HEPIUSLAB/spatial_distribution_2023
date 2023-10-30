clear; close all;
% Purpose: To capture pixelwise spatial distribution evolution over time.
% This script was written for use with blood flow data stored as DICOMs
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% 

%% Step 1: Initialize data to analyze

root_path = 'D:\Data\';

% generate list of all autoreg rats
a = dir(root_path);
rootexpsets = {a(contains({a.name}, 'Rat SCAR') | contains({a.name}, 'Rat Vasc Reac') | contains({a.name}, 'propofol')).name};

%% Step 2: Identify position changes through experiment and establish injury location

for rat = 2:length(rootexpsets)

    % check if the data is stored in DICOMs
    try
        load(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"))
        disp(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"))
    catch
        continue
    end
    
    % create stream of entire experiment
    vid = cat(4,listing.image);
    
    % label which frames contain a position change or other status change
    % (injury). Includes first frame as a status change
    frames = user_record_frames(vid);
    frames = [1;frames];
    frames = sort(frames);

    % draw location of injury in each position/state. Starts from back so
    % user can reference injury image.
    Injuries = draw_injuries(vid, frames');
    
    % store injury location in structure
    for ff = 1:length(frames)-1
        [listing(frames(ff):frames(ff+1)).InjuryROI] = deal(Injuries(:,:,ff));
    end    
    [listing(frames(end):end).InjuryROI] = deal(Injuries(:,:,end));


    % save data    
    save(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"), 'listing','-append');
end


%% Step 3: User defines SC borders

for rat = 2:length(rootexpsets)

    % check if the data is stored in DICOMs
    try
        load(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"))
        disp(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"))
    catch
        continue
    end
    
    % create stream of entire experiment
    vid = cat(4,listing.image);

    injROIs = cat(3, listing.InjuryROI);

    frames = [1; find(diff(injROIs(1,1,:), 1,3)~=0)+1];
%     frames = 1;
    
    disp('Draw dorsal border')
    dorsal = draw_borders(vid, frames');
    disp('Draw ventral border')
    ventral = draw_borders(vid, frames');
    
    spineROI = ROI_from_manual(vid, dorsal, ventral, frames);
    
        
    % store injury location in structure
    for ff = 1:length(frames)-1
        [listing(frames(ff):frames(ff+1)).SpineROI] = deal(spineROI(:,:,ff));
    end    
    [listing(frames(end):end).SpineROI] = deal(spineROI(:,:,end));


    % save data    
    save(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"), 'listing','-append');
end

%% Step 3: Analyze blood flow parameter through space

for rat = 2:length(rootexpsets)

    % check if the data is stored in DICOMs
    try
        load(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"))
        disp(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"))
    catch
        continue
    end    

    % find velocity maps in analysis folders
    vel_dir = fullfile(root_path, rootexpsets{rat}, "\US\", listing(1).AcquisitionDate, "Velocity Maps");
    maps = dir(vel_dir);
    maps = {maps(contains({maps.name}, 'Map')).name};    

    % check and/or generate save directory
    save_dir = fullfile(vel_dir,'..','Spatial Flow');
    if ~isdir(save_dir)
        mkdir(save_dir)
    %%% below commented section is for if analysis paradigm is changed    
%     else
%         movefile(save_dir, fullfile(save_dir,'..','Spatial Flow old'));
%         mkdir(save_dir)
    end    

    % for status update in command window
    lineLength = 0;
    
    % set impossible ROI for future creation of distance image
    lastROI = zeros(2);
    
%     disttest = zeros([size(velocity, [1,2]),length(maps)]);
    % go through each velocity map
    for mm = 1:length(maps)

        % check for prior analysis
        if isfile(fullfile(save_dir, maps{mm}))
            continue
        end
        
        % update command window
        fprintf(repmat('\b',1,lineLength));
        lineLength = fprintf('Current Map: %d, Total Maps: %d\n', mm, length(maps));
        
        % load file
        load(fullfile(vel_dir, maps{mm}))

        % injury and spatial resolution variables
        injROI = listing(strcmp({listing.name}, scanner_parameters.curr_dataset(1:5))).InjuryROI;
        delxy = listing(strcmp({listing.name}, scanner_parameters.curr_dataset(1:5))).delxy;
        spineROI = listing(strcmp({listing.name}, scanner_parameters.curr_dataset(1:5))).SpineROI;

        % unanalyzable maps will not have delxy
        if isempty(delxy)
            continue
        end
        
        % establish distance axis and resolution
        dist = -12:0.1:12;
        dist_flo = zeros([size(velocity, 3), length(dist)]);
        
        save(fullfile(save_dir, maps{mm}), "dist");
        
        % if a new ROI is encountered, regenerate distance image
        if injROI ~= lastROI
            
            % ROI that excludes off-target signal/noise and filters out
            % major arteries (ASA and PSA)
%             spineROI = remove_ASA(velocity);
            lastROI = injROI;

            % generate image of mxn where the value of each entry is
            % distance from injury ROI
            distIm = flow_through_distance(injROI, velocity(:,:,1), delxy);
            distIm(~spineROI) = nan;
            
            % calculates total pixels at each distance bin for calculation
            % of 'I' parameter
            for dd = 1:length(dist)
                d = dist(dd);
                totalpx(dd) = sum(distIm >=d & distIm <d+1, "all");
            end

            im2rite = 0.5*spineROI+(velocity(:,:,end));
            im2rite(isnan(im2rite)) = 0;
            imwrite(im2rite,fullfile(save_dir, '..',sprintf('%s.png', scanner_parameters.curr_dataset)));
        
        end
        
        
        % framewise calculation, could probably be done at matrix level
        for ff = 1:size(velocity, 3)
            
            % isolation frame
            US = velocity(:,:,ff);
            US(~spineROI) = nan;
%             imshow(US,[])
            
            % measure flow and distance of each non-zero pixel in ROI
            distance = distIm(US>0);    
            flo = US(US>0);
            
            % initialize binned flow
            avflo = nan(size(dist));
            
            % calculate average 'I' for each distance bin and store
            for dd = 1:length(dist)
                d = dist(dd);                
                avflo(dd) = mean(flo(distance >=d & distance < d+1))*length(flo(distance >=d & distance < d+1))/totalpx(dd);             
            end
            dist_flo(ff,:)  = avflo;

        end
        
        % save in folder
        save(fullfile(save_dir, maps{mm}), "dist_flo", '-append');
        
    end
end

%% Step 4: Link files from previous step through time

% this file contains all necessary metadata for each experiment
load('surgery_data.mat')

% structure to save each rat in one file
all_rats = struct();
analysis_date = datetime("today");
all_save = sprintf('%d%02d%02d All Autoreg Rats.mat', year(analysis_date)-2000,  month(analysis_date), day(analysis_date));
save(fullfile(root_path,'..', '\Figures\Spatial Distribution\Time Plotting', all_save), "all_rats");

for rat = 2:length(rootexpsets)
    
    % generate save name    
    save_name = sprintf('%d%02d%02d Spatial Flow.mat', year(analysis_date)-2000,  month(analysis_date), day(analysis_date));
    
    % check if the data is stored in DICOMs
    try
        load(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"))
        disp(fullfile(root_path, rootexpsets{rat}, "\US\DICOM listing.mat"))
    catch
        continue
    end
    
    % load all files with data
    folder = fullfile(root_path, rootexpsets{rat}, "\US\", listing(1).AcquisitionDate, "Spatial Flow");
    
    [t_USpre, USpre, t_USpost, USpost, inj_time] = time_stitching(folder, surgery_data, 'dist_flo');
    
    % store into parent structure (with significant downsampling
    all_rats(rat).Date = listing(1).folder(9:end-3);
    all_rats(rat).InjuryTime = inj_time;
    all_rats(rat).TimePre = t_USpre(1:100:end);
    all_rats(rat).TimePost = t_USpost(1:100:end);
    all_rats(rat).FlowPre = USpre(1:100:end,:);
    all_rats(rat).FlowPost = USpost(1:100:end,:);

    surf([all_rats(rat).FlowPre; all_rats(rat).FlowPost],'EdgeColor','flat');
    
    % save variables into according files (including original sizes)
    disp('Saving...')
    save(fullfile(folder,'..','..','..', save_name), 't_USpre', 'USpre', 't_USpost', 'USpost', 'inj_time', 'listing', 'dist', 'save_name');
    save(fullfile(root_path,'..', '\Figures\Spatial Distribution\Time Plotting', all_save), "all_rats", 'dist', '-append');
    
end


%% Step 5: Calculate injury and penumbra location and levels


load('D:\Figures\Spatial Distribution\Time Plotting\230522 All Autoreg Rats.mat');
%%
figure(3)
[pre_values, post_values] = find_injury_regions_in_time(all_rats, dist, true);


%% Fig. 7B,C: umbra, penumbra, normal over time.   Also average and downsample for Prism

figure(6)
set(gcf, 'Position', [188,675,516,211])


% plotting parameters/rats/vars
% times = -20:.1:-5;
times = 0:.1:30;
rats = [2 3 4 10 12];
vars2plot = [2, 3, 5];  % cranial penumbra/distal
% vars2plot = [2, 4, 6];  % caudal penumbra/distal

% colormap decision
paru = (parula(4));
colors = {'r', paru(2,:), paru(2,:), [234,139,0]/255,  [234,139,0]/255};
ls = {'--', '-', '-', ':', ':'};

% initialize other variables (choose pre or post)
% varnames = fieldnames(pre_values);
varnames = fieldnames(post_values);
post_combined = struct();

for vn = 1:length(varnames)
    
    [mean_flo, mean_time, time_var, flo_var] = deal(nan(length(times), length(rats)));
    
    % average for each rat
    for rat = 1:length(rats)        
        
        % for pre plots (supplement)
%         time_samps = cat(1, pre_values(rats(rat)).Time);
%         flo_samps = cat(1, pre_values(rats(rat)).(varnames{vn}));
        
        % for post plots (main)
        time_samps = cat(1, post_values(rats(rat)).Time);
        flo_samps = cat(1, post_values(rats(rat)).(varnames{vn}));
        
        % average at each time (easiest way to write the code)
        for tt = 1:length(times)-1
            inds = find(time_samps>=times(tt) & time_samps<times(tt+1));
            
            mean_time(tt,rat) = mean(time_samps(inds), 'omitnan');
            mean_flo(tt,rat) = mean(flo_samps(inds), 'omitnan');
            
            time_var(tt,rat) = sqrt(var(time_samps(inds), 'omitnan'));
            flo_var(tt,rat) = sqrt(var(flo_samps(inds), 'omitnan'));        
        end    
    end
    
    % store averages
    post_combined.Time = mean_time;
    post_combined.TimeVariance = time_var;
    post_combined.(varnames{vn}) = mean_flo;
    post_combined.([varnames{vn}, 'Variance']) = flo_var;

    % plot if indicated above
    if ismember(vn, vars2plot)
        disp(varnames{vn});
        
        % mean among rats
        mean_time = mean(post_combined.Time, 2, 'omitnan');
        mean_flo = mean(post_combined.(varnames{vn}), 2, 'omitnan');
        
        % variance among rats
        flo_var = sqrt(var(post_combined.(varnames{vn}), 0, 2, 'omitnan'));
        
        flo_var= flo_var(~isnan(flo_var));
        mean_flo = mean_flo(~isnan(mean_flo));
        mean_time = mean_time(~isnan(mean_time));
    
        fill([mean_time;flipud(mean_time)],[mean_flo-flo_var;flipud(mean_flo+flo_var)],[.9 .9 .9], ...
            'linestyle','none', 'FaceColor', colors{vn-1}, 'FaceAlpha', '0.3');
        hold on
        plot(mean_time,mean_flo, 'Color', colors{vn-1}, 'LineWidth', 2, 'LineStyle', ls{vn-1})
    
        % need to commment areas in find_injury_regions for this to work
        errorbar(-2.6+vn*0.3, mean([pre_values(rats).(varnames{vn})]),std([pre_values(rats).(varnames{vn})]), ...
            'Marker', '.', 'MarkerSize', 14, 'Color', colors{vn-1}, 'LineWidth', 2)
    end
end


% pre supplement
% legend({'','Umbra','','Penumbra', '','Distal'},'Location', 'best')
% xlim([-20 -5])
% xlabel('Time Until Injury (min)');

% post main
legend({'','Umbra','','','Penumbra','', '','Distal'},'Location', 'best')
xlabel('Time After Injury (min)');
% xlim([-3 30])
xlim([-3 15])


ylabel({'Area-Adjusted','Velocity Index (AU)'})
ylim([0 0.7])

set(gca,'FontWeight','bold')
set(gca,'LineWidth',1.5)
hold off

%% Save figure
exportgraphics(figure(6), "C:\Users\Denis\OneDrive - Johns Hopkins\Lab Personal\230824 Spatial Distribution Paper\Figures\7. Change over time\231019 Cranial post.pdf",'ContentType','vector');

%% Fig. 7A: Representative 3D surface plot


figure(2)
set(gcf, 'Position', [182,77,819*0.7,350*0.7]);
rat = 3;

% mesh grid for smoothing of plot (for better visualization)
[xk,yk] = meshgrid(-5:5);
K = exp(-(xk.^2 + (yk.^2)/3)/5); K = K/sum(K(:));


% pre-injury reference
time = [all_rats(rat).TimePre];
flow = conv2([all_rats(rat).FlowPre],K,'same');
[X,Y] = meshgrid(dist, (time-min(time))/60);

pr = plot3(X(1,:), zeros(1,size(flow,2))-1, mean(flow,1), 'LineWidth', 3, 'Color',[0.3,0.3,0.3], 'LineStyle',':');
hold on


% post-injury plot
time = [all_rats(rat).TimePost];
flow = conv2([all_rats(rat).FlowPost],K,'same');  % smoothed for reflection artifact

% easier referencing of parameters
umbra_loc = post_values(rat).UmbraLoc;
umbra_val = post_values(rat).InjuryHeight;
cran_loc = post_values(rat).CranLoc;
caud_loc = post_values(rat).CaudLoc;
transition = post_values(rat).Transition;


% flow throughout all distance
[X,Y] = meshgrid(dist, (time-min(time))/60);
Spost = surf(X(1:10:end,:), Y(1:10:end,:), flow(1:10:end,:));
colormap autumn
alpha(Spost, 0.5)
        
% set correct shading
shading(gca,'interp')

% separate surface for umbra
flowumb = flow;
for uu = 1:length(transition)
%     flowumb(uu, 1:transition(uu,1)) = nan;
%     flowumb(uu,transition(uu,2):end) = nan;
    
    if ~isnan(flow(uu,umbra_loc(uu)))
        flowumb(uu,find(flowumb(uu,:)>0.025+flow(uu,umbra_loc(uu)))) = nan;
    %     flow(uu,umbra_loc(uu))
    else
        flowumb(uu,:) = nan;
    end

end

% plot umbra with deep color
umb = surf(X(1:10:end,:),Y(1:10:end,:),flowumb(1:10:end,:), 'EdgeColor','flat', 'FaceColor','r');
dummyumb = surf(nan(2), nan(2), nan(2), 'EdgeColor','flat', 'FaceColor','r');
alpha(umb, 1)

% plot different locations along surface
pumb1 = plot3(dist(cran_loc), post_values(rat).Time, post_values(rat).CranPenHeight(1:length(cran_loc)), 'LineWidth',4, 'Color',paru(2,:));
pumb2 = plot3(dist(caud_loc), post_values(rat).Time, post_values(rat).CaudPenHeight(1:length(caud_loc)), 'LineWidth',4, 'Color',paru(2,:));
dummypumb = plot3(nan, nan, nan, 'LineWidth',4, 'Color',paru(2,:));
%         plot3(dist(transition), post_values(rat).Time, [post_values(rat).CranPenHeight(1:length(umbra_loc))/2, ...
%             post_values(rat).CaudPenHeight(1:length(umbra_loc))/2, ], 'LineWidth',4, 'Color','y');
%         plot3(dist(umbra_loc), post_values(rat).Time, post_values(rat).InjuryHeight(1:length(umbra_loc)), 'LineWidth',4, 'Color',paru(1,:));        

hold off
% title(all_rats(rat).Date)

xlabel('Distance from Injury (mm)')
ylabel({'Time since','Injury (min)'})
zlabel({'Velocity','Index (AU)'})

% legend({'Baseline Flow','','','','Penumbra Peaks','','Umbra Borders','','Umbra'}, 'Location', 'best')
lege = legend({'Baseline Flow','','','Umbra','','','Penumbra Peaks'}, 'Location', 'northwest','AutoUpdate','off');

% axis labels and position
ylh = get(gca,'ylabel');
set(ylh, 'Rotation',-48, 'Position',[-10,0,0])
xlh = get(gca,'xlabel');
xlp = get(xlh, 'Position');
set(xlh, 'Rotation',2);

xlim([-8 6])
zlim([0 0.9])
ylim([-1 15])

set(gca,'FontWeight','bold')
set(gca,'LineWidth',1.5)
%     set(gca, 'Color', 'k')

set(gca, 'CameraPositionMode', 'manual')
% set(gca, 'CameraPosition',[-22.6274 -290.3636 4.1734]);
set(gca, 'CameraPosition', [-17.928366945994735,-75.422746354860140,4.233415396181620])
set(gca, 'CameraTargetMode','manual');
% set(gca, 'CameraTarget',[-1 20.9317 0.4000]);
set(gca, 'CameraTarget',[-1,7.500000000000000,0.300000000000000]);
set(figure(2),'defaultLegendAutoUpdate','off');

%% Save suface plot as rasterized and everything else as vectorized

BitmapRender(gca, [pr, dummypumb, dummyumb]);
exportgraphics(figure(2), "C:\Users\Denis\OneDrive - Johns Hopkins\Lab Personal\230824 Spatial Distribution Paper\Figures\7. Change over time\231023 3D.pdf",'ContentType','vector');


