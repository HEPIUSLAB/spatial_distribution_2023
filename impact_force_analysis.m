clear; close all;
% Purpose: Generate space-flow plots of NCUS parameters wrt impact force
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% 



%% Section 1: Generate data structure from NCUS, calculate desired parameters

% these rats are the ones with adequate image quality
% load('good_CEUS_data.mat');
load('force_impact_rats.mat');
save_path = "C:\Users\Denis\OneDrive - Johns Hopkins\";

% find data location
data_path = 'D:\Data';
% data_path = 'C:\Users\Denis\OneDrive for Business\Data';
datasets = dir(data_path);

% extract datasets from data folder
datasets = {datasets(contains({datasets.name}, 'Rat SCBF')).name};

[S_FM, S_NCpx] = get_NC_data(data_path, datasets, good_data);

analysis_date = datetime("today");
save_name = sprintf('%d%02d%02d Spatial Data.mat', year(analysis_date)-2000,  month(analysis_date), day(analysis_date));
save(fullfile(save_path,'Figures','Spatial Distribution', 'Injury Severity', save_name), 'S_FM', 'S_NCpx','-v7.3');

%% Section 1.5: Draw Spine Borders

for rat = 13:length(S_NCpx)
    fprintf('%d out of %d images\n', rat, length(S_NCpx))
    disp('Draw dorsal border')
    dorsal = draw_borders(S_NCpx(rat).velim);
    disp('Draw ventral border')
    ventral = draw_borders(S_NCpx(rat).velim);
    
    spineROI = ROI_from_manual(S_NCpx(rat).velim, dorsal, ventral);
    S_NCpx(rat).SpineROI = spineROI;
end


disp('Saving...');
save(fullfile(save_path,'Figures','Spatial Distribution', 'Injury Severity', save_name), 'S_NCpx','-append');

%% Step 1.2: Check ROI
figure

for rat = 1:length(S_NCpx)
    subplot(ceil(length(S_NCpx)/4),4,rat)
    imshow((S_NCpx(rat).SpineROI).*S_NCpx(rat).velim,[])
end

%% Load most recent file

S_FM = load("C:\Users\Denis\OneDrive - Johns Hopkins\Figures\Spatial Distribution\Injury Severity\230525 Spatial Data.mat").S_FM;
load("C:\Users\Denis\OneDrive - Johns Hopkins\Figures\Spatial Distribution\Injury Severity\230526 SMI averaged px data.mat")
% load("C:\Users\Denis\OneDrive - Johns Hopkins\Figures\Spatial Distribution\Injury Severity\230526 ADF px data.mat")

%% Section 2: Sort by force
[~,sortforce] = sort([S_FM.ImpactForce]);


if sum(~(sortforce==(1:length(S_FM))))
    S_NCpost = S_NCpx(1:2:end);
    S_NCpre = S_NCpx(2:2:end);

    S_FM = S_FM(sortforce);
    S_NCpost = S_NCpost(sortforce);
    S_NCpre = S_NCpre(sortforce);

    clear S_NCpx;

    disp('Sorted')

else
    disp('I already sorted you goon')
end


%% Section 3: Calculate pixelwise analysis (with plots)

figure(1)
S_NCpre = flow_pixelwise(S_NCpre, S_FM, 'pre', true);
figure(2)
S_NCpost = flow_pixelwise(S_NCpost, S_FM, 'post',true);


%% Export 230105 Rat2 as representative rat (for Fig. 5)

rat = 17;
% without ROI
rawflo = S_NCpost(rat).velim;
imwrite(rawflo, ...
    fullfile("C:\Users\Denis\OneDrive - Johns Hopkins\230124 Spatial Distribution Paper\Figures\5. BF Patterns (pretty)", sprintf('%s vmap raw.png', S_FM(rat).Rat)));

% with ROI
rawflo(~S_NCpost(rat).SpineROI) = 0;
imwrite(rawflo, ...
    fullfile("C:\Users\Denis\OneDrive - Johns Hopkins\230124 Spatial Distribution Paper\Figures\5. BF Patterns (pretty)", sprintf('%s vmap trimmed.png', S_FM(rat).Rat)));

%% Export 230105 Rat2 as representative rat (for Fig. 3)
rat = 17;
distIm = round(S_NCpost(rat).DistanceImage, 1);
editROI = S_NCpost(rat).SpineROI;
distIm(~editROI) = 0;
distIm = abs(distIm);
fprintf('Max Dist in Bar: %02g mm\n',max(distIm(:)));
distIm = distIm./max(distIm(:));
distIm(334:373, 150:645) = repmat((0:(1/500):0.99), 40, 1);
editROI(334:373, 150:645) = 1;
distIm = cat(3, 0.7*(1-distIm).*(editROI), 0.5*distIm, distIm);
distIm(repmat(logical(imdilate(S_NCpost(rat).InjuryImage, strel('disk', 5))), 1, 1, 3)) = 0.6;

imshow(distIm)
imwrite(distIm, ...
    fullfile("C:\Users\Denis\OneDrive - Johns Hopkins\230124 Spatial Distribution Paper\Figures\2. NCUS Analysis", sprintf('%s distance trimmed draft2.png', S_FM(rat).Rat)));

%% Section 4: Calculate values for injury regions (for Fig. 5/6)

% includes ouput of post values for Fig. 5

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

Sratio_pre = Sval_pre;
Sratio_post = Sval_post;

% set distance binning
binmm = 0.1;
dist = -12:binmm:12;

% Calculate distance distribution for each rat
for rat = 1:length(S_NCpost)

    % load unbinned data
    predist = S_NCpre(rat).DistanceNCUS;
    postdist = S_NCpost(rat).DistanceNCUS;    
    preflo = S_NCpre(rat).NCUSFlow;
    postflo = S_NCpost(rat).NCUSFlow;

    plottitle = S_FM(rat).ImpactForce;
    
    distImpre = S_NCpre(rat).DistanceImage;
    distImpost = S_NCpost(rat).DistanceImage;


    % initialize storage variables
    avpre = nan(size(dist));
    avpost = nan(size(dist));
    

    for dd = 1:length(dist)
        d = dist(dd);

        % V parameter (Vel Index)
%         avpre(dd) = mean(preflo(predist >=d & predist < d+1));
%         avpost(dd) = mean(postflo(postdist >=d & postdist < d+1));

        % A parameter (Area)
%         avpre(dd) = length(preflo(predist >=d & predist < d+1));
%         avpost(dd) = length(postflo(postdist >=d & postdist < d+1));
        
        % I parameter (Area-Adjusted Vel Index)
        avpre(dd) = mean(preflo(predist >=d & predist < d+1))*length(preflo(predist >=d & predist < d+1))/sum(distImpre >=d & distImpre <d+1, "all");
        avpost(dd) = mean(postflo(postdist >=d & postdist < d+1))*length(postflo(postdist >=d & postdist < d+1))/sum(distImpost >=d & distImpost <d+1, "all");        
    end
    
    %%%% FOR FIGURE 5: Save data for 230105 Rat2 as representative rat %%%%
%     if rat==17
%         writematrix([dist; avpre; avpost]', ...
%             fullfile("C:\Users\Denis\OneDrive - Johns Hopkins\230124 Spatial Distribution Paper\Figures\5. BF Patterns (pretty)", sprintf('%s prepostflo.xls', S_FM(rat).Rat)));
%     end
    
    % values for pre and post separately
    [pre_values, post_values] = find_injury_regions(avpre, avpost, 0, dist,0.03);
    lims_val = [0 0.8];

    % store values
    Sval_pre(rat) = pre_values;
    Sval_post(rat) = post_values;

    
    
    
    
    % plot resulting labels
    figure(6);
    subplot(length(S_NCpost)/3,3,rat);
    plot(dist, avpost)
%     plot(dist, ratio)
    hold on
    xline(dist(post_values.UmbraLoc), 'Color','r')
    xline(dist(post_values.CaudLoc), 'Color','g')
    xline(dist(post_values.CranLoc), 'Color','g')
    xline(dist(post_values.Transition), 'Color',"#D95319")
    title(plottitle)
    ylim(lims_val)
%     ylim(lims_ratio)
    xlim([-12 12])
    hold off
end


%% Section 5: sort rats into bins and save variables for output to Prism (Fig. 6)

numvars = 12;

% bin into Mild, Moderate, Severe
% although moderate injury was intended for 175 kDyn, actual mean was ~200 kDyn
force_bins = [100, 200, 250];
[~,I] = min(abs([S_FM.ImpactForce]'-force_bins), [], 2);

% find mode for matrix formatting
[~,val_l] = mode(I);
values = nan(val_l, length(force_bins), numvars*2);


for force = 1:length(force_bins)

    % calculate pre and post values
    for vr = 1:2
        if vr == 1
            rats = Sval_pre(I==force);
        elseif vr == 2
            rats = Sval_post(I==force);
            
            %%% same unused ration %%%
%             rats = Sratio_post(I==force);
        end
        
        
        % need to store each rat individually
        for rat = 1:length(rats)
            
            % Umbra
            values(rat, force, 1+numvars*(vr-1)) = [rats(rat).InjuryHeight];
            vartitles{1,1} = 'Injury Height';

            % not useful, change from post to post
            values(rat, force, 2+numvars*(vr-1)) = [rats(rat).InjuryHeight]-[rats(rat).InjuryHeight];
            vartitles{2,1} = 'Injury Height Change';
            
            % height of each penumbra
            values(rat, force, 3+numvars*(vr-1)) = [rats(rat).CranPenHeight];
            vartitles{3,1} = 'Cranial Penumbra Height';
            values(rat, force, 4+numvars*(vr-1)) = [rats(rat).CaudPenHeight];
            vartitles{4,1} = 'Caudal Penumbra Height';
            
            % Injury widths (only valid for post-injury)
            if vr ~= 1
                tempvar = vertcat(rats(rat).Transition);
                values(rat, force, 5+numvars*(vr-1)) = binmm*(tempvar(:,2)-tempvar(:,1));
                vartitles{5,1} = 'Half Maximum Width';
                values(rat, force, 6+numvars*(vr-1)) = binmm*([rats(rat).CaudLoc]-[rats(rat).CranLoc]);
                vartitles{6,1} = 'Full Maximum Width';

                % Slope of umbra
                values(rat, force, 11+numvars*(vr-1)) = -0.5*([rats(rat).CranPenHeight]-[rats(rat).InjuryHeight])./(binmm*([rats(rat).Transition(1)]-[rats(rat).UmbraLoc]));
                vartitles{11,1} = 'Cranial Slope';
                values(rat, force, 12+numvars*(vr-1)) = 0.5*([rats(rat).CaudPenHeight]-[rats(rat).InjuryHeight])./(binmm*([rats(rat).Transition(2)]-[rats(rat).UmbraLoc]));
                vartitles{12,1} = 'Caudal Slope';
            end
            
            % Penumbra prominence of each side
            values(rat, force, 7+numvars*(vr-1)) = [rats(rat).CranPenHeight]-[rats(rat).InjuryHeight];
            vartitles{7,1} = 'Cranial Penumbra Difference';
            values(rat, force, 8+numvars*(vr-1)) = [rats(rat).CaudPenHeight]-[rats(rat).InjuryHeight];
            vartitles{8,1} = 'Caudal Penumbra Difference';
            
            % Distal (assumed healthy) tissue
            values(rat, force, 9+numvars*(vr-1)) = [rats(rat).CranNormHeight];
            vartitles{9,1} = 'Cranial Healthy Flow';
            values(rat, force, 10+numvars*(vr-1)) = [rats(rat).CaudNormHeight];
            vartitles{10,1} = 'Caudal Healthy Flow';

            
        end
    end
end

%% Seciton 5B: Organizing data structure for Prism plots (Fig. 6)

% this section is done for easy copy and pasting into the Prism file

S_bars = struct();
floratio = {'Pre', 'Flow', 'Ratio'};

figure(4)
for vv = 13:size(values, 3)

    %%% plot for early comparison %%%
    subplot(5,5,vv);
    scatter(force_bins, values(:,:,vv));
    try title(vartitles{mod(vv-1,numvars)+1})
    end
    
    % the weird format is due to how I plot in Prism
    entry = nan(5, 24*3);
    flatpre = values(:,:,mod(vv-1,10)+1);
    flatpre = flatpre(~isnan(flatpre));
    
    entry(1,49:end) = [flatpre' nan(1,24-length(flatpre))];
    entry(3:end, 1:10) = values(:,:,mod(vv-1,numvars)+1)';
    entry(3:end, 25:34) = values(:,:,vv)';
    
    % save chosen variable
    S_bars.([erase(vartitles{mod(vv-1,numvars)+1},[" ", "%"]), floratio{ceil(vv/numvars)}]) = entry;
%     S_bars.([erase(vartitles{mod(vv-1,10)+1},[" ", "%"]), floratio{ceil(vv/10)}]) = [values(:,:,mod(vv-1,10)+1); values(:,:,vv)]';
end

% save("C:\Users\Denis\OneDrive - Johns Hopkins\230124 Spatial Distribution Paper\Figures\6. Spatial bar graphs\230606 barvars2.mat", 'S_bars');

%% Section 6: calculated error of injury label

% discrepancy between umbra location and injury drawn
injury_label_error = dist([Sval_post.UmbraLoc]);

% value of error
error_mean = mean(injury_label_error);
error_std = std(injury_label_error);
[~, perr] = ttest((injury_label_error));

% magnitude of error
error_mag = mean(abs(injury_label_error));
error_mag_std = std(abs(injury_label_error));
[~, pmag] = ttest(abs(injury_label_error));

fprintf('Average Error: %f +- %f; p = %f \n', error_mean, error_std, perr)
fprintf('Average Error Magnitude: %f +- %f; p = %f\n', error_mag, error_mag_std, pmag)

