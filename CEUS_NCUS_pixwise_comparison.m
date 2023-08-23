clear; close all;
% Purpose: Generate space-flow plots with CEUS and NCUS side-by-side
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% 



%% Section 1: Generate data structure from NCUS, calculate desired parameters

% these rats are the ones with adequate image quality
load('good_CEUS_data.mat');

% find data location
data_path = 'D:\Data';
% data_path = 'C:\Users\Denis\OneDrive for Business\Data';
datasets = dir(data_path);

% figure path to save
figure_path = 'C:\Users\Denis\OneDrive - Johns Hopkins\Figures';

% extract datasets from data folder
datasets = {datasets(contains({datasets.name}, 'Rat SCBF and CEUS')).name};

% compile all the non-contrast data (both FlowMorph and pixelwise)
[S_FM, S_NCUS] = get_NC_data(data_path, datasets, good_data);

% save data structure
analysis_date = datetime("today");
save_name = sprintf('%d%02d%02d Spatial Data.mat', year(analysis_date)-2000,  month(analysis_date), day(analysis_date));
% save(fullfile(figure_path,'Spatial Distribution', 'Contrast Comparison' , save_name), 'S_FM','-v7.3');


%% Section 2: Compile CEUS (and pixelwise NCUS) results from good rats into 1 structure

% load good rats again
load('good_CEUS_data.mat');

% Initialize structure for CEUS measurements
S_pixelwise = struct('Name',[],'Image',[],'ImageFilt',[]);
% S_pixelwise = struct();

% again, specify location of data
data_path = 'D:\Data';
% data_path = 'C:\Users\Denis\OneDrive for Business\Data';
datasets = dir(data_path);
datasets = {datasets(contains({datasets.name}, 'Rat SCBF and CEUS')).name};

% for each experiment date
for ii = 1:length(datasets)
    dset = datasets{ii};
    
    % find good rats
    goodrats = good_data(strcmp({good_data.Date},dset(1:6))).Rats;
    
    % load CEUS structure
    S_CEUS_im = load(fullfile(data_path, dset, 'Spatial Distribution',sprintf('%s 1 px CEUS.mat',dset(1:6)))).S_CEUS_im;
    
    % load each good rat into CEUS structure
    for jj = 1:length(goodrats)
        rat = goodrats{jj};
        S_pixelwise = [S_pixelwise S_CEUS_im(contains({S_CEUS_im.Name},rat))];
    end
end

% trim empty beginning
S_pixelwise = S_pixelwise(2:end);

% also renames several fields for specificity
S_pixelwise = cell2struct([struct2cell(S_pixelwise);struct2cell(S_NCUS)],[{'Name'};{'CEUS_Image'}; {'CEUS_ImageFilt'};fieldnames(S_NCUS)]);

% save(fullfile(figure_path,'Spatial Distribution', 'Contrast Comparison' , save_name), 'S_pixelwise', '-append');

%%

% NOTE!!!! INSTEAD OF RUNNING THIS AGAIN, I TRANSFERRED THE SPINEROI
% VARIABLE FROM THE IMPACT FORCE ANALYSIS RESULT, INSTEAD OF DOING THE
% MORPH DETERMINATION


%% Section 3: Pixelwise processing of CEUS and NCUS images

data_path = 'D:\Data';

for ii = 5%1:length(S_pixelwise)
    plottitle = S_pixelwise(ii).Name(1:11);

    % BS way to find files
    foldertemp = S_pixelwise(ii).Name;
    foldertemp = fullfile(data_path,sprintf('%s Rat SCBF and CEUS', foldertemp(1:6)), foldertemp(1:strfind(foldertemp, ' CEUS')-1));
    
    % need to find ROI from preprocessing folder, unfortunately
    listing = load(fullfile(foldertemp,'DICOM listing.mat')).listing;
    ROI_idx = load(fullfile(foldertemp,listing(1).AcquisitionDate, 'PreProcessing',sprintf('PreProcessing - %s.mat', listing(end).name))).ROI_idx;
    ROI_wind = fliplr([min(ROI_idx, [], 1)-20; max(ROI_idx, [], 1)+20]);
    
    % extract CEUS and NCUS images
    CEUSim = S_pixelwise(ii).CEUS_ImageFilt;
    NCUSim = S_pixelwise(ii).velim;

    % morph processing to remove ASA and PSV
%     spineROI = imerode(imclose(imopen(CEUSim>0.02,strel('disk',3)),strel('disk',60)), strel('disk',20)); % also since I copied from impact analysis
    spineROI = S_pixelwise(ii).SpineROI;

%     vess_bounds = imdilate(CEUS_center, strel('disk',5))~=CEUS_center;
    
    % define injury as line drawn using draw_injury.m
%     [injLocY,injLocX] = bresenham(listing(1).InjuryROI(1),listing(1).InjuryROI(3),listing(1).InjuryROI(2),listing(1).InjuryROI(4));
    delxy = listing(find([listing.Mode]=="SMI" & [listing.NumberOfFrames]>1,1)).delxy;
    
    % distance function per pixel
%     injIm = zeros(size(CEUSim));
%     injLocRound = round([injLocX,injLocY]-ROI_wind(1,:));

    %
%     injIm(sub2ind(size(injIm),injLocRound(:,1),injLocRound(:,2))) = 1;
%     distIm = bwdist(injIm);
%     distIm(:,1:mean(injLocRound(:,2))) = -distIm(:,1:mean(injLocRound(:,2)));
    distIm = flow_through_distance(listing(1).InjuryROI-fliplr(ROI_wind(1,:)), NCUSim, delxy);


    % no processing of pixels not in ROI
    CEUSim(~logical(spineROI)) = NaN;
    NCUSim(~logical(spineROI)) = NaN;
    distIm(~logical(spineROI)) = NaN;

        
    % distance of each pixel within ROI
    distCEUS = distIm(CEUSim>0);
    distNCUS = distIm(NCUSim>0);
    
    % measured flow of each pixel in ROI
    CEUSflo = CEUSim(CEUSim>0);
    NCUSflo = NCUSim(NCUSim>0);

%     % pixelwise correlation
%     figure(4)
%     hold on
%     CEUSpix = reshape(CEUSim,[],1);
%     NCUSpix = reshape(NCUSim,[],1);
%     scatter(CEUSpix, NCUSpix)
%     [rho, pval] = corr(CEUSpix(CEUSpix>0 & NCUSpix>0),NCUSpix(CEUSpix>0 & NCUSpix>0),'rows','complete','type','Spearman');
%     fprintf('%s rho: %.02f; p: %d\n', plottitle, rho, pval)
    
    % store variables back in structure
%     S_pixelwise(ii).SpineROI = spineROI; % since I copied from impact analysis
    S_pixelwise(ii).DistanceImage = distIm;
    S_pixelwise(ii).DistanceCEUS = distCEUS;
    S_pixelwise(ii).FlowCEUS = CEUSflo;
    S_pixelwise(ii).DistanceNCUS = distNCUS;
    S_pixelwise(ii).NCUSFlow = NCUSflo;
    
    % process pre and post for visual confirmation
    [distCEUS, Isort] = sort(distCEUS);
    CEUSflo = CEUSflo(Isort);
    
    % rough smooth pre and post for visual inspection
    filt_size = 1000;
    avfilt = ones(filt_size,1)/filt_size;
    distCEUS = conv(distCEUS, avfilt, 'valid');
    CEUSflo = conv(CEUSflo, avfilt, 'valid');
    
    % plot quant
    figure(2)
    subplot(length(S_pixelwise)/2,2,ii);
    plot(distCEUS, CEUSflo, 'LineStyle', 'none', 'Marker', '.')
    ylim([0, 0.06])
    title(S_pixelwise(ii).Name(1:16));
    
    % display ROI
    figure(3)
    subplot(length(S_pixelwise)/2,2,ii);
%     centerpoint(injLocRound(1):(injLocRound(1)+10),injLocRound(2):(injLocRound(2)+10))=0;
%     imshow(isnan(CEUSim)+injIm,[0,0.1])
    imshow(~isnan(CEUSim)+NCUSim,[0 3])
    title(S_pixelwise(ii).Name(1:16));
end

% save(fullfile(figure_path,'Spatial Distribution', 'Contrast Comparison', save_name), 'S_pixelwise', '-append');

%% Laod data for further processing

load("C:\Users\Denis\OneDrive - Johns Hopkins\Figures\Spatial Distribution\Contrast Comparison\230609 Spatial Data.mat");



%% Section 5: Calculate correlations and draft plot for Fig. 4E

rhos = [];
ps = [];
ns = [];

X = [];
Y = [];

disp('Correlations:')

% uncomment for single rat plot
% rat2plot = '230117 Rat2';
% BU_S_pixelwise = S_pixelwise;
% S_pixelwise = S_pixelwise(contains({S_pixelwise.Name}, rat2plot));

% uncomment for entire dataset
try
    S_pixelwise = BU_S_pixelwise;
end


% Initialize structures for storing information
Sce_post  = struct('Time',[], ...
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

Sce_pre = struct('InjuryHeight',[], ...
                    'CranPenHeight',[], ...
                    'CaudPenHeight',[], ...
                    'CranNormHeight',[], ...
                    'CaudNormHeight',[]);

Snc_pre = Sce_pre;
Snc_post = Sce_post;


for rat = [2:2:length(S_pixelwise)]
    plottitle = S_pixelwise(rat).Name(1:11);
    
    % raw NCUS distance and flo variables
    predist = S_pixelwise(rat).DistanceNCUS;
    postdist = S_pixelwise(rat-1).DistanceNCUS;    
    preflo = S_pixelwise(rat).NCUSFlow;
    postflo = S_pixelwise(rat-1).NCUSFlow;

    distImpre = S_pixelwise(rat).DistanceImage;
    distImpost = S_pixelwise(rat-1).DistanceImage;

    % binned variables
    dist = -12:0.1:12;
    avpre = nan(size(dist));
    avpost = nan(size(dist));

    dist_edges = [-7.5 7.5]';

    [~, dist_inds] = find(dist==dist_edges);

    % calculate NCUS parameters by distance bin
    for ii = 1:length(dist)
        d = dist(ii);

        % V parameter (Vel Index)
%         avpre(ii) = mean(preflo(predist >=d & predist < d+1));
%         avpost(ii) = mean(postflo(postdist >=d & postdist < d+1));

        % A parameter (Area)
%         avpre(ii) = length(preflo(predist >=d & predist < d+1));
%         avpost(ii) = length(postflo(postdist >=d & postdist < d+1));
        
        % I parameter (Area-Adjusted Vel Index)
        avpre(ii) = mean(preflo(predist >=d & predist < d+1))*length(preflo(predist >=d & predist < d+1))/sum(distImpre >=d & distImpre <d+1, "all");
        avpost(ii) = mean(postflo(postdist >=d & postdist < d+1))*length(postflo(postdist >=d & postdist < d+1))/sum(distImpost >=d & distImpost <d+1, "all");
        
        % Old I parameter (incorrect)
%         avpre(ii) = mean(preflo(predist >=d & predist < d+1))*length(preflo(predist >=d & predist < d+1));
%         avpost(ii) = mean(postflo(postdist >=d & postdist < d+1))*length(postflo(postdist >=d & postdist < d+1));
    end
    
    % Export this variable for Prism plotting. (Fig. 4C,D)
%     NCUSvar = (avpost-avpre)./avpre;
%     NCUSvar = avpost;
    NCUSvar = avpost(dist_inds(1):dist_inds(2));
    NCUSvar(isnan(NCUSvar)) = 0;


%     % values for pre and post separately
    [pre_values, post_values] = find_injury_regions(avpre, avpost, 0, dist,0.03);
    lims_val = [0 0.8];

    % store values
    Snc_pre(rat/2) = pre_values;
    Snc_post(rat/2) = post_values;


    % raw CEUS distance and flo variables
    predist = S_pixelwise(rat).DistanceCEUS;
    postdist = S_pixelwise(rat-1).DistanceCEUS;    
    preflo = S_pixelwise(rat).FlowCEUS;
    postflo = S_pixelwise(rat-1).FlowCEUS;
    
    % binned variables
    avpre = nan(size(dist));
    avpost = zeros(size(dist));

    % calculate CEUS parameters by distance bin
    for ii = 1:length(dist)
        d = dist(ii);

        % V parameter (Vel Index)
%         avpre(ii) = mean(preflo(predist >=d & predist < d+1));
%         avpost(ii) = mean(postflo(postdist >=d & postdist < d+1));

        % A parameter (Area)
%         avpre(ii) = length(preflo(predist >=d & predist < d+1));
%         avpost(ii) = length(postflo(postdist >=d & postdist < d+1));
        
        % I parameter (Area-Adjusted Vel Index)
        avpre(ii) = mean(preflo(predist >=d & predist < d+1))*length(preflo(predist >=d & predist < d+1))/sum(distImpre >=d & distImpre <d+1, "all");
        avpost(ii) = median(postflo(postdist >=d & postdist < d+1))*length(postflo(postdist >=d & postdist < d+1))/sum(distImpost >=d & distImpost <d+1, "all");

%         disp(sum(distImpost >=d & distImpost <d+1, "all"));
        
        % Old I parameter (incorrect)
%         avpre(ii) = mean(preflo(predist >=d & predist < d+1))*length(preflo(predist >=d & predist < d+1));
%         avpost(ii) = mean(postflo(postdist >=d & postdist < d+1))*length(postflo(postdist >=d & postdist < d+1));
    end
    
    
%     % values for pre and post separately
    [pre_values, post_values] = find_injury_regions(avpre, avpost, 0, dist,0.003);
    lims_val = [0 0.8];

    % store values
    Sce_pre(rat/2) = pre_values;
    Sce_post(rat/2) = post_values;


    % export this variable for Prism plotting. (Fig. 4C,D)
%     CEUSvar = (avpost-avpre)./avpre;
%     CEUSvar = avpost;
    CEUSvar = avpost(dist_inds(1):dist_inds(2));
    CEUSvar(isnan(CEUSvar)) = 0;
    
    
    Y = [Y CEUSvar];
    X = [X NCUSvar];
    
    % Calculate correlations and fit linear models
    [rho, pval] = corr(CEUSvar',NCUSvar','rows','complete','type','Spearman');
    
    rhos = [rhos; rho];
    ps = [ps; pval];
    ns = [ns; length(CEUSvar)];

    mdl = fitlm(NCUSvar',CEUSvar');

    fprintf('%s rho: %.02f; p: %d; R^2: %02f\n', plottitle, rho, pval, mdl.Rsquared.Ordinary)

    
    %%%%%% Sample plotting of linear models
%     figure(5)
%     subplot(ceil(length(S_pixelwise)/4),2,rat/2);
%     imshow(distImpost, []);

    if rat == 6
        figure(2)
        plot(mdl, 'Marker', 'none', 'LineWidth', 1.5)
        hold on
        
        ylabel({'Area-Adjusted','Decay Constant (1/s)'})
        xlabel({'Area-Adjusted Velocity Index (AU)'})
        xlim([0 0.6])
        ylim([0 0.03])
        
        hFIT = findobj(gca,'DisplayName','Fit');
        [hFIT.LineWidth]=deal(2);
        [hFIT.Color]=deal([0.2 0.2 0.2 0.7]);
        
        hBounds = [findobj(gca,'DisplayName','') findobj(gca,'DisplayName','Confidence bounds')];
        [hBounds.LineWidth]=deal(1);
        [hBounds.LineStyle]=deal('--');
        [hBounds.Color]=deal([1 0 0 0.6]);

        scatter(NCUSvar, CEUSvar, 'Marker','.', 'MarkerEdgeColor', [0.4 0.4 0.4])

%         legend({'','Fitted Model', 'Conf. Limits','','Points'})
        legend off
        title('')
        set(gca, 'LineWidth', 1.5)
        set(gca,'FontWeight','bold')
        set(figure(2), 'Position', [360,738,229,184]);
        hold off
    end

    figure(6);
%     scatter(CEUSvar,NCUSvar)
    plot(mdl, 'Marker', 'none', 'LineWidth', 1.5)
%     plot(mdl)
    title(plottitle)
    ylabel({'Area-Adjusted','Decay Constant (1/s)'})
    xlabel({'Area-Adjusted Velocity Index (AU)'})
    xlim([0 0.6])
    ylim([0 0.05])
    hold on


    %%%%%%% Shows sample Fig. 4D traces
    figure(1)
    subplot(3,1,1), plot(dist(dist_inds(1):dist_inds(2)), CEUSvar)
    xline(dist(Sce_post(rat/2).UmbraLoc));
    xline(dist(Sce_post(rat/2).CranLoc));
    xline(dist(Sce_post(rat/2).CaudLoc));
    ylabel({'Area-Adjusted','Decay Constant (px/s)'})
    title('Contrast-Enhanced Flow Measurement')

    subplot(3,1,2), plot(dist(dist_inds(1):dist_inds(2)), (NCUSvar))
    xline(dist(Snc_post(rat/2).UmbraLoc));
    xline(dist(Snc_post(rat/2).CranLoc));
    xline(dist(Snc_post(rat/2).CaudLoc));
    xlabel('Distance from Injury Epicenter (mm)')
    ylabel({'Area-Adjusted','Velocity Index (AU)'})
    title('Non-Contrast Flow Measurement')

    subplot(3,1,3), plot(dist(dist_inds(1):dist_inds(2)), log(NCUSvar./CEUSvar))
    xlabel('Distance from Injury Epicenter (mm)')
    ylabel({'Log Ratio'})
    title('Log of Ratio between the two')
    
    ylim([0 5]);
    pause(3)

end

figure(6)
hold off

hFIT = findobj(gca,'DisplayName','Fit');
[hFIT.LineWidth]=deal(2);
[hFIT.Color]=deal([0.2 0.2 0.2 0.7]);

hBounds = [findobj(gca,'DisplayName',''); findobj(gca,'DisplayName','Confidence bounds')];
[hBounds.LineWidth]=deal(1);
[hBounds.LineStyle]=deal('--');
[hBounds.Color]=deal([1 0 0 0.6]);

title('')
legend({'','Fitted Model', 'Conf. Limits'})
set(gca, 'LineWidth', 1.5)
set(gca,'FontWeight','bold')
set(figure(6), 'Position', [360,738,229,184]);

% set(figure(1), 'Position', [360,660,309,262]);
%% Save Figures

exportgraphics(figure(2), "C:\Users\Denis\OneDrive - Johns Hopkins\230124 Spatial Distribution Paper\Figures\4. CEUS vs non-CEUS\230713 (230117 Rat2) corr.pdf",'ContentType','vector');
exportgraphics(figure(6), "C:\Users\Denis\OneDrive - Johns Hopkins\230124 Spatial Distribution Paper\Figures\4. CEUS vs non-CEUS\230713 All rats corr.pdf",'ContentType','vector');


%% Section 6: Calculation of rho and p from multiple measurements

% all rats
rhos2use = rhos;
ns2use = ns;

% excluding flat relationship
% rhos2use = rhos([1:4 7]);
% ns2use = ns([1:4, 7]);

% z-transformation of rhos and calculation of average
av_z = sum(atanh(rhos2use).*ns2use)./sum(ns2use);
av_rho = tanh(av_z);
av_p = 1-normcdf(av_z*sqrt(sum(ns2use)));

fprintf('Total rho: %.02f; p: %d\n', av_rho, av_p)



%% Figure 4A, B: representative images of NCUS and CEUS

% figure path to save
figure_path = 'C:\Users\Denis\OneDrive - Johns Hopkins\Figures';

imwrite(12*S_pixelwise(1).CEUS_ImageFilt(25:end-50,50:end-20), fullfile(figure_path,'Spatial Distribution', 'Contrast Comparison', '220117 Rat2 CEUS.png'),'png');
imwrite(S_pixelwise(1).velim(25:end-50,50:end-20), fullfile(figure_path,'Spatial Distribution', 'Contrast Comparison', '220117 Rat2 NCUS.png'),'png');



%% Sample plot: Plot chosen pixelwise results


% X = [];
% 
% for rat = 2:2:length(S_pixelwise)
% %     predist = S_pixelwise(rat).DistanceNCUS;
% %     postdist = S_pixelwise(rat-1).DistanceNCUS;
% %     
% %     preflo = S_pixelwise(rat).NCUSFlow;
% %     postflo = S_pixelwise(rat-1).NCUSFlow;
% 
%     predist = S_pixelwise(rat).DistanceCEUS;
%     postdist = S_pixelwise(rat-1).DistanceCEUS;
%     
%     preflo = S_pixelwise(rat).FlowCEUS;
%     postflo = S_pixelwise(rat-1).FlowCEUS;
% 
%     plottitle = S_pixelwise(rat).Name(1:11);
%     
%     dist = -12:0.1:12;
%     avpre = nan(size(dist));
%     avpost = nan(size(dist));
% 
%     for ii = 1:length(dist)
%         d = dist(ii);
% %         avpre(ii) = mean(preflo(predist >=d & predist < d+1));
% %         avpost(ii) = mean(postflo(postdist >=d & postdist < d+1));
% %         avpre(ii) = length(preflo(predist >=d & predist < d+1));
% %         avpost(ii) = length(postflo(postdist >=d & postdist < d+1));
% 
%         avpre(ii) = mean(preflo(predist >=d & predist < d+1))*length(preflo(predist >=d & predist < d+1));
%         avpost(ii) = mean(postflo(postdist >=d & postdist < d+1))*length(postflo(postdist >=d & postdist < d+1));
%     end
%     
%     X = [X, (avpost-avpre)./avpre];
% 
%     figure(4);
%     subplot(length(S_pixelwise)/2,1,rat/2);
% %     plot(dist, (avpost-avpre)./avpre)
%     plot(dist, (avpost))
%     title(plottitle)
% %     ylim([0 3000])
%     xlim([-10 10])
% end


