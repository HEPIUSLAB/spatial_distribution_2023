clear; close all;
% Purpose: Run this script to generate all required images/
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% 


%% Figure 2

% replace with the location of data folder
data_drive = 'D:\Data\';

% load chosen rat (230119 Rat2)
load(fullfile(data_drive,'230119 Rat SCBF and CEUS\Spatial Distribution\230119 Rat2 Post CEUS 1 px.mat'));
load(fullfile(data_drive,'230119 Rat SCBF and CEUS\Spatial Distribution\230119 1 px CEUS.mat'));

% choose save location
save_path = 'C:\Users\Denis\OneDrive - Johns Hopkins\230124 Spatial Distribution Paper\Figures\3. CEUS Quantification';

% load images and normalize for better visibility
unfilt_image = S_CEUS_im(3).Image/max(S_CEUS_im(3).ImageFilt,[],'all');
repr_image = S_CEUS_im(3).ImageFilt/max(S_CEUS_im(3).ImageFilt,[],'all');

% image enhancement for visibility and saving
% imwrite(unfilt_image,fullfile(save_path, 'unfilteredCEUS.png'));
% imwrite(repr_image,fullfile(save_path, 'filteredCEUS.png'));

% individual pixel traces
trace1 = RawFlow(198,627,:);
trace2 = RawFlow(196,620,:);

% run pixel time-course through image generating function (plots if size of
% trace 1 is 1x1xn
TIC_to_image(t, trace1);

% figure aesthetics
set(gcf, 'Position', [228,533,455,195])
set(gca,'FontWeight','bold')
set(gca, 'LineWidth', 1.5)
hline = findobj(gcf, 'type', 'line');
set(hline(1),'LineWidth',2)
set(hline(2),'LineWidth',2)
set(hline(2),'Color','blue')
set(hline(3),'Color',[0 0.5 0.8])
set(hline(1),'LineStyle','--')
legend({'Unfiltered Signal', 'Filtered Signal', 'Fitted Model'})
ylabel('Intensity (AU)')
xlabel('Time (s)')

% save raw CHI images with original colormap
for tt = 0:20:t(end)
    [~, ind] = min(abs(t-tt));
    cmap = permute(colormap(copper), [1,3,2]);
    rawIm = round(255*RawFlow(:,:,ind))+1;
    rawIm(1:158, 1000:end) = 1;
    rawIm = reshape(cmap(rawIm,:), [size(rawIm),3]);
    
    imwrite(rawIm,fullfile(save_path, 'CHI images', sprintf('CHI %ds.png', tt)));
end

%% Figure 4

% open this script for more granularity
CEUS_NCUS_pixwise_comparison

%% Figure 3,5,6 

% open this script for more granularity
impact_force_analysis

%% Figure 7

% open this script for more granularity
time_evolution
