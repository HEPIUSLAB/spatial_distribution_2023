% Function: movie_saving_function
%
% Purpose: Create .avi movie file comparing scanner screenshots, calculated
% velocity map, and the quantitative plot
%
% Input parameters:
%   scanner_parameters
%       scanner_parameters.curr_colorbarbounds: double
%       scanner_parameters.curr_dataset: char
%       scanner_parameters.curr_img_mode: char
%       scanner_parameters.curr_scanner: char
%       scanner_parameters.date: char
%
% Output parameters: None
%
% Created by: Kelley Kempski (kkempski@jhmi.edu)
% Last edited: 1/26/2022 (Denis Routkevitch droutke1@jhmi.edu)
%       Added looping and different folder selection mode

function movie_saving_function(scanner_parameters)
curr_dataset = scanner_parameters.curr_dataset; % Extract current dataset
curr_colorbarbounds = scanner_parameters.curr_colorbarbounds; % Extract current colorbar bounds
curr_img_mode = scanner_parameters.curr_img_mode; % Extract current imaging modality
date = scanner_parameters.date; % Extract surgery date
cmap = scanner_parameters.cmap; % Extract colormap
cmap_centercol = squeeze(cmap(:,9,:));

% forward compatibility for selecting different folder
if isfield(scanner_parameters, 'folder')
    base_folder = scanner_parameters.folder;
else 
    base_folder = '..';
end

foldername = fullfile(base_folder,date,'Images',sprintf('%s',curr_dataset)); % Define current dataset foldername

frame_rate = 4; % Set frame rate for video file

% Define filename for re-ordered velocity map
filename = fullfile(base_folder,date,'Velocity Maps',sprintf('Velocity Map - ordered - %s.mat',curr_dataset));
file1 = load(filename); % Load velocity map

number_frames = file1.number_frames; % Extract number of frames
velocity_ordered = file1.velocity_ordered; % Extract velocity map
velocity_ordered_frame1 = velocity_ordered(:,:,1); % Isolate the first frame of the velocity map matrix

[row_idx,col_idx] = find(velocity_ordered_frame1); % Locate the flow in the velocities
row_range = min(row_idx)-50:max(row_idx)+50; % Define range of rows around the ROI for zoomed in display
col_range = min(col_idx)-50:max(col_idx)+50; % Define range of columns around the ROI for zoomed in display

filename = fullfile(base_folder,date,'Quantitative Data',sprintf('Matlab Quant - %s.mat',curr_dataset)); % Define filename for quantitative data
file2 = load(filename); % Load quantitative data

writerObj1 = VideoWriter(fullfile(base_folder,date,'Videos',sprintf('%s.avi',curr_dataset))); % Define video filename
writerObj1.FrameRate = frame_rate; % Set frame rate for the video
img = cell(number_frames,1); % Initialize cell for video frames
open(writerObj1); % Open video writer object

listing = dir(foldername); % Create directory for files in the current dataset folder

if strcmp(curr_img_mode,'CDI') % If CDI, continue    
    % Extract velocity and intensity values
    I_blu = file2.I_blu; 
    I_red = file2.I_red;
    v_red = file2.v_red;
    v_blu = file2.v_blu;

elseif strcmp(curr_img_mode,'ADF') % If ADF, continue    
    % Extract velocity and intensity values
    I_blu = file2.I_blu;
    I_red = file2.I_red;
    v_red = file2.v_red;
    v_blu = file2.v_blu;    

elseif strcmp(curr_img_mode,'SMI') % If SMI, continue
    % Extract velocity and intensity values
    I_color = file2.I_color;
    v_color = file2.v_color;

elseif strcmp(curr_img_mode,'PDI') % If PDI, continue
    % Extract velocity and intensity values
    I_color = file2.I_color;
    v_color = file2.v_color;
end

f = waitbar(number_frames,{sprintf('Initializing video - %s',curr_dataset),[sprintf('%0.1f',0),'%']}); % Initialize waitbar
for x = 1:number_frames % Loop through each frame
    waitbar(x/number_frames,f,{sprintf('Saving video - %s',curr_dataset),[sprintf('%0.1f',(x/number_frames)*100),'%']}) % Update wait bar
    tmp_name = listing(x+2).name; % Create temporary filename for current frame
    img_filename = fullfile(foldername,tmp_name); % Define filename for the current frame

    curr_name1 = tmp_name(1:end-4); % Extract frame number from jpg filename
    str_check = zeros(size(curr_name1));
    for y = 1:length(str_check)
        tmp = str2num(curr_name1(y));
        if isempty(tmp)
            str_check(y) = 0;
        else
            str_check(y) = 1;
        end
    end
    idx = find(str_check==0);
    curr_frame_number = str2double(curr_name1(idx(end)+1:end));

    curr_img = imread(img_filename); % Load current screenshot
    img{curr_frame_number} = curr_img(row_range,col_range,:); % Add zoomed in current screenshot to the img cell
end
close(f); % Close waitbar

f1 = figure('visible','off');set(gcf,'position',[229 558 2043 420],'color','w'); % Initialize figure parameters
f = waitbar(number_frames,{sprintf('Saving video - %s',curr_dataset),[sprintf('%0.1f',0),'%']}); % Initialize waitbar
for y = 1:number_frames % Loop through each frame
    waitbar(y/number_frames,f,{sprintf('Saving video - %s',curr_dataset),[sprintf('%0.1f',(y/number_frames)*100),'%']}) % Update wait bar

    subplot(1,3,1)
    imagesc(img{y}); axis image; axis off; % Display zoomed in screenshot
    
    subplot(1,3,2)
    curr_frame = squeeze(velocity_ordered(row_range,col_range,y)); % Isolate zoomed in velocity map
    imagesc(curr_frame,[-curr_colorbarbounds curr_colorbarbounds]); % Display zoomed in velocity map
    colormap(cmap_centercol); colorbar; axis image; axis off;
    set(gca,'fontsize',18);
    
    subplot(1,3,3)
    if strcmp(curr_img_mode,'CDI') || strcmp(curr_img_mode,'ADF') % If CDI or ADF, continue
        plot(1:number_frames,v_red,'r','Linewidth',1.5); hold on % Display quantitative plot
        plot(1:number_frames,abs(v_blu),'b','linewidth',1.5);

    elseif strcmp(curr_img_mode,'SMI') || strcmp(curr_img_mode,'PDI') % If SMI or PDI, continue
        plot(1:number_frames,v_color,'color',[0.5 0 0.5],'Linewidth',1.5); hold on; % Display quantitative plot
    end

    xline(y,'Linewidth',1.5,'color','k'); hold off % Display vertical line showing which part of the quantitative 
                                                   % plot corresponds to the currently displayed screenshot and 
                                                   % velocity map
    xlabel('Frame'); ylabel('Velocity (cm/s)')
    set(gca,'fontsize',18);
    frame = getframe(f1);
    writeVideo(writerObj1,frame); % Write current figure to the video
end

close(writerObj1); % Close video object
close(f); % Close waitbar
end