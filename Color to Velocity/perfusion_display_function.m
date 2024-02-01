% Function: perfusion_display_function
%
% Purpose: Calculate mean velocity and mean intensity for each frame and
% save this data
%
% Input parameters:
%   scanner_parameters: struct
%       scanner_parameters.curr_colorbarbounds: double
%       scanner_parameters.curr_dataset: char
%       scanner_parameters.curr_img_mode: char
%       scanner_parameters.curr_scanner: char
%       scanner_parameters.date: char
%   ROI_idx: double
%
% Output parameters: None
%
% Created by: Kelley Kempski (kkempski@jhmi.edu)
% Last edited: 1/26/2022 (Denis Routkevitch droutke1@jhmi.edu)
%       Added looping and different folder selection mode
% Edit: 2/7/2022 added compatibility with DICOM splitting and new image
%       naming system
% Edit: 2/15/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%       added ROI_idx input for correct I calculations. Added velocity
%       input for faster runtime and easier folder-change compatibility
% Edit: 10/12/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%       changed ROI_idx processing pipeline to be correct

function perfusion_display_function(scanner_parameters, ROI_idx, velocity)
curr_dataset = scanner_parameters.curr_dataset; % Extract current dataset name
curr_img_mode = scanner_parameters.curr_img_mode; % Extract current imaging modality
date = scanner_parameters.date; % Extract surgery date

% forward compatibility for selecting different folder
if isfield(scanner_parameters, 'folder')
    base_folder = scanner_parameters.folder;
else 
    base_folder = '..';
end

% Define filename for the re-ordered velocity map
number_frames = size(velocity,3); % Extract the total number of frames


if strcmp(curr_img_mode,'CDI') || strcmp(curr_img_mode,'ADF') % If CDI or ADF, continue
    % Initialize velocity and intensity arrays for positive (red) and
    % negative (blue) directional flow
    v_red = zeros(1,number_frames); 
    v_blu = zeros(1,number_frames);
    I_red = zeros(1,number_frames);
    I_blu = zeros(1,number_frames);
elseif strcmp(curr_img_mode,'SMI') || strcmp(curr_img_mode,'PDI') % If SMI or PDI, continue
    % Initialize velocity and intensity arrays for nondirectional flow
    I_color = zeros(1,number_frames);
    v_color = zeros(1,number_frames);
end

% convert inds into linear format
ROI_inds = sub2ind(size(velocity(:,:,1)), ROI_idx(:,2), ROI_idx(:,1));

% f = waitbar(number_frames,{sprintf('Saving Quantitative Data - %s',curr_dataset),[sprintf('%0.1f',0),'%']}); % Initialize waitbar
for frame = 1:number_frames % Loop through each frame
%     waitbar(frame/number_frames,f,{sprintf('Saving Quantitative Data - %s',curr_dataset),[sprintf('%0.1f',(frame/number_frames)*100),'%']}); % Update waitbar
    curr_frame = squeeze(velocity(:,:,frame)); % Isolate current frame of velocity map
    
    if strcmp(curr_img_mode,'CDI') || strcmp(curr_img_mode,'ADF') % If CDI or ADF, continue
        % edited the area function to be correct
        A_red = sum(curr_frame(ROI_inds)>0); % Area of all red colored pixels
        A_blu = sum(curr_frame(ROI_inds)<0); % Area of all blue colored pixels

        A_ROI = size(ROI_idx,1); % Area of the ROI
        
        % extract only colored pixels
        curr_frame = curr_frame(ROI_inds);

        if A_red > 0 % If any red pixels are present in the current frame, continue
            v_red(frame) = mean(curr_frame(curr_frame>0)); % [cm/s] Mean of all red pixels
        end

        if A_blu > 0 % If any blue pixels are present in the current frame, continue
            v_blu(frame) = mean(curr_frame(curr_frame<0)); % [cm/s] Mean of all blue pixels
        end

        I_red(frame) = A_red*v_red(frame)/A_ROI; % [cm/s] Mean red flow intensity of the ROI
        I_blu(frame) = A_blu*v_blu(frame)/A_ROI; % [cm/s] Mean blue flow intensity of the ROI
    
    elseif strcmp(curr_img_mode,'SMI') || strcmp(curr_img_mode,'PDI') % If SMI or PDI, continue
        % edited the area function to be correct
        A_color = sum(curr_frame(ROI_inds)~=0); % Area of all colored pixels

        A_ROI = size(ROI_idx,1); % Area of the ROI

        if A_color > 0 % If any colored pixels are present in the current frame
            curr_frame = abs(curr_frame); % Convert all values to positive since nondirectional flow
            % edited for better calculation
            v_color(frame) = mean(nonzeros(curr_frame(ROI_inds))); % [cm/s] Mean of all colored pixels
        end

        I_color(frame) = A_color*v_color(frame)/A_ROI; % [cm/s] Mean red flow intensity of the ROI
    end
end

% Save mean velocity and intensity
if strcmp(curr_img_mode,'CDI') || strcmp(curr_img_mode,'ADF') 
    save(fullfile(base_folder,date,'Quantitative Data',sprintf('Matlab Quant - %s',curr_dataset)),'v_red','v_blu','I_red','I_blu','curr_img_mode');
elseif strcmp(curr_img_mode,'SMI') || strcmp(curr_img_mode,'PDI')
    save(fullfile(base_folder,date,'Quantitative Data',sprintf('Matlab Quant - %s',curr_dataset)),'v_color','I_color','curr_img_mode');
end
% close(f); % Close waitbar
end