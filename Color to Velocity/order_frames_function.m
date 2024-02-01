% Function: order_frames_function
%
% Purpose: Re-order frames in the velocity map. The directory orders frames
% 1, 10, 11, 12..., 2, 20, 21, ... This script reorders to 1, 2, 3...
%
% Input parameters:
%   base_folder: char (original folder to locate images)
%   velocity: double (matrix of velocity map - m x n x r, where m and n are
%   the rows and columns, respectively, and r is the number of frames)
%   number_frames: double
%   curr_dataset: double (number of dataset, e.g., 4)
%
% Output parameters:
%   velocity_ordered: double (matrix of velocity map- m x n x r)
%
% Created by: Kelley Kempski (kkempski@jhmi.edu)
% Last edited: 1/26/2022 (Denis Routkevitch droutke1@jhmi.edu)
%       Added looping and different folder selection mode
%       Added compatibility with frame downsampling

function velocity_ordered = order_frames_function(base_folder,velocity,number_frames,curr_dataset)
a = dir(fullfile(base_folder,'Images',curr_dataset)); % Create directory of files contained in the dataset folder
b = zeros(number_frames,1); % Initialize matrix for correct index numbers
f = waitbar(number_frames,{sprintf('Re-ordering Frames - %s',curr_dataset),[sprintf('%0.1f',0) '%']}); % Initialize waitbar
for i = 1:length(b) % Loop through each frame
    waitbar(i/number_frames,f,{sprintf('Re-ordering Frames - %s',curr_dataset),[sprintf('%0.1f',(i/number_frames)*100) '%']}); % Update waitbar
    curr_name = a(i+2).name; % Extract the name of the current file
    curr_name1 = curr_name(1:end-4); % Extract frame number from jpg filename
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
    number = str2double(curr_name1(idx(end)+1:end));
    b(i) = number; % Assign the correct frame number to the 'b' vector
end

b_sort = sort(b); % added for downsampling compatibility

velocity_ordered = zeros(size(velocity)); % Initialize re-ordered velocity matrix
for i = 1:length(b) % Loop through each frame
    velocity_ordered(:,:,find(b_sort == b(i))) = velocity(:,:,i); % Assign the correct velocity map in the correct frame location
end
close(f); % Close waitbar
end
