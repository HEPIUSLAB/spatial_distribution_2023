% Function: check_data_path
%
% Purpose: to replace some file paths in an old analysis file
%
% Input parameters: 
%       files: cell array of strings
%       vess_files: cell array of strings
%       imfiles: cell array of strings
%       root_path: string
%
% Output parameters:
%       files: cell array of strings
%       vess_files: cell array of strings
%       imfiles: cell array of strings
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function [files,vess_files,imfiles] = check_data_path(files,vess_files,imfiles, root_path)

    if ~exist(files{1}, 'file')    
    
        for ff = 1:length(files)
            old_file_path = files{ff};
            branch_path = old_file_path(strfind(old_file_path,'Data')+4:end);
            files{ff} = fullfile(root_path, branch_path);
            files{ff} = strrep(files{ff},'Rat SCBF AH','221208 Rat SCBF Compiled');
        end
    
        for ff = 1:length(vess_files)
            old_file_path = vess_files{ff};
            branch_path = old_file_path(strfind(old_file_path,'Data')+4:end);
            vess_files{ff} = fullfile(root_path, branch_path);
            vess_files{ff} = strrep(vess_files{ff},'Rat SCBF AH','221208 Rat SCBF Compiled');
        end
    
        for ff = 1:length(imfiles)
            old_file_path = imfiles{ff};
            branch_path = old_file_path(strfind(old_file_path,'Data')+4:end);
            imfiles{ff} = fullfile(root_path, branch_path);
            imfiles{ff} = strrep(imfiles{ff},'Rat SCBF AH','221208 Rat SCBF Compiled');
        end
    end
end