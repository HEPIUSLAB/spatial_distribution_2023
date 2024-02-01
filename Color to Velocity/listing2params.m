% Function: struct2params
%
% Purpose: Translate information from listing structure (from DR code) to
% scanner_parameters structure (for KK code)
%
% Input parameters:
%   listing: struct (1 entry from listing in outer scripts)
%
% Output parameters:
%   scanner_parameters: struct (5 fields)
%       scanner_parameters.curr_colorbarbounds: double
%       scanner_parameters.curr_dataset: char
%       scanner_parameters.curr_img_mode: char
%       scanner_parameters.curr_scanner: char
%       scanner_parameters.date: char
%
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% Last edited: 1/25/2021

function scanner_parameters = listing2params(listing)

modelist = {'MFI', 'CDI', 'PDI', 'SMI', 'ADF', 'BMODE', 'SWE', 'CHI'};
mode = listing.Mode;

whichmode = strcmp(modelist, mode);

% Choose colormap or throw error if wrong type of modality
if whichmode(1)
    % placeholder
elseif whichmode(2)
    cmap = load('CDI Color Map.mat');
    cmap = cmap.cdi_map;
elseif whichmode(3)
    cmap = load('PDI Color Map.mat');
    cmap = cmap.pdi_map;
elseif whichmode(4)
    cmap = load('SMI Color Map.mat');
    cmap = cmap.smi_map;
elseif whichmode(5)
    cmap = load('ADF Color Map.mat');
    cmap = cmap.adf_map;
elseif whichmode(6)
    throw(MException('funcInput:nonBloodFlow', ['This DICOM is B-mode and ' ...
        'not suitable for this analysis']));
elseif whichmode(7)
    throw(MException('funcInput:nonBloodFlow', ['This DICOM is SWE and ' ...
        'not suitable for this analysis']));
elseif whichmode(8)
    cmap = load('CHI Color Map.mat');
    cmap = cmap.chi_map;
else
    throw(MException('funcInput:nonBloodFlow', ['This DICOM is spectral and ' ...
        'not suitable for this analysis']));
end



scanner_parameters = struct;
scanner_parameters.curr_colorbarbounds = listing.bounds;
scanner_parameters.curr_dataset = listing.name;
scanner_parameters.curr_img_mode = mode;
scanner_parameters.curr_scanner = 'Canon';
scanner_parameters.date = listing.AcquisitionDate;
scanner_parameters.folder = listing.folder;
scanner_parameters.cmap = cmap;


end