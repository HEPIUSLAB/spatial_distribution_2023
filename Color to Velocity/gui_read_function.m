% Function: gui_read_funtion
%
% Purpose: Extract information from the GUI
%
% Input parameters:
%   A: scanner_parameters_selection
%
% Output parameters:
%   scanner_parameters: struct (5 fields)
%       scanner_parameters.curr_colorbarbounds: double
%       scanner_parameters.curr_dataset: char
%       scanner_parameters.curr_img_mode: char
%       scanner_parameters.curr_scanner: char
%       scanner_parameters.date: char
%
% Created by: Kelley Kempski (kkempski@jhmi.edu)
% Lasted edited: 12/2/2021

function scanner_parameters = gui_read_function(A)
flag = true;
while flag % This loop prevents the code from continuing until the GUI is closed
    if ~A.UIFigure.Visible
        flag = false;
    end
    pause(0.25);
end

clear flag

MFI = A.MFIButton.Value; % Extract the value of the MFI button (on/off)
CDI = A.CDIButton.Value; % Extract the value of the CDI button (on/off)
PDI = A.PDIButton.Value; % Extract the value of the PDI button (on/off)
SMI = A.SMIButton.Value; % Extract the value of the SMI button (on/off)
ADF = A.ADFButton.Value; % Extract the value of the PDI button (on/off)

% The imaging modality whose button is 'on' is set as the current imaging
% modality (curr_img_mode)
if MFI 
    curr_img_mode = 'MFI';
elseif CDI
    curr_img_mode = 'CDI';
    cmap = load('CDI Color Map.mat');
    cmap = cmap.cdi_map;
elseif PDI
    curr_img_mode = 'PDI';
    cmap = load('PDI Color Map.mat');
    cmap = cmap.pdi_map;
elseif SMI
    curr_img_mode = 'SMI';
    cmap = load('SMI Color Map.mat');
    cmap = cmap.smi_map;
elseif ADF
    curr_img_mode = 'ADF';
    cmap = load('ADF Color Map.mat');
    cmap = cmap.adf_map;
end

clear MFI SMI CDI PDI ADF

Canon = A.CanonButton.Value; % Extract the value of the Canon button (on/off)
Philips = A.PhilipsButton.Value; % Extract the value of the Philips button (on/off)

% The scanner whose button is 'on' is set as the current scanner
% (curr_scanner)
if Canon
    curr_scanner = 'Canon';
elseif Philips
    curr_scanner = 'Philips';
end

clear Canon Philips

% Extract the value of the colorbar bounds
curr_colorbarbounds = str2double(A.ColorbarBoundscmsEditField.Value); 

% Extract the name of the current dataset
curr_dataset = A.DatasetNameEditField.Value;

% Extract the date from the GUI and convert to a new format (MM-DD-YYYY)
date = char(A.SurgeryDateDatePicker.Value);
hyphens = find(date=='-');
day = str2double(date(1:hyphens(1)-1));
month_name = date(hyphens(1)+1:hyphens(2)-1); 
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
month = find(contains(months,month_name));
year = str2double(date(hyphens(2)+1:end));
date = sprintf('%d-%d-%d',month,day,year);

% Add all relevant info to scanner_parameters structure
scanner_parameters = struct;
scanner_parameters.curr_colorbarbounds = curr_colorbarbounds;
scanner_parameters.curr_dataset = curr_dataset;
scanner_parameters.curr_img_mode = curr_img_mode;
scanner_parameters.curr_scanner = curr_scanner;
scanner_parameters.date = date;
scanner_parameters.cmap = cmap;
end