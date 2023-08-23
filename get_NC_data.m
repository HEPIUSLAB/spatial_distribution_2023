% Function: get_NC_data
%
% Purpose: 
%
% Input parameters:
%   data_path: string
%   datasets: cell array of strings
%   good_data: struct of rats to analyze
%
% Output parameters:
%   S_FM: struct with rat single-vessel data
%   S_NCUS: struct with pixelwise analysis
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function [S_FM, S_NCUS] = get_NC_data(data_path, datasets, good_data)
    
    % initialize data structures
    S_FM = struct();
    S_NCUS = struct('skel',[],'skelorig',[],'velim',[],'numvess',[],'image',[],'ROI',[]);
    ratnum = 1;
    
    % run through each experiment date
    for ii = 1:length(datasets)
        dset = datasets{ii};
        
        % find the appropriate rats
        goodrats = good_data(strcmp({good_data.Date},dset(1:6))).Rats;
        
        % for each rat
        for jj = 1:length(goodrats)
            rat = goodrats{jj};
            
            % load impact data
            try
                impact_force = import_impact(fullfile(data_path, dset, 'Impacts', sprintf('%s %s.exp',dset(1:6),rat)));
            catch
                impact_force = import_impact(fullfile(data_path, dset, 'Impact Data', sprintf('%s.exp',rat)));
            end
            
            % load vessel links file
            a = dir(fullfile(data_path, dset,'Spatial Distribution'));
            analysis_file = {a(contains({a.name}, rat) & contains({a.name}, 'links')).name};
            if isempty(analysis_file)
                analysis_file = {a(contains({a.name}, rat(1:6))).name};
                analysis_file = analysis_file(1);
            end
            plotname = analysis_file{1};
            disp(plotname)
            S_path = fullfile(data_path, dset,'Spatial Distribution', plotname);
            load(S_path)
            plotname = plotname(strfind(plotname,'202'):strfind(plotname,'202')+8);
        
            % find number of separate vessels
            [bins, binsizes] = conncomp(graph(Adj_revise));
            numbins = length(binsizes);
            
            % edit file paths if on different computer
            [files,vess_files,imfiles] = check_data_path(files,vess_files,imfiles, data_path);
    
            
            % load injury location and pixel size parameters
            listing = load(fullfile(files{1},'..','..','..','DICOM listing.mat')).listing;
            [injLocY,injLocX] = bresenham(listing(1).InjuryROI(1),listing(1).InjuryROI(3),listing(1).InjuryROI(2),listing(1).InjuryROI(4));
            injLoc_pre = [injLocX,injLocY];
            delxy = listing(find([listing.Mode]=="SMI" & [listing.NumberOfFrames]>1,1)).delxy;
            
            % load this for future use
            ROI_idx = load(fullfile(files{1},'..','..','PreProcessing',sprintf('PreProcessing - %s.mat', listing(end-1).name))).ROI_idx;
            ROI_wind_pre = fliplr([min(ROI_idx, [], 1)-20; max(ROI_idx, [], 1)+20]);
                    
            % load this variable to find amount of vessels in pre
            load(vess_files{1});
            lastpre = length(vessels);

                        
            
            
                    
            % merge vessels by connected component in post
            load(vess_files{2});

            % load injury location and pixel size parameters
            listing = load(fullfile(files{2},'..','..','..','DICOM listing.mat')).listing;
            [injLocY,injLocX] = bresenham(listing(1).InjuryROI(1),listing(1).InjuryROI(3),listing(1).InjuryROI(2),listing(1).InjuryROI(4));
            injLoc_post = [injLocX,injLocY];
            delxy = listing(find([listing.Mode]=="SMI" & [listing.NumberOfFrames]>1,1)).delxy;
            
            % load this for future use
            ROI_idx = load(fullfile(files{2},'..','..','PreProcessing',sprintf('PreProcessing - %s.mat', listing(end-1).name))).ROI_idx;
            ROI_wind_post = fliplr([min(ROI_idx, [], 1)-20; max(ROI_idx, [], 1)+20]);

            
            
            
            % load data into structure
            S_FM(ratnum).Rat = sprintf('%s%s',plotname,rat);
            S_FM(ratnum).ImpactForce = impact_force;
            S_FM(ratnum).delxy = delxy;
            S_FM(ratnum).InjuryLocationPre = injLoc_pre;
            S_FM(ratnum).ROIWindowPre = ROI_wind_pre;
            S_FM(ratnum).InjuryLocationPost = injLoc_post;
            S_FM(ratnum).ROIWindowPost = ROI_wind_post;
            
            
            ratnum = ratnum+1;
            
            % store pixelwise data in format compatible with Section 2
            try
                S_NCUS = [S_NCUS flip(S)];
            catch
                S = backfill_S(S, files, vess_files);
                save(S_path, 'S', '-append');
                S_NCUS = [S_NCUS flip(S)];
            end    
        end
    
    end
    
    % trim empty beginning
    S_NCUS = S_NCUS(2:end);

end