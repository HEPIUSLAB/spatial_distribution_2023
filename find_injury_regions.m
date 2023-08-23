% Function: find_injury_regions
%
% Purpose: to calculate injury regions from flow surfaces
%
% Input parameters: 
%       all_rats: struct
%       dist: double vector
%
% Output parameters:
%       pre_values: struct
%       post_values: struct
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function [pre_values, post_values] = find_injury_regions(avpre, avpost, time, dist, prom, filt_size)

    if nargin == 4
        prom = 0.05;
        filt_size = 40;
    elseif nargin == 5
        filt_size = 40;
    end
        
    % initialize pre- and post-injury structure
    post_values = struct('Time',[], ...
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

    pre_values = struct('InjuryHeight',[], ...
                        'CranPenHeight',[], ...
                        'CaudPenHeight',[], ...
                        'CranNormHeight',[], ...
                        'CaudNormHeight',[]);
    

    
    
    % reset temporary storage variables
    [penumbra, transition] = deal(nan(length(time),2));
    umbra_loc = nan(length(time),1);
    injBin = zeros(size(avpost));
        
    
    % unfortunately the for-loop was easiest way to write this
    for tt = 1:size(avpost, 1)

        % flow at each time point for spatial analysis
%         sample_flow = (avpost(tt,:)-avpre)./avpre;
        sample_flow = avpost(tt,:);
        sample_flow(isnan(sample_flow)) = 0;
%         findpeaks(sample_flow, dist);
        
        % locate peaks to determine penumbra
        [~, penumbra_locs] = findpeaks(sample_flow, 'MinPeakProminence',prom, 'MaxPeakWidth',200);
        peak_dist = dist(penumbra_locs);
        
        % choose only the two closest peaks on either side of 0. Skip any
        % points that fail
        try penumbra_locs = penumbra_locs(peak_dist==min(peak_dist(peak_dist>0)) | peak_dist==max(peak_dist(peak_dist<0)));
        catch penumbra_locs = [penumbra_locs(end), length(dist)];
        end

        if isempty(penumbra_locs)
            continue
        end
        
        % determine injury epicenter as minimum point between penumbras
        % I2 is in case of non singular min
        [~, injury_ind] = min(sample_flow(penumbra_locs(1):penumbra_locs(2)));
        [~, I2] = min(fliplr(sample_flow(penumbra_locs(1):penumbra_locs(2))));
        
        % average injury center if non-singular
        if penumbra_locs(1)+injury_ind-1~=penumbra_locs(2)-I2+1
            injury_ind = round(mean([penumbra_locs(1)+injury_ind-1 penumbra_locs(2)-I2+1]));
        else
            injury_ind = penumbra_locs(1)+injury_ind-1;
        end
        
        % binary image of penumbra outline (unused)
        injBin(tt,penumbra_locs) = 1;
        
        % store all the indices
        penumbra(tt,:) = (penumbra_locs);
        umbra_loc(tt) = (injury_ind);
    end
    
    % median filtering to remove sudden jumps in location due to noise/
    % other peaks around
    if length(time)>10
        penumbra = medfilt1(penumbra, filt_size, 'omitnan', 'truncate');
        umbra_loc = medfilt1(umbra_loc, filt_size, 'omitnan', 'truncate');
    end
    
    if sum(isnan(penumbra), "all")>0
        return;
    end
    
    % in some injuries, the penumbra was not immediately apparent, but
    % developed later. This labels the penumbra at times that it is not
    % apparent
    step = find(abs(diff(penumbra(:,1)))>0.5, 1,'last');
    penumbra(1:step,1) = penumbra(step+1,1);
    step = find(abs(diff(penumbra(:,2)))>0.5, 1,'last');
    penumbra(1:step,2) = penumbra(step+1,2);
    
    % rounding for indices
    umbra_loc = round(umbra_loc);
    cran_loc = round(penumbra(:,1));
    caud_loc = round(penumbra(:,2));    
    
    % finds umbra and penumbra height at filtered locations 
    cranpen  = avpost(sub2ind(size(avpost),(1:length(cran_loc))',cran_loc));
    caudpen = avpost(sub2ind(size(avpost),(1:length(caud_loc))',caud_loc));
    umbra = avpost(sub2ind(size(avpost),(1:length(umbra_loc))',umbra_loc));
    
    % all nan values are actually zero, they are there due to earlier part
    % of analysis
    umbra(isnan(umbra)) = 0;
    caudpen(isnan(caudpen)) = 0;
    cranpen(isnan(cranpen)) = 0;
    
    % skip non-injury rats
    if isempty(umbra)
        return;
    end

    % compute umbra-penumra border via  half-height difference
    % alas, for loop was fastest way to write the code
    for tt = 1:length(umbra)

        % flow at each time point for spatial analysis
%         sample_flow = (avpost(tt,:)-avpre)./avpre;
        sample_flow = avpost(tt,:);
        sample_flow(isnan(sample_flow)) = 0;
        
        % set at half-max, can theoretically be any fraction
        frac = 0.5;
        
        % calculate half-max values and find indices
        inj_thresh_cran = abs(umbra(tt) - cranpen(tt)) * frac + umbra(tt);
        inj_thresh_caud = abs(umbra(tt) - caudpen(tt)) * frac + umbra(tt);
        [~, border_cran] = min(abs(sample_flow(cran_loc(tt):umbra_loc(tt))-inj_thresh_cran));
        [~, border_caud] = min(abs(sample_flow(umbra_loc(tt):caud_loc(tt))-inj_thresh_caud));
        
        border_cran = border_cran+cran_loc(tt)-1;
        border_caud = border_caud+umbra_loc(tt)-1;

        % store indices
        transition(tt,:) = ([border_cran, border_caud]);
    end
    
    % same median filtering as before
    if length(time)>10
        transition = round(medfilt1(transition, filt_size, 'includenan', 'truncate'));
    end
    
    % no step this time, not necessary as it has already been done for
    % penumbra
%     step = find(abs(diff(transition(:,1)))>0.5, 1,'last');
%     transition(1:step,1) = transition(step+1,1);
% 
%     step = find(abs(diff(transition(:,2)))>0.5, 1,'last');
%     transition(1:step,2) = transition(step+1,2);

    % find normal blood flow as average of rest of tissue
    [crannorm, caudnorm] = deal(nan(size(cran_loc)));
    
    % another for loop for easy writing
    for tp = 1:length(cran_loc)
        crannorm(tp)  = mean(avpost(tp,1:cran_loc(tp)),'omitnan');
        caudnorm(tp) = mean(avpost(tp,caud_loc(tp):end), 'omitnan');
    end

    % repeat for average pre values
    [caudnormpre, crannormpre] = deal(nan(size(avpre,1),1));
    for tp = 1:size(avpre,1)
        crannormpre(tp)  = mean(avpre(tp,1:cran_loc(1)),'omitnan');
        caudnormpre(tp) = mean(avpre(tp,caud_loc(1):end), 'omitnan');
    end


    % store values in structure
    post_values.Time = time'/60;
    post_values.InjuryHeight = umbra;
    post_values.CranPenHeight = cranpen;
    post_values.CaudPenHeight = caudpen;    

    post_values.CranNormHeight = crannorm;
    post_values.CaudNormHeight = caudnorm;

    post_values.CranPenWidth = (cran_loc-transition(1:length(cran_loc),1))*(dist(2)-dist(1));
    post_values.CaudPenWidth = (caud_loc-transition(1:length(caud_loc),2))*(dist(2)-dist(1));

    post_values.UmbraLoc = umbra_loc;
    post_values.CranLoc = cran_loc;
    post_values.CaudLoc = caud_loc;
    post_values.Transition = transition;
    
    
    % pre values are averaged over time
    pre_values.InjuryHeight = mean(avpre(:,umbra_loc(1)),1);
    pre_values.CranPenHeight = mean(avpre(:,cran_loc(1)),1);
    pre_values.CaudPenHeight = mean(avpre(:,caud_loc(1)),1);
    pre_values.CranNormHeight = mean(crannormpre);
    pre_values.CaudNormHeight = mean(caudnormpre);
    
    
    

end