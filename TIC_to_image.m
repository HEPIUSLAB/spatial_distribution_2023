% Function: TIC_to_image
%
% Purpose: 
%
% Input parameters:
%   RawFlow: double array (m x n x n_frames)
%   time: double array (n_frames x 1)
%
% Output parameters:
%   newflow: double array (m x n)
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function newflow = TIC_to_image(t, RawFlow)
    
    % average filter
    filt_size = 32;
    avfilt = reshape(ones(filt_size,1)/filt_size,1,1,[]);    
    FlowSmooth = convn(RawFlow, avfilt, 'valid');
    
    newflow = zeros(size(FlowSmooth,[1,2]));
    
    % for each "pixel" in TIC image
    JL = size(FlowSmooth,2);
    lineLength = 0;
    for ii = 1:size(FlowSmooth,1)

        % progress display
        if mod(ii,10)==0
            lineLength = 0; % comment if not parfor
            fprintf(repmat('\b',1,lineLength));
            lineLength = fprintf('%d / %d \n', ii, size(FlowSmooth,1));
        end

        for jj = 1:JL

            % fit first-order exponoential model from max value of TIC
            flow_smooth = squeeze(FlowSmooth(ii,jj,:));
            [M,I] = max(flow_smooth);
            
            % filtering 
            if M>0.015 && I<floor(length(flow_smooth)/2)
                f2 = fit(t(I:end-filt_size+1)', flow_smooth(I:end),'exp1');
                coeffs = coeffvalues(f2);
                
                % save decay constant as pixel flow value
                if coeffs(1)>0.015 && coeffs(2) < 0
                    newflow(ii,jj) = -coeffs(2);
                end
            end
        end
    end
    
    % display result
    if size(RawFlow, [1,2])==[1,1]
        raw_flow = squeeze(RawFlow);
        plot(t, squeeze(RawFlow));
        hold on
        plot(t(1:(end-length(avfilt)+1)), flow_smooth);
        plot(f2);
        hold off

        save('CEUS_fitting.mat', 't', 'raw_flow', 'flow_smooth', 'f2');
    end

    
end