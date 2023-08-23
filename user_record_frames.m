% Function: user_record_frames
%
% Purpose: simple GUI to mark where the probe view switches
%
% Input parameters: 
%       
%
% Output parameters:
%       
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function critFrames = user_record_frames(vid)
    
    fig = uifigure("Position",[1500 100 650 475]);
    g = uigridlayout(fig);
    g.RowHeight = {'1x','fit'};
    g.ColumnWidth = {'1x','10x','1x'};
    
    ax = uiaxes(g);
    ax.Layout.Row = 1;
    ax.Layout.Column = [1 length(g.ColumnWidth)];
    
    % cg = uigauge(g);
    % cg.Layout.Row = 1;
    % cg.Layout.Column = [1 3];
    
    if length(size(vid))==3
        imshow(vid(:,:,1),'Parent',ax);
    elseif length(size(vid))==4
        imshow(vid(:,:,:,1),'Parent',ax);
    else
        error('Input is not a video');
    end


    critFrames = [];
    
    sld = uislider(g, ...
        "ValueChangingFcn",@(src,event)updateFrame(src,event,vid,ax));
    sld.Layout.Row = 2;
    sld.Layout.Column = 2;
    sld.Limits = [1, size(vid, length(size(vid)))];

    record_btn = uibutton(g,"Text","Record", ...
        "ButtonPushedFcn", @(src,event) recordButtonPushed(src,event,round(sld.Value)));
    record_btn.Layout.Row = 2;
    record_btn.Layout.Column = 1;


    close_btn = uibutton(g,"Text","Done", ...
        "ButtonPushedFcn", @(src,event) closeButtonPushed(src,event));
    close_btn.Layout.Row = 2;
    close_btn.Layout.Column = 3;
    
    waitfor(fig);

    function updateFrame(src,event,vid,ax)
        if length(size(vid))==3
            imshow(vid(:,:,round(event.Value)),'Parent',ax);
        else
            imshow(vid(:,:,:,round(event.Value)),'Parent',ax);
        end
        
    end
    
    function recordButtonPushed(src,event,value)
        critFrames = [critFrames; value];
    end

    function closeButtonPushed(src,event)
        closereq();
    end
end
