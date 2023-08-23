% Function: BitmapRender
%
% Purpose: found this function on the internet (https://stackoverflow.com/questions/44286749/bitmap-render-part-of-plot-during-vector-graphics-export-in-matlab)
% it essentially rasterizes all aspects of a figure except what you ask it
% not to.
% (Example in time_evolution)
%
% Input parameters: 
%       Axes: axes of plot
%       KeepObjects: objects on plot (to keep)_
%       RelativePosition: figure it out (I've ignored so far)
%       Draft: figure it out (I've ignored so far)
%       Key: figure it out (I've ignored so far)
%
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function BitmapRender(Axes, KeepObjects, RelativePosition, Draft, Key)
    
    if nargin < 2
        KeepObjects = [];
    end
    if nargin < 3
        RelativePosition = [0 0 1 1];
    end
    if nargin < 4
        Draft = false;
    end
    if nargin < 5
        Key = '';
    end
    
    Figure = Axes.Parent;
    FigureInnerWH = Figure.InnerPosition([3 4 3 4]);
    PixelPosition = round(RelativePosition .* FigureInnerWH);
    
    if isempty(Key)
        OverlayAxes = axes(Figure, 'Units', 'Normalized', 'Position', PixelPosition ./ FigureInnerWH);
        if Draft
            OverlayAxes.Box = 'on';
            OverlayAxes.Color = 'none';
            OverlayAxes.XTick = [];
            OverlayAxes.YTick = [];
            OverlayAxes.HitTest = 'off';
        else
            uistack(OverlayAxes, 'bottom');
            OverlayAxes.Visible = 'off';
        end
        setappdata(Figure, 'BitmapRenderOriginalVisibility', get(Axes.Children, 'Visible'));
    
        Axes.CLimMode = 'manual';
        Axes.XLimMode = 'manual';
        Axes.YLimMode = 'manual';
        Axes.ZLimMode = 'manual';
    
        hManager = uigetmodemanager(Figure);
        [hManager.WindowListenerHandles.Enabled] = deal(false);
        set(Figure, 'KeyPressFcn', @(f, e) BitmapRender(gca, KeepObjects, RelativePosition, Draft, e.Key));
    elseif strcmpi(Key, 'space')
        OverlayAxes = findobj(Figure, 'Tag', 'BitmapRenderOverlayAxes');
        delete(get(OverlayAxes, 'Children'));
        OriginalVisibility = getappdata(Figure, 'BitmapRenderOriginalVisibility');
        [Axes.Children.Visible] = deal(OriginalVisibility{:});
    else
        return;
    end
    
    if Draft
        return;
    end
    
    axpos = get(Axes, 'Position');
    Axes.Visible = 'off';
    set(Axes, 'Position',axpos);
    
    KeepObjectsVisibility = get(KeepObjects, 'Visible');
    [KeepObjects.Visible] = deal('off');
    
    drawnow;
    
    print('-r1200', Figure,'temp','-dpng');
    Frame = imread('temp.png');
    delete('temp.png');
    % Frame = getframe(Figure, PixelPosition);
    
    [Axes.Children.Visible] = deal('off');
    Axes.Visible = 'on';
    Axes.Color = 'none';
    if numel(KeepObjects) == 1
        KeepObjects.Visible = KeepObjectsVisibility;
    else
        [KeepObjects.Visible] = deal(KeepObjectsVisibility{:});
    end
    
    % Image = imagesc(OverlayAxes, Frame.cdata);
    Image = imagesc(OverlayAxes, Frame);
    uistack(Image, 'bottom');
    OverlayAxes.Tag = 'BitmapRenderOverlayAxes';
    OverlayAxes.Visible = 'off';

end