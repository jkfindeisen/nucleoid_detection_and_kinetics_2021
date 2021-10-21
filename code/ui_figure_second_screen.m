function f = ui_figure_second_screen(varargin)
% Creates new figure and shows on second screen by default, if there is
% one.
%
% See: https://www.mathworks.com/matlabcentral/answers/16663-is-it-possible-to-viewing-the-figure-window-on-second-display
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

monitor_positions = get(0, 'MonitorPositions');
if size(monitor_positions, 1) > 1
    % At least two monitors, put on second
    
    % order can be reversed
    [~, idx] = sort(monitor_positions(:, 1));
    monitor_positions = monitor_positions(idx, :);
    
    f = figure(varargin{:});
    
    % shift
    units = f.Units;
    f.Units = 'pixels';
    f.OuterPosition = [monitor_positions(2, 1), monitor_positions(2, 2) + monitor_positions(2, 4) - f.OuterPosition(4), f.OuterPosition(3:4)];
    f.Units = units;
else
    % no second monitor return empty
    f = [];
end

end