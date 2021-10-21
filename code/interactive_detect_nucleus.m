function [nucleus, p, quit] = interactive_detect_nucleus(imgs, p, name)
% Interactively detects the nucleus in a cell image with some parameters
%
% Keyboard instructions:
%   escape  close the figure, continue with processing
%   q       close the figure, indicates that an abort is desired
%   return  switch between different views of the raw (if available), cycles
%   arrow down/up       decreases/increases threshold (larger, smaller area)
%   arrow left/right    more smoothing/less smoothing (more coherent area)
%   a/s     larger or smaller minimal area size (in pixel)
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

assert(nargin == 3);

quit = false;
F = 1.1;
A = 1.5;
counter = 1;
fig = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'Name', name);

    function update()
        img = imgs{counter};
        nucleus = img_detect_nucleus(img, p);
        B = bwboundaries(nucleus);
        hold off;
        im = imagesc(img);
        im.Parent.YDir = 'normal';
        colormap(hot);
        axis image;
        hold on;
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'w')
        end
        title(sprintf('smooth %.0f, threshold %.2f, area %.0f', p.smooth, p.threshold, p.area));
    end

    function key(~, event)
        switch event.Key
            case 'escape'
                uiresume();
            case 'q'
                quit = true;
                uiresume();
            case 'return'
                % switch image
                counter = mod(counter, numel(imgs)) + 1;
                update();
            case 'downarrow'
                p.threshold = p.threshold / F;
                update();
            case 'uparrow'
                p.threshold = p.threshold * F;
                update();
            case 'leftarrow'
                p.smooth = p.smooth / F;
                update();
            case 'rightarrow'
                p.smooth = p.smooth * F;
                update;
            case 'a'
                p.area = p.area * A;
                update();
            case 's'
                p.area = p.area / A;
                update();
        end
    end

fig.KeyPressFcn = @key;
update();
uiwait(fig);
close(fig);

end