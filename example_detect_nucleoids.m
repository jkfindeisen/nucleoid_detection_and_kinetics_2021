function example_detect_nucleoids()
% This script shows exemplarily how DNA and EdU spots are detected in the
% STED images. Some exemplary measurements are loaded. The user is given
% the possibility to change an initial threshold and interactively sees
% how many spots would be selected. After manual threshold selection is
% done, Gaussian peak functions are fit to each identified spot position to
% refine the positions.
%
% Requires that the nucleus is determined before.
%
% Part of "The TFAM to mtDNA ratio defines inner-cellular nucleoid
% populations with distinct activity levels"
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

% load some example data
dna = imread('data/HDFa_EdU-incubation-18h_ROI9_DNA_AlexaFluor594_20nm-px.tiff');
edu = imread('data/HDFa_EdU-incubation-18h_ROI9_EdU_StarRed_20nm-px.tiff');
nucleus = imread('data/HDFa_EdU-incubation-18h_ROI9_nucleus-mask_20nm-pixelsize.tiff');
nucleus = nucleus > 0;
px = 20e-9; % 20nm pixel size

%% threshold DNA spot detection and fit them
fprintf('threshold DNA spots\n');

threshold = 0.1;

% thresholding
[im, d, object] = prepare_image(dna, nucleus, px);
[output, quit] = interactive_threshold(im, d, object, px, threshold, 'DNA example image');

% fitting spots
output.spots = fit_spots(output.positions(:, 1), output.positions(:, 2), dna, px, 100e-9);

% save output
save('data/dna_spots.mat', 'output');

%% threshold EdU spot detection and fit them
fprintf('threshold EdU spots\n');

threshold = 0.2;

% thresholding
[im, d, object] = prepare_image(edu, nucleus, px);
[output, quit] = interactive_threshold(im, d, object, px, threshold, 'EdU example image');

% fitting spots
output.spots = fit_spots(output.positions(:, 1), output.positions(:, 2), edu, px, 100e-9);

% output would now be saved
save('data/edu_spots.mat', 'output');


end

function spots = fit_spots(xi, yi, im, px, fwhm)
fprintf(' fit spots\n');
fwhmp = fwhm / px;
spots = fit_nucleoids(xi, yi, im, fwhmp);
spots(:, [5, 10, 11]) = spots(:, [5, 10, 11]) * px;
spots = [xi, yi, spots]; % add xi, yi to the beginning
end

function [im, data, object] = prepare_image(im, mask, px)

im = double(im);

% dilate mask a bit more
mask = imdilate(mask, ones(ceil(0.5e-6/px)));

% add border (include in mask)
B = ceil(0.5e-6/px);
[x,y] = ndgrid(1:size(mask,1),1:size(mask,2));
mask = mask | x <= B | x > size(mask,1)-B | y <= B | y > size(mask,2)-B;

% invert mask
object = ~mask;

% smooth image to be able to find peaks
% im_sm = img_smooth_mask(im, object, 0.08e-6 / px);
im_sm = img_smooth_mask(im, object, 0.12e-6 / px);

% smooth image to estimate background
im_bg = img_smooth_mask(im, object, 0.6e-6 / px);

% subtract both and set negative and in mask to zero
data = im_sm - 0.75 * im_bg;
data(~object) = 0;
data(data<0) = 0;

% suppress intensity in ~object of image to 1/6th
h = [max(im(~object)), max(im(object))];
im(~object) = im(~object) / h(1) * h(2) / 6;

% shrink object further by 2 px
object = imerode(object, ones(5));

end

function [output, quit] = interactive_threshold(im, data, object, px, t, fig_name)
% image, data, object, pixel size, initial threshold

% boundaries of object
object_boundaries = bwboundaries(object);

% min, max of image at least 50
R = [0, max(max(data(:)), 50)];

% step size
dt = 0.01;

% caxis_limits (start with min, max)
C = [0, 0.7 * max(im(:))];
dc = C(2) / 10;

% show figure (maximized)
fig = figure('units','normalized','outerposition',[0 0 1 1], 'Name', fig_name);
fig2 = ui_figure_second_screen();

% default parameters
show = true;
quit = false;
xi = [];
yi = [];

    function update()
        
        figure(fig);
        hold off;
        
        % display image
        h = imagesc(im);
        axis image;
        colormap(hot);
        caxis(C);
        colorbar();
        hold on;
        
        % this is just to maximize the size of the image without loosing
        % the scaling
        pos = h.Parent.Position;
        a = 0.02;
        dpos = min([pos(1)-a, 1-a-pos(1)-pos(3), pos(2) - a, 1-a-pos(2)-pos(4)]) .* pos(3:4) / max(pos(3:4));
        pos = pos + [-dpos(1), -dpos(2), 2*dpos(1), 2*dpos(2)];
        h.Parent.Position = pos;
        
        % display object boundaries
        for i = 1 : length(object_boundaries)
            b = object_boundaries{i};
            plot(b(:, 2), b(:, 1), 'm');
        end
        
        % apply threshold
        T = t*R(2)+(1-t)*R(1);
        
        % specific things should be shown
        % local max
        [idx, vi, xi, yi] = omex_local_max(data, 'max', 8, T);
        m = object(idx);
        xi = xi(m);
        yi = yi(m);
        vi = vi(m);
        
        if ~isempty(fig2)
            % show histogram in figure2
            figure(fig2);
            histogram(vi, R(1):diff(R)/20:R(2));
            figure(fig);
        end
        
        if show
            % show circles around positions
            plot(yi, xi, 'o', 'MarkerSize', round(150e-9/px), 'Color', [0, 0.4, 0]);
        end
        
        % update title
        title(sprintf('t=%.2f, n=%d', t, length(xi)));
        
        
    end

    function adjust(~, event)
        
        switch event.Key
            case 'escape' % finished with this figure
                uiresume();
            case 'q' % quit application
                fprintf('Will quit. Nothing will be saved.\n');
                quit = true;
                uiresume();
            case 'uparrow' % less positions, higher threshold
                if ~isempty(event.Modifier) && strcmp(event.Modifier{1}, 'shift')
                    t = min(1, t + 10 * dt);
                else
                    t = min(1, t + dt);
                end
                update();
            case 'downarrow' % more positions, lower threshold
                if ~isempty(event.Modifier) && strcmp(event.Modifier{1}, 'shift')
                    t = max(dt, t - 10 * dt);
                else
                    t = max(dt, t - dt);
                end
                update();
            case 'v' % less saturation
                C(2) = C(2) + dc;
                update();
            case 'b' % more saturation,
                C(2) = max(C(1) + dc, C(2) - dc);
                update();
            case 'c' % toggle show
                show = ~show;
                update();
        end
    end

update();
fig.KeyPressFcn = @adjust;
uiwait(fig);
close(fig);
if ~isempty(fig2)
    close(fig2);
end

output = struct('threshold', t, 'R', R, 'positions', [xi, yi]);

end
