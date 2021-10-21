function [p, model] = fit_2d_gaussian_rotated(img, center_guess, width_guess, x, y, w, opt)
% Fits an image with a rotated, 2d, assymetric gaussian peak and returns
% the best parameters.
%
% Order of parameters = [background, amplitude, center-x, center-y,
% fwhm-small, fwhm-large, orientation or large axis with x-axis in rad]
%
% See also: lsqcurvefit
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

assert(nargin >= 3, 'Not enough arguments!')

assert(numel(center_guess) == 2);
assert(numel(width_guess) == 2);

img = double(img);

dims = size(img);
if nargin < 5
    [x, y] = ndgrid(1 : dims(1), 1 : dims(2));
end
if nargin < 6 % default weights all 1
    w = ones(size(img));
end
if nargin < 7
    opt = optimset('Display', 'off');
end

sqrt_w = sqrt(w);
xdata.x = x;
xdata.y = y;
xdata.sqrt_w = sqrt_w;
img_w = img .* sqrt_w;

bg = max(min(img(:)), 0);
br = max(img(:)) - bg;

S = 4;
lb = [0, 0, min(x(:)), min(y(:)), mean(width_guess) * [1,1] / S, -Inf];
ub = [Inf, Inf, max(y(:)), max(y(:)), mean(width_guess) * [1,1] * S, Inf];
p0 = [bg, br, center_guess, width_guess, 0];
best_chisq = Inf;
for i = 0 : pi / 4 : pi
    p0(end) = i;
    [p, chisq] = lsqcurvefit(@fit_func, p0, xdata, img_w, lb, ub, opt);
    if chisq < best_chisq
        best_p = p;
    end
end
p = best_p;

% exchange small and large axis
if p(5) > p(6)
    p(5:6) = p([6,5]);
    p(7) = p(7) + pi / 2;
end
p(7) = mod(p(7), pi);

model = fit_func(p, xdata);

end

% fit function model
function f = fit_func(p, xdata)

bg = p(1);
br = p(2);
center = p(3:4);
widths = p(5:6);
orient = p(7);

x = xdata.x;
y = xdata.y;

% center and rotate grid
xc = x - center(1);
yc = y - center(2);

co = cos(orient);
so = sin(orient);

xr = co * xc + so * yc;
yr = -so * xc + co * yc;

% compute 2d gaussian
f = bg + br * power(2., - (xr.^2 / (widths(1) / 2)^2 + yr.^2 / (widths(2) / 2)^2));

% multiply with weight
f = f .* xdata.sqrt_w;
end