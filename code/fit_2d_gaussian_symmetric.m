function [params, model, chisq] = fit_2d_gaussian_symmetric(img, center_guess, width_guess, x, y, opt)
% Fits an image with a 2d symetric gaussian peak and returns the best parameters.
%
% Order of parameters = [background, amplitude, center-x, center-y, fwhm]
%
% See also: lsqcurvefit
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

assert(nargin >= 1, 'Not enough arguments!')
img = double(img);

dims = size(img);
if nargin < 5
    [x, y] = ndgrid(1 : dims(1), 1 : dims(2));
end

if nargin < 6
    opt = optimset('Display', 'off');
end

xdata.x = x;
xdata.y = y;

if nargin < 3
    width_guess = mean(dims) / 3;
end
assert(numel(width_guess) == 1);

if nargin < 2
    center_guess = (dims + 1) / 2;
end

% guess for background and brightness
bg = max(min(img(:)), 0);
br = max(img(:)) - bg;

S = 4;
p0 = [bg, br, center_guess, width_guess];
lb = [0, 0, min(x(:)), min(y(:)), width_guess / S];
ub = [Inf, Inf, max(y(:)), max(y(:)), width_guess * S];
params = lsqcurvefit(@fit_func, p0, xdata, img, lb, ub, opt);

model = fit_func(params, xdata);
chisq = sum((model(:)-img(:)).^2 ./ max(model(:), 1e-6));

end

% fit function model
function f = fit_func(p, xdata)

bg = p(1);
br = p(2);
center = p(3:4);
width = p(5);

x = xdata.x;
y = xdata.y;

% center and rotate grid
xc = x - center(1);
yc = y - center(2);

% compute 2d gaussian
f = bg + br * power(2., - (xc.^2 + yc.^2) / (width / 2)^2);
end