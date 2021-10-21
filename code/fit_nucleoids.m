function r = fit_nucleoids(xi, yi, data, fwhmp)
% does a symmetric and an asymmetric fit
% symmetric fit output order: [background, amplitude, center-x, center-y, fwhm]
% assymetric fit output order: [background, amplitude, center-x, center-y, fwhm-small, fwhm-large, orientation]
% both outputs one after another
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

assert(nargin == 4);

% prepare fit
B = ceil(1.2*max(fwhmp));
g = -B:B;
[xj, yj] = ndgrid(g, g);
w = ones(size(xj));

% do not process those far away
dims = size(data);
ignore = xi - B < 1 | xi + B > dims(1) | yi - B < 1 | yi + B > dims(2);

% loop over local max and fit
n = length(xi);
r = zeros(n, 12);
% opt = optimset('Display', 'iter', 'FunctionTolerance', 1e-6);
opt = optimoptions('lsqcurvefit', 'Display', 'off', 'TolFun', 1e-6);
for i = 1 : n
    if ignore(i)
        continue;
    end
    xk = xi(i);
    yk = yi(i);
    cut = data(g+xk,g+yk); % take raw data here
    [fp, ~] = fit_2d_gaussian_symmetric(cut, [0, 0], fwhmp, xj, yj, opt);
    fp(3) = fp(3) + xk;
    fp(4) = fp(4) + yk;
    [fp2, ~] = fit_2d_gaussian_rotated(cut, [0, 0], fwhmp * [0.9, 1.1], xj, yj, w, opt);
    fp2(3) = fp2(3) + xk;
    fp2(4) = fp2(4) + yk;
    r(i, :) = [fp, fp2];
end

end