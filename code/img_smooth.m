function [smooth, kernel_max, kernel] = img_smooth(img, s, type)
% Just blurs an image with a gaussian of size s in pixels. If length of s
% is two than a variable x and y size is used. Works only on 2D images.
%
% Syntax
%   function [smooth, kernel_max] = img_smooth(img, s, type)
%
% Comment: 2 important properties: No shift! Total sum stays constant!
%
% type == 0 -> gaussian with FWHM s (=default)
% type == 1 -> use disk of radius s
% type == 2 -> use square of length 2*s+1
% type == 3 -> use ring
% type == 4 -> use rotated gaussian with FWHM s(1), s(2) and rotation angle
% s(3) (rad)
%
% See also: convolution, img_fourier_grid
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

if nargin < 2
    error('Must provide arguments!');
end

if nargin >= 3
    mode = type;
else
    mode = 0;
end

if length(s) == 1
    s = [s, s];
end

% not working with int32
img = double(img);

% compute grid with (0, 0) at origin
[x, y] = img_fourier_grid(size(img));

% compute smoother
switch mode
    case 0
        kernel = power(2., - (x.^2 / (s(1) / 2)^2 + y.^2 / (s(2) / 2)^2));
        kernel = kernel / sum(kernel(:)); % normalization
    case 1
        kernel = (x.^2+y.^2) <= s(1)^2;
    case 2
        kernel = abs(x) <= s(1) & abs(y) <= s(2);
    case 3
        r2 = x.^2 / (s(1) / 2)^2 + y.^2 / (s(2) / 2)^2;
        kernel = r2 .* power(2., -r2);
        kernel = kernel / sum(kernel(:));
    case 4
        cp = cos(s(3));
        sp = sin(s(3));
        xr = cp * x + sp * y;
        yr = -sp * x + cp * y;
        kernel = power(2., - (xr.^2 / (s(1) / 2)^2 + yr.^2 / (s(2) / 2)^2));
        kernel = kernel / sum(kernel(:)); % normalization        
    otherwise
        error('Unknown type!');
end

% convolve
smooth = real(ifft2(fft2(img) .* fft2(kernel)));
kernel_max = max(kernel(:));

end