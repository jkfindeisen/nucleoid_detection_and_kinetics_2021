function sm = img_smooth_mask(im, msk, fwhmp, mode)
% Smoothes 2D image on the pixels marked by a mask with a 2D Gaussian with
% a FWHM given in pixels, mode is {'conv2', 'fft'}
%
% See also: conv2
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

assert(nargin >= 3);
if nargin < 4
    mode = 'fft';
end

switch mode
    case 'conv2'
        % smoothing kernel
        L = ceil(3 * fwhmp);
        [x, y] = ndgrid(-L:L,-L:L);
        k = exp(-4*log(2)*(x.^2+y.^2)/fwhmp^2);
        k = k / sum(k(:));
        
        % convolve two times
        im_sm = conv2(im .* msk, k, 'same');
        msk_sm = conv2(msk, k, 'same');
    case 'fft'
        % compute grid with (0, 0) at origin
        [x, y] = img_fourier_grid(size(im));
        % compute gaussian kernel
        k = exp(-4*log(2)*(x.^2+y.^2)/fwhmp^2);
        k = k / sum(k(:));
        % convolutions
        otf = fft2(k);
        im_sm = real(ifft2(fft2(im .* msk) .* otf));
        msk_sm = real(ifft2(fft2(msk) .* otf));
    otherwise
        error('Unknown mode');
end

% assemble output
T = 0.1;
m = msk_sm > T;
sm = zeros(size(im));
sm(m) = im_sm(m) ./ msk_sm(m);

end