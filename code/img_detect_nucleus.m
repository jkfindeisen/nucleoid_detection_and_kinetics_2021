function [m, a] = img_detect_nucleus(d, p)
% Using some heuristics, detect the nucleus in an image or any large bright
% segments with a certain number of pixels inside.
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

assert(nargin == 2);

% smooth
P = p.pad;
dp = padarray(d, [P, P], 'symmetric');
ds = img_smooth(dp, p.smooth);
ds = ds(1+P:end-P,1+P:end-P);

% binarize
T = p.threshold;
m = ds > max(ds(:)) * T + min(ds(:)) * (1 - T);

% fill holes
m = bwfill(m, 'holes');

% only keep the largest with more than X pixels
S = p.area;
m = bwareafilt(m, [S, Inf]);

% convex hull
m = bwconvhull(m, 'objects');

% dilate slightly (to really get everything)
m = imdilate(m, ones(5));

% get area
a = bwarea(m);

end