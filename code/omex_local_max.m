function [idx, vi, xi, yi] = omex_local_max(stack, type, connection, threshold)
% Finds local minima, maxima in a 2D data stack above or under a threshold. And
% returns the positions and values sorted by descending value.
%
% Syntax
%   function [idx, vi, xi, yi] = omex_local_max(stack, type, connection,
%   threshold)
%
% Input parameters
%   stack   2D data stack
%   type    'min'/'max'(=default) for either minima or maxima
%   connection  4/8(=default) for either only direct neighbours or
%               diagonals also
%   threshold   upper/lower threshold for minima/maxima (default: Inf/-Inf)
%
% Output parameters
%   idx     positions of minima/maxima
%   vi      values at these positions
%   xi      x position (convenience)
%   yi      y position (convenience)
%
% See also: circshift
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

%% error checks
assert(nargin >= 1, 'Not enough arguments!');

dims = size(stack);
assert(length(dims) == 2, 'Parameter stack must be 2D matrix!');

if nargin < 2
    type = 'max';
end

if nargin < 3
    connection = 8;
end

if nargin < 4
    threshold = -Inf;
end

switch type
    case 'max'
    case 'min'
        % if it were minima, mirror
        stack = -stack;
    otherwise
        error('Unknown type!');
end

%% search for maxima only from here
msk = stack > threshold ...
    & stack > circshift(stack, [1, 0]) & stack > circshift(stack, [-1, 0]) ...
    & stack > circshift(stack, [0, 1]) & stack > circshift(stack, [0, -1]);
% add diagonal connection if wished
if connection == 8
    msk = msk & stack > circshift(stack, [1, 1]) & stack > circshift(stack, [-1, 1]) ...
        & stack > circshift(stack, [1, -1]) & stack > circshift(stack, [-1, -1]);
end
idx = find(msk);
[xi, yi] = ind2sub(dims, idx);
vi = stack(idx);

% sort
[vi, h]  = sort(vi, 'descend');
idx = idx(h);
xi = xi(h);
yi = yi(h);

% if it were minima mirror again
if strcmp(type, 'min')
    vi = -vi;
end

end