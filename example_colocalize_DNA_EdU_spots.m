function example_colocalize_DNA_EdU_spots()
% This script examplarily shows the combination of the fitted DNA and EdU
% spot detection to determine how many nucleoids were EdU-positive, i.e.
% how many of the detected DNA spots also had an EdU spot close by. The
% decision is done by matching DNA and EdU spot pairs and then choosing all
% pairs that have a mutual distance below a certain threshold.
%
% Requires that the nucleoids are determined before.
%
% Part of "The TFAM to mtDNA ratio defines inner-cellular nucleoid
% populations with distinct activity levels"
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

nucleus = imread('data/HDFa_EdU-incubation-18h_ROI9_nucleus-mask_20nm-pixelsize.tiff');
nucleus = nucleus > 0;

AB = 15e-9; % antibody size
px = 20e-9; % 20nm pixel size
T = 0.1e-6; % minimal distance between next
TR = 1;

dna = load('data/dna_spots.mat');
dna = dna.output.spots;

edu = load('data/edu_spots.mat');
edu = edu.output.spots;

% get those out which are on nucleus (because we sometimes modified the nucleus mask afterwards)
ix = in_nucleus(edu(:, 5:6), nucleus);
edu(ix, :) = [];
ix = in_nucleus(dna(:, 5:6), nucleus);
dna(ix, :) = [];

%% reduce by nearest neighbors too close
ix = reduce_events(edu(:, [5,6]), edu(:, 7), T / px);
edu = edu(ix, :);

ix = reduce_events(dna(:, [5,6]), dna(:, 7), T / px);
dna = dna(ix, :);

n_edu = size(edu, 1);
n_dna = size(dna, 1);

%% sort those edu points out where there is no dna point
d = sqrt((edu(:, 5) - dna(:, 5).').^2+(edu(:, 6) - dna(:, 6).').^2)*px; % distance edu - dna
R1 = repmat(edu(:, 7), [1, n_dna]); % size of edu
R2 = repmat(dna(:, 7).', [n_edu, 1]); % size of dna

% subtract antibody-shell size
R1 = max(0, R1 - AB);
R2 = max(0, R2 - AB);

% relative distance
dr = d ./ (R1 + R2);

% edu positive
edu_positive = any(dr < TR, 1); % very important here, otherwise Matlab changes the direction for 0 or 1 entries in first dimension!
n_edu_positive = sum(edu_positive);

fprintf('%d DNA spots, %d EdU spots, %d DNA spots with a single EdU spot within the valid distance\n', n_dna, n_edu, n_edu_positive);

end

function ix = in_nucleus(pos, nucleus)
% checks if positions are in a nucleus mask

pos = round(pos);
ix = sub2ind(size(nucleus), pos(:, 1), pos(:, 2));
ix = nucleus(ix);

end

function idx = reduce_events(spots, fwhm, t)
% reduce those which are too close to each other

while true
    
    nn = omex_nearest_neighbour(spots);
    
    if ~any(nn(:, 1) < t)
        break;
    end
    
    [~, i1] = min(nn(:, 1));
    i2 = nn(i1, 2);
    if fwhm(i1) < fwhm(i2)
        i = i2;
    else
        i = i1;
    end
    
    spots(i, :) = [Inf, Inf];
end

idx = nn(:, 1) ~= Inf;

end