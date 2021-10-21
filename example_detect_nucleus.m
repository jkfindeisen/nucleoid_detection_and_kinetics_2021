function example_detect_nucleus()
% This script shows exemplarily how the cell nucleus was detected in 2-color
% confocal and STED nucleoid measurements. Both STED channels are important
% and given to an interactive nucleus detection, i.e. one can set
% thresholds but typically the initial values are already quite good. The
% resulting mask is stored as a binary image.
%
% Part of "The TFAM to mtDNA ratio defines inner-cellular nucleoid
% populations with distinct activity levels"
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

initialize();

% load some example data
dna = imread('data/HDFa_EdU-incubation-18h_ROI9_DNA_AlexaFluor594_20nm-px.tiff');
edu = imread('data/HDFa_EdU-incubation-18h_ROI9_EdU_StarRed_20nm-px.tiff');

% nucleus detection parameters
params = struct('pad', 20, 'smooth', 30, 'threshold', 0.3, 'area', 30000);

% interactive nucleus detection (returns updated parameters)
[nucleus, params] = interactive_detect_nucleus({dna+edu, dna, edu}, params, 'Example data');

% write nucleus
imwrite(uint8(255*nucleus), 'data/HDFa_EdU-incubation-18h_ROI9_nucleus-mask_20nm-pixelsize.tiff');

end