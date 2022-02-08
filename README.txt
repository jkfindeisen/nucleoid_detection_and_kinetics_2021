Software demonstration set for

"The TFAM to mtDNA ratio defines inner-cellular nucleoid populations with distinct activity levels"
https://doi.org/10.1016/j.celrep.2021.110000

by Christian Brüser, Jan Keller-Findeisen, and Stefan Jakobs

Department of NanoBiophotonics, Research Group Mitochondrial Structure and Dynamics, Max Planck Institute for Biophysical Chemistry; Clinic of Neurology, High Resolution Microscopy of the Cell, University Medical Center Göttingen; Department of NanoBiophotonics, Max Planck Institute for Biophysical Chemistry; Fraunhofer Institute for Translational Medicine and Pharmacology ITMP, Translational Neuroinflammation and Automated Microscopy
Göttingen, Germany

Software description
====================

The software is written in Matlab and tested on Matlab 2020.

example_colocalize_DNA_EdU_spots.m 	- combines estimated DNA and EdU spots to determine the fraction of EdU
					  active spots
example_detect_nucleoids.m 		- identifies DNA and EdU spots in STED images
example_detect_nucleus.m 		- identifies the nucleus of STED images
fit_timeseries_edu_positive_nucleoids.m - fits the determined fraction of EdU containing spots with suitable 					  model functions

folder code 		- additional software
folder data 		- example data

Usage instructions and example data set
=======================================

Please run the example_xxx.m scripts to see exemplary steps of the analysis procedure focused on a dataset
of a HDFa cell incubated with EdU for 18 hours and measured for mitochonrdrial DNA as well as EdU in STED mode.

Please run script fit_timeseries_edu_positive_nucleoids.m to reproduce values given in Supplementary Table 1 and Fig. 1f in the manuscript.

Run initialize.m once to add all relevant folders to the Matlab path.

License
=======

The code is MIT licensed (see LICENSE.txt)
