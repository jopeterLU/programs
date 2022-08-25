
function out = Func_LoadImages(fname,fpath)
%% Script info

% Used to find fluorescent vesicles in single microscope images. 
% In: fname, fpath

% Author: Alexandra Andersson
% Date: 2018-01-26

%% File path

fpath = [fpath fname];

%% Load single image, e.g. test image

I_org = double(imread(fpath));

out = I_org;





