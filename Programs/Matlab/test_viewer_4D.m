clc
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure');

% === Parameters ==========================================================

fname = '/home/ljp/Science/Projects/ImageAnalysis/3DCellTracking/Data/emb1-C2-1-50.tif';

nL = 54;

% -------------------------------------------------------------------------

Info = imfinfo(fname);
nf = numel(Info);
w = Info(1).Width;
h = Info(1).Height;

% =========================================================================





View.stack4d

