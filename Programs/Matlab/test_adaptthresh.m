clc

% === Parameters ==========================================================

fname = '/home/ljp/Science/Projects/ImageAnalysis/3DCellTracking/Data/emb1-C2-1-50.tif';

% =========================================================================

Img = imread(fname, 1);

% Img = imgaussfilt(Img, [1 1]*1);
% Img = medfilt2(Img, [1 1]*3);


T = adaptthresh(Img, 0.3, ...
    'NeighborhoodSize', [1 1]*7, ...
    'ForegroundPolarity', 'bright', ...
    'Statistic', 'gaussian');

BW = imbinarize(Img, T);

BW = bwpropfilt(BW, 'Area', [50 Inf]);
BW = imclose(BW, strel('disk', 1));

% --- Display -------------------------------------------------------------

figure(1)
clf
hold on

% imshow(BW)
imshowpair(Img, BW, 'montage')

% caxis auto
colorbar