clc

% === Parameters ==========================================================

fname = '/home/ljp/Science/Projects/ImageAnalysis/3DCellTracking/Data/emb1-C2-1-50.tif';

% =========================================================================

% --- Load

Raw = double(imread(fname, 40));

% --- Pre-process --------------------------------------------------------- 

Img = Raw;

% --- Noise filtering

Img = imnlmfilt(Img, 'DegreeOfSmoothing', 20);

% --- Morphological recontruction

Tmp = imerode(Img, strel('disk', 5));
Res = imreconstruct(Tmp, Img);

% T = adaptthresh(Img, 0.3, ...
%     'NeighborhoodSize', [1 1]*7, ...
%     'ForegroundPolarity', 'bright', ...
%     'Statistic', 'gaussian');
% 
% BW = imbinarize(Img, T);
% 
% BW = bwpropfilt(BW, 'Area', [50 Inf]);
% BW = imclose(BW, strel('disk', 1));

% --- Display -------------------------------------------------------------

figure(1)
clf
hold on

% imshow(BW)
imshowpair(Raw, Img, 'montage')
% imshowpair(Img, BW, 'montage')

axis tight

% caxis auto