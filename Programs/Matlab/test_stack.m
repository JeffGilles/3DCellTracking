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

if ~exist('Stack', 'var')
    
    fprintf('Building stack ...');
    tic
    
    Stack = NaN(h, w, nL);
    
    for i = 1:nL
        
        Stack(:,:,i) = double(imread(fname, i))/255;
        
    end
    
    fprintf(' %.02f sec\n', toc);
    
end

% --- Display -------------------------------------------------------------

figure(1)
clf
hold on

R = squeeze(Stack(10,:,:))';

imshow(R)

% caxis auto
colorbar

axis on image

daspect([3 1 1])