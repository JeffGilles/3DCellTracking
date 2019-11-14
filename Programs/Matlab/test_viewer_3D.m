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

% xyz = bsxfun(@plus, [246 ; 246 ; 17], 0:5:20)';
xyz = [221 276 27 ; ...
    285 110 27 ; ...
    228 24 27 ; ...
    188 256 18 ; ...
    256 179 33 ; ...
    256 42 51];

View.stack3d(Stack, 'points', xyz, 'zfactor', 3);
