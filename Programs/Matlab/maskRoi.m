clc
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure')

% === Parameters ==========================================================

tSlices = 2700;
nL = 54;
nTime = 5;
fDir = '/Users/jean-francoisgilles/Documents/Data/2020/Pauline/datatest/';

runC1 = 'MAX_emb1-C1-1-50';
runC2 = 'MAX_emb1-C2-1-50';
run = 'emb1-C2-1-50';

fname1 = ['/Users/jean-francoisgilles/Documents/Data/2020/Pauline/datatest/' runC1 '.tif'];
fname2 = ['/Users/jean-francoisgilles/Documents/Data/2020/Pauline/datatest/' runC2 '.tif'];
PIVDir = ['/Users/jean-francoisgilles/Documents/Data/2020/Pauline/datatest/' run '/'];

gridSize = 25;
extSize = 45;
factor = 3;

min_Int_xy = 0.05;
min_Int_z = 0.08;
min_corr = 0.65;

% Display
df = 5;

force = true;

% -------------------------------------------------------------------------

Info = imfinfo(fname1);
nf = numel(Info);
w = Info(1).Width;
h = Info(1).Height;


% =========================================================================
% --- Grid ----------------------------------------------------------------

% --- 1D coordinates

xg = 1:gridSize:w;
yg = 1:gridSize:h;
zg = round(1:gridSize/factor:nL);

% --- 2D grid

%[X_,Y_] = meshgrid(xg, yg);
[X_,Y_] = ndgrid(xg, yg);
%grid = [X_(:), Y_(:)];

% --- 3D grid

% [X_, Y_, Z_] = meshgrid(xg, yg, zg);
% Xg = X_(:);
% Yg = Y_(:);
% Zg = Z_(:);

%--init
Raw1 =  NaN(h, w, nTime); %Channel 1
Raw2 =  NaN(h, w, nTime); %Channel 2

for ti=1:nTime
    fprintf('=== %i =======================================\n\n', ti);
    
    % --- Raw stack -----------------------------------------------------------
    
    if ~exist('Raw1', 'var') || force
        
        fprintf('Loading raw images ...');
        tic
        %for i = 1:nL
            Raw1(:,:,ti) = double(imread(fname1, ti));
%             Raw1(:,:,2) = double(imread(fname1, 2))/255;
        %end
        
        fprintf(' %.02f sec\n', toc);
        
    end

end

s = length(X_(:));

% X26 = zeros(s(1), nTime);
% Y26 = zeros(s(1), nTime);
piv = zeros(s, 2, nTime);

for t = 1:nTime
    T1 = load([PIVDir num2str(t, '%03i') '.mat']);
    indZ = T1.Zg==26; %keep size
%     indZ = find(T1.Zg==26);
    piv(:,1,t) = T1.Dx(indZ);
    piv(:,2,t) = T1.Dy(indZ);
end

% overlap points with data

View.stack2DTime(Raw1, grid, piv)

% drawpolygon with green channel contours

roi = drawpolygon('Color','r')

roiEvo = zeros(length(roi.Position), 2, nTime);
% interpolate (interp2)

%resh = reshape(piv(:,1,1),size(X_));

for t = 1:nTime
%     lx = interp1(unique(X_), unique(Y_), unique(piv(:, 1, t)));
    gridX = reshape(piv(:, 1, t),size(X_));
    gridY = reshape(piv(:, 2, t),size(Y_));
    
    gridX(isnan(gridX)) = 0;
    gridY(isnan(gridY)) = 0;
    
    gy = griddedInterpolant(X_, Y_, gridX);
    gx = griddedInterpolant(X_, Y_, gridY);
    
    
% rx = griddedInterpolant(reshape(piv(:,1,1), size(X_)), roi.Position(:,1))

% gdx = zeros(size(roi.Position(:,1)));
% gdy = zeros(size(roi.Position(:,2)));
% gx = zeros(size(roi.Position(:,1)));
% for pos = 1:size(roi.Position(:,1))
%     roi.Position(pos,:)
%     indx = knnsearch(X_(:),roi.Position(pos,1));
%     indy = knnsearch(Y_(:),roi.Position(pos,2));
%     gx(pos, 1) = X_(indx);
%     [X_(indx), Y_(indy)]
%     gdx(pos, 1) = piv(indx,1,1);
%     gdy(pos, 1) = piv(indy,2,1);
%     [gdx(pos, 1), gdy(pos, 1)]
% end
    
    
    for p = 1:length(roi.Position)
        if t == 1
            roiEvo(p, 1, 1) = roi.Position(p, 1);
            roiEvo(p, 2, 1) = roi.Position(p, 2);
        else
            roiEvo(p, 1, t) = roiEvo(p, 1, t-1) + 10*gx(roiEvo(p, 1, t-1), roiEvo(p, 2, t-1));
            roiEvo(p, 2, t) = roiEvo(p, 2, t-1) + 10*gy(roiEvo(p, 1, t-1), roiEvo(p, 2, t-1));
%             roiEvo(p, 1, t) = roiEvo(p, 1, t-1) + lx(roiEvo(p, 1, t-1), roiEvo(p, 2, t-1))
%             roiEvo(p, 2, t) = roiEvo(p, 2, t-1);
        end
    end
    
    
end


View.stack2DTimeroi(Raw1, grid, piv, roiEvo)
%roi2 = interp2(size(X_),size(Y_), reshape(piv(:,1,1),size(X_)), roi.Position(:,1), roi.Position(:,2));

% make the evolution with PIV

% adjust...
