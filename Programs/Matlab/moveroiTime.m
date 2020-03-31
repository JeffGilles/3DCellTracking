%function [dx, dy, dz] = moveRoi(roi1, roi2, varargin)

% makes the interpolation between Roi1 & Roi2 with the number of time
% between the two positions

%% ------------ Init parameters -------------------------------------------

nL = 54;
nT = 100;
fDir = '/Users/jean-francoisgilles/Documents/Data/2020/Pauline/';

run = 'MAX_180503_emb2-C1_z54t208';
runShape = '180503_emb2-C1_z54t208';

PIVDir = [fDir '180503_emb2-C2_z54t208/'];

fname1 = ['/Users/jean-francoisgilles/Documents/Data/2020/Pauline/' run '.tif'];
fname2 = ['/Users/jean-francoisgilles/Documents/Data/2020/Pauline/' runShape '.tif'];

factor = 3;

gridSize = 25;
extSize = 45;

min_Int_xy = 0.05;
min_Int_z = 0.08;
min_corr = 0.65;

% Display
df = 5;

Info = imfinfo(fname1);
nf = numel(Info);
w = Info(1).Width;
h = Info(1).Height;

force = false;

% =========================================================================

% --- Grid ----------------------------------------------------------------

% --- 1D coordinates

xg = 1:gridSize:w;
yg = 1:gridSize:h;
% zg = round(1:gridSize/factor:nL);
%zg = 100;

% --- 2D slices

[X_, Y_] = meshgrid(xg, yg);
Gxy = [X_(:) Y_(:)];

Nxy = size(Gxy,1);

%% ------- Load images -------------

if ~exist('raw', 'var') || force
    raw =  NaN(h, w,nL, nT); %Channel 1
    count = 1;
    for ti=1:nT
        for zi=1:nL
            raw(:,:, zi, ti) = double(imread(fname2, count))/255;
            count = count+1;
        end
    end
    raw2D =  NaN(h, w, nT);
    for ti=1:nT
            raw2D(:,:, ti) = double(imread(fname1, ti))/255;
    end
end


if ~exist('piv', 'var') || force
    piv = zeros(length(X_(:)), 2, nT);
    for t = 1:nT
        T1 = load([PIVDir num2str(t, '%03i') '.mat']);
        indZ = find(T1.Zg==26);
        piv(:,2,t) = T1.Dx(indZ);
        piv(:,1,t) = T1.Dy(indZ);
    end
end


% ----------- ROI ------------------

po1 = [296 338; 234 330; 214 278; 260 246; 342 222; 400 246; 416 310; 344 362];
%x,y
po2 = [250 316; 186 296; 186 242; 259 188; 332 168; 400 170; 446 220; 456 270; 372 306];

if length(po1) >= length(po2)
    po2 = Functions.reshapePolygon(po1, po2);
else
    po1 = Functions.reshapePolygon(po1, po2);
end

% --------- MOVE point -------------------

roiEvoTemp = zeros(length(po1), 2, nT);%(x, y, z)
roiEvo = zeros(length(po1), 2, nT);
%evoLinear = zeros(length(po1), 2, nL);

dx = zeros(length(Gxy), nT);%(x, y, z)
dy = zeros(length(Gxy), nT);%(x, y, z)

dy(:, :) = piv(:, 1, :);
dx = piv(:, 2, :);

for ti = 1:nT
    gridX = reshape(dx(:, ti), size(X_(:, :)));
    gridY = reshape(dy(:, ti), size(Y_(:, :)));
    
    gridX(isnan(gridX)) = 0;
    gridY(isnan(gridY)) = 0;
    
    gx = griddedInterpolant(X_(:, :)', Y_(:, :)', gridX');
    gy = griddedInterpolant(X_(:, :)', Y_(:, :)', gridY');
    
    for p = 1:length(po1)
        if ti == 1
            roiEvoTemp(p, :, 1) = po1(p, :);
        else
            roiEvoTemp(p, 1, ti) = roiEvoTemp(p, 1, ti-1) + gx(roiEvoTemp(p, 1, ti-1), roiEvoTemp(p, 2, ti-1));%x
            roiEvoTemp(p, 2, ti) = roiEvoTemp(p, 2, ti-1) + gy(roiEvoTemp(p, 1, ti-1), roiEvoTemp(p, 2, ti-1));%y
        end
    end
end

% ---- Result move with drawing indication

dt = po2-po1;%from roi
D = roiEvoTemp(:, :, end)-roiEvoTemp(:, :, 1);

for t = 1:nT
    for p = 1:length(po1)
        
        roiEvo(p, 1, t) = po1(p,1)+(roiEvoTemp(p, 1, t)-roiEvoTemp(p, 1, 1))*dt(p, 1)/D(p, 1);
        roiEvo(p, 2, t) = po1(p,2)+(roiEvoTemp(p, 2, t)-roiEvoTemp(p, 2, 1))*dt(p, 2)/D(p, 2);
            %evoLinear (:, :, z) = (po2(:, :)-po1(:, :))*z/nL;
        if (roiEvo(p, :, t) - roiEvo(p, :, 1)) > dt(p, :)   % Linear correction
            roiEvo(p, :, t) = po1(p,:)+((po2(p, :)-po1(p, :)))*t/nT;
        end
    end
end

%% ---------- LAST DISPLAY ----------

hold on;

o1 = 'PIV-ROI';
%o2 = 'out2';

View.stack2DTimeroi(raw2D*255, o1, Gxy, piv, roiEvo);


