%#ok<*AGROW>
clc
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure')

% === Parameters ==========================================================

tSlices = 2700;
nSlices = 23;
nTime = 97;
fDir = '/Users/jean-francoisgilles/Documents/Data/';
width = 512;
height = 512;
xy = 0.4805308;
xz = 2.0141621;
xzratio = xy/xz; %sizeZ(um)/sizeXY(um)

% Segment

sigmaG = 3; %gaussian
dilsize = 15; %Dilate factor
thr = 2; %Detection threshold %2=optimal value
%rad = 5;

% Tracking
fnumel = [5 Inf];

force = false;
savefile = true;

% -------------------------------------------------------------------------

% DS = dataSource;
% imDir = [DS.root study filesep run filesep 'Images' filesep];
% fDir = [DS.root study filesep run filesep 'Files' filesep];

% =========================================================================

% --- Preparation ---------------------------------------------------------

% D = dir(imDir);
% D([D.isdir]) = [];
% n = numel(D);
% ext = D(1).name(end-3:end);

% lImg = @(i) double(imread([imDir 'frame_' num2str(i, '%06i') ext]))/255;

if ~exist('stack', 'var') || force %if it already exists, don't do this
    stack = zeros(height, width, nSlices, nTime);
    count = 1;
    for t=1:nTime
        for z=1:nSlices
            %stack(:, :, z, t) = imread([fDir '181017_Pos001_with eyes-C2_z23t97.tif'], count);
            stack(:, :, z, t) = imread([fDir 'emb1-C2-1-50.tif'], count);
            count = count+1;
        end
    end
end

% --- Localization --------------------------------------------------------

% if ~exist('X', 'var')  % || force
%     
%     fprintf('Image processing .');
%     tic
%     
%     T = [];
%     X = [];
%     Y = [];
%     F = [];
%     
%     for i = 1:n
%         
%         % Process image
%         B = IP.(fIP)(lImg(i) - Bkg);
%         
%         % Prepare tracking inputs
%         np = numel(B);
%         pos = reshape([B.pos], [2 np])';
%         
%         T = [T ; i*ones(np,1)];
%         X = [X ; pos(:,1)];
%         Y = [Y ; pos(:,2)];
%         F = [F ; [B.fluo]'];
%         
%         if ~mod(i, round(n/10)), fprintf('.'); end
%         
%     end
%     
%     fprintf(' %.02f sec\n', toc);
%     
% end

if ~exist('points', 'var') || force
    d = zeros(height, width, nSlices, nTime);
    points = cell(nTime, 1);
    intensites = cell(nTime, 1);
%     pl = zeros(15, nTime); %decrease related to number points detected
%     pl = zeros(nTime, 1);
    for t=1:nTime
        g = imgaussfilt3(stack(:, :, :, t), round(sigmaG*[1 1 1/xzratio]));
        d(:, :, :, 1) = imdilate(g, strel('cuboid', round(dilsize*[1 1 1/xzratio])))==g;
        %d(:, :, :, 1) = imdilate(g, strel('sphere',7))==g;
        A = find(d);
%         for q=0:14 % graph of decrease nb points
%             A(g(A)<q) = []; %max points from g
%             pl(q+1, t) = size(A, 1);
%         end
        A(g(A)<thr) = [];
%         pl(t) = size(A, 1);
        [I,J,K] = ind2sub(size(d), A);  %localization from indexes
        nb = size(I);
        points{t} = [J*xy, I*xy, ((nSlices+1)-K)*xz]; %cell {}
        intensites{t} = g(A);
    end
end


% --- Tracking ------------------------------------------------------------

if ~exist('Tr', 'var') || force
    
    % --- Stats
    
%     Tracking.pdf_distance_n1([T X Y Z]);
    
    % --- Definitions
    
    Tr = Tracking.Tracker;
    
    Tr.parameter('position', 'max', 10, 'norm', 10);
    Tr.parameter('intensity', 'max', 20, 'norm', 0.01);
    
%   Tr.parameter('position', 'hard', 'n1', 'max', 40, 'norm', 10);
%   Tr.parameter('intensity', 'active', false);
    
    % --- Tracking
    
    for i = 1:nTime
    
%         Idx = T==i;
%         Tr.set('position', [X(Idx) Y(Idx) Z(Idx)]);
%         Tr.set('intensity', F(Idx));
%         Idx = T==i;
        Tr.set('position', [points{i}(:, 1) points{i}(:, 2) points{i}(:, 3)]);
        Tr.set('intensity', [intensites{i}]);
        
        
        Tr.match('method', 'fast', 'verbose', false);
    
    end

    % --- Assemble

    Tr.assemble('method', 'fast', 'max', 10, 'norm', 1);

    % --- Filtering
        
    Tr.filter('numel', fnumel);
    
end



% --- Save ----------------------------------------------------------------

if savefile
    
    save([fDir 'tracking.mat'], 'Tr');
    
end

% ---- Interpolation ------------------------------------------------------


f = scatteredInterpolant(linspace(0, 512), linspace(0, 512), ...
        Tr.traj(k).t, [Tr.traj(k).position(:,1); ...
        Tr.traj(k).position(:,2); ...
        Tr.traj(k).position(:,3)]);
        


% === Display =============================================================

Tr.disp

figure(1)
clf
hold on
for k = 1:numel(Tr.traj)
    plot3(Tr.traj(k).position(:,1), ...
        Tr.traj(k).position(:,2), ...
        Tr.traj(k).t, '-');
%     size(Tr.traj(k).position(:,1))
end
% 
% for k = 2:numel(Tr.traj)
% %     plot3(Tr.traj(k).position(:,1), ...
% %         Tr.traj(k).position(:,2), ...
% %         Tr.traj(k).t, '-');
%     quiver3(Tr.traj(k-1).position(:,1), ...
%         Tr.traj(k-1).position(:,2), ...
%         Tr.traj(k-1).t, ...
%         Tr.traj(k).position(:,1), ...
%         Tr.traj(k).position(:,2), ...
%         Tr.traj(k).t, '-');
% end

axis ij image on
daspect([1 1 0.3])
view(40,20)

box on