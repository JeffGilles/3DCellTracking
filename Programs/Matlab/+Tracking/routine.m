%#ok<*AGROW>
clc
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure')

% === Parameters ==========================================================

study = '190424 D1_16b';
run = 'P2t';
fIP = 'LowMag';

% Tracking
fnumel = [5 Inf];

force = false;
savefile = true;

% -------------------------------------------------------------------------

DS = dataSource;
imDir = [DS.root study filesep run filesep 'Images' filesep];
fDir = [DS.root study filesep run filesep 'Files' filesep];

% =========================================================================

% --- Preparation ---------------------------------------------------------

D = dir(imDir);
D([D.isdir]) = [];
n = numel(D);
ext = D(1).name(end-3:end);

lImg = @(i) double(imread([imDir 'frame_' num2str(i, '%06i') ext]))/255;

% --- Background ----------------------------------------------------------

if ~exist('Bkg', 'var') % || force

    fprintf('Computing background ...');
    tic
    
    Bkg = IP.Bkg.(fIP)(lImg(1));

    fprintf(' %.02f sec\n', toc);

end

% --- Localization --------------------------------------------------------

if ~exist('X', 'var')  % || force
    
    fprintf('Image processing .');
    tic
    
    T = [];
    X = [];
    Y = [];
    F = [];
    
    for i = 1:n
        
        % Process image
        B = IP.(fIP)(lImg(i) - Bkg);
        
        % Prepare tracking inputs
        np = numel(B);
        pos = reshape([B.pos], [2 np])';
                        
        T = [T ; i*ones(np,1)];
        X = [X ; pos(:,1)];
        Y = [Y ; pos(:,2)];
        F = [F ; [B.fluo]'];
        
        if ~mod(i, round(n/10)), fprintf('.'); end
        
    end
    
    fprintf(' %.02f sec\n', toc);
    
end

% --- Tracking ------------------------------------------------------------

if ~exist('Tr', 'var') || force
    
    % --- Stats
    
    % Tracking.pdf_distance_n1([T X Y]);
    
    % --- Definitions
    
    Tr = Tracking.Tracker;

    Tr.parameter('position', 'max', 30, 'norm', 10);
    Tr.parameter('intensity', 'max', 0.2, 'norm', 0.01);

%   Tr.parameter('position', 'hard', 'n1', 'max', 40, 'norm', 10);
%   Tr.parameter('intensity', 'active', false);
    
    % --- Tracking
    
    for i = 1:n
    
        Idx = T==i;
        Tr.set('position', [X(Idx) Y(Idx)]);
        Tr.set('intensity', F(Idx));
        
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

% === Display =============================================================

Tr.disp

figure(1)
clf
hold on

for k = 1:numel(Tr.traj)
    plot3(Tr.traj(k).position(:,1), ...
        Tr.traj(k).position(:,2), ...
        Tr.traj(k).t, '-');   
end

axis ij image on
daspect([1 1 0.3])
view(40,20)

box on