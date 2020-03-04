clc
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure');

% === Parameters ==========================================================

run = '180503_emb1';

fname = ['/home/ljp/Science/Projects/ImageAnalysis/3DCellTracking/Data/' run '.tif'];
outDir = ['/home/ljp/Science/Projects/ImageAnalysis/3DCellTracking/Data/Files/' run '/'];

nL = 54;
factor = 3;

gridSize = 25;
extSize = 45;

min_Int_xy = 0.05;
min_Int_z = 0.08;
min_corr = 0.65;

% Display
df = 5;

force = true;

% -------------------------------------------------------------------------

Info = imfinfo(fname);
nf = numel(Info);
w = Info(1).Width;
h = Info(1).Height;

nf

return

% =========================================================================

% --- Grid ----------------------------------------------------------------

% --- 1D coordinates

xg = 1:gridSize:w;
yg = 1:gridSize:h;
zg = round(1:gridSize/factor:nL);

% --- 2D slices

[X_, Y_] = meshgrid(xg, yg);
Gxy = [X_(:) Y_(:)];

[X_, Z_] = meshgrid(xg, zg);
Gxz = [X_(:) Z_(:)];

Nxy = size(Gxy,1);
Nxz = size(Gxz,1);

% --- 3D grid

[X_, Y_, Z_] = meshgrid(xg, yg, zg);
Xg = X_(:);
Yg = Y_(:);
Zg = Z_(:);

for ti = 1%:50
    
    fprintf('=== %i =======================================\n\n', ti);
    
    % --- Raw stack -----------------------------------------------------------
    
    if ~exist('Raw1', 'var') || force
        
        fprintf('Loading raw images ...');
        tic
        
        Raw1 =  NaN(h, w, nL);
        Raw2 =  NaN(h, w, nL);
        
        for i = 1:nL
            Raw1(:,:,i) = double(imread(fname, nL*(ti-1)+i))/255;
            Raw2(:,:,i) = double(imread(fname, nL*(ti)+i))/255;
        end
        
        fprintf(' %.02f sec\n', toc);
        
    end
    
    % --- XY displacements ----------------------------------------------------
    
    % --- All layers
    
    if ~exist('Dx', 'var') || force
        
        fprintf('Computing xy displacements ...');
        tic
        
        dx = zeros(Nxy, nL);
        dy = zeros(Nxy, nL);
        Cxy = zeros(Nxy, nL);
        
        for l = 1:nL
            
            for k = 1:Nxy
                
                i1 = max(Gxy(k,2)-extSize, 1);
                i2 = min(Gxy(k,2)+extSize, h);
                j1 = max(Gxy(k,1)-extSize, 1);
                j2 = min(Gxy(k,1)+extSize, w);
                
                Sub1 = Raw1(i1:i2, j1:j2, l);
                Sub2 = Raw2(i1:i2, j1:j2, l);
                
                if mean(Sub1(:))<min_Int_xy || mean(Sub2(:))<min_Int_xy
                    continue
                end
                
                [dx(k,l), dy(k,l), Cxy(k,l)] = fcorr(Sub1, Sub2);
                
            end
            
        end
        
        Cxy(Cxy<min_corr) = 0;
        
        fprintf(' %.02f sec\n', toc);
        
        fprintf('Averaging xy displacements ...');
        tic
        
        % --- z-grid averaging
        
        Dx = [];
        Dy = [];
        
        for k = zg
            
            k1 = max(round(k-extSize/factor/2), 1);
            k2 = min(round(k+extSize/factor/2), nL);
            
            dx_ = sum(dx(:,k1:k2).*Cxy(:,k1:k2),2)./sum(Cxy(:,k1:k2),2);
            dy_ = sum(dy(:,k1:k2).*Cxy(:,k1:k2),2)./sum(Cxy(:,k1:k2),2);
            
            Dx = [Dx ; dx_];
            Dy = [Dy ; dy_];
            
        end
        
        fprintf(' %.02f sec\n', toc);
        
        
    end
    
    % --- Z displacements -----------------------------------------------------
    
    % --- All layers
    
    if ~exist('Dz', 'var') || force
        
        fprintf('Computing z displacements ...');
        tic
        
        dx2 = zeros(Nxz, h);
        dz = zeros(Nxz, h);
        Cz = zeros(Nxz, h);
        
        for l = 1:h
            
            for k = 1:Nxz
                
                i1 = max(Gxz(k,2)-extSize/factor, 1);
                i2 = min(Gxz(k,2)+extSize/factor, nL);
                j1 = max(Gxz(k,1)-extSize, 1);
                j2 = min(Gxz(k,1)+extSize, w);
                
                Sub1 = permute(Raw1(l, j1:j2, i1:i2), [3 2 1]);
                Sub2 = permute(Raw2(l, j1:j2, i1:i2), [3 2 1]);
                
                if mean(Sub1(:))<min_Int_z || mean(Sub2(:))<min_Int_z
                    continue
                end
                
                [dx2(k,l), dz(k,l), Cz(k,l)] = fcorr(Sub1, Sub2);
                
            end
            
        end
        
        Cz(Cz<min_corr) = 0;
        
        fprintf(' %.02f sec\n', toc);
        
        fprintf('Averaging z displacements ...');
        tic
        
        % --- z-grid averaging
        
        Dz = [];
        
        for k = yg
            
            k1 = max(round(k-extSize/2), 1);
            k2 = min(round(k+extSize/2), h);
            
            dz_ = sum(dz(:,k1:k2).*Cz(:,k1:k2),2)./sum(Cz(:,k1:k2),2)*factor;
            Dz = [Dz ; dz_];
            
        end
        
        % Reorganize
        tmp = permute(reshape(Dz, [numel(zg) numel(xg) numel(yg)]), [3 2 1]);
        Dz = tmp(:);
        
        fprintf(' %.02f sec\n\n', toc);
        
    end
    
    % --- Save ------------------------------------------------------------
    
%     save([outDir num2str(ti, '%03i') '.mat'], 'Xg', 'Yg', 'Zg', ...
%         'Dx', 'Dy', 'Dz');
    
end

% --- Display -------------------------------------------------------------

% Regularization of NaNs
Dx(isnan(Dx)) = 0;
Dy(isnan(Dy)) = 0;
Dz(isnan(Dz)) = 0;

clf
hold on

imshow(Raw1(:,:,end));

q = quiver3(Xg, Yg, nL+1-Zg, Dx*df, Dy*df, Dz*df/factor, 0);

colormap(flipud(hot));
View.colorquiver(q)

colormap(gray);

% box on
grid on

axis([1 h 1 w 0 nL], 'xy', 'on');
daspect([1 1 1/factor])

view(-100, 35)

% % % subplot(1,2,1)
% % % hold on
% % %
% % % % l = 25;
% % %
% % % imshowpair(Raw1(:,:,l), Raw2(:,:,l))
% % % quiver(Gxy(:,1), Gxy(:,2), dx(:,l)*sf, dy(:,l)*sf, 0, 'color', 'y');
% % % axis([1 h 1 w], 'xy', 'on');
% % %
% % % subplot(1,2,2)
% % % hold on
% % %
% % % C = reshape(c(:,l), [numel(yg) numel(xg)]);
% % % C(C<min_corr) = 0;
% % %
% % % imshow(C)
% % %
% % % axis xy on tight
% % % caxis([0 1])
% % % colorbar
