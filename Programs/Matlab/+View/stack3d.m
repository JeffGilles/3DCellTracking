function stack3d(varargin)

% ----------------------------
%       3D STACK VIEWER
% ----------------------------

% === Parameters ==========================================================

p = inputParser;
addRequired(p, 'Stack', @isnumeric);
addParameter(p, 'zfactor', 1, @isnumeric);
addParameter(p, 'points', NaN(0,3), @isnumeric);
addParameter(p, 'ndist', 10, @isnumeric);

parse(p, varargin{:});
Stack = p.Results.Stack;
zf = p.Results.zfactor;
P = p.Results.points;
ndist = p.Results.ndist;

% Visualization
space = 5;
pdiam = 8;
pointer = struct('ext',  5:20, 'color', [1 0 1]);

% === Initialization ======================================================

Stack = flip(Stack,3);
[Ny, Nx, Nz] = size(Stack);

x = Nx/2;
y = Ny/2;
z = Nz/2;

block = false;

% --- Display

figure(1)
clf
% hold on

% Axis
ax = axes('units', 'pixels', 'Position', [0 0 1 1]);

% Title
tl = uicontrol('style', 'text', 'units','pixel', 'position', [0 0 1 1]);

% --- Controls

% Sliders

sx = uicontrol('style','slider', 'position', [0 0 1 1], ...
    'min', 1, 'max', Nx, 'value', x, 'SliderStep', [1 1]/(Nx+1));

sy = uicontrol('style','slider', 'position', [0 0 1 1], ...
    'min', 1, 'max', Ny, 'value', y, 'SliderStep', [1 1]/(Ny+1));

sz = uicontrol('style','slider', 'position', [0 0 1 1], ...
    'min', 1, 'max', Nz, 'value', z, 'SliderStep', [1 1]/(Nz+1));

addlistener(sx, 'Value', 'PostSet', @updateImage);
addlistener(sy, 'Value', 'PostSet', @updateImage);
addlistener(sz, 'Value', 'PostSet', @updateImage);

% Control size callback
set(1, 'ResizeFcn', @updateControlSize);
updateControlSize();

updateImage();

% === Controls ============================================================

    function updatePosition(varargin)
        
        event = varargin{2};
        i = round(event.IntersectionPoint(2));
        j = round(event.IntersectionPoint(1));
        
        if i<=zf*Nz                 % XZ image
            
            x = min(j, Nx);
            z = round(i/zf);
            
        elseif j>Nx+space           % YZ image
            
            y = max(1, i - zf*Nz - space);
            z = Nz + 1 - round((j - Nx - space)/zf);
            
        elseif i>zf*Nz+space && j<Nx   % XY image
            
            x = j;
            y = i - zf*Nz - space;
            
        end
        
        block = true;
        sx.Value = x;
        sy.Value = y;
        sz.Value = z;
        block = false;
        
        updateImage();
    end

    function updateControlSize(varargin)
        
        % Get figure size
        tmp = get(1, 'Outerposition');
        w = tmp(3);
        h = tmp(4);
        
        % Set widgets size
        
        sx.Position = [35 10 w-80 20];
        sy.Position = [10 35 20 h-45];
        sz.Position = [w-30 35 20 h-45];
        
        ax.Position = [75 75 w-100 h-150];
        tl.Position = [w/2-100 h-70 200 20];
        
    end

% === Image ===============================================================

    function updateImage(varargin)
        
        if block, return; end
        
        % --- Get sliders values
        
        x = round(get(sx, 'Value'));
        y = round(get(sy, 'Value'));
        z = round(get(sz, 'Value'));
        
        % --- Image
        
        Ixy = Stack(:,:,z);
        Ixz = squeeze(Stack(y,:,:))';
        Iyz = flip(squeeze(Stack(:,x,:)),2);
        
        if zf>1
            Ixz = imresize(Ixz, [zf*Nz Nx]);
            Iyz = imresize(Iyz, [Ny zf*Nz]);
        end
        
        Img = ones(Ny+zf*Nz+space, Nx+zf*Nz+space)*240/255;
        
        Img(1:zf*Nz, 1:Nx) = Ixz;
        Img(zf*Nz+space+1:end, 1:Nx) = Ixy;
        Img(zf*Nz+space+1:end, Nx+space+1:end) = Iyz;
        
        % --- Display -----------------------------------------------------
        
        % --- RGB channels
        
        R = Img;
        G = Img;
        B = Img;
        
        % --- Data points
        
        pR = zeros(size(Img));
        %         pG = zeros(size(Img));
        %         pB = zeros(size(Img));
        
        % XY
        Pxy = P(abs(P(:,3)-z)<=ndist/zf, :);
        vxy = 1 - abs(Pxy(:,3)-z)/(ndist/zf+1);
        pR(sub2ind(size(pR), zf*Nz+space+Pxy(:,2), Pxy(:,1))) = vxy;
        
        % XZ
        Pxz = P(abs(P(:,2)-y)<=ndist, :);
        vxz = 1 - abs(Pxz(:,2)-y)/(ndist+1);
        pR(sub2ind(size(pR), zf*Pxz(:,3), Pxz(:,1))) = vxz;
        
        % YZ
        Pyz = P(abs(P(:,1)-x)<=ndist, :);
        vyz = 1 - abs(Pyz(:,1)-x)/(ndist+1);
        pR(sub2ind(size(pR), zf*Nz+space+Pyz(:,2), Nx+space+zf*(Nz+1-Pyz(:,3)))) = vyz;
        
        % Overlay
        Mask = imdilate(pR, strel('disk', round(pdiam/2)));
        R = R + Mask;
        G(Mask>0) = 0;
        B(Mask>0) = 0;
        
        % --- Pointer
        
        pij = [zf*Nz+space+y x ; zf*z x ; zf*Nz+space+y Nx+space+zf*(Nz+1-z)];
        
        % Left
        il = repmat(pij(:,1), [1 numel(pointer.ext)]);
        jl = max(pij(:,2) - pointer.ext, 1);
        
        % Right
        ir = repmat(pij(:,1), [1 numel(pointer.ext)]);
        jr = min(pij(:,2) + pointer.ext, size(Img,2));
        
        % Up
        iu = min(pij(:,1) + pointer.ext, size(Img,1));
        ju = repmat(pij(:,2), [1 numel(pointer.ext)]);
        
        % Down
        id = max(pij(:,1) - pointer.ext, 1);
        jd = repmat(pij(:,2), [1 numel(pointer.ext)]);
        
        % Color
        R(sub2ind(size(Img), [il(:) ; ir(:) ; iu(:) ; id(:)], ...
            [jl(:) ; jr(:) ; ju(:) ; jd(:)])) = pointer.color(1);
        G(sub2ind(size(Img), [il(:) ; ir(:) ; iu(:) ; id(:)], ...
            [jl(:) ; jr(:) ; ju(:) ; jd(:)])) = pointer.color(2);
        B(sub2ind(size(Img), [il(:) ; ir(:) ; iu(:) ; id(:)], ...
            [jl(:) ; jr(:) ; ju(:) ; jd(:)])) = pointer.color(3);
        
        cla
        hold on
        
        % --- Image
        
        h = imshow(cat(3, R, G, B));
        set(h, 'ButtonDownFcn', @updatePosition);
        
        axis xy tight
        
        tl.String = ['x = ' num2str(x) ', y = ' num2str(y) ', z = ' num2str(z)];
        
    end

end