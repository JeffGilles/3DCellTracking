function stack2DTime(varargin)

% ----------------------------
%    STACK 2D + TIME VIEWER
% ----------------------------

% === Parameters ==========================================================

clc
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure');

p = inputParser;
addRequired(p, 'Stack', @isnumeric);
addRequired(p, 'grid', @isnumeric);
addRequired(p, 'piv', @isnumeric);

parse(p, varargin{:});
Stack = p.Results.Stack;
grid = p.Results.grid;
piv = p.Results.piv;

% Visualization
% pointer = struct('ext',  5:20, 'color', [1 0 1]);

% === Initialization ======================================================

% --- Stack dimensions

Stack = flip(Stack,3);
[Ny, Nx, Nt] = size(Stack);

% --- piv dimensions


% --- Shared variables

Img = NaN;

% --- Display

figure(1)
clf
% hold on

% Axis
ax = axes('units', 'pixels', 'Position', [0 0 1 1]);

% Title
tl = uicontrol('style', 'text', 'units','pixel', 'position', [0 0 1 1]);

% --- Controls

st = uicontrol('style','slider', 'position', [0 0 1 1], ...
    'min', 1, 'max', Nt, 'value', 1, 'SliderStep', [1 1]/(Nt+1));

% sz = uicontrol('style','slider', 'position', [0 0 1 1], ...
%      'min', 1, 'max', Nt, 'value', 1, 'SliderStep', [1 1]/(Nt+1));

% Control size callback
set(1, 'ResizeFcn', @updateControlSize);
updateControlSize();

addlistener(st, 'Value', 'PostSet', @updateImage);
% addlistener(sz, 'Value', 'PostSet', @updateImage);

updateImage();

% === Controls ============================================================

    function updateControlSize(varargin)
       
        % Get figure size       
        tmp = get(1, 'Outerposition');
        w = tmp(3);
        h = tmp(4);
        
        % Set widgets size
        st.Position = [35 10 w-45 20];
%         sz.Position = [10 35 20 h-45];
        ax.Position = [75 75 w-100 h-150];
        tl.Position = [w/2-100 h-70 200 20];
        
    end

% === Image ===============================================================

    function updateImage(varargin)
        
        % --- Get sliders values
        
%         zi = round(get(sz, 'Value'));
        ti = round(get(st, 'Value'));
        
        % --- Image
        
        Img = double(Stack(:,:, ti));
        
        % --- Display
        
        cla
        
        imshow(Img);
        
        colormap(flipud(hot));
        
        axis on xy tight
        caxis([0 255])
        colorbar
        
%         tl.String = ['z = ' num2str(zi) ' ; t = ' num2str(ti) ];
        tl.String = ['t = ' num2str(ti) ];
        
        hold on;
        
        quiver(grid(:,1), grid(:,2), piv(:,1,ti), piv(:,2,ti));
        
        hold off;
        
    end

end