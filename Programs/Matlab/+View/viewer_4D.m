function viewer_4D(Data, points)

% ------------------------
%       VIEWER 4D
% ------------------------

clc
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure');

% === Initialization ======================================================

% --- Load demo data

%Data = load('mri');

% --- Stack dimensions

[Ny, Nx, Nz, Nt] = size(Data);

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

sz = uicontrol('style','slider', 'position', [0 0 1 1], ...
     'min', 1, 'max', Nz, 'value', 1, 'SliderStep', [1 1]/(Nz+1));

% Control size callback
set(1, 'ResizeFcn', @updateControlSize);
updateControlSize();

addlistener(st, 'Value', 'PostSet', @updateImage);
addlistener(sz, 'Value', 'PostSet', @updateImage);

updateImage();

% === Controls ============================================================

    function updateControlSize(varargin)
       
        % Get figure size       
        tmp = get(1, 'Outerposition');
        w = tmp(3);
        h = tmp(4);
        
        % Set widgets size
        st.Position = [35 10 w-45 20];
        sz.Position = [10 35 20 h-45];
        ax.Position = [75 75 w-100 h-150];
        tl.Position = [w/2-100 h-70 200 20];
        
    end

% === Image ===============================================================

    function updateImage(varargin)
        
        % --- Get sliders values
        
        zi = round(get(sz, 'Value'));
        ti = round(get(st, 'Value'));
        
        % --- Image
        
        Img = double(Data(:,:, zi, ti));
        
        % --- Display
        
        cla
        
        imshow(Img);
        
        axis on xy tight
        caxis([0 255])
        colorbar
        
        tl.String = ['z = ' num2str(zi) ' ; t = ' num2str(ti) ];
        positionP = cell2mat(points(ti));
        
        ind = find(positionP(:, 3) >= (Nz+1-zi)-1 & positionP(:, 3) <= (Nz+1-zi)+1);
        
        %diam = find(positionP(:, 3) >= (Nz+1-zi));
%         if positionP(:, 3) == (Nz+1-zi)
%             diam = 40;
%         end
        hold on;
        
        
        scatter(positionP(ind, 1), positionP(ind, 2), 30, 'filled', 'red');
        
    end

end