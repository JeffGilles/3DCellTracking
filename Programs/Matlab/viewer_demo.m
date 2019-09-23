function viewer_demo()

% ----------------------------
%       VIEWER DEMO
% ----------------------------

clc
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure');

% === Initialization ======================================================

% --- Load demo data

Data = load('mri');

% --- Stack dimensions

Nx = Data.siz(1);
Ny = Data.siz(2);
Nz = Data.siz(3);
Nt = 100;

% --- Shared variables

Img = NaN;

% --- Display

figure(1)
clf
% hold on

% Axis
ax = axes('units', 'pixels', 'Position', [0 0 1 1]);

% Title
ti = uicontrol('style', 'text', 'units','pixel', 'position', [0 0 1 1]);

% --- Controls

st = uicontrol('style','slider', 'position', [0 0 1 1], ...
    'min', 1, 'max', Nt, 'value', 1, 'SliderStep', [1 1]/(Nz+1));

sz = uicontrol('style','slider', 'position', [0 0 1 1], ...
     'min', 1, 'max', Nz, 'value', 1, 'SliderStep', [1 1]/(Nt+1));

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
        ti.Position = [100 h-75 200 20];
        ti.String = 'ok';
        
    end

% === Image ===============================================================

    function updateImage(varargin)
        
        % --- Get tranlation and rotation values
        
        zi = round(get(sz, 'Value'));
        ti = round(get(st, 'Value'));
                
        % --- Move model
        Img = Image.rotate(M.state(1).img, param.alpha, 'center', M.state(1).center);
        Img = Image.translate(Img, param.tx, param.ty);
        
        % --- Locate objects
        findObject()
        
        % --- Display
        updateDisplay();
        
    end

% === Object ==============================================================

    function findObject()
        
        tic
        
        % Remove small objects
        CC = bwconncomp(Img>=0.5);
        [~, mi] = max(cellfun(@numel, CC.PixelIdxList));
        Tmp = Img;
        for j = setdiff(1:numel(CC.PixelIdxList), mi)
            Tmp(CC.PixelIdxList{j}) = 0;
        end
        
        % Define state and get contours
        S = State(Tmp);
        S.getContour();
                        
        % Get best match
       [Ci, E, D2] = Cr.fiton(S.contour, ref);
                   
        % Update computation time
        set(ptime, 'string', [ '[' num2str(param.tx) ...
            ', ' num2str(param.ty) ...
            ', ' num2str(param.alpha) ...
            '] | Computation time: ' num2str(toc*1000, '%.02f') ' ms']);
        
    end

% === Display =============================================================

    function updateDisplay()
        
        figure(1)
        cla
                   
        if gom>2
            
            imshow(cat(3,Img, zeros(size(Img)), zeros(size(Img))));
            if ~isnumeric(Ci), Ci.plot('y-'); end
                        
            param
            
        else
            
            imshow(Img);
            if ~isnumeric(Ci), Ci.plot('r.-'); end
            
            
        end
        
        axis xy tight
        
    end

end