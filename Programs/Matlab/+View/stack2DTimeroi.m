function stack2DTimeroi(varargin)

% ----------------------------
%    STACK 2D + TIME VIEWER
% ----------------------------

% === Parameters ==========================================================

%clc
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure');

p = inputParser;
addRequired(p, 'Stack', @isnumeric);
addRequired(p, 'output', @ischar);
addRequired(p, 'grid', @isnumeric);
addRequired(p, 'piv', @isnumeric);
addOptional(p, 'roi', @isnumeric);


parse(p, varargin{:});
Stack = p.Results.Stack;
grid = p.Results.grid;
piv = p.Results.piv;
roi = p.Results.roi;
output = p.Results.output;

% Visualization
% pointer = struct('ext',  5:20, 'color', [1 0 1]);

% === Initialization ======================================================

% --- Stack dimensions

%Stack = flip(Stack,3);
[Ny, Nx, Nt] = size(Stack);

% --- Shared variables

Img = NaN;

% --- Display

figure(1)
clf
%hold on

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

% Set up the movie.
fDir = '/Users/jean-francoisgilles/Documents/Data/2020/Pauline/';
writerObj = VideoWriter([fDir output '.avi']); % Name it.
writerObj.FrameRate = 5; % How many frames per second.
open(writerObj);


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
        
        for ti=1:Nt
            % --- Get sliders values

    %         zi = round(get(sz, 'Value'));
            %ti = round(get(st, 'Value'));

            % --- Image

            Img = double(Stack(:,:, ti));

            % --- Display

            cla

            imshow(Img);

            colormap(flipud(hot));

            %axis on xy tight
            caxis([0 255])
            colorbar

    %         tl.String = ['z = ' num2str(zi) ' ; t = ' num2str(ti) ];
            %tl.String = ['t' num2str(ti) ];

            hold on;
            
            txtnumber = {'t' ti};
            ylimits = ylim;
            ymax = ylimits(2);
            spacing = ymax/20; %for different lines
            text(10, ymax-spacing*1, txtnumber); % *2 for line 2
            
            quiver(grid(:,2), grid(:,1), piv(:,2,ti), piv(:,1,ti));%s!! sens grille si pb affichage
            
%             x = roi(:, 1, ti);
%             y = roi(:, 2, ti);
%             x(end+1) = x(1);
%             y(end+1) = y(1);
%             plot(x, y, 'g.-', 'MarkerSize', 20);
            drawpolygon('Position', roi(:, :, ti));
            
%             if ti==1
%                 roi(:, :, ti)
%             end
            
            hold off;

            %axis([1 size(Img,2) 1 size(Img,1)]);

            pause(0.1);
            %updateImage();

            frame = getframe(1); % 'gcf' can handle if you zoom in to take a movie.
            writeVideo(writerObj, frame);
            
        end
        
    end

% ===== Movie =============================================================


%updateImage();

    

end