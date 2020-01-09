clc

tSlices = 2700;
nSlices = 54;
nTime = 5;
path = '/Users/jean-francoisgilles/Downloads/';
width = 512;
height = 512;
xzratio = 2.0141621/0.4805308; %sizeZ(um)/sizeXY(um)

sigmaG = 3; %gaussian
dilsize = 15; %Dilate factor
thr = 2; %Detection threshold %2=optimal value
%rad = 5;

if ~exist('stack', 'var') %if it already exists, don't do this
    stack = zeros(height, width, nSlices, nTime);
    count = 1;
    for t=1:nTime
        for z=1:nSlices
            stack(:, :, z, t) = imread([path 'emb1-C2-1-50.tif'], count);
            count = count+1;
        end
    end
end

%%

if ~exist('d', 'var')
    d = zeros(height, width, nSlices, nTime);
    points = cell(nTime,1);
%     pl = zeros(15, nTime); %decrease related to number points detected
    pl = zeros(nTime, 1);
    for t=1:nTime
        g = imgaussfilt3(stack(:, :, :, t), round(sigmaG*[1 1 1/xzratio]));
        %mean(g, 'all')
        d(:, :, :, 1) = imdilate(g, strel('cuboid', round(dilsize*[1 1 1/xzratio])))==g;
        %d(:, :, :, 1) = imdilate(g, strel('sphere',7))==g;
        A = find(d);
        
%         for q=0:14 % graph of decrease nb points
%             A(g(A)<q) = []; %max points from g
%             pl(q+1, t) = size(A, 1);
%         end
        A(g(A)<thr) = [];
        pl(t) = size(A, 1);
        [I,J,K] = ind2sub(size(d), A);  %localization from indexes
        points{t} = [J, I, (nSlices+1)-K]; %cell {}
    end
end

% clf
% hist(stack(A),50)
% 
% return

%%
%viewer_3D([path 'emb1-C2-1-50.tif'], points{1}); 

% figure  %plot graph decrease
% plot(0:14, pl)
% set(gca, 'YScale', 'log')
View.viewer_4D(stack, points);

%figure
%plot(0:nTime-1, pl)

% for t=1:nTime
%     size(points{t}(:, 1))
% end
%save([path 'centers.mat'], 'points');