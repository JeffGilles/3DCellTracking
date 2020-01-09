
clc
fDir = '/Users/jean-francoisgilles/Documents/Data/';
im = imread([fDir 'eyes-test.tif'], 1);

%Gaussian + substract background 
threshold = 4; 
gaus = imgaussfilt(im, 2);
im2 = uint8(gaus > threshold ); %need to be converted
gaus = gaus .* im2; 

bw = imbinarize(gaus,'adaptive'); %adaptative = bradley's thresholding

% Medial axis skeleton
sk = bwmorph(bw, 'skel', 'inf');

%distance map
dm = bwdist(~bw);

% Medial-Axis Transform
MAT = sk .* dm;

%Shape thinning result
skthin = bwmorph(bw, 'thin', 'inf');
resTh = skthin .* dm;

%show
% multi = cat(2, resTh, MAT);
% montage(multi)

%imshow(bw);
% hold on;
% contour(skthin, 'red'); 
% hold off;

regions = detectMSERFeatures(gaus);
%figure; imshow(gaus); hold on;
%plot(regions,'showPixelList',true,'showEllipses',false);
figure; imshow(gaus); 
hold on;
plot(regions);

title('Thinning result and MAT of ~bw')

set(gca, {'YDir'}, {'reverse'});
axis on xy tight
%caxis([0 10])
%colorbar
%colormap default
