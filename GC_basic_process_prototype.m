%% Prototyping:
% Goal: understand image seg in matlab better
%  - use to experiment with different methods
clear all; close all; clc
%% Load Image
imageName = 'best_images/c3_344b_2203K_4.9atm_0012.tif';
img = imread(sprintf('data/%s',imageName));

%% Convert to Grayscale if needed
if size(img,3) == 3
    img = rgb2gray(img);
end

%% Preprocess Image
img = imadjust(img); % Contrast enhancement
img_filt = medfilt2(img, [3,3]); % Noise reduction

%% Segment Agglomerates using Watershed
bw = imbinarize(img_filt, 'adaptive', 'Sensitivity', 0.5);
bw = imfill(bw, 'holes'); % Fill holes within particles
bw = bwareaopen(bw, 50); % Remove small noise

% Morphological Closing to connect broken regions
% Reduces fragmentation of individual particles.
se = strel('disk', 5);
bw = imclose(bw, se);

% Distance Transform and Watershed
% watershed helps separate overlapping agglomerates.
D = -bwdist(~bw);
D(~bw) = -Inf;
L = watershed(D);
bw(L == 0) = 0;

% Remove cropped objects
% Edge Removal (imclearborder): Excludes cropped agglomerates.
ebw = imclearborder(bw);

%% Label and Measure Agglomerates
stats = regionprops(ebw, 'Area', 'Centroid', 'Perimeter', 'EquivDiameter');
numAgglomerates = length(stats);

%% Compute Projected Area and Mobility Diameter
projectedAreas = [stats.Area]; % In pixels
mobilityDiameters = 2 * sqrt(projectedAreas / pi); % Area-equivalent diameter

%% Estimate Primary Particle Diameters
% Use Euclidean Distance Transform to identify particle cores
distMap = bwdist(~ebw);

% Apply Hough Circle Transform for primary particle detection
% The Hough Transform may be detecting too many small circles, 
% possibly noise or incorrectly identifying overlapping particles
radiiRange = [5 25]; % Adjust based on expected particle sizes
sensitivity = 0.9;
[centers, radii] = imfindcircles(distMap, radiiRange, 'Sensitivity', sensitivity, 'EdgeThreshold', 0.1);

%% Display Results
figure;
imshow(img);
hold on;

% Overlay detected agglomerates
for i = 1:numAgglomerates
    text(stats(i).Centroid(1), stats(i).Centroid(2), num2str(i), 'Color', 'red');
end

% Overlay detected primary particles
viscircles(centers, radii, 'EdgeColor', 'blue');
hold off;

disp('Projected Areas (pixels):');
disp(projectedAreas);
disp('Mobility Diameters (pixels):');
disp(mobilityDiameters);
disp('Estimated Primary Particle Diameters (pixels):');
disp(radii);
