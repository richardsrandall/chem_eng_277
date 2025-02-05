clc
clear

% To Dos
% 1) Aggs contains all needed data, save to .mat or .csv 

%% Try scale bar first
%manual = imread('data/single_test_images/test_image_1.tif');
%tools.ui_scale_bar(manual)
%hardcode_scale = 0.8081;


%% Load image and attempt to find scale in footer:
%[Imgs, imgs, pixsizes] = tools.load_imgs('data/single_test_images', 1);
[Imgs, imgs, pixsizes] = tools.load_imgs('data/test_folder_1');
%imgs = {Imgs.cropped};                              % copy variables locally

disp("Pixel sizes found: ")
disp(pixsizes)
disp("  ")

if any(isnan(pixsizes(:)))
    disp("Scale detection failed on: ")
    for i=1:length(pixsizes)
        if isnan(pixsizes(i))
            disp(Imgs(i).fname)
            %open(strcat('data/test_folder_1/',Imgs(i).fname))
        end
    end
end

%% Try manually adding scale bar

% Manually assign for now
%Imgs.pixsizes = {hardcode_scale};

%% Return values
pixsizes = [Imgs.pixsize];                          % pixel size for each image
fname = {Imgs.fname};

%% Do Kmeans, return segmented image
% Notes: didnt work great on 346a
imgs_binary = agg.seg_kmeans(imgs, pixsizes);
%% Agglomerate analysis
Aggs = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname);

%% Particle scale analysis
Aggs = pp.pcm(Aggs); % apply pair correlation method

%% Save the data
tools.write_excel(Aggs, strcat('processed/test_1/kmeans/process_results.xlsx'));

tools.imwrite_agg(Aggs, 'processed/test_1/kmeans')