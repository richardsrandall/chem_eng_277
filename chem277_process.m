clc
clear

% To Dos
% 1) Aggs contains all needed data, save to .mat or .csv 

%% Try scale bar first
manual = imread('Data\346a 2247K 4.5atm\c3_346a_2247K_4.5atm_0016.tif');
tools.ui_scale_bar(manual)


%% Load image and attempt to find scale in footer:
Imgs = tools.load_imgs('Data\346a 2247K 4.5atm');
imgs = {Imgs.cropped};                              % copy variables locally

%% Try manually adding scale bar
% Lets us crop scale bar and click on left/right to specify length
% scale = 0.6157; % nm/px
scale_346a_0016 = 0.6146; % nm/ox

% Manually assign for now
Imgs.pixsize = scale_346a_0016;

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
tools.write_excel(Aggs, strcat('processed/5pCH4_2247K/kmeans/process_results.xlsx'));

tools.imwrite_agg(Aggs, 'processed/5pCH4_2247K/kmeans')