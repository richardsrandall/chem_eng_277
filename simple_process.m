clc
clear


Imgs = tools.load_imgs('TEM_images');
imgs = {Imgs.cropped}; % copy variables locally
pixsizes = [Imgs.pixsize]; % pixel size for each image
fname = {Imgs.fname};

imgs_binary = agg.seg_kmeans(imgs, pixsizes);
Aggs = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname);
Aggs = pp.pcm(Aggs); % apply pair correlation method
tools.write_excel(Aggs, strcat('processed/kmeans/process_results.xlsx'));

tools.imwrite_agg(Aggs, 'processed/kmeans')