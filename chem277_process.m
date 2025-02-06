clc
clear

%% Manually detect scale bars where needed
%which_image = "data/all_pyrolysis_data/c3_346a_2247K_4.5atm_0023.tif";
%which_image = "data/all_pyrolysis_data/c3_347a_2132K_4.2atm_0001.tif";
%which_image = "data/all_pyrolysis_data/c3_348b_1875K_4.5atm_0001.tif";
%which_image = "data/all_pyrolysis_data/c3_344b_2203K_4.9atm_0005.tif";
%which_image = "data/all_pyrolysis_data/c3_344b_2203K_4.9atm_0026.tif";

%disp("Output from manual UI scale bar identification:")
%disp(tools.ui_scale_bar(imread(which_image)));

%% Store results from manual scale bar detection in a map
% These values are the outputs of tools.ui_scale_bar
% They're hardcoded so we don't have to manually do the rescale every time.
% Use the code cell above to do the re-scale if desired.
override_scale = containers.Map();
override_scale("c3_346a_2247K_4.5atm_0023.tif")=2.3680;
override_scale("c3_347a_2132K_4.2atm_0001.tif")=0.3060;
override_scale("c3_348b_1875K_4.5atm_0001.tif")=2.3686;
override_scale("c3_344b_2203K_4.9atm_0005.tif")=0.3055;
override_scale("c3_344b_2203K_4.9atm_0026.tif")=0.3054;

%% Auto-detect scale bars on all images except those in the above map.

%[Imgs, imgs, pixsizes] = tools.load_imgs('data/single_test_images', 1);
%[Imgs, imgs, pixsizes] = tools.load_imgs('data/test_folder_1');
[Imgs, imgs, pixsizes] = tools.load_imgs('data/all_pyrolysis_data');
disp("Autodetection of scale bars finished running successfully.")
%imgs = {Imgs.cropped};                              % copy variables locally

%disp("Pixel sizes found: ")
%disp(pixsizes)
%disp("  ")

%% Override the detection results where necessary
for i=1:length(pixsizes)
    if override_scale.isKey(Imgs(i).fname)
        % Override the pixel size and report that we did it
        pixsizes(i)=override_scale(Imgs(i).fname);
        disp("Overrode scale of "+string(Imgs(i).fname)+" to "+string(override_scale(Imgs(i).fname)))
        % Label the validation image so we know we already took care of it
        RGB = imread("validation/"+string(Imgs(i).fname));
        RGB = insertText(RGB,[50,100],"OVERRIDDEN","FontSize",50,"BoxColor","red");
        imwrite(RGB,("validation/"+string(Imgs(i).fname)));
    end
end

% Report any remaining failures
if any(isnan(pixsizes(:)))
    disp("Scale detection failed on: ")
    for i=1:length(pixsizes)
        if isnan(pixsizes(i))
            fn = Imgs(i).fname;
            disp(fn)
        end
    end
else
    disp("Scale detection (or override) succeeded on all images.")
end

%% Store the names and pixel sizes
fnames=string({Imgs.fname});
table_to_store = table(fnames',pixsizes');
table_to_store.Properties.VariableNames = ["file_name","pixel_size"];
writetable(table_to_store,"processed/cache_pixel_sizes_latest_run.csv");
disp("Wrote pixel sizes to .csv as a backup.")

%% Bail since I'm only testing the scale bar finding for now and I always
% instinctively hit run-all instead of run-section :P 
return %GET RID OF THIS LINE TO PROCEED WITH PROCESSING!!!!

%% Return values
%pixsizes = [Imgs.pixsize];   % Shouldn't need this since it's stored
%before.
fname = {Imgs.fname};

%% Do Kmeans, return segmented image
% Notes: didnt work great on 346a
imgs_binary = agg.seg_kmeans(imgs, pixsizes);
%% Agglomerate analysis
Aggs = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname);

%% Particle scale analysis
Aggs = pp.pcm(Aggs); % apply pair correlation method

%% Save the data
tools.write_excel(Aggs, strcat('processed/all_pyrolysis/kmeans/process_results.xlsx'));

tools.imwrite_agg(Aggs, 'processed/all_pyrolysis/kmeans')