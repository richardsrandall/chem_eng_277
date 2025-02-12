clc
clear

%% User Controls

% loading/savingdata
data_Dir_name           = 'all_pyrolysis_data'; % where to pull data from
%data_Dir_name = 'cleaned_pyrolysis_data_1-5x_threshold';
                                         % OR: 'all_pyrolysis_data', 'test_folder_1'
                                         % also where to write Imgs and Excel results

% Writing New Data
cache_pixel_sizes       = true;          % write csv of pixel sizes of latest run into processed/ directory
checkObjDet_diagnostics = true;          % Save side-by-side black/white and greyscale imaged as a diagnostic
writeImgsFolder         = 'all_pyrolysis_data'; % save img results under processed/
writeExcelFolder        = 'all_pyrolysis_data'; % save excel results under processed/


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

%% Try denoising beforehand for specific noisy images
% % Basic, will take some more work
% imgName = 'c3_346a_2247K_4.5atm_0009.tif';
% in_path = sprintf("data/best_images/%s",imgName);
% out_path   = 'data/best_images/denoised/c3_346a_2247K_4.5atm_0009.tif';
% GC_denoise_preprocess(in_path, out_path)


%% Load & Auto-detect scale bars on all images except those in the above map.
location = sprintf('data/%s',data_Dir_name);

[Imgs, imgs, pixsizes] = tools.load_imgs(location); % find ones that work best
%[Imgs, imgs, pixsizes] = tools.load_imgs('data/single_test_images', 1);

% Imgs = 1 x nImgs struct with 6 fields
% (1) & (2) - filename & directory
% (3) Raw image data (1024x1024 px before processing)
% (4) OCR text data (optical character recognition output)
% (5) Cropped image data - some ROI has been cropped. still 1024x1024 px
% (6) pixsize (nm/pixel) for each image
% imgs = 1 x nImgs cellArr. each cell contains the raw image data
% pixsizes = 1 x nImgs double. each entry containing nm/px calibration

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
if cache_pixel_sizes
    fnames=string({Imgs.fname});
    table_to_store = table(fnames',pixsizes');
    table_to_store.Properties.VariableNames = ["file_name","pixel_size"];
    writetable(table_to_store,"processed/cache_pixel_sizes_latest_run.csv");
    disp("Wrote pixel sizes to .csv as a backup.")
end

%% Bail since I'm only testing the scale bar finding for now and I always
% instinctively hit run-all instead of run-section :P 
%return %GET RID OF THIS LINE TO PROCEED WITH PROCESSING!!!!

%% Return values
%pixsizes = [Imgs.pixsize];   % Shouldn't need this since it's stored
%before.
fname = {Imgs.fname};

%% Do Kmeans, return segmented image
% Notes: didnt work great on 346a
imgs_binary = agg.seg_kmeans(imgs, pixsizes);

%% Mess with imgs_binary
%if checkObjDet_diagnostics
%     for i = 1:length(pixsizes)
%        img = (cat(2,(cell2mat(imgs(i))),255*cell2mat(imgs_binary(i))));
%        imwrite(img,("check_object_detection/"+string(Imgs(i).fname)))
%     end
%     disp("Saved side-by-side black/white and greyscale imaged as a diagnostic.")
% end

%% Agglomerate analysis
Aggs = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname);

%% Particle scale analysis
Aggs = pp.pcm(Aggs); % apply pair correlation method

%% Save the data in processed folder
%tools.write_excel(Aggs, strcat(sprintf('processed/%s/kmeans/process_results.xlsx',data_Dir_name)));
tools.imwrite_agg(Aggs, sprintf('processed/%s/kmeans',data_Dir_name))
%close all

%% Analyze aggregates and mess with imgs_binary
agg_fnames = ({Aggs.fname});            % filename associated with each aggregate
disp("Saving diagnostic images...")

% Loop over each image, write a side-by-side actual and binary. Label with
% Green or Red box. Red if no aggregates found by algorithm
for i = 1:length(pixsizes)
    if ismember(Imgs(i).fname,agg_fnames)                       % if no, then current image is NOT in agg_fnames
        loc = sprintf("processed/%s/kmeans/",data_Dir_name);
        og_img = imread((location+"/"+string(Imgs(i).fname)));
        %imwrite(og_img,("data/cleaned_pyrolysis_data_1-5x_threshold/"+string(Imgs(i).fname))); % This is our way of only including images that have aggregates
        if ~isfile(loc + string(Imgs(i).fname))
            disp("FAILURE ON: "+string(Imgs(i).fname))
            continue
        end
        final_img = imread(loc + string(Imgs(i).fname));
        img = 255*cell2mat(imgs_binary(i));
        RGB = insertText(img,[50,100],"     ","FontSize",50,"BoxColor","green");
        size_original = size(img);
        size_final = size(final_img);
        final_img = final_img(:,round(size_final(2)*0.5-size_final(1)*0.5):round(size_final(2)*0.5+size_final(1)*0.5),:);
        final_img_resized = imresize(final_img,[1024,1024]);
        RGB = cat(2,final_img_resized,RGB);
        RGB = insertText(RGB,[50,100],"AGGREGATES FOUND","FontSize",50,"BoxColor","green");
    else
        img = (cat(2,(cell2mat(imgs(i))),255*cell2mat(imgs_binary(i))));    % combine image and resulting binary image
        RGB = insertText(img,[50,100],"NO AGGREGATES","FontSize",50,"BoxColor","red");
    end
    imwrite(RGB,("check_object_detection/"+string(Imgs(i).fname)))
end
disp("Saved annotated side-by-side black/white and greyscale images that contain aggregates.")
