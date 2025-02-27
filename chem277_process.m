clc
clear
close all

%% User Controls

% Enter the name of the folder we full data from AND write it to
%data_dir_name           = 'validation_data';
%data_dir_name           = 'small_validation_data';
data_dir_name            = 'all_pyrolysis_data';

% Writing New Data
cache_pixel_sizes       = true;          % write csv of pixel sizes of latest run into processed/ directory
checkObjDet_diagnostics = true;          % Save side-by-side black/white and greyscale imaged as a diagnostic

%% Empty out the validation image folders
delete('validate_scale_bar_detection/*')
delete('validate_object_detection/*')

%% Manually detect scale bars where needed
% You only need to do this once, then the results are stored in the dict in
% the next section.
%which_image = "raw_tem_data/all_pyrolysis_data/c3_346a_2247K_4.5atm_0023.tif";
%which_image = "raw_tem_data/all_pyrolysis_data/c3_347a_2132K_4.2atm_0001.tif";
%which_image = "raw_tem_data/all_pyrolysis_data/c3_348b_1875K_4.5atm_0001.tif";
%which_image = "raw_tem_data/all_pyrolysis_data/c3_344b_2203K_4.9atm_0005.tif";
%which_image = "raw_tem_data/all_pyrolysis_data/c3_344b_2203K_4.9atm_0026.tif";

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

%% Load & Auto-detect scale bars on all images except those in the above map.
location = sprintf('raw_tem_data/%s',data_dir_name);
delete('validate_scale_bar_detection/*') % Empty out the validation contents
[Imgs, imgs, pixsizes] = tools.load_imgs(location);

% Imgs = 1 x nImgs struct with 6 fields
% (1) & (2) - filename & directory
% (3) Raw image data (1024x1024 px before processing)
% (4) OCR text data (optical character recognition output)
% (5) Cropped image data - some ROI has been cropped. still 1024x1024 px
% (6) pixsize (nm/pixel) for each image
% imgs = 1 x nImgs cellArr. each cell contains the raw image data
% pixsizes = 1 x nImgs double. each entry containing nm/px calibration

disp("Autodetection of scale bars finished running successfully.")
%disp("Pixel sizes found: ")
%disp(pixsizes)

%% Override the detection results where necessary
for i=1:length(pixsizes)
    if override_scale.isKey(Imgs(i).fname)
        % Override the pixel size and report that we did it
        pixsizes(i)=override_scale(Imgs(i).fname);
        disp("Overrode scale of "+string(Imgs(i).fname)+" to "+string(override_scale(Imgs(i).fname)))
        % Label the scale bar validation image so we know we already took care of it
        RGB = imread("validate_scale_bar_detection/"+string(Imgs(i).fname));
        RGB = insertText(RGB,[50,100],"OVERRIDDEN","FontSize",50,"BoxColor","red");
        imwrite(RGB,("validate_scale_bar_detection/"+string(Imgs(i).fname)));
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
    writetable(table_to_store,strcat(sprintf('processed_data/%s/cache_pixel_sizes.csv',data_dir_name)));
    disp("Wrote pixel sizes to .csv as a backup.")
end
fname = {Imgs.fname};

%% Forcibly resize the images to get rid of the scale bar (!)
for ii=1:length(Imgs)
    raw = Imgs(ii).raw;
    s = size(raw);
    scale_bar_height = 150;
    Imgs(ii).raw = raw(1:s(1)-scale_bar_height,:);
    Imgs(ii).cropped = raw(1:s(1)-scale_bar_height,:);
end
imgs = {Imgs.cropped};
imshow(cell2mat(imgs(1)))


%% Do Kmeans, return segmented image
imgs_binary = agg.seg_kmeans(imgs, pixsizes);

%% Agglomerate analysis
Aggs = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname);

%% Particle scale analysis
%Aggs = pp.pcm(Aggs); % apply pair correlation method - DEPRECATED

delete('validate_primary_particle_detection/*')
Aggs = pp.hough_kook2(Aggs);

%% Manually aggregate the Hough-Kook primary particles into a spreadsheet...
output = strings(0,3);
for q = 1:length(Aggs)
    new_pp = Aggs(q).Pp_kook.dp;
    fname = Aggs(q).fname;
    where_K = strfind(fname,"K_");
    where_K = where_K(1);
    temp = str2double(fname(where_K-4:where_K-1));
    if length(new_pp)>0
        %stack_me_int = [new_pp, temp*ones(length(new_pp),1)];
        %pp_diams = [pp_diams; stack_me_int];
        dp = string(new_pp);
        temp_strs = string(temp*ones(length(new_pp),1));
        file_names = fname+strings(length(new_pp),1);
        disp(temp_strs);
        disp(file_names);
        disp(dp);
        stack_me_str = [temp_strs, dp, file_names];
        output = [output; stack_me_str];
    end
end

%% Save the data in processed folder

delete(sprintf('processed_data/%s/kmeans_imgs/*',data_dir_name))
tools.write_excel(Aggs, strcat(sprintf('processed_data/%s/kmeans_results.xlsx',data_dir_name)));
tools.imwrite_agg(Aggs, sprintf('processed_data/%s/kmeans_imgs/',data_dir_name))
% Putting a 'close all' here causes it to beach-ball and crash for some
% reason.

%% Save the Hough-Kook primaries...
writematrix(output,strcat(sprintf('processed_data/%s/all_primary_particles.xlsx',data_dir_name)));

%% Analyze aggregates and mess with imgs_binary
agg_fnames = ({Aggs.fname});            % filename associated with each aggregate
disp("Saving diagnostic images...")
delete('validate_object_detection/*') % Empty out the validation contents

% Loop over each image, write a side-by-side actual and binary. Label with
% Green or Red box. Red if no aggregates found by algorithm
for i = 1:length(pixsizes)
    if ismember(Imgs(i).fname,agg_fnames)                       % if no, then current image is NOT in agg_fnames
        loc = sprintf("processed_data/%s/kmeans_imgs/",data_dir_name);
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
        final_img_resized = imresize(final_img,[1024-scale_bar_height,1024]);
        RGB = cat(2,final_img_resized,RGB);
        RGB = insertText(RGB,[50,100],"AGGREGATES FOUND","FontSize",50,"BoxColor","green");
    else
        img = (cat(2,(cell2mat(imgs(i))),255*cell2mat(imgs_binary(i))));    % combine image and resulting binary image
        RGB = insertText(img,[50,100],"NO AGGREGATES","FontSize",50,"BoxColor","red");
    end
    imwrite(RGB,("validate_object_detection/"+string(Imgs(i).fname)))
end
disp("Saved annotated side-by-side black/white and greyscale images that contain aggregates.")

close all