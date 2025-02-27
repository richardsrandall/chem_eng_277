
% Quick script to facilitate manually validating aggregate data
% It'll sequentially show every image in the target directory, 
% then focus the command window so you can just type 'y' or 'n' 
% based on whether the image is acceptable (yes/no). The filenames 
% and y/n get saved in a .csv file in a location that you specify.

% Best to downsize the Matlab window before running
% Use ctrl-c to bail; the 'stop' button doesn't work for some reason.
target_dir = 'validate_object_detection/';
save_file = 'processed_data/all_pyrolysis_data/manual_validation.csv';

imgs = dir(target_dir);
output = strings(length(imgs),2);

for i = 1:length(imgs)
    fname = (imgs(i).name);
    output(i,1) = fname;
    output(i,2) = "n";
    try
        imshow(imread(strcat([target_dir,fname])));
        commandwindow;
        clc;
        result = "";
        result = input("Input: ","s");
        output(i,2) = result;
    catch
    end
end

%%

writematrix(output,save_file);
disp("Saved results.")

close all