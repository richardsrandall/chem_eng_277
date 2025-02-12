% This script copies one data folder to a different one, ignoring
% files that are listed in a manually-organized text document created
% by manually clicking through images in 'check_object_detection'.

clc

dir_to_clone = "data/cleaned_pyrolysis_data_1-5x_threshold/";
target_dir = "data/best_raw_images/";
rejection_txt_doc = "data/bad_images_v1.txt";

txt_lines = strsplit(fileread(rejection_txt_doc),"\n");

problem_tags = ["scale","macroscopic","duplicate","tiny"];
indices = zeros(1,length(problem_tags));
for i=(1:length(problem_tags))
    mask=cellfun(@any,strfind(txt_lines,problem_tags(i)));
    indices(i) = (find(mask==1,1));
end
indices = sort(indices);
% Comment out any that you don't want to exclude
tags_to_exclude = [...
                    "Have issues with the scale bar:",...
                    "Are duplicates:",...
                    "Have spurious tiny aggregates:",...
                    "Misidentify macroscopic aggregates:"...
                    ];
files_to_exclude = [];
disp("Flagging images to exclude for the following reasons (may include duplicates):")
for i=(1:length(problem_tags))
    index = indices(i);
    issue = txt_lines(index);
    %disp(" ")
    %disp(issue)
    if i==length(problem_tags)
        problem_fns = txt_lines(index+1:(length(txt_lines)-1));
    else
        problem_fns = txt_lines(index+1:(indices(i+1)-1));
    end
    %disp(problem_fns);
    files_to_exclude = cat(1,[files_to_exclude, problem_fns]);
    if ismember(issue, tags_to_exclude)
        disp(string(issue)+" "+string(length(problem_fns)));
    end
end

disp("")
disp("Beginning file copying...")

files = dir(dir_to_clone);
rejected_files = 0;
kept_files = 0;
for i=1:length(fns)
    file = fns(i);
    fn = file.name;
    exclude_me = 0;
    try
        img = imread(dir_to_clone+fn);
    catch
        continue
    end
    for j=1:length(files_to_exclude)
        if contains(files_to_exclude(j),fn)
            exclude_me = 1;
        end
    end
    if exclude_me == 0
        kept_files = kept_files+1;
        imwrite(img,target_dir+fn);
    else
        rejected_files = rejected_files+1;
    end
end

disp("...done. Processed "+string(kept_files+rejected_files)+" total files.")
disp("Copied "+string(kept_files)+" files.")
disp("Excluded "+string(rejected_files)+" files.")
disp("")