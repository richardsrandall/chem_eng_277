
---------- introduction ------------

This is a readme for Richard and Gibson's ChemE 277 final project. Since we're just turning this in for a class and not releasing or publishing it, we're not adding a license for our own modifications. The MIT license for the original atems package is included in 'README_atems_original.md" and the "LICENSE" file.

Our project and the atems package are implemented in MATLAB. Richard has been running MATLAB R2020b on an M1 Mac and Gibson has been running MATLAB R2024b on a PC, so any versions in between and either OS will probably work fine. To run our code, you also need to install some MATLAB toolboxes. The atems documentation identifies the following MATLAB toolboxes as dependencies: the curve fitting toolbox, the financial toolbox, the image toolbox, the optimization toolbox, and the computer vision toolbox. We are not sure that we actually used the financial toolbox. 

---------- basic procedure to run the code ------------

The general flow of data through the code is: 
1) Data begins as .tif image files in a folder in raw_tem_data
2) The 'chem277_process.m' script identifies the image scales, finds aggregates, and finds primary particles, then writes the results into the 'processed_data' folder. Three types of diagnostic images are written in the three folders named 'validate_...'. 
3) Using the results written in 'processed_data', as well as user inputs like lists of rejected images that also live in 'processed_data', a postprocessing script is run to generate histograms. 'chem277_postprocess_validation.m' is used to generate plots that are overlaid with the 'ground truth' histogram that results from manual particle measurements in imageJ. 'chem277_postprocess_all_data.m' just generates histograms without ground truths, and so can be used with any data set, not just the ones for which we have ground-truth measurements. The histograms are stored in the 'postprocessed_data' folder.

It is important to note that 'raw_tem_data', 'processed_data', and 'postprocessed_data' all need to have folders with the same name in which data from a given set of image is stored. The way we've set it up, each of these folders has three subfolders: 'all_pyrolysis_data', 'validation_data', and 'small_validation_data'. This way, the results from processing one batch of data won't be overwritten if we decide to process a different batch (without any need to manually copy-paste data to preseve it). However, to save memory, the three 'validation_...' image folders don't have subfolders; they just contain the diagnostic images from whatever batch of data was processed most recently using 'chem277_process.m'. 

In 'chem277_process.m', 'chem277_postprocess_all_data.m', and 'chem277_postprocess_validation.m', you can set the batch of data you're working on by setting the value of the 'data_dir_name' variable to 'small_validation_data', 'validation_data', or 'all_pyrolysis_data'. Note that you can't run 'chem277_postprocess_validation.m' on 'all_pyrolysis_data'. Then, you should just need to hit 'run' for each script to do its job. The diagnostic image folders are already populated, so you can also just look at them to see representative samples of what they look like after processing a batch of images.

The most important outputs of the particle detection routine are stored in subfolders of 'processed_data' in the files 'cache_pixel_sizes.xlsx', 'all_primary_particles.xlsx', and 'kmeans_results.xlsx'. However, looking at the diagnostic images and final outputs is a much easier way to tell how the primary particle and aggregate detection routines are behaving.

---------- some more tools ------------

Two utility scripts are also included. 'empty_output_folders.m' just cleans out the diagnostic image ('validation_...') folders. We never really needed to use this -- putting those folders in the gitignore was sufficient to speed up our pushed to Github -- but it's there if you need it. The other script is 'manual_image_checker.m', which we used to make the list of images whose results are ignored due to errors in aggregate detection. You'll need to fill in the path of a save file, then scrunch the MATLAB window into the side/bottom of your screen so that it doesn't obstruct the default location where images appear. Each time an image flashes up, type 'y' or 'n' and then return in the MATLAB console to approve or reject an image. At the end, a list of the file names of images you rejected will be stored in the location that you specified.

---------- notes on postprocessing ------------

There are a couple of settings in the postprocessing scripts that you can mess with, though we've already got them set to values that produce pretty good results. We originally collected TEM images of particles from pyrolysis at 6 different temperatures: 1875K, 1984K, 2132K, 2203K, 2247K, and 2322K. However, since the carbon at that temperature does not form distinct particles, the data at 1875K aren't really usable. Additionally, after applying manual rejections and filtering for the less problematic image scales, 1984K and 2322K had too few acceptable particles/aggregates detected to make good histograms. That's why we've left it with 2132K, 2203K, and 2247K as the default temperatures that get histograms.

You can change the criteria for accepting images based on their scale  by changing the line where 'pixsize_excluded_imgs' is set, down in the guts of the postprocessing scripts. Currently, only images with pixel sizes of more than 1nm (i.e., fairly zoomed-out images) are kept, since more zoomed-in images tend to misidentify TEM noise as tiny primary particles.

You can also mess with 'dp_min', the smallest particle diameter (in nm) that it allows to be detected. This value is set at the top of the postprocessing scripts. It defaults to 9.