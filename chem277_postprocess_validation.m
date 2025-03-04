% Post-Process atems results -- with validation from manual data

% This code reads the atems-processed data, reads the corresponding manually-processed data, 
% and plots/saves histograms that let us compare them.
% It only works on the 'validation data' folder, since those images are the
% ones for which Gibson manually recorded the primary particle sizes.

clear all; close all; clc
%% (1) User inputs

% For analysis
filterBy_dpMin          = true;    % Remove all entries below dp_min
dp_min                  =  9;      % nm
temperatures_to_show = [2132, 2203, 2247];
% 1984K, 2322K, and 1875K have too little usable data after manual
% filtering, so we're excluding them here. You can add 1984 to the list to
% see what some messy data look like.

% For Binning
% Manual Binning ---------------------- (other options 'fd' or 'scotts')
bin_method              = 'manual';  % 'manual' 'fd' 'scotts' binning - see notes
n_edges                 = 20;
edge_min                = dp_min - 1;
edge_max                = 40;
preset_edges            = linspace(edge_min, edge_max, n_edges); % goes into histcounts() to return counts and edges

% For plotting
show_kept_images = false; %Quickly shows all images whose particles/aggregates are being analyzed, to make sure exclusions worked
save_histograms = true;
weight_by_agg_area = false;
overlap = true;
plots_to_keep_open = []; % One can specify plots by the figure number - e.g., only keep Figs 1 and 2 open.
open_all_plots = true; % Overrides plots_to_keep_open
save_plots_mode = false; % Closes all figures at the end if all you wanted to do was save the plots

% File locations
data_dir_name = 'validation_data';
save_subfolder='default';
save_dir_name = strcat([data_dir_name,'/',save_subfolder,'/']);
folderregex = sprintf('postprocessed_data/%s/*',save_dir_name);
delete(folderregex)

%% (2) Load all Data. Initialize Variables.

% LOAD MANUAL (by hand) DATA ----------------------------------------------
% Data shape: Each column is a different temperature, rows are dp vals
% Date shape: toal of 5 temperatures (1984, 2132, 2203, 2247, 2322K) (1875K has no data)
filename         = '5pCH4_0pH2_midP_T5_dPfromImageJ.csv'; 
full_path_manual = strcat(sprintf('manual_measurements/%s',filename));
data_hand        = readmatrix(full_path_manual);
% Setup by Hand variables
T5s_hand = data_hand(1,2:end);            % Temperatures are the headers - remove 1875 K
data_hand = data_hand(2:end,2:end);       % Remove column 1 because 1875K had no values 
num_conditions = size(data_hand, 2);           % (5 temperatures) Number of columns (each corresponding to a different temperature experiment)
colors = lines(num_conditions);           % Define colors for each temperature - MATLAB's 'lines' colormap for distinct colors
T5s_hand;

% LOAD ATEMS DATA ---------------------------------------------------
% Data Shape: table, 1st col are filnames containing T5b data
% Every row is a different segmented agglomerate with unique average dp
% (one dp estimate per agglomerate)
full_path_atems = strcat(sprintf('processed_data/%s/kmeans_results.xlsx',data_dir_name));
data_atems = readtable(full_path_atems);
% Setup atems variables
T5s_all_atems = zeros(height(data_atems), 1);   % Extract temperature from 'fname'
P5s_all_atems = zeros(height(data_atems), 1);
for i = 1:height(data_atems)
    fname_str = data_atems.fname{i};    
    T5_str = fname_str(9:12);
    P5_str = fname_str(15:17);
    T5s_all_atems(i) = str2double(T5_str);    
    P5s_all_atems(i) = str2double(P5_str);
    %disp(T5_str)
end
data_atems.Temperature = T5s_all_atems;       % Assign these values to the data_atems table so we can sort by T5
data_atems.Pressure    = P5s_all_atems;
T5s_unique_atems = unique(T5s_all_atems);
exclude_1875K           = true;    % MUST BE TRUE FOR NOW. Set to true to exclude 1875K data, 
if exclude_1875K                             % Option to exclude 1875K
    T5s_unique_atems = T5s_unique_atems(T5s_unique_atems  ~= 1875);
    data_atems = data_atems(data_atems.Temperature ~= 1875, :);
end
T5s_unique_atems;
original_data_atems = data_atems;

% Define range of expected dp values for plotting
x_range     = linspace(2, 50, 200);         % Uniform x range dp for all plots
log_xrange      = log(x_range);             % Log of particle diameters for plotting  
dlog_dp     = diff(log_xrange);             % Bin width in log-space
dlog_dp     = [dlog_dp, dlog_dp(end)]';     % Ensure same length (and column)

% Initialize results storage
% By Hand Data (5 temperatures) --------------------------------------
geo_means_hand      = zeros(1, num_conditions);
geo_stds_hand       = zeros(1, num_conditions);
pdf_logNorm_hand       = zeros(length(x_range), num_conditions);     % prob density function fit
dp_vals_hand        = nan(length(data_hand(:,1)), num_conditions);  % Store dp's
N_pp_tot_hand       = zeros(1,num_conditions);                      % total dp measurments per T5
% By Hand Data (5 temperatures) --------------------------------------
geo_means_atems     = zeros(1, 5);
geo_stds_atems      = zeros(1, 5);
pdf_logNorm_atems      = zeros(length(x_range), 5);      % prob density function fit
dp_vals_atems       = nan(height(data_atems(:,1)), 5);  % for  dp's
N_agg_tot_atems     = zeros(1,5);                       %  total # of aggs id'd per T5

%%

% LOAD PRIMARY PARTICLE DATA
full_path_atems = strcat(sprintf('processed_data/%s/all_primary_particles.xlsx',data_dir_name));
data_atems = readtable(full_path_atems);
data_atems.Properties.VariableNames = ["temp","diameter","fname"];
data_atems.dp = str2double(data_atems.diameter);
data_atems.area = data_atems.dp*0+1; % Normalizing this will do nothing
data_atems.Temperature = str2double(data_atems.temp);

% Load the manual override data
full_path_manual_exclusions = sprintf("processed_data/%s/manual_validation.csv",data_dir_name);
opts = detectImportOptions(full_path_manual_exclusions,'Delimiter',',');
manual_exclusions = readtable(full_path_manual_exclusions,opts);
manual_exclusions.Properties.VariableNames = ["fname","status"];
manually_excluded_imgs = manual_exclusions(manual_exclusions.status ~= "y", :).fname;

% Filter by pixsize
pixsize_excluded_imgs = original_data_atems(original_data_atems.pixsize>1.0, :).fname;

%Apply exclusions to the primary particle data and print some
%diagnostics...
fprintf("Images with aggregates: %d\n",length(unique(data_atems.fname)));
data_atems = data_atems(~ismember(data_atems.fname, manually_excluded_imgs),:);
fprintf("Images after manual exclusions: %d\n",length(unique(data_atems.fname)));
data_atems = data_atems(~ismember(data_atems.fname, pixsize_excluded_imgs),:);
fprintf("Images after pixel size exclusions: %d\n",length(unique(data_atems.fname)));

%Apply exclusions to the aggregate data
original_data_atems = original_data_atems(~ismember(original_data_atems.fname, manually_excluded_imgs),:);

% Filter based on image scale and/or manual overrides
%data_atems = data_atems(data_atems.Temperature ~= 1875, :);
kept_fnames = unique(data_atems.fname);
for i=1:length(kept_fnames)
    fn = string(kept_fnames(i));
    full_path_pic = sprintf("processed_data/%s/kmeans_imgs/%s",data_dir_name,fn);
    if show_kept_images
        imshow(imread(full_path_pic));
    end
end
pause(0.5)
close all


%% Loop over each temperature, compute and plot things (histograms, log-normal fits)

close all
% Loop over each temperature
for i = 1:num_conditions
    % By Hand data ------------------------------------------
    T5_hand = T5s_hand(i);
    if ~ismember(T5_hand, temperatures_to_show)
        continue
    end
    temp_data_hand = data_hand(:, i);               % Extract column of dp data
    dp_temp_hand = temp_data_hand(~isnan(temp_data_hand)); 
    N_pp_tot_hand(i) = length(dp_temp_hand);        % number of pp measured
    % ATEMS Data -------------------------------------------
    T5_atems = T5s_unique_atems(i);
    temp_data_atems = data_atems(data_atems.Temperature == T5_atems, :);
    if (length(temp_data_atems.Temperature))==0 % Too few data to plot this one!
        fprintf("Skipping plots of Temp = %dK because there are too few PP's/aggregates.\n",T5_atems);
        continue
    else
        fprintf("Plotting Temp = %dK...\n",T5_atems);
    end
    
    % FILTERING ATEMS DATA
    % (e.g., dp limits or filename exclusions)
    % FILTER 1: pp less than XX  nm are most definitely noise
    if filterBy_dpMin
        dp_min = 7; % nm
        temp_data_atems = temp_data_atems(temp_data_atems.dp > dp_min, :);
        dp_temp_atems = temp_data_atems.dp;
        N_agg_tot_atems(i) = length(dp_temp_atems);     % number of aggs Identified by atems for this T5
    else
        dp_temp_atems = temp_data_atems.dp;
        N_agg_tot_atems(i) = length(dp_temp_atems);   
    end
    
    % (a) Geometric Mean & SD
    geo_means_hand(i) = geomean(dp_temp_hand);          % does same thing as exp(mean(log(dp)))
    geo_stds_hand(i) = exp(std(log(dp_temp_hand)));
    geo_means_atems(i) = exp(mean(log(dp_temp_atems)));
    geo_stds_atems(i) = exp(std(log(dp_temp_atems)));
    
    % (b) Log-Normal Distribution Fitting (to plot vs histogram of data)
    % Let X = dp values for a give T5, mu = mean of (natural) logarithm of dp
    % Ln(dp) should be normally distributed if dp is lognormal
    % BY HAND DATA ----------------------------------------------------
    dp_log_hand  = log(dp_temp_hand); mu_hand = mean(dp_log_hand);
    pd_hand      = fitdist(dp_temp_hand, 'Lognormal');         % RETURNS MU, SIGMA params for lognormal dist 
    mu_fit_hand  = pd_hand.mu;
    pdf_logNorm_hand(:,i) = pdf(pd_hand, x_range);             % produce pdf over range of dp vals prespecified
    % ATEMS DATA -------------------------------------------------------
    dp_log_atems  = log(dp_temp_atems); mu_atems = mean(dp_log_atems);
    pd_atems      = fitdist(dp_temp_atems, 'Lognormal');         % RETURNS MU, SIGMA params for lognormal dist 
    mu_fit_atems  = pd_atems.mu;
    pdf_logNorm_atems(:,i) = pdf(pd_atems, x_range);             % produce pdf over range of dp vals prespecified

    % (c) Plotting Individual histograms and lognormal fit
    % Notes:
    % - should control the bin size, for consistency across temperature
    % Data by Hand ------------------------------------------------------
    figure("Name",sprintf("PRIMARY PARTICLES: %dK",T5_hand));
    if ~overlap
        subplot(2,1,1);
    end
    [counts, edges] = histcounts(dp_temp_hand, preset_edges);
    %figure("Name",sprintf("Hist.Hand: %dK",T5_hand));
    histogram(dp_temp_hand, edges,'FaceColor',"#0072BD",'Normalization', 'pdf','FaceAlpha',0.6); % Histogram of data
    hold on;
    plot(x_range, pdf_logNorm_hand(:,i), 'b-', 'LineWidth', 2);              % Plot fitted log-normal curve  
    %title(['Manual: ', num2str(T5_hand),' K. N_pp: ', num2str(N_pp_tot_hand(i))]);
    xlabel('dp [nm]');
    ylabel('Probability Density');
    if ~overlap
        lg_hand = legend('Manual Data', 'Log-Normal Hand Fit');
    end
    % Data ATEMS (orange) ---------------------------------------------------------
    % use same bin size/edges
    [counts, edges] = histcounts(dp_temp_atems, preset_edges);
    if ~overlap
        hold off;
        subplot(2,1,2);
        hold on;
    end
    %figure("Name",sprintf("Hist.ATEMS: %dK",T5_atems));
    which_hist_data = dp_temp_atems;
    histogram(which_hist_data, edges,'FaceColor',	"#D95319",'Normalization', 'pdf','FaceAlpha',0.6); % Histogram of data
    text(0.02,0.98,...
        sprintf('Number of ATEMS  primary particles: %d\nNumber of manual primary particles: %d',length(which_hist_data), length(dp_temp_hand)),...
        'Units','Normalized','HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    plot(x_range, pdf_logNorm_atems(:,i), 'r-', 'LineWidth', 2);              % Plot fitted log-normal curve  
    %title(['ATEMS: ', num2str(T5_atems),' K. N_agg: ', num2str(N_agg_tot_atems(i))]);
    xlabel('Primary Particle Diameter [nm]');
    ylabel('Probability Density');
    if ~overlap
        lg_atems = legend('ATEMS Data', 'Log-Normal ATEMS Fit');
    else
        lg = legend('Manual Data', 'Log-Normal Hand Fit','ATEMS Data', 'Log-Normal ATEMS Fit');
    end
    hold off;
    
    figure("Name",sprintf("AGGREGATES: %dK",T5_hand));
    agg_diams = (original_data_atems(original_data_atems.Temperature==T5_hand,:).da);
    log_agg_diams  = log(agg_diams); mu_agg_diams = mean(log_agg_diams);
    fit_agg_diams      = fitdist(agg_diams, 'Lognormal');
    histogram(agg_diams,'FaceColor',"#CF9FFF",'Normalization', 'pdf','FaceAlpha',0.6);
    hold on;
    text(0.02,0.98,...
        sprintf('Number of ATEMS  aggregates: %d',length(agg_diams)),...
        'Units','Normalized','HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    xvs = linspace(10,max(agg_diams),200);
    plot(xvs, pdf(fit_agg_diams, xvs), 'm-', 'LineWidth', 2);
    legend('Aggregate Size Data', 'Log-Normal Fit');
    xlabel('Aggregate Diameter [nm]');
    ylabel('Probability Density');
    hold off;

end

% Close all but a certain figure to avoid manually x'ing out of ones you
% don't want. Intended for debugging. Must disable to save figures.
all_open_figures = findobj(0, 'type', 'figure');
if ~save_plots_mode
    if ~open_all_plots
        delete(setdiff(all_open_figures, plots_to_keep_open)); 
    end
end

%% Save figures

if save_histograms
    % Specify the folder where the figures will be saved
    figoutputFolder = sprintf('postprocessed_data/%s/',save_dir_name);
    delete(folderregex)
    if ~exist(figoutputFolder, 'dir')
        mkdir(figoutputFolder);
    end
    % Get all open figure handles
    figHandles = findall(0, 'Type', 'figure');
    % Save each figure as PNG or SVG
    for i = 1:length(figHandles)
        figHandle = figHandles(i);
        figName   =  figHandle.Name;   % Access the figure name
    
        % Sanitize the name for use as a file name
        % Remove invalid characters for file names (like :, /, etc.)
        sanitizedFigName = regexprep(figName, '[:\/\\\*\?\"<>\|]', '_');
        
        % Generate file paths
        pngFile = fullfile(figoutputFolder, [sanitizedFigName '.png']);
        svgFile = fullfile(figoutputFolder, [sanitizedFigName '.svg']);
        
        % Save the figure with trimmed whitespace
        try
            exportgraphics(figHandle, pngFile, 'Resolution', 300); % High-res PNG
        catch
        end
        % saveas(figHandle, svgFile);
        % exportgraphics(figHandle, svgFile, 'ContentType', 'vector'); % Scalable SVG

    end
    disp('All figures saved successfully.');
    
end

if save_plots_mode
    if ~open_all_plots
        close all
    end
end


