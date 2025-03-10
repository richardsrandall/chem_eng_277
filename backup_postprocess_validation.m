% Post-Process atems results

% This code reads the atems-processed data, reads the corresponding manually-processed data, 
% and plots/saves histograms that let us compare them.

clear all; close all; clc
%% (1) User inptus

% For analysis
exclude_1875K           = true;    % MUST BE TRUE FOR NOW. Set to true to exclude 1875K data, 
filterBy_dpMin          = true;    % Remove all entries below dp_min
dp_min                  =  9;      % nm

% For Binning
% Manual Binning ---------------------- (other options 'fd' or 'scotts')
bin_method              = 'manual';  % 'manual' 'fd' 'scotts' binning - see notes
n_edges                 = 20;
edge_min                = dp_min - 1;
edge_max                = 40;
preset_edges            = linspace(edge_min, edge_max, n_edges); % goes into histcounts() to return counts and edges

% For plotting
show_kept_images = true;
save_histograms = true;
weight_by_agg_area = false;
overlap = true;
open_all_plots = false;
plots_to_keep_open = [1,2,3,4];
save_plots_mode = false;

% File locations
%data_dir_name = 'small_validation_data';
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
    % disp(T5_str)
end
data_atems.Temperature = T5s_all_atems;       % Assign these values to the data_atems table so we can sort by T5
data_atems.Pressure    = P5s_all_atems;
T5s_unique_atems = unique(T5s_all_atems);
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

%% OVERRIDE THE ATEMS DATA
if true
    full_path_atems = strcat(sprintf('processed_data/%s/all_primary_particles.xlsx',data_dir_name));
    data_atems = readtable(full_path_atems);
    data_atems.Properties.VariableNames = ["temp","diameter","fname"];
    data_atems.dp = str2double(data_atems.diameter);
    data_atems.area = data_atems.dp*0+1; % Normalizing this will do nothing
    data_atems.Temperature = str2double(data_atems.temp);
end

% Load the manual override data
full_path_manual_exclusions = sprintf("processed_data/%s/manual_validation.csv",data_dir_name);
opts = detectImportOptions(full_path_manual_exclusions,'Delimiter',',');
manual_exclusions = readtable(full_path_manual_exclusions,opts);
manual_exclusions.Properties.VariableNames = ["fname","status"];
manually_excluded_imgs = manual_exclusions(manual_exclusions.status ~= "y", :).fname;

% Filter by pixsize
pixsize_excluded_imgs = original_data_atems(original_data_atems.pixsize>1.0, :).fname;

%disp(excluded_imgs);
fprintf("Images with aggregates: %d\n",length(unique(data_atems.fname)));
data_atems = data_atems(~ismember(data_atems.fname, manually_excluded_imgs),:);
fprintf("Images after manual exclusions: %d\n",length(unique(data_atems.fname)));
data_atems = data_atems(~ismember(data_atems.fname, pixsize_excluded_imgs),:);
%original_data_atems = original_data_atems(~ismember(original_data_atems.fname, pixsize_excluded_imgs),:);
fprintf("Images after pixel size exclusions: %d\n",length(unique(data_atems.fname)));

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


%% Loop over each T5, compute things
% - a histogram
% - fit a lognormal, 
% - plot both
% - compute geometric avg dp and std dev
% - store the number of agglomerates or particles measured 

% Prepare figures for overlay later
%figure_combined_dp = figure; hold on;
%figure_combined_Rg = figure; hold on;

% Loop over each temperature
for i = 1:num_conditions
    % By Hand data ------------------------------------------
    T5_hand = T5s_hand(i);
    temp_data_hand = data_hand(:, i);               % Extract column of dp data
    dp_temp_hand = temp_data_hand(~isnan(temp_data_hand)); 
    N_pp_tot_hand(i) = length(dp_temp_hand);        % number of pp measured
    % ATEMS Data -------------------------------------------
    T5_atems = T5s_unique_atems(i);
    temp_data_atems = data_atems(data_atems.Temperature == T5_atems, :);
   
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
    % WEIGHT BY AREA
    weights = round(temp_data_atems.area/100);
    weighted_hist_data = [];
    for q=1:length(weights)
        add_me = dp_temp_atems(q)*ones(1,weights(q));
        weighted_hist_data = [add_me, weighted_hist_data];
    end
    %disp(T5_atems+", "+length(weighted_hist_data))
    
    
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
    
    % % (c) Compute the standard format lognormal distribution
    % % The first pdf_vals provides values per unit of diameter. 
    % % But, in log-space, 
    % % we want values per unit logarithm of diameter, so we divide by dlog_dp.
    % dN_dlog_dp(:,i) = pdf_vals(:,i) ./ dlog_dp; % Convert PDF to log-scale format
    % % Compute dN/dlog(dp) / Ntot
    % % Normalize approach 1 - by integral over domain
    % dN_dlog_dp_norm(:,i) = dN_dlog_dp(:,i) / trapz(dp_log, dN_dlog_dp(:,i)); % Normalize by total area under curve
    % % Normalize approach 2 - Count - based normalization
    % % NOTE: This approach does not integrate to zero . . . .
    % % can be useful if you want to compare different datasets that may have varying total particle counts.
    % N_tot = N_tots(i);
    % dN_dlog_dp_norm2(:,i) = dN_dlog_dp(:,i) / N_tot;
    % % NOTE: norm2 APPROACH does not integrate to zero . . 
    % % ORRRR NEW APPROACH
    % % COMPUTE PSDF with dN/dlog(dp)/Ntot vs dp as axes
    % % Try different edges 
    % nedges = 15;
    % edges = linspace(5, 40, nedges);
    % [counts, edges] = histcounts(dp_temp, edges);
    % bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    % % Convert to dN/dlog(dp)/Ntot
    % bin_widths_log = log10(edges(2:end)) - log10(edges(1:end-1));
    % PSDF = counts ./ bin_widths_log / sum(counts);

    % PLOTTING below this line ----------------------------------------

    % (d) Plotting Individual histograms and lognormal fit
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
    if weight_by_agg_area
        which_hist_data = weighted_hist_data;
    else
        which_hist_data = dp_temp_atems;
    end
    histogram(which_hist_data, edges,'FaceColor',	"#D95319",'Normalization', 'pdf','FaceAlpha',0.6); % Histogram of data
    plot(x_range, pdf_logNorm_atems(:,i), 'r-', 'LineWidth', 2);              % Plot fitted log-normal curve  
    %title(['ATEMS: ', num2str(T5_atems),' K. N_agg: ', num2str(N_agg_tot_atems(i))]);
    xlabel('dp [nm]');
    ylabel('Probability Density');
    if ~overlap
        lg_atems = legend('ATEMS Data', 'Log-Normal ATEMS Fit');
    else
        lg = legend('Manual Data', 'Log-Normal Hand Fit','ATEMS Data', 'Log-Normal ATEMS Fit');
    end
    hold off;
    
    figure("Name",sprintf("AGGREGATES: %dK",T5_hand));
    histogram(original_data_atems(original_data_atems.Temperature==T5_hand,:).da);

    % % PLOT PSDF with dN/dlog(dp)/Ntot vs dp as axes
    % figure(figure_combined_dp);
    % plot(bin_centers, PSDF, 'DisplayName', [num2str(temp), 'K']);
    % set(gca, 'XScale', 'log');
    % xlabel('dp (nm)');
    % ylabel('dN/dlog(dp)/Ntot');
    % title('PSDF - dp for All Temperatures');
    % legend show;
    % 
    % % Store other data I want for plotting later
    % dp_vals(1:length(temp_data_hand),i)  = temp_data_hand;      % for subplots of historgrams

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


%% HISTOGRAMS NOTES
% Questions
%   - Am I working with particle sizes spanning orders of magnitude? → Use log bins and log x-axis.
%   - Do I want a smooth lognormal curve to overlay on a log histogram? → Use logspace() for the curve.
%   - Should I log-transform my particle size data before calling pdf()? → No. Pass raw dp values into pdf().
% BINNING:
% - 3 Main options for choosing bins
%   - 1) Manual: edges = linspace(min(dp), max(dp), 15); % 15 bins on linear scale
%   - 2) Freedman-Diaconis ('fd') (automatic/t):  [counts, edges] = histcounts(dp, 'BinMethod', 'fd');
%   - 3) Scotts Rule: good if data norm dist: [counts, edges] = histcounts(dp, 'BinMethod', 'scott');
% - Data Range Narrow (~10-50 nm) 
%       - Linear binning is often fine, but lognormal fit is still reasonable.
%- Date Range wide (1 - 200 nm)
%       - On a linear x-axis, bins are uniform, and small particles may dominate.
%       - On a logarithmic x-axis, bins cover orders of magnitude and better represent small + large particles together
% NORMALIZING: Can normalize by count, pdf, or probability (rel frequency)
%   - pdf normalization = area under curve is 1 

% DATA NOTES
% BY HAND DATA
% - generally narrow range, so fits lognormal and normal distributions okay
% - using lognormal, but linear binning 

% Difference for PSDF
% - special type of histogram

%  Notes on Lognormal Distrbutions
% - x = dp, log(x) can often be described as normally distributed
% - y-axis: dN/dlog(dp)/Ntot
%   - dN = num particles in small group range
%   - dlog(dp) = Increment in the logarithm of particle diameter.
%   - Ntot = total # particles
%   - dN/dlog(dp) = particles per logarithmic interval
%   - normalized by Ntot
% - Plotting:
%   - Answers: How are particles distributed across logarithmic size intervals?" 
%       rather than: How many particles are between 1 nm and 2 nm?
%       When size spans orders of magnitude, option one is best.

% For project
% -  Added folder containing the manual analysis of individual dp vs T
% -  Fitting lognormal distribution to data uses maxlikelihood esimtator
% -  What we might need is atems to spit out all the dp measurements per
%   agglomerate, instead of just 1, can you make it do this?

% Currently:
% - loads atems and validation data
% - organizes it and loops through Temperatures, does some stats
% - plots Histogram lognormal fit for validation data and atems data
% - saves histograms to manual_measurements/histograms
% Next steps:
% - get atems to return MORE dp data (so dp for each particle in a given
%   agglomerate
% - filter atems data 
%   - by the surface area of the particle, 
%   - by the bad_images text file 
% - Compute the standard PSDF for these things,  dN/dlog(dp)/Ntot vs dp 
