% Post-Process atems results

% Goal
% - Load processed results 
% - For each temperature (5-6 total), 
%   - Extract agglomerate data: L, W, aspect ratio, num_px, da, area, Rg
%                               perim, circularity, zbar_opt, dp_pcm
%   - FILTER: (many metrics)
%   - CALCULATE: geometric avg and standard deviation of dp 
%   - PLOT: area, dp and Rg histograms
%   - FIT: lognormal distirbutions
%   - Report the above


% VERSIONS
% V2 (Feb 17 2024)
%   - runs on the 143 best images, but produces some really wonky PSDFs
%   - going to try a new version of post processing the results


% Notes on HISTOGRAMS
% Can normalize by count, pdf, or probability (rel frequency)
%   - pdf normalization = area under curve is 1 
% On a linear x-axis, bins are uniform, and small particles may dominate.
% On a logarithmic x-axis, bins cover orders of magnitude and better represent small + large particles together

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

clear all; close all; clc
%% User inptus

% For analysis
exclude_1875K           = true;    % Set to true to exclude 1875K data

% For plotting
load_additional_data    = false;    % Set to true to load manual dataset


%% Run - define range of dp values for plotting

% Define range for plotting
x_range     = linspace(3, 50, 100);     % Uniform x range dp for all plots
log_dp      = log(x_range);             % Log of particle diameters for plotting  
dlog_dp     = diff(log_dp);             % Bin width in log-space
dlog_dp     = [dlog_dp, dlog_dp(end)]'; % Ensure same length (and column)




%% Load data
data_Dir_name = 'cleaned_pyrolysis_data_1-5x_threshold';
full_path = strcat(sprintf('processed/%s/kmeans/process_results.xlsx',data_Dir_name));
data = readtable(full_path);

% Extract temperature from 'fname'
% Assuming 'data' is your table and 'fname' is a cell array of strings
T5s_all = zeros(height(data), 1);
P5s_all = zeros(height(data), 1);
for i = 1:height(data)
    fname_str = data.fname{i}; % Convert from cell to string
    T5_str = fname_str(9:12);
    P5_str = fname_str(15:17);
    T5s_all(i) = str2double(T5_str);
    P5s_all(i) = str2double(P5_str);
end
% Assign these values to the data table so we can sort by T5
data.Temperature = T5s_all;
data.Pressure    = P5s_all;

% Extract the relevant columns
dp = data.dp_pcm; % Primary particle diameter in nm
Rg = data.Rg;     % Radius of gyration

% Initialize storage for temperature-wise analysis
T5s_unique = unique(T5s_all);

% Option to exclude 1875K
if exclude_1875K
    T5s_unique(T5s_unique == 1875) = [];
end

% Prepare figures for overlay later
figure_combined_dp = figure; hold on;
figure_combined_Rg = figure; hold on;

% Loop over each temperature
results = struct();
for i = 1:length(T5s_unique)
    temp = T5s_unique(i);
    temp_data = data(data.Temperature == temp, :);

    % Apply any custom filters (e.g., dp limits or filename exclusions)
    % Example: temp_data = temp_data(temp_data.dp_pcm > 5 & temp_data.dp_pcm < 100, :);

    % Extract dp and Rg for this temperature
    dp_temp = temp_data.dp_pcm;
    Rg_temp = temp_data.Rg;

    % Remove NaNs (sometimes image analysis outputs them)
    dp_temp = dp_temp(~isnan(dp_temp));
    Rg_temp = Rg_temp(~isnan(Rg_temp));

    % Calculate geometric mean and standard deviation for dp
    geo_mean_dp = exp(mean(log(dp_temp)));
    geo_std_dp = exp(std(log(dp_temp)));

    % Lognormal fit for dp
    [mu_dp, sigma_dp] = lognfit(dp_temp);

    % Lognormal fit for Rg
    [mu_Rg, sigma_Rg] = lognfit(Rg_temp);

    % Store results
    results(i).Temperature = temp;
    results(i).geo_mean_dp = geo_mean_dp;
    results(i).geo_std_dp = geo_std_dp;
    results(i).mu_dp = mu_dp;
    results(i).sigma_dp = sigma_dp;
    results(i).mu_Rg = mu_Rg;
    results(i).sigma_Rg = sigma_Rg;

    % 

    % Plot histograms for dp (Linear + Lognormal)
    figure;
    subplot(1, 2, 1);
    histogram(dp_temp, 'Normalization', 'pdf');         % area integrates to 1
    xlabel('dp (nm)');
    ylabel('Probability Density');
    title(['dp Histogram at ', num2str(temp), 'K']);
    subplot(1, 2, 2);
    edges = logspace(log10(min(dp_temp)), log10(max(dp_temp)), 20);
    histogram(dp_temp, edges, 'Normalization', 'pdf');
    set(gca, 'XScale', 'log');
    xlabel('dp (nm)');
    ylabel('Probability Density');
    title(['dp Log-Scale at ', num2str(temp), 'K']);

    % Plot histograms for Rg (Linear + Lognormal)
    figure;
    subplot(1, 2, 1);
    histogram(Rg_temp, 'Normalization', 'pdf');
    xlabel('Rg');
    ylabel('Probability Density');
    title(['Rg Histogram at ', num2str(temp), 'K']);
    subplot(1, 2, 2);
    edges = logspace(log10(min(Rg_temp)), log10(max(Rg_temp)), 20);
    histogram(Rg_temp, edges, 'Normalization', 'pdf');
    set(gca, 'XScale', 'log');
    xlabel('Rg');
    ylabel('Probability Density');
    title(['Rg Log-Scale at ', num2str(temp), 'K']);

    % Plot Particle Size Distribution Function (PSDF) for dp
    % fd : 	The Freedman-Diaconis rule is less sensitive to outliers in the
    % data, might be more suitable for data with heavy-tailed distributions. 
    % It uses a bin width of 2*iqr(X(:))*numel(X)^(-1/3).
    % [counts, edges] = histcounts(dp_temp, 'BinMethod', 'fd');
    % bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

    % Try different edges 
    nvals = 20;
    edges = logspace(log10(min(dp_temp(dp_temp > 0))), log10(max(dp_temp)), nvals);
    [counts, edges] = histcounts(dp_temp, edges);
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

    % Convert to dN/dlog(dp)/Ntot
    bin_widths_log = log10(edges(2:end)) - log10(edges(1:end-1));
    PSDF = counts ./ bin_widths_log / sum(counts);

    figure(figure_combined_dp);
    plot(bin_centers, PSDF, 'DisplayName', [num2str(temp), 'K']);
    set(gca, 'XScale', 'log');
    xlabel('dp (nm)');
    ylabel('dN/dlog(dp)/Ntot');
    title('PSDF - dp for All Temperatures');
    legend show;

    % Repeat for Rg PSDF
    [counts_Rg, edges_Rg] = histcounts(Rg_temp, 'BinMethod', 'fd');
    bin_centers_Rg = (edges_Rg(1:end-1) + edges_Rg(2:end)) / 2;
    bin_widths_log_Rg = log10(edges_Rg(2:end)) - log10(edges_Rg(1:end-1));
    PSDF_Rg = counts_Rg ./ bin_widths_log_Rg / sum(counts_Rg);

    figure(figure_combined_Rg);
    plot(bin_centers_Rg, PSDF_Rg, 'DisplayName', [num2str(temp), 'K']);
    set(gca, 'XScale', 'log');
    xlabel('Rg');
    ylabel('dN/dlog(Rg)/Ntot');
    title('PSDF - Rg for All Temperatures');
    legend show;
end

%% Optional: Load additional manual data and overlay it
if load_additional_data
    manual_data = readtable('manual_dataset.xlsx'); % Adjust filename

    manual_temps = unique(manual_data.Temperature);
    for i = 1:length(manual_temps)
        temp = manual_temps(i);
        dp_manual = manual_data.dp_pcm(manual_data.Temperature == temp);

        [counts, edges] = histcounts(dp_manual, 'BinMethod', 'fd');
        bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
        bin_widths_log = log10(edges(2:end)) - log10(edges(1:end-1));
        PSDF_manual = counts ./ bin_widths_log / sum(counts);

        figure(figure_combined_dp);
        plot(bin_centers, PSDF_manual, '--', 'DisplayName', ['Manual ', num2str(temp), 'K']);
        legend show;
    end
end