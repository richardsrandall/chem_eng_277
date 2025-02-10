function GC_denoise_preprocess(input_path, output_path)
    % PREPROCESS_IMAGE: Denoise and enhance borders of an input image
    %
    % Usage:
    %   preprocess_image(input_path, output_path)
    %
    % Inputs:
    %   input_path  - Path to the input noisy image file
    %   output_path - Path where the processed image will be saved
    %
    % This function applies a Gaussian filter and a median filter for denoising,
    % enhances edges using sharpening, adjusts contrast, displays results, and saves the processed image.    


    %% Load Image
    img = imread(input_path);
    
    %% Convert to Grayscale if Needed
    if size(img,3) == 3
        img = rgb2gray(img);
    end
    
    %% Denoising
    img_filtered = imgaussfilt(img, 1.5);  % Gaussian filter (adjust sigma)
    img_filtered = medfilt2(img_filtered, [3,3]); % Median filter
    
    %% Edge Enhancement
    img_edges = imsharpen(img_filtered, 'Radius', 2, 'Amount', 1.5);
    
    %% Contrast Adjustment
    img_enhanced = imadjust(img_edges);
    
    %% Display Results
    figure;
    subplot(1,3,1), imshow(img), title('Original');
    subplot(1,3,2), imshow(img_filtered), title('Denoised');
    subplot(1,3,3), imshow(img_enhanced), title('Enhanced Borders');
    
    %% Save Processed Image
    imwrite(img_enhanced, output_path);
    
    disp("DeNoise Preprocessing complete. Saved as: " + output_path);
end
