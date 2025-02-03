
% SEG_SLIDER_ORIG  A method to use a slider GUI to segment images.
%  
%  This represents a largely manual technique originally developed by 
%  Dastanpour et al. (2016). The function enacts a GUI-based 
%  method with a slider for adaptive, semi-automatic thresholding 
%  of the image (*adaptive* in that small sections of the 
%  image can be cropped and assigned individually-selected 
%  thresholds). This is done in several steps: 
%  
%  1. The user is first prompted to **crop** around a single
%   aggregate. This allows the user to zoom in on the image 
%  for the remaining steps. 
%  
%  2. The user uses a lasso tool to draw around the aggregate. 
%  The excluded regions are used for **background subtraction** 
%  in the cropped region of the image. 
%  
%  3. Gaussian blurring is then performed on the image to reduce 
%  the noise in the output binary image. Then, the user is prompted 
%  with a <strong>slider</strong> that is used to manually adjust the level 
%  of the threshold in the cropped region of the image. Very dark regions 
%  denotes sections above the selected threshold. The optimal threshold 
%  normally occurs when parts of the surrounding region are also considered 
%  above the threshold (i.e., show up as black) but do not connect to the 
%  main aggregate (see Step 4 in the figure below depicting the window 
%  progression). 
%  
%  4. The user is prompted to **select** which regions are aggregate, 
%  ignoring any white regions that may be above the threshold but 
%  are not part of the aggregate. 
%  
%  5. Finally, the user will be prompted to **check** 
%  whether the segmentation was successful by referring to an image with  
%  the resultant binary overlaid on the original image.
%  
%  The progression of windows is shown in docs/manual_progress.png. 
%  
%  <strong>NOTE</strong>: The mostly manual nature of this approach will resulting 
%  in variability and subjectiveness between users. However, the human 
%  input often greatly improves the quality of the segmentations and, 
%  while more time-intensive, can act as a reference in considering the 
%  appropriateness of the other segmentation methods. 
%  
%  <strong>NOTE</strong>: Several sub-functions are included within this file. 
%  
%  <strong>NOTE</strong>: We note that this code saw important bug updates since 
%  the original code by Dastanpour et al. (2016). This includes fixing 
%  how the original code would repeatedly apply a Gaussian filter 
%  every time the user interacted with the slider in the GUI (which may 
%  cause some backward compatibility issues), a reduction in the use of 
%  global variables, memory savings, and other performance improvements. 
%  
%  ------------------------------------------------------------------------
%  
%  [IMGS_BINARY] = agg.seg_slider_orig(IMGS) applies the slider method to the
%  images speified in IMGS, an Imgs data structure. IMGS_BINARY is a binary
%  mask resulting from the procedure.
%  
%  [IMGS_BINARY] = agg.seg_slider_orig({IMGS}) applies the slider method to the
%  images speified in IMGS, a cell of cropped images.
%  
%  [IMG_BINARY] = agg.seg_slider_orig(IMG) applies the slider method to the
%  single, copped image given by IMG.
%  
%  [IMG_BINARY] = agg.seg_slider_orig(...,IMGS_BINARY) adds the options for a
%  pre-classified binary image, which allows for modification of an
%  existing binary mask. 
%  
%  [IMG_BINARY] = agg.seg_slider_orig(...,IMGS_BINARY,F_CROP) adds a flag
%  specifying whether of not to run the image crop tool to focus in on
%  small part of the image. By default, F_CROP = 1 and does use the crop
%  tools. Note, IMGS_BINARY can be empty.
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2019-10-11 (modified)
%   Ramin Dastanpour, 2016-02 (based on, github.com/rdastanpour)
%   Developed at the University of British Columbia

function [imgs_binary] = seg_slider_orig(imgs, imgs_binary, f_crop) 


%== Parse input ==========================================================%
% Use common inputs parser, noting that pixel size is not
% used as an input to this function.
[imgs, ~, n] = agg.parse_inputs(imgs, []);

if ~exist('f_crop','var'); f_crop = []; end
if isempty(f_crop); f_crop = 1; end

% Initial cellular array of image binaries, if image 
% binary is not provided. Empty if not provided.
if ~exist('imgs_binary','var'); imgs_binary = []; end
if isempty(imgs_binary); imgs_binary{n} = []; end
if ~iscell(imgs_binary); imgs_binary = {imgs_binary}; end

% Check to make sure there is an equal number of binary 
% and imags inputs, if supplied. If not supplied, imgs_binary{n} = [];
% ensures appropriate length.
if length(imgs_binary) ~= n
    error('Size mismatch between imgs and imgs_binary');
end
%=========================================================================%


f0 = figure; % initialize a new figure for UI
f0.WindowState = 'maximized'; % maximize the figure window


for kk=1:n
    
    if n>1 % if more than one image, output text indicating image number
        disp(['[== IMAGE ',num2str(kk), ' OF ', ...
            num2str(length(imgs)), ' ============================]']);
    end
    
    img = imgs{kk}; % image for this iteration
    
    % intialize binary for this iteration
    if isempty(imgs_binary{kk}); img_binary0 = zeros(size(img));
        
    % use partially binarized image, covert to logical, if necessary
    else img_binary0 = logical(imgs_binary{kk}); end
    
    
%== CORE FUNCTION ========================================================%
    
    moreaggs = 1;
    while moreaggs==1
        img_binary = []; % declare nested variable (allows GUI feedback)

        %== STEP 1: Crop image ===========================================%
        if f_crop
            figure(f0); clf; % intialize plot of image with initial binary overlaid
            tools.imshow_binary(img, img_binary0);
            
            uiwait(msgbox('Please crop the image around missing region.'));
            [~, rect] = imcrop; % user crops image
            rect = round(rect);
            
            if or(rect(3)==0, rect(4)==0)
                uiwait(msgbox('Rectangle undefined, try again.'));
                continue;
            end
            
            % Get pixel values from original image instead of binary overlay.
            inds1 = rect(2):(rect(2) + rect(4) - 1);
            inds2 = rect(1):(rect(1) + rect(3) - 1);
            img_cropped = img(inds1, inds2);
            
        else
            img_cropped = img; % originally bypassed in Kook code
            rect = [];
        end
        
        
        %== STEP 2: Image refinment ======================================%
        %-- Step 1-1: Apply Lasso tool -----------------------------------%
        img_binary = lasso_fnc(img_cropped);

        %-- Step 1-2: Refining background brightness ---------------------%
        img_refined = background_fnc(img_binary, img_cropped);



        %== STEP 3: Thresholding =========================================%
        figure(f0); clf;

        hax = axes('Units','Pixels');
        tools.imshow(img_refined);
        title('Applying threshold...');
        
        %-- Add a slider uicontrol ---------------------------------------%
        level = graythresh(img_refined); % Otsu thresholding
        hst = uicontrol('Style', 'slider',...
            'Min', 0.4-level, 'Max', 1-level, 'Value', 0.6-level,...
            'Position', [20 390 150 15],...
            'Callback', {@thresh_slider,hax,img_refined,img_binary});
        get(hst,'value');
        
        % add a text uicontrol to label the slider
        htxt = uicontrol('Style','text',...
            'Position', [20 370 150 15],...
            'String','Threshold level');

        %-- Pause program while user changes the threshold level ---------%
        h = uicontrol('Position',[20 335 150 25],'String','Finished',...
            'Callback','uiresume(gcbf)');
        message = sprintf(['Move the slider to the right or left to change ', ...
            'threshold level\nWhen finished, click on ''Finished'' continue.']);
        uiwait(msgbox(message));
        disp('Waiting for the user to apply the threshold to the image.');
        uiwait(gcf);
        
        delete(h); % delete ui components used for thresholding
        delete(hst);
        delete(htxt);
        disp('Thresholding is applied.');
        
        
        %== STEP 4: Select particles and format output ===================%
        uiwait(msgbox(['Please selects (left click) particles satisfactorily ', ...
            'detected; and press enter or double click.']));
        img_binary = bwselect(img_binary,8);
        
        
        %-- Check if result is satisfactory ------------------------------%
        figure(f0); clf;
        tools.imshow_binary(img_cropped, img_binary);
        choice2 = questdlg(['Satisfied with aggregate detection? ', ...
            'If not, try drawing an edge around the aggregate manually...'], ...
            'Agg detection','Yes','No','Yes');
        if strcmp(choice2,'No'); continue; end
            % if 'No', then go back to crop without incorporating binarys


        %-- Subsitute rectangle back into orignal image ------------------%
        inds1 = rect(2):(rect(2) + rect(4) - 1);
        inds2 = rect(1):(rect(1) + rect(3) - 1);
        img_binary0(inds1,inds2) = ...
            or(img_binary0(inds1,inds2), img_binary);
        
        
        %-- Save a temporary copy of image on each iteration -------------%
        %   This is done in case of an error during segmentation.
        if ~isfolder('temp'); mkdir('temp'); end
        imwrite(img_binary0, ['temp/slider_',num2str(kk),'.tif']);
        

        %-- Query user ---------------------------------------------------%
        figure(f0); clf;
        tools.imshow_binary(img, img_binary0);

        choice = questdlg('Are there any particles poorly or not detected?',...
            'Missing particles','Yes','No','No');
        if strcmp(choice,'Yes')
            moreaggs=1;
        else
            moreaggs=0;
        end
        
        
        
    end
%=========================================================================%
    
    imgs_binary{kk} = img_binary0; % copy binary to cellular array
    
    if n>1 % if more than one image, output text
        disp('[== Complete. ==============================]');
        disp(' ');
        disp(' ');
    end
end

close(f0); % close figure used during segmentation
delete('temp/slider_*.tif'); % delete temporary files upon success

% If a single image, cell arrays are unnecessary.
% Extract and just output images. 
if length(imgs)==1
    imgs_binary = imgs_binary{1};
end






%===================%
%== SUB-FUNCTIONS ==%
%===================%

%=========================================================================%
%== BACKGROUND_FNC =======================================================%
% Smooths out background using curve fitting
% Author:    Ramin Dastanpour, Steven N. Rogak, Last updated in Feb. 2016
% Modified:  Timothy Sipkens, 2019-07-16
%
% Notes:
%   This function smoothens background brightness, specially on the edges of
%   the image where intensity (brightness) has a curved planar distribution.
%   This improves thresholding in the following steps of image processing
function img_refined = background_fnc(img_binary,img_cropped)

nagg = nnz(img_binary); % pixels within the aggregate
ntot = numel(img_cropped); % pixels within the whole cropped image 
nbg = ntot-nagg; % pixels in the backgound of the aggregate


%-- Computing average background intensity -------------------------------%
burned_img = img_cropped;
burned_img(img_binary) = 0;
mean_bg =  mean(mean(burned_img))*ntot/nbg;


%-- Replace aggregate pixels' with intensity from the background ---------%
img_bg = img_cropped;
img_bg(img_binary) = mean_bg;


%-- Fit a curved surface into Filled_img data ----------------------------%
[x_d,y_d] = meshgrid(1:size(img_bg,2),1:size(img_bg,1));
xdata = {x_d,y_d};
fun = @(c,xdata) c(1).*xdata{1}.^2+c(2).*xdata{2}.^2+c(3).*xdata{1}.*xdata{2}+...
    c(4).*xdata{1}+c(5).*xdata{2}+c(6);

c_start = [0 0 0 0 0 mean_bg];
options = optimset('MaxFunEvals',1000);
options = optimset(options,'MaxIter',1000); 
[c] = lsqcurvefit(fun,c_start,xdata,double(img_bg),[],[],options);


%-- Build the fitted surface ---------------------------------------------%
img_bg_fit = zeros(size(img_bg));
for ii = 1:size(img_bg,1)
    for jj = 1:size(img_bg,2)
        img_bg_fit(ii,jj) = ...
            c(1)*ii^2+c(2)*jj^2+c(3)*ii*jj+c(4)*ii+c(5)*jj+c(6);
    end
end

%-- Refine Cropped_img, using fitted surface -----------------------------%
img_refined = mean_bg+double(img_cropped)-img_bg_fit;
img_refined = uint8(img_refined);

end




%=========================================================================%
%== LASSO_FNC ============================================================%
%   This function allows user to draw an approximate boundary around the
%   particle. Region of interest (ROI))
%   Author:   Ramin Dastanpour & Steven N. Rogak
%             Developed at the University of British Columbia
%   Modified: Yiling Kang, 2018-05-10
%             Timothy Sipkens
% 
%   Yiling Kang updates/QOL changes:
%   - Asks user if their lasso selection is correct before applying the
%     data
%   - QOL - User will not have to restart program if they mess up the lasso
function img_mask = lasso_fnc(img_in)

%-- Displaying cropped image ---------------------------------------------%
clf;
tools.imshow(img_in);
title('Applying lasso tool...');

%-- Freehand drawing. Selecting region of interest (ROI) -----------------%
drawing_correct = 0; % this variable is used to check if the user drew the lasso correctly
while drawing_correct==0 
    message = sprintf('Please draw an approximate boundary around the aggregate.\nLeft click and hold to begin drawing.\nLift mouse button to finish');
    uiwait(msgbox(message));
    
    %-- Draw around the aggregate (tool depends on Matlab version) -------%
    [~,vdate] = version; vdate = year(vdate); % get Matlab version year
    if vdate>2018
        fh = drawfreehand(); % alternate for MATLAB 2019+
    else
        fh = imfreehand(); % for older versions of Matlab
    end
    finished_check = questdlg('Are you satisfied with your drawing?','Lasso Complete?','Yes','No','No');
    
    if strcmp(finished_check, 'Yes') % if user is happy with their selection...
        drawing_correct = 1;
        
    else % if user would like to redo their selection...
        delete(fh);
    end     
end


%-- Create a binary masked image from the ROI object ---------------------%
img_mask = fh.createMask();


end




%=========================================================================%
%== THRESH_SLIDER ========================================================%
%   Thresholding the image using a slider GUI
%   Function to be used with the pair correlation method (PCM) package
%   Author:   Ramin Dastanpour & Steven N. Rogak, 2016-02
%             Developed at the University of British Columbia
%   Modified: Timothy Sipkens
function thresh_slider(hObj,~,hax,img_in,img_binary0)

%-- Average filter -------------------------------------------------------%
hav = fspecial('average');
img_mod = imfilter(img_in, hav);


%-- Median ---------------------------------------------------------------%
% Examines a neighborhood of WxW matrix, takes and makes the centre of that
% matrix the median of the original neighborhood
W = 5;
for ii=1:6 % repeatedly apply median filter, which will result in artifacts on edges
    img_mod = medfilt2(img_mod, [W W], 'symmetric');
end
% NOTE: The loop is intended to imitate the increasing amounts of 
% median filter that is applied each time the slider button is clicked
% in the original code. This was a bug in the previous software. 


%-- Binary image via threshold value -------------------------------------%
adj = get(hObj, 'Value');
level = graythresh(img_mod); % default threshold is Otsu
level = level + adj;
img_binary1 = imbinarize(img_mod, level);


% Binary image via dilation, which
% reduces initial noise and fills initial gaps.
img_binary2 = imdilate(~img_binary1, strel('square',1));


%-- Refining binary image. Before refinig, thresholding causes some ------%
%   Errors, initiating from edges, grows towards the aggregate. In
%   this section, external boundary, or background region, is utilized to
%   eliminate detection errors in the background region.
img_binary3 = 0 .* img_binary2;
img_binary3(img_binary0) = img_binary2(img_binary0);
img_binary = logical(img_binary3);

% Impose the binary on the cropped image for display to user.
% This will adjust as the threshold is updated.
img_toshow = double(img_mod) .* (double(~img_binary)+1) ./ 2;

axes(hax);
tools.imshow(img_toshow);
title('Applying threshold...');

end

end
