
% PARSE_INPUTS  A common input parser for agg.* functions. 
%  For example, extracts and formats necessary information from Imgs
%  structure and outputs a common format. 
%  
%  [IMGS,PIXSIZES] = agg.parse_inputs(IMGS,[]) extracts all of the
%  necessary information from IMGS, if it is a Imgs data structure.
%  Returns the cropped images, IMGS, as a cell and the pixel sizes,
%  PIXSIZES, as an appropriately lengthed vector. 
%  
%  [IMGS,PIXSIZES] = agg.parse_inputs({IMGS},PIXSIZES) ensures the
%  appropriate format for the pixel sizes and image cell. Should result in
%  minimal changes. 
%  
%  [IMGS,PIXSIZES] = agg.parse_inputs({IMGS},[]) ignores the pixsizer
%  argument (under the assumption it is not required for the calling
%  method (e.g., seg_slider). PIXSIZES will be returned as an empty vector. 
%  
%  [IMGS,PIXSIZES] = agg.parse_inputs(IMG,PIXSIZE) acts on a single image
%  and pixel size, wrapping the IMG input in a cell and double checking the
%  value of the pixel size. As with preceeding option, PIXSIZE can be
%  empty, in which case it will be ignored and an empty vector returned. . 
%  
%  [IMG,PIXSIZE,N] = agg.parse_inputs(...) also adds output for the number 
%  of images in the set, N. This is a simple calculation of N =
%  length(IMGS). 
%  
%  AUTHOR: Timothy Sipkens, 2021-01-29

function [imgs, pixsizes, n] = parse_inputs(imgs, pixsizes)

%== OPTION 1: An Imgs structure is provided ==============================%
if isstruct(imgs)
    Imgs = imgs;
    imgs = {Imgs.cropped};  % convert input images to a cell array
    pixsizes = [Imgs.pixsize];
    
    n = length(imgs);  % also return the number of images in the set
    
    return;  % done, so return
end



%== OPTION 2: Cropped images are provided directly =======================%
%   No action on imgs, unless input is not a cell, then...
if ~iscell(imgs)
    imgs = {imgs};  % convert to a cell for future handling
end
n = length(imgs);  % also return the number of images in the set


% Handle PIXSIZES variable.
if ~exist('pixsizes', 'var'); pixsizes = []; end % make sure pixsizes exists

% See if there are any entries missing on non-sensical.
% Return a warning about missing pixel sizes and replace with ones.
if or(any(isnan(pixsizes)), any(pixsizes==0))
    warning(['Some pixel sizes were missing and assigned a value of 1 nm/px. ', ...
        'This may cause errors in classification and should be checked.']);
    disp(' ');
    
    if isempty(pixsizes); pixsizes = ones(size(imgs));
    else; pixsizes(pixsizes==0) = 1; pixsizes(isnan(pixsizes)) = 1; end
end

% Extend pxiel sizes if input is a scalar. 
% Output will be the same lengths as imgs.
if length(pixsizes) == 1
    pixsizes = pixsizes .* ones(size(imgs));
    
elseif and(~isempty(pixsizes), length(pixsizes) ~= length(imgs))
    error('IMGS and PIXSIZES size mismatch.')
    
end



end

