
% MANUAL Allows for manual primary particle sizing on an array of aggregates.
%  
%  ------------------------------------------------------------------------
%  
%  NOTE:  While this function can be used to run through an entire 
%  aggregate structure, it is often safer to loop through the Aggs 
%  structure in the parent code, which then saves the result of each 
%  aggregate to the workspace before moving to the next aggregate. 
%  Specifically, code akin to:
%  
%  for ii=1:length(Aggs)
%   Aggs = pp.manual(Aggs, ii);
%  end
%  
%  ------------------------------------------------------------------------
%  
%  MODIFIED:  Timothy Sipkens, 2019-07-23
%  
%   Pieces of this code are adapted from code by Ramin Dastanpour, 
%   Hugo Tjong, Arka Soewono from the University of British Columbia, 
%   Vanouver, BC, Canada. 


function [Aggs, Pp, dp] = manual(Aggs, idx)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('idx','var'); idx = []; end
if isempty(idx); idx = 1:length(Aggs); end
    % if ind was not specified, analyze all of the aggregates
%-------------------------------------------------------------------------%


disp('Performing manual analysis:');

%-- Check whether the data folder is available ---------------------------%
if exist('data','dir') ~= 7 % 7 if exist parameter is a directory
    mkdir('data') % make output folder
end

f0 = figure; % figure handle used during manual sizing
f0.WindowState = 'maximized';


%== Process image ========================================================%
tools.textbar(0);
Pp(length(idx)) = struct(); % re-initialize data structure
for ll = idx % run loop as many times as aggregates selected
    
    idx0 = [Aggs.img_id]==Aggs(ll).img_id; %  index of first aggregate in image
    idx_agg = 1:length(Aggs);
    idx_agg = idx_agg(idx0);
    a1 = idx_agg(1);
    
    pixsize = Aggs(ll).pixsize; % copy pixel size locally
    img_cropped = imcrop(Aggs(a1).image, Aggs(ll).rect);
    
    %== Step 3: Analyzing each aggregate =================================%
    f_finished = 0;
    jj = 0; % intialize particle counter
    
    figure(f0); % plot aggregate
    clf;
    imagesc(img_cropped);
    colormap gray; axis image off;
    hold on;
    
    uiwait(msgbox('Please select two points on the image that correspond to the length of the primary particle',...
        ['Process Stage: Length of primary particle ' num2str(jj)...
        '/' num2str(jj)],'help'));
    
    while f_finished == 0
        
        jj = jj+1;
        
        % prompt user to draw first line
        [x,y] = ginput(2);
        Pp(ll).length(jj,1) = pixsize*sqrt((x(2)-x(1))^2+(y(2) - y(1))^2);
        line([x(1),x(2)],[y(1),y(2)], 'linewidth', 3);
        
        % prompt user to draw second line
        [a,b] = ginput(2);
        Pp(ll).width(jj,1) = pixsize*sqrt((a(2)-a(1))^2+(b(2) - b(1))^2);
        line([a(1),a(2)],[b(1),b(2)],'Color', 'r', 'linewidth', 3);
        
        %-- Save center of the primary particle --------------------------%
        Pp(ll).centers(jj,:) = find_centers(x,y,a,b);
        Pp(ll).radii(jj,:) = (sqrt((a(2)-a(1))^2 + (b(2)-b(1))^2)+...
        	sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2)) / 4;
            % takes an average over drawn lines (given in pixels)
        Pp(ll).dp = 2 .* pixsize .* Pp(ll).radii; % particle diameter (given in nm)
        
        %-- Check if there are more primary particles --------------------%
        choice = questdlg('Do you want to analyze another primary particle ?',...
        'Continue?', 'Yes', 'No', 'Yes');
        if strcmp(choice,'Yes')
        	f_finished = 0;
        else
        	f_finished = 1;
        end
        
    end
    
    % Allow for refinement of circles by
    % using handles and prompting the user.
    Pp(ll) = tools.refine_circles(img_cropped, Pp(ll));
    
    commandwindow; % return focus to Matlab window
    
    %== Save results =====================================================%
    %   Format output and autobackup data ------------------------%
    disp(' Saving temporary data ...');
    Aggs(ll).Pp_manual = Pp(ll); % copy Pp data structure into Aggs
    Aggs(ll).dp_manual = median(Pp(ll).dp);
    Aggs(ll).dp = Aggs(ll).dp_manual;
    save(['temp',filesep,'Pp_manual.mat'],'Pp'); % backup Pp
    tools.textdone(2);
    
    disp(' Overall progress:');
    tools.textbar(0);
    tools.textbar(find(ll==idx) / length(idx));
    disp(' ');
end

close(f0); % close existing figure
delete(['temp',filesep,'Pp_manual.mat']); % delete temporary datas
dp = [Aggs.dp];

disp('Complete.');
disp(' ');

end




%== FIND_CENTERS =========================================================%
%   Computes the intersection of the two lines to get particle center.
%   Author: Yeshun (Samuel) Ma, Timothy Sipkens, 2019-07-12
%
%   Notes:
%   x,y,a,b are single column vectors of two entries denoting the x or y
%   coordinate at the points the user had clicked on the image during
%   primary particle sizing.
%
%   Utilizes derived formulae to compute the intersection of the linear
%   functions composed of the provided parameters.  Returns center: a
%   two column vector corresponding to the x and y coordinates of the
%   primary particle
function centers = find_centers(x,y,a,b)

tol = 1e-10; % used to prevent division by zero

%-- Initialize xy parameters ---------------------------------------------%
x1 = x(1);
x2 = x(2);
if x1==x2; x1 = x1+tol; end % prevents division by zero
y1 = y(1);
y2 = y(2);
if y1==y2; y1 = y1+tol; end % prevents division by zero

%-- Initialize ab parameters ---------------------------------------------%
a1 = a(1);
a2 = a(2);
if a1==a2; a1 = a1+tol; end % prevents division by zero
b1 = b(1);
b2 = b(2);
if b1==b2; b1 = b1+tol; end % prevents division by zero

%-- Compute slopes -------------------------------------------------------%
m = (y1-y2)./(x1-x2);    % For xy line (vectorized)
n = (b1-b2)./(a1-a2);    % For ab line
if any(m==n); m = m+tol; end % prevents division by zero

%-- Assign appropriate fields and return ---------------------------------%
centers(1,:) = (b1-y1-n.*a1+m.*x1)./(m-n); % x-coordinate
centers(2,:) = (centers(1,:)-a1).*n+b1; % y-coordinate

end




