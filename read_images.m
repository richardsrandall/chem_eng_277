Imgs = tools.load_imgs('stanford');
imgs = {Imgs.cropped}; % copy variables locally
pixsizes = [Imgs.pixsize]; % pixel size for each image
fname = {Imgs.fname};