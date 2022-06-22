%% Movie animation for thesis
% Use this file to animate your png files from the model. Change the ourdir
% and runID accoridingly. To save the movie file in a specific folder, add
% the path to the VideoWriter(...). 
% Script generates a folder called Movies and adds the files in a subfolder
% with the runID name.
% PAel June 2022

outdir      = '/plot/out/';
runID       = '2D_Ta8_interm_N200';

path        = strcat(outdir,runID);

%add to path
addpath(path);

writerObj = VideoWriter( '/plot/out/Movies/2D_Ta8_interm_N200/Sgr.avi');
writerObj.FrameRate = 5;
open(writerObj);
for n = 1 : 54
  filename = sprintf('2D_Ta8_interm_N200_sgr_%d.png', n);
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
end
if ~isfolder([outdir,'/','Movies','/',runID])
        mkdir([outdir,'/','Movies','/',runID]);
end
close(writerObj);