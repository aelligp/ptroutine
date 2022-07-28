%% Movie animation for thesis
% Use this file to animate your png files from the model. Change the ourdir
% and runID accoridingly. To save the movie file in a specific folder, add
% the path to the VideoWriter(...). 
% Script generates a folder called Movies and adds the files in a subfolder
% with the runID name.
% PAel June 2022

srcdir      = '../plot/out/';
outdir      = '/plot/out/';
runID       = '2D_Ta4_bas_N200';


path        = strcat(srcdir,runID);

%add to path
addpath(path);
if ~isfolder([outdir,'/','Movies','/',runID])
        mkdir([outdir,'/','Movies','/',runID]);
end


writerObj = VideoWriter( '/plot/out/Movies/2D_Ta4_bas_N200/vep.avi');
writerObj.FrameRate = 10;
open(writerObj);
for n = 1 : 356
  filename = sprintf('2D_Ta4_bas_N200_vep_%d.png', n);
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
end

close(writerObj);