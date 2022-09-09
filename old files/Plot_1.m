
function [outdir,runID,path,parfile,contfile] = Plot_1(runID)

    outdir      = '../Cluster/out/';  
    %runID       = '1D_Ta4_bas';
    
    path        = strcat(outdir,runID);
    addpath(path);
    addpath('./src/');

    %Files we need to plot 
    parfile = [path ,'/', runID, '_par.mat']; % parameter file
    contfile   = [path, '/', runID, '_0.mat']; %change to continuum matfile

end