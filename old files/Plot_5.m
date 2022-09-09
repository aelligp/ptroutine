
function [outdir,path,parfile,contfile] = Plot_5(runID)

    outdir      = '../Cluster/out/';  
    %runID       = '1D_Ta8_bas';
    
    path        = strcat(outdir,runID);
    addpath(path);
    addpath('../src/');

    %Files we need to plot 
    parfile = [path ,'/', runID, '_par.mat']; % parameter file
    contfile   = [path, '/', runID, '_80.mat']; %change to last matfile
end