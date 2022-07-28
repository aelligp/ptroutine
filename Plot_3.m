
function [outdir,path,parfile,contfile] = Plot_3(runID)

    outdir      = '../Cluster/out/';  
%     runID       = '1D_Ta4_bas';
    
    path        = strcat(outdir,runID);
    addpath(path);
    addpath('../src/');

    %Files we need to plot 
    parfile = [path ,'/', runID, '_par.mat']; % parameter file
    contfile   = [path, '/', runID, '_40.mat']; %change to last matfile
end