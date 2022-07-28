
function [outdir,path,parfile,contfile] = Plot_8(runID)

    outdir      = '../Cluster/out/';  
%     RunID       = '1D_Ta8_rhy';
    
    path        = strcat(outdir,runID);
    addpath(path);
    addpath('../src/');

    %Files we need to plot 
    parfile = [path ,'/', runID, '_par.mat']; % parameter file
    contfile   = [path, '/', runID, '_80.mat']; %change to last matfile
end