
function [outdir,path,parfile,contfile] = Plot_7(runID)

    outdir      = '../Cluster/out/';  
%     RunID       = '1D_Ta8_rhy';
    
    path        = strcat(outdir,runID);
    addpath(path);
    addpath('../src/');

    %Files we need to plot 
    parfile = [path ,'/', RunID, '_par.mat']; % parameter file
    contfile   = [path, '/', RunID, '_80.mat']; %change to last matfile
end