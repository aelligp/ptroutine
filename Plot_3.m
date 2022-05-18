
function [outdir,RunID,path,parfile,contfile] = Plot_3()

    outdir      = '../out/';  
    RunID       = '1D_anh_rhy_TT4_Ta8_070_wall_4H2O';
    
    path        = strcat(outdir,RunID);
    addpath(path);
    addpath('../src/');

    %Files we need to plot 
    parfile = [path ,'/', RunID, '_par.mat']; % parameter file
    contfile   = [path, '/', RunID, '_cont.mat']; %change to last matfile
end