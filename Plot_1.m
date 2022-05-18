
function [outdir,RunID,path,parfile,contfile] = Plot_1()

    outdir      = '../out/';  
    RunID       = '1D_anh_bas_TT4_ta8_050_wall_4H2O';
    
    path        = strcat(outdir,RunID);
    addpath(path);
    addpath('../src/');

    %Files we need to plot 
    parfile = [path ,'/', RunID, '_par.mat']; % parameter file
    contfile   = [path, '/', RunID, '_cont.mat']; %change to continuum matfile

end