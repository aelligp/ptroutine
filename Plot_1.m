
function [outdir,RunID,path,parfile,contfile] = Plot_1()

    outdir      = '../out/0D_simulations/';  
    RunID       = '0D_anh_Tau_T4<<Tau_a8_bas049_wall_4H2O';
    
    path        = strcat(outdir,RunID);
    addpath(path);
%     addpath('../src/');

    %Files we need to plot 
    parfile = [path ,'/', RunID, '_par.mat']; % parameter file
    contfile   = [path, '/', RunID, '_cont.mat']; %change to continuum matfile

end