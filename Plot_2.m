
function [outdir,RunID,path,parfile,contfile] = Plot_2()

    outdir      = '../out/0D_simulations/';  
    RunID       = '0D_anh_Tau_T4<<Tau_a8_interm060_wall_4H2O';
    
    path        = strcat(outdir,RunID);
    addpath(path);
    %addpath('../src/');

    parfile = [path ,'/', RunID, '_par.mat']; % parameter file
    contfile   = [path, '/', RunID, '_cont.mat']; %change to last matfile
% 
%     if exist(contfile,'file'); load(contfile); end
%     if exist(parfile,'file'); load(parfile); end
% 
%     % define Nx and Nz for output file to recognize if 0D, 1D or 2D plots
%     X         = -h/2:h:L+h/2;
%     Z         = -h/2:h:D+h/2;
%     Nx = length(X);
%     Nz = length(Z);
% 
%     %Run output file to plot figures
%     output_plots;
% %     legend(runID);
% %     hold  on
end