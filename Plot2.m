
function [] = Plot2

    outdir      = '../out/';  
    runID       = '0D_assml_anh_Krafla_wet_wall4H2O';
    
    path        = strcat(outdir,runID);
    addpath(path);
    addpath('../src/');
%     filePattern = fullfile(path, '*.mat');
%     dir (filePattern)

%     
    %for i = 1:
     
  %  Num = sscanf(cat(2, :), '%d');
    parfile = [path ,'/', runID, '_par.mat']; % parameter file
    lastfile   = [path, '/', runID, '_cont.mat']; %change to lat matfile

    if exist(lastfile,'file'); load(lastfile); end
    if exist(parfile,'file'); load(parfile); end

    % define Nx and Nz for output file to recognize if 0D, 1D or 2D plots
    X         = -h/2:h:L+h/2;
    Z         = -h/2:h:D+h/2;
    Nx = length(X);
    Nz = length(Z);

    %Run output file to plot figures
    run('output_plots.m');
%     legend(runID);
%     hold  on
end