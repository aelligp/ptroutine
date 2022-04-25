
function [] = Plot1

    outdir      = '../out/';  
    runID       = '0D_assml_sat_Krafla_4H20_wet_wall';
    
    path        = strcat(outdir,runID);
    addpath(path);
    addpath('../src/');

    %Files we need to plot 
    parfile = [path ,'/', runID, '_par.mat']; % parameter file
    contfile   = [path, '/', runID, '_cont.mat']; %change to lat matfile

    if exist(contfile,'file'); load(contfile); end
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