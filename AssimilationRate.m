%% Assimilation assessment simulations for 2D (and possibly 1D)
% use this script to analyse the assimilation of SI, C, V and assimilation
% in general
% PAel July 2022

% source directory and loading of variables

i = 150; % specify which time step you would want to plot, important compare each sim @ same time step

for n = 1:6
    sourcedir   = '../Cluster/';
    outdir      = 'out/';

    switch n 
        case 1
            runID       = '2D_Ta4_bas_N200';
        case 2
            runID       = '2D_Ta8_bas_N200';
        case 3
            runID       = '2D_Ta4_interm_N200';
        case 4
            runID       = '2D_Ta8_interm_N200';
        case 5
            runID       = '2D_Ta4_rhy_N200';
        case 6
            runID       = '2D_Ta8_rhy_N200';
    end
            path        = strcat(sourcedir,outdir,runID);
            src         = strcat(sourcedir,'src');

            %add to path
            addpath(path);
            addpath(src);

            %Files we need to plot
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            % ocean = load('ocean.mat'); %ocean colormap

            name = runID;
            contfile = ([runID '_' num2str(i) '.mat']);

            %loading workspace
            if exist(parfile,'file'); load(parfile); end
            if exist(contfile,'file')
                load(contfile,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
            end


% Initiate figure
% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
LW = {'LineWidth',1};


LS = {'LineStyle','-','--','-.',':'};
CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95], [0.20 0.25 0.75],[0.35 0.70 0.50]};

fh(1) = figure(1); 
plot(hist.time/hr,hist.RaSI,'-',CL{[1,i]},LW{:}); axis xy tight; box on; hold on
hold on 

fh(2) = figure(2);
plot(hist.time/hr,hist.Ra,'-',CL{[1,i]},LW{:}); axis xy tight; box on; hold on
hold on 


clear
end