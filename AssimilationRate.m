%% Assimilation assessment simulations for 2D (and possibly 1D)
% use this script to analyse the assimilation of SI, C, V and assimilation
% in general
% PAel July 2022

close all; clear; clc;

% source directory and loading of variables

for i = 1:1
    sourcedir   = '../Cluster/';
    outdir      = 'out/';
    n = 80; % specify which time step you would want to plot, important compare each sim @ same time step

    switch i 
        case 1
            runID       = '1D_Ta4_bas';
            txt        = 'Ta4 basaltic';
        case 2
            runID       = '2D_Ta8_bas_N200';
            txt        = 'Ta8 basaltic';
        case 3
            runID       = '2D_Ta4_bas_6wt_N200';
            txt        = 'Ta4 6wt basaltic';
        case 4
            runID       = '2D_Ta8_bas_6wt_N200';
            txt        = 'Ta8 6wt basaltic';
        case 5
            runID       = '2D_Ta4_interm_N200';
            txt        = 'Ta4 intermediate';
        case 6
            runID       = '2D_Ta8_interm_N200';
            txt        = 'Ta8 intermediate';
        case 7
            runID       = '2D_Ta4_interm_6wt_N200';
            txt        = 'Ta4 6wt intermediate';
        case 8
            runID       = '2D_Ta8_interm_6wt_N200';
            txt        = 'Ta8 6wt intermediate';
        case 9
            runID       = '2D_Ta4_rhy_N200';
            txt        = 'Ta4 rhyolitic';
        case 10
            runID       = '2D_Ta8_rhy_N200';
            txt        = 'Ta8 rhyolitic';
        case 11
            runID       = '2D_Ta4_rhy_6wt_N200';
            txt        = 'Ta4 6wt rhyolitic';  
%         case 12
%             runID       = '2D_Ta8_rhy_6wt_N200';
%             txt        = 'Ta8 6wt rhyolitic';
    end

            path        = strcat(sourcedir,outdir,runID);
            src         = strcat(sourcedir,'src');

            %add to path
            addpath(path);
            addpath(src);

            %Files we need to plot
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            % ocean = load('ocean.mat'); %ocean colormap

            

            contfile = ([runID '_' num2str(n) '.mat']);

            %loading workspace
            if exist(parfile,'file'); load(parfile); end
            if exist(contfile,'file')
                load(contfile,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
            end
xbar = [2.15];
ybar = [100];


% Initiate figure
% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
LW = {'LineWidth',1};


LS = {'LineStyle','-','--','-.',':'};
%CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};

fh(1) = figure(1); 
plot(hist.time/hr, si0+hist.RaSI,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('stable isotope assimilation',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('Stable isotope [\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});

legend
fh(2) = figure(2);
plot(hist.time/hr, hist.Ra,'-',LW{:},'DisplayName',txt);  box on; hold on %axis xy tight;
title('assimilation rate',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('dimensionless',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend
fh(3) = figure(3); 
plot(hist.time/hr, c0+hist.RaC,'-',LW{:},'DisplayName',txt); box on; hold on %axis xy tight;
title('assimilation of composition',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('composition [wt \%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend
fh(4) = figure(4); 
plot(hist.time/hr, v0+hist.RaV,'-',LW{:},'DisplayName',txt);  box on; hold on %axis xy tight;
title('assimilation of volatiles',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('volatiles [wt \%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

fh(6) = figure(6); 
plot(hist.time/hr, hist.ct,'-',LW{:},'DisplayName',txt);  box on; hold on %axis xy tight;
title('assimilation of CT',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('CT [wt \%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

fh(5) = figure(5); 
yyaxis left
plot(si0+hist.RaSI, (hist.RaSI/10.6).*-100,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
%bar(xbar,ybar, 2.3)
title('stable isotope assimilation',TX{:},FS{:}); xlabel('$\delta ^18 O$ [\permil]',TX{:},FS{:}); ylabel('Assimilation SI [\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
xlim([-10, 7.5]), ylim([0,100])
yyaxis right
plot(si0+hist.RaSI, hist.time/hr, '-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('stable isotope assimilation',TX{:},FS{:}); ylabel('time [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

fh(7) = figure(7);
plot(hist.time/hr, hist.si(:,1), '-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('stable isotope assimilation',TX{:},FS{:}); ylabel('time [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend


clear 
end