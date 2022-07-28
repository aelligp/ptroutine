%% Nakhla: Magma, mush and pluton assessment simulations for 2D
% use this script to analyse the assimilation of SI, C, V and assimilation
% in general
% PAel July 2022

close all; clear; clc;

% source directory and loading of variables

for i = 1:2
    sourcedir   = '../Cluster/';
    outdir      = 'out/';
    n = 330; % specify which time step you would want to plot, important compare each sim @ same time step

    switch i 
        case 1
            runID       = '2D_Ta4_bas_N200';
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
        case 12
            runID       = '2D_Ta8_rhy_6wt_N200';
            txt        = 'Ta8 6wt rhyolitic';
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



% calculate necessary variables for analyses
X         = -h/2:h:L+h/2;
Z         = -h/2:h:D+h/2;
Nx        = length(X);
Nz        = length(Z);
[XX,ZZ]   = meshgrid(X,Z);
Xfc       = (X(1:end-1)+X(2:end))./2;
Zfc       = (Z(1:end-1)+Z(2:end))./2;
[XXu,ZZu] = meshgrid(Xfc,Z);
[XXw,ZZw] = meshgrid(X,Zfc);
rhom = rhom0 .* (1 - aTm.*(T-perT-273.15) - gCm.*(cm-(perCx+perCm)/2));
rhox = rhox0 .* (1 - aTx.*(T-perT-273.15) - gCx.*(cx-(perCx+perCm)/2));

rhyoliticmelt = mean(chi(ZZ>0.01 & ZZ<0.3)<=0.15) & mean(C(ZZ>0.01 & ZZ<0.5)>0.68) %& hist.Cmagma>=0.70)
if rhyoliticmelt == true
sprintf('rhyolitic melt layer')
else sprintf('no rhyolitic melt layer')
end
% Initiate figure
% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
LW = {'LineWidth',1};


LS = {'LineStyle','-','--','-.',':'};
%CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};

%magma

fh(1) = figure(1); 
plot(hist.time/hr, hist.Fmagma,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('fraction of magma ($\mu>$0.55)',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('Percentage [$\%$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

fh(2) = figure(2);
plot(hist.time/hr, hist.Cmagma,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('composition of magma ($\mu>$0.55)',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('composition [wt \%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

fh(3) = figure(3); 
plot(hist.time/hr, hist.Tmagma-273.15,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on;
title('temperature of magma ($\mu>$0.55)',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('temperature [$^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

%mush

fh(4) = figure(4); 
plot(hist.time/hr, hist.Fmush,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('fraction of mush (0.15$<\mu<$0.55)',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('fraction [\%/100]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

fh(5) = figure(5); 
plot(hist.time/hr, hist.Cmush,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('composition of mush (0.15$<\mu<$0.55)',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('composition [wt \%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

fh(6) = figure(6);
plot(hist.time/hr, hist.Tmush-273.15,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('temperature of mush (0.15$<\mu<$0.55)',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('temperature [$^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

%pluton

fh(7) = figure(7); 
plot(hist.time/hr, hist.Fpluton,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('fraction of pluton ($\mu<$0.15)',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('fraction [\%/100]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

fh(8) = figure(8); 
plot(hist.time/hr, hist.Cpluton,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('composition of pluton ($\mu<$0.15)',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('composition [wt \%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

fh(9) = figure(9); 
plot(hist.time/hr, hist.Tpluton-273.15,'-',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
title('temperature of pluton ($\mu<$0.15)',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('temperature [$^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
legend

clear 
end