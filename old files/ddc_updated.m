%% Double Diffusive Convection Analyses
% Use this script to analyse the different factors which influence the
% layered convection in the nakhla model
% different density of T, c and x
% Also applicable to differences in vertical convection speed of different
% parameter sets

clc; clear all; close all;

for j = 1:3
    if j == 1
        %define environment and load files
        runID   = '2D_Ta4_interm_N200';
        txt ='intermediate \tau 4h';
    elseif j == 2
        runID   = '2D_Ta8_interm_N200';
        txt ='intermediate \tau 8h';
    else
        runID   = '2D_Ta8_rhy_N200';
        txt ='felsic \tau 8h';
    end
    outdir  = '../Cluster/out/';
    path    = strcat(outdir,runID);
    for i = 80
        addpath(path);



        parfile =  [path ,'/', runID, '_par.mat'];
        contfile=  [path ,'/', runID, '_' num2str(i) '.mat']; %set to desired .mat file number
        if exist(parfile,'file'); load(parfile); end
        if exist(contfile,'file'); load(contfile,'U','W','P','Pt','f','x','m', ...
                'phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf', ...
                'IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt', ...
                'dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx', ...
                'rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt', ...
                'time','step','hist','VolSrc','wf','wx');
        end
%%
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
        rhom = rhom0 .* (1 - aTm.*(T-perT) - gCm.*(cm-(perCx+perCm)/2));
        rhox = rhox0 .* (1 - aTx.*(T-perT) - gCx.*(cx-(perCx+perCm)/2));

                %% Density difference analyses

     
%         % mixture density difference rhobar = effect on temp + compo + xtals
%         drhobar = chibot.*(-rhox0.*aTx.*(Tbot-Ttop))+(mubot.*(-rhom0.*aTm.*(Tbot-Ttop))) + ...
%             (chibot.*(-rhox0.*gCx.*(cxbot-cxtop)))+(mubot.*(-rhom0.*gCm.*(cmbot-cmtop))) + ...
%             (rhoxtop -rhomtop).*(chibot-chitop)
% 
% 
%         rhobar = rhobot -rhotop % model mixture density


%Initial temp, compo, x

chi0 = 0.1308;

% difference in temp and compo per layer    
        dTemp = mean(mean(T))-T0;
        dcompo = mean(mean((cm)))-c0;
        dchi = mean(mean(chi))-chi0;
        drhoX = rhox- rhox0;
        %dX = 
        drhoT = (chi.*(-rhox0.*aTx.*(T)))+(mu.*(-rhom0.*aTm.*(T)))


        drhoC = (chi.*(-rhox0.*gCx.*(cx)))+(mu.*(-rhom0.*gCm.*(cm))) %compositional difference


        %  drhoX = (rhox -rhom).*(chi)

        drho = drhoT+drhoC+drhoX;


        % Xtal density
    


        % Crystall settling
%         R_x = (drhoX.*dx^2)./

% thermal diffusivity
        kT = kTm./(mean(mean(rho)).*Cp)

        %settling
        %R_x = drhx

 %Rayleigh numbers
        Ra_T = (aTm.*g0.*dTemp.*D^3)/(kT.*mean(mean(eta./rho)))
        Ra_c = (gCm.*g0.*dcompo.*D^3)/(kT.*mean(mean(eta./rho)))
        Ra_x = (gCx.*g0.*dchi.*D^3)/(kT.*mean(mean(eta./rho)))
      %  R_assml = (dw.*mean(mean((eta)))./(tau_a.*g0^2))
        
        R_rhotc = Ra_c/Ra_T
        R_rhotx = Ra_x/Ra_T
%         R_rhosum = R_rhoc + R_rhox;
        
% Prandtl number
%         Pr = mean(mean(eta./rho))/aTm

        % if drhobar <0 sprintf('unstable')
        % else sprintf('stable')
        % end


        % % Initiate figure
        % % prepare for plotting
        TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
        TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
        UN = {'Units','Centimeters'};
        LW = {'LineWidth',1};
        %
        LS = {'LineStyle','-','--','-.',':'};
        copper = colormap(copper(3));

        CL = {'Color', copper(j,:)};
        %

    end
    fh(1) = figure(1);
    scatter(drhoT, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
    %xlim([49.5, 50.5]);
    ylim([2, 10])
    title('Characteristic crystal settling ',TX{:},FS{:}); xlabel('$\rho_T$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    legend

    fh(2) = figure(2);
    scatter(drhoC, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
    %xlim([49.5, 50.5]);
    ylim([2, 10])
    title('Characteristic crystal settling ',TX{:},FS{:}); xlabel('$\rho_c$',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    legend

    fh(3) = figure(3);
    scatter(drhoX, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
    %xlim([49.5, 50.5]);
    ylim([2, 10])
    title('Characteristic crystal settling ',TX{:},FS{:}); xlabel('$\rho_\chi$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    legend
% 
    fh(4) = figure(4);
    scatter(drhobar, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
    %xlim([49.5, 50.5]);
    ylim([2, 10])
    title('Characteristic crystal settling ',TX{:},FS{:}); xlabel('$\rho_\chi$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    legend
% 
%     fh(4) = figure(4);
%     scatter(Ra, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); hold on;
%     scatter(Ra_c, tau_a./3600,'d', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
%     ylim([2, 10])
%     title('Rayleigh Number of T and \bar{c}',TX{:},FS{:}); xlabel('$Rayleigh number$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%     legend
    
    fh(5) = figure(5);
    scatter(R_rhotop, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','top layer'); box on; hold on;
    scatter(R_rhobot, tau_a./3600,'d', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','Bottom layer'); hold on;
    
    ylim([2, 10])
    title('Rayleigh Number of T and c',TX{:},FS{:}); xlabel('$R_\rho$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    legend
    
    fh(6) = figure(6);
    scatter(Pr,Ra_top, 'o', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T top layer'); box on; hold on;set(gca,'xscale','log',TL{:},TS{:});
    scatter(Pr,Ra_bot, 's', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T bottom layer'); hold on;set(gca,'yscale','log',TL{:},TS{:});
    scatter(Pr,Ra_ctop,'d',LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','c top layer'); hold on;
    scatter(Pr,Ra_cbot,'^', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','c bottom layer'); hold on;
    xlim([7.5e7, 12e7])
    title('Prandtl vs Rayleigh Number of T',TX{:},FS{:}); xlabel('Rayleigh number ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); 
    legend 


    fh(7) = figure(7);
    scatter(R_rhotc,tau_a./hr, 'o', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T top layer'); box on; hold on;set(gca,'xscale','log',TL{:},TS{:});
    %scatter(Pr,Ra_bot, 's', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T bottom layer'); hold on;set(gca,'yscale','log',TL{:},TS{:});
    %xlim([7.5e7, 12e7])
    title('Prandtl vs Rayleigh Number of T',TX{:},FS{:}); xlabel('Rayleigh number ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); 
    legend 

    fh(8) = figure(8);
    scatter(R_rhotc,tau_a./hr,'o', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T top layer'); box on; hold on;set(gca,'xscale','log',TL{:},TS{:});
    
    title('Prandtl vs Rayleigh Number of T',TX{:},FS{:}); xlabel('Rayleigh number ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); 
    legend 
%     fh(7) = figure(7);
% 
%     %  xlim([2, 10])
%     title('Prandtl vs Rayleigh Number of T',TX{:},FS{:}); xlabel('Prandtl',TX{:},FS{:}); ylabel('Rayleigh number',TX{:},FS{:}); set(gca,'yscale','log',TL{:},TS{:});
%     legend

end
