%% Double Diffusive Convection Analyses
% Use this script to analyse the different factors which influence the
% layered convection in the nakhla model
% different density of T, c and x
% Also applicable to differences in vertical convection speed of different
% parameter sets

clc; clear all; close all;

for j = 1:6
    switch j
        case 1
            runID       = '2D_Ta4_bas_N200';
            txt         = '\tau_a = 4h basaltic';
            sh          = 'o';
        case 2
            runID       = '2D_Ta8_bas_N200';
            txt         = '\tau_a = 8h basaltic';
            sh          = 'o';
        case 3
            runID       = '2D_Ta4_interm_N200';
            txt         = '\tau_a = 4h intermediate';
            sh          = 'd';
        case 4
            runID       = '2D_Ta8_interm_N200';
            txt         = '\tau_a = 8h intermediate';
            sh          = 'd';
        case 5
            runID       = '2D_Ta4_rhy_N200';
            txt         = '\tau_a = 4h rhyolitic';
            sh          = 's';
        case 6
            runID       = '2D_Ta8_rhy_N200';
            txt         = '\tau_a = 8h rhyolitic';
            sh          = 'd';

        case 7 
            runID       = '2D_Ta8_interm_6wt_N200';
            txt         ='\tau_a = 8h intermediate 6wt';
    end
    outdir  = '../Cluster/out/';
    path    = strcat(outdir,runID);
    for i = 0
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


        %dT, dC, dX
        dT = T0 - Twall;
        dC = c0 - cwall;
        dX = 0.1; 
       
        % drho, eta, K_t
        drho = (rhox0+rhom0)./2; 
        eta  = (etam0 + etax0)./2;
        kT = kTm./((drho).*Cp);

        % Cx
        Cx = dx^2./eta; 

        %specific velocities
        U_diff = kT./D;
        U_settle = drho.*g0.*Cx;
        U_assml = dw./tau_a;

        U_conv_T = (aTm.*drho.*g0.*dT.*D^2)./(eta./drho); 
        U_conv_C = (gCm.*drho.*g0.*dC.*D^2)./(eta./drho); 
        U_conv_X = (gCx.*drho.*g0.*dX.*D^2)./(eta./drho); 

        sum_U_conv = U_conv_T + U_conv_C + U_conv_X;


        %Rayleigh numbers
        Ra_T = (aTm.*g0.*dT.*D^3)/(kT.*(eta./drho));
        Ra_c = (gCm.*g0.*dC.*D^3)/(kT.*(eta./drho));
        Ra_x = (gCx.*g0.*dX.*D^3)/(kT.*(eta./drho));
    
        %Buoyancy ratios
        R_rho_cT  = U_conv_C./sum_U_conv;
        R_rho_x,T = U_conv_X./sum_U_conv;
        R_assml   = U_assml./sum_U_conv;
        R_settle  = U_settle./sum_U_conv;

%         % difference in temp and compo per layer
%         dTemp = mean(mean(T))-T0;
%         dcompo = (mean(mean(c))-c0);
%         dchi = mean(mean(chi))-chi0;

% 
%         drhoT = aTm.*dTemp.*rhom0
% 
% 
%         drhoC = gCm.*(mean(mean(c))-c0).*rhom0 %compositional difference
% 
% 
%         drhox = (rhox0 -rhom0).*dchi
% 
%         drho = drhoT+drhoC+drhox
% 

        % Xtal density



        % Crystall settling
        %         R_x = (drhoX.*dx^2)./

        % thermal diffusivity
       
        %settling
        %R_x = drhx


       
        for ii = 1:3
            if ii ==1
                dx = 1e-3;
                txt1 = '1mm';
            elseif ii ==2
                dx = 0.5e-3;
                txt1 = '0.5 mm';
            else
                dx = 1.5e-3;
                txt1 = '1.5 mm';
            end

            R_x = abs((drhox.*dx^2)./(drho.*D^2));
            R_assml = abs((dw.*mean(mean((eta)))./(tau_a.*drho.*D^2)));

            R_rhotc = Ra_c/Ra_T
            R_rhotx = Ra_x/Ra_T
            R_rhosum = R_rhotc + R_rhotx;

        
            % % Initiate figure
            % % prepare for plotting
            TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
            TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
            UN = {'Units','Centimeters'};
            LW = {'LineWidth',1};
            %
            LS = {'LineStyle','-','--','-.',':'};
            copper = colormap(copper(6));

            CL = {'Color', copper(j,:)};
            OL = {'o','s', 'd'}; 
            %

            %     fh(1) = figure(1);
            %     scatter(drhoT, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
            %     %xlim([49.5, 50.5]);
            %     ylim([2, 10])
            %     title('Characteristic crystal settling ',TX{:},FS{:}); xlabel('$\rho_T$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            %     legend
            %
            %     fh(2) = figure(2);
            %     scatter(drhoC, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
            %     %xlim([49.5, 50.5]);
            %     ylim([2, 10])
            %     title('Characteristic crystal settling ',TX{:},FS{:}); xlabel('$\rho_c$',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            %     legend
            %
            %     fh(3) = figure(3);
            %     scatter(drhoX, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
            %     %xlim([49.5, 50.5]);
            %     ylim([2, 10])
            %     title('Characteristic crystal settling ',TX{:},FS{:}); xlabel('$\rho_\chi$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            %     legend
            % %
            %     fh(4) = figure(4);
            %     scatter(drhobar, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
            %     %xlim([49.5, 50.5]);
            %     ylim([2, 10])
            %     title('Characteristic crystal settling ',TX{:},FS{:}); xlabel('$\rho_\chi$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            %     legend
            % %
            % %     fh(4) = figure(4);
            % %     scatter(Ra, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); hold on;
            % %     scatter(Ra_c, tau_a./3600,'d', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); axis xy tight; box on; hold on;
            % %     ylim([2, 10])
            % %     title('Rayleigh Number of T and \bar{c}',TX{:},FS{:}); xlabel('$Rayleigh number$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            % %     legend
            %
            %     fh(5) = figure(5);
            %     scatter(R_rhotop, tau_a./3600, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','top layer'); box on; hold on;
            %     scatter(R_rhobot, tau_a./3600,'d', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','Bottom layer'); hold on;
            %
            %     ylim([2, 10])
            %     title('Rayleigh Number of T and c',TX{:},FS{:}); xlabel('$R_\rho$ ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            %     legend
            %
            %     fh(6) = figure(6);
            %     scatter(Pr,Ra_top, 'o', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T top layer'); box on; hold on;set(gca,'xscale','log',TL{:},TS{:});
            %     scatter(Pr,Ra_bot, 's', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T bottom layer'); hold on;set(gca,'yscale','log',TL{:},TS{:});
            %     scatter(Pr,Ra_ctop,'d',LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','c top layer'); hold on;
            %     scatter(Pr,Ra_cbot,'^', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','c bottom layer'); hold on;
            %     xlim([7.5e7, 12e7])
            %     title('Prandtl vs Rayleigh Number of T',TX{:},FS{:}); xlabel('Rayleigh number ',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:});
            %     legend


            fh(1) = figure(1);
            scatter(R_rhotc,tau_a./hr, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,'xscale','log',TL{:},TS{:});
            %scatter(Pr,Ra_bot, 's', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T bottom layer'); hold on;set(gca,'yscale','log',TL{:},TS{:});
            ylim([2, 10])
            title('Buoyancy ratio of $\bar{c}$ and T',TX{:},FS{:}); xlabel('$R_{\rho,c}$',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:});
%             legend

            fh(2) = figure(2);
            scatter(R_rhotx,tau_a./hr,sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
            ylim([2, 10])
            title('Buoyancy ratio of $\chi$ and T',TX{:},FS{:}); xlabel('$R_{\rho,x}$',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:});
%             legend

            fh(3) = figure(3);
            scatter(R_rhosum,R_assml,sh,LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,'yscale','log',TL{:},TS{:});
%             OL{ii}, 
            %scatter(Pr,Ra_bot, 's', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T bottom layer'); hold on;set(gca,'yscale','log',TL{:},TS{:});
            %xlim([7.5e7, 12e7])
            title('Assimilation',TX{:},FS{:}); xlabel('$\Sigma R_{\rho}$ ',TX{:},FS{:}); ylabel('$R_{assml}$',TX{:},FS{:});
%             legend

            fh(4) = figure(4);
            scatter(R_rhosum,R_x,OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt1); box on; hold on;set(gca,'yscale','log',TL{:},TS{:});

            title('Crystal settling',TX{:},FS{:}); xlabel('$\Sigma R_{\rho}$',TX{:},FS{:}); ylabel('$R_{settle}$',TX{:},FS{:});
            if j == 1
            legend('DDC','Plume', 'no special behavior')
            end
            %     fh(7) = figure(7);
            %
            %     %  xlim([2, 10])
            %     title('Prandtl vs Rayleigh Number of T',TX{:},FS{:}); xlabel('Prandtl',TX{:},FS{:}); ylabel('Rayleigh number',TX{:},FS{:}); set(gca,'yscale','log',TL{:},TS{:});
            %     legend
        end
    end
end
