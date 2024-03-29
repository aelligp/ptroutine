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
            %dT, dC, dX
            dT = T0 - Twall; % temperature difference
            dC = cwall - c0; % wall rock composition - initial magma compo
            dX = 0.1; % approximate crystallinity difference

            % drho, K_t, Cx
            drho = rhox0-rhom0;
            kT   = kTm./(rhom0.*Cp);
            Cx   = dx^2./etam0;
            
            %density differences
            drho_T   = aTm.*dT.*rhom0;
            drho_c   = gCm.*dC.*rhom0;
            drho_X   =      dX.*drho;

            sum_drho = drho_T + drho_c + drho_X;
            
            %specific velocities
            U_diff   = kT./D;
            U_settle = (rhox0-rhom0).*g0.*Cx;
            U_assml  = dw./tau_a;
                        

            U_conv_T = (aTm.*dT.*rhom0.*g0.*D^2)    ./etam0;
            U_conv_C = (gCm.*dC.*rhom0.*g0.*D^2)    ./etam0;
            U_conv_X = (     dX.*(rhox0-rhom0).*D^2)./etam0;

            sum_U_conv = U_conv_T + U_conv_C + U_conv_X;


            %Rayleigh numbers
            Ra_T = (aTm.*dT.*rhom0.*g0.*D^3)/(kT.*etam0);
            Ra_c = (gCm.*dC.*rhom0.*g0.*D^3)/(kT.*etam0);
            Ra_x = (dX.*(rhox0-rhom0).*g0.*D^3)/(kT.*etam0);

            R_settle  = ((rhox0-rhom0).*dx^2)./(sum_drho.*D^2);
            R_assml   = (dw.*etam0)./(sum_drho.*tau_a.*D^2);
            

            R_rho_cT  = Ra_c./Ra_T;
            R_rho_xT  = Ra_x./Ra_T;


            %             R_x = abs((drhox.*dx^2)./(drho.*D^2));
            %             R_assml = abs((dw.*mean(mean((eta)))./(tau_a.*drho.*D^2)));
            %
            %             R_rhotc = Ra_c/Ra_T
            %             R_rhotx = Ra_x/Ra_T
            %             R_rhosum = R_rhotc + R_rhotx;


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

            if ii ==1
            fh(1) = figure(1);
            scatter(Ra_xT,R_rho_cT, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
            title('Rayleigh number of T and assimilation speed',TX{:},FS{:}); xlabel('$Ra_T$',TX{:},FS{:}); ylabel('$R_{assml}$',TX{:},FS{:});
            legend
            fh(2) = figure(2);
            scatter(R_assml,R_settle, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
            title('Rayleigh number of $\bar{c}$ and assimilation speed',TX{:},FS{:}); xlabel('$Ra_T$',TX{:},FS{:}); ylabel('$R_{assml}$',TX{:},FS{:});
            fh(3) = figure(3);
            end

%             if ii == 1
%             %Assimilation
%             fh(1) = figure(1);
%             scatter(Ra_T,R_assml, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Rayleigh number of T and assimilation speed',TX{:},FS{:}); xlabel('$Ra_T$',TX{:},FS{:}); ylabel('$R_{assml}$',TX{:},FS{:});
%             legend
%             fh(2) = figure(2);
%             scatter(Ra_c,R_assml, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Rayleigh number of $\bar{c}$ and assimilation speed',TX{:},FS{:}); xlabel('$Ra_T$',TX{:},FS{:}); ylabel('$R_{assml}$',TX{:},FS{:});
%             fh(3) = figure(3);
%             scatter(Ra_x,R_assml, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Rayleigh number of $\chi$ and assimilation speed',TX{:},FS{:}); xlabel('$Ra_T$',TX{:},FS{:}); ylabel('$R_{assml}$',TX{:},FS{:});
%             
%                 %Buoyancy ratio
%             fh(4) = figure(4);
%             scatter(R_rho_cT,R_assml, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Buoyancy ratio of $\bar{c}$ and assimilation speed',TX{:},FS{:}); xlabel('$Ra_\rho,c$',TX{:},FS{:}); ylabel('$R_{assml}$',TX{:},FS{:});
%             
%             fh(5) = figure(5);
%             scatter(R_rho_xT,R_assml, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Buoyancy ratio of $\chi$ and assimilation speed',TX{:},FS{:}); xlabel('$R_\rho,x$',TX{:},FS{:}); ylabel('$R_{assml}$',TX{:},FS{:});
% 
%             end
% 
%             % Settling
%             fh(6) = figure(6);
%             scatter(Ra_T,R_settle, OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Rayleigh number of T and settling speed',TX{:},FS{:}); xlabel('$Ra_T$',TX{:},FS{:}); ylabel('$R_{settle}$',TX{:},FS{:});
%             legend
%             fh(7) = figure(7);
%             scatter(Ra_c,R_settle, OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Rayleigh number of $\bar{c}$ and settling speed',TX{:},FS{:}); xlabel('$Ra_\bar{c}$',TX{:},FS{:}); ylabel('$R_{settle}$',TX{:},FS{:});
%             fh(8) = figure(8);
%             scatter(Ra_x,R_settle, OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Rayleigh number of $\chi$ and settling speed',TX{:},FS{:}); xlabel('$Ra_\chi$',TX{:},FS{:}); ylabel('$R_{settle}$',TX{:},FS{:});
%             
%             %Buoyancy ratios
%             fh(9) = figure(9);
%             scatter(R_rho_cT,R_settle, OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Buoyancy ratio of $\bar{c}$ and settling speed',TX{:},FS{:}); xlabel('$R_\rho,c$',TX{:},FS{:}); ylabel('$R_{settle}$',TX{:},FS{:});
%             
%             fh(10) = figure(10);
%             scatter(R_rho_xT,R_settle, OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Buoyancy ratio of $\chi$ and settling speed',TX{:},FS{:}); xlabel('$Ra_\rho,x$',TX{:},FS{:}); ylabel('$R_{settle}$',TX{:},FS{:});

%             fh(7) = figure(7);
%             scatter(R_rho_cT,tau_a./hr, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Buoyancy ratio of $\bar{c}$ and T',TX{:},FS{:}); xlabel('$R_{\rho,c}$',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:});
%             fh(8) = figure(8);
%             scatter(R_rho_xT,tau_a./hr, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             title('Buoyancy ratio of $\chi$ and T',TX{:},FS{:}); xlabel('$R_{\rho,x}$',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:});

            
            %scatter(Pr,Ra_bot, 's', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T bottom layer'); hold on;set(gca,'yscale','log',TL{:},TS{:});
%             ylim([2, 10])
%             title('Buoyancy ratio of $\bar{c}$ and T',TX{:},FS{:}); xlabel('$R_{\rho,c}$',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:});
            %             legend

%             fh(4) = figure(4);
%             scatter(R_rho_cT,R_settle,OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%             fh(5) = figure(5);
%             scatter(R_rho_xT,R_settle, OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt);box on; hold on;set(gca,TL{:},TS{:});
%             fh(6) = figure(6);
%             scatter(R_rho_sum,R_settle, OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
%           
            
            %ylim([2, 10])


           % sgtitle('R_settle',TX{:},FS{:}); %xlabel('$R_{\rho,x}$',TX{:},FS{:}); ylabel('$\tau_a$ [hr]',TX{:},FS{:});
            %             legend

%             fh(3) = figure(3);
%             scatter(R_rhosum,R_assml,sh,LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,'yscale','log',TL{:},TS{:});
%             %             OL{ii},
%             %scatter(Pr,Ra_bot, 's', LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName','T bottom layer'); hold on;set(gca,'yscale','log',TL{:},TS{:});
%             %xlim([7.5e7, 12e7])
%             title('Assimilation',TX{:},FS{:}); xlabel('$\Sigma R_{\rho}$ ',TX{:},FS{:}); ylabel('$R_{assml}$',TX{:},FS{:});
%             %             legend
% 
%             fh(4) = figure(4);
%             scatter(R_rhosum,R_x,OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt1); box on; hold on;set(gca,'yscale','log',TL{:},TS{:});
% 
%             title('Crystal settling',TX{:},FS{:}); xlabel('$\Sigma R_{\rho}$',TX{:},FS{:}); ylabel('$R_{settle}$',TX{:},FS{:});
%             if j == 1
%                 legend('DDC','Plume', 'no special behavior')
%             end
%             %     fh(7) = figure(7);
%             %
            %     %  xlim([2, 10])
            %     title('Prandtl vs Rayleigh Number of T',TX{:},FS{:}); xlabel('Prandtl',TX{:},FS{:}); ylabel('Rayleigh number',TX{:},FS{:}); set(gca,'yscale','log',TL{:},TS{:});
            %     legend
        end
    end
end
