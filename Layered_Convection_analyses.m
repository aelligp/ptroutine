%% Double Diffusive Convection Analyses
% Use this script to analyse the different factors which influence the
% layered convection in the nakhla model
% PAEL Aug. 2022

clc; clear all; close all;

for j = 1:6
    switch j
        case 1
            runID       = '2D_Ta4_bas_N200';
            txt         = '\tau_a = 4h basaltic';
            sh          = 'o'; %define shape according to conv. behaviour
        case 2
            runID       = '2D_Ta8_bas_N200';
            txt         = '\tau_a = 8h basaltic';
            sh          = 'o'; %define shape according to conv. behaviour
        case 3
            runID       = '2D_Ta4_interm_N200';
            txt         = '\tau_a = 4h intermediate';
            sh          = 'd'; %define shape according to conv. behaviour
        case 4
            runID       = '2D_Ta8_interm_N200';
            txt         = '\tau_a = 8h intermediate';
            sh          = 'd'; %define shape according to conv. behaviour
        case 5
            runID       = '2D_Ta4_rhy_N200';
            txt         = '\tau_a = 4h rhyolitic';
            sh          = 's'; %define shape according to conv. behaviour
        case 6
            runID       = '2D_Ta8_rhy_N200';
            txt         = '\tau_a = 8h rhyolitic';
            sh          = 'd'; %define shape according to conv. behaviour
    end
    %define direction and path 
    outdir  = '../Cluster/out/';
    path    = strcat(outdir,runID);
    %Define desired output .mat file
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

        % crystal diameter
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

            % K_t, Cx
            kT   = kTm./(rhom0.*Cp); %thermal diffusitivity
            Cx   = dx^2./etam0;

            %density differences
            drho_T   = aTm.*dT.*rhom0;
            drho_c   = gCm.*dC.*rhom0;
            drho_X   =      dX.*(rhox0-rhom0);

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

            %Speed ratio's for settling and assimilation
            R_settle  = ((rhox0-rhom0).*dx^2)./(sum_drho.*D^2);
            R_assml   = (dw.*etam0)./(sum_drho.*tau_a.*D^2);

            %Buoyancy ratio
            R_rho_cT  = Ra_c./Ra_T; %Buoyancy ratio after (Hansen&Yuen,1990)
            R_rho_xT  = Ra_x./Ra_T; %Buoyancy ratio to evaluate crystallinity


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
                scatter(R_rho_xT,R_rho_cT, sh, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
                title('Buoyancy ratio of $\chi$ and $\bar{c}$',TX{:},FS{:}); xlabel('$R_{\rho,x}$',TX{:},FS{:}); ylabel('$R_{\rho,c}$',TX{:},FS{:});
                legend
            end

            fh(3) = figure(3);
            scatter(R_assml,R_settle, OL{ii}, LW{:},'MarkerEdgeColor',copper(j,:),'MarkerFaceColor', copper(j,:),'DisplayName',txt); box on; hold on;set(gca,TL{:},TS{:});
            title('Speed ratio of assimilation and settling',TX{:},FS{:}); xlabel('$R_{assml}$',TX{:},FS{:}); ylabel('$R_{settle}$',TX{:},FS{:});

        end
    end
end
