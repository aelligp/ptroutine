%% Main plotting routine %%
% This routine enables you to plot multiple simulations together to see the
% differences. It uses the output_plots.m file to plot the different 0D/1D
% plots. IF desired you can also change the line color to change with every
% run rather than the style --> change to LC{i}
% PAel June 2022

close all; clear all; clc;

outdir      = '../Cluster/out/';
runID       = '2D_Ta4_bas_N200';

for i = 1:5
    switch i
        case  1
            %Files we need to plot
            %             outdir      = '../Cluster/out/';
            %             runID       = '2D_Ta4_bas';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_100.mat']; %change to continuum matfile
            n = 1;
            txt = 't = 100';
        case  2
            outdir      = '../Cluster/out/';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_200.mat']; %change to continuum matfile
            txt = 't = 200';
        case  3
            outdir      = '../Cluster/out/';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_300.mat']; %change to continuum matfile
            txt = 't = 300';
        case  4
            outdir      = '../Cluster/out/';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_400.mat']; %change to continuum matfile                txt = 't = 12000';
            txt = 't = 400';
        case  5
            outdir      = '../Cluster/out/';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_420.mat']; %change to continuum matfile                txt = 't = 16000';
            txt = 't = 500';
    end

    %for 0D fh1-3 and for 1D fh1-4
    %     [Nx,Nz,fh1,fh2,fh3,fh4] = output_plots_new(name,txt,n,runID,i,parfile, contfile);
    %     hold on;
    %     legend('time step 1','basaltic, time step 40','intermediate, time step 40','rhyolitic, time step 40', ...
    %          'Interpreter','latex','location','best');
    %         %'basaltic, time step 80', 'intermediate, time step 80','rhyolitic, time step 80', ...

    % outdir      = '../Cluster/out/';
    % path        = strcat(outdir,runID);
    addpath(path);


    %loading workspace
    if exist(parfile,'file'); load(parfile); end
    if exist(contfile,'file')
        load(contfile,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
    end

    addpath('../src/')
    load ocean %ocean colormap

    % define Nx and Nz for output file to recognize if 0D, 1D



    X         = -h/2:h:L+h/2;
    Z         = -h/2:h:D+h/2;
    Nx = length(X);
    Nz = length(Z);

    % minimum cutoff phase, component fractions
    TINY     =  1e-16;

    % update phase densities
    rhom = rhom0 .* (1 - aTm.*(T-perT-273.15) - gCm.*(cm-(perCx+perCm)/2));
    rhox = rhox0 .* (1 - aTx.*(T-perT-273.15) - gCx.*(cx-(perCx+perCm)/2));
    rhof = rhof0 .* (1 - aTf.*(T-perT-273.15) + bPf.*(Pt-Ptop ));

    % update effective viscosity
    etam  = etam0 .* exp(Em./(8.3145.*T)-Em./(8.3145.*(perT+273.15))) ...
        .* Fmc.^((cm-(perCx+perCm)/2)./(cphs1-cphs0)) ...
        .* Fmv.^(vm./0.01);                                          % variable melt viscosity
    etaf  = etaf0.* ones(size(f));                                             % constant fluid viscosity
    etax  = etax0.* ones(size(x));

    % get geochemical phase compositions
    itm  = it./(m + x.*KIT); itx = it./(m./KIT + x);
    ctm  = ct./(m + x.*KCT); ctx = ct./(m./KCT + x);
    ripm = rip./(m + x.*KRIP); ripx = rip./(m./KRIP + x);
    ridm = rid./(m + x.*KRID); ridx = rid./(m./KRID + x);

    %Velocity field
    Div_V  =  0.*P;  Div_rhoV = 0.*P;  Div_rhoVo = Div_rhoV;

    % update velocity divergence
    Div_V(2:end-1,2:end-1) = ddz(W(:,2:end-1),h) ...                           % get velocity divergence
        + ddx(U(2:end-1,:),h);
    Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
    Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

    %decay
    dcy_rip = rho.*rip./HLRIP.*log(2);
    dcy_rid = rho.*rid./HLRID.*log(2);


    %necessary for 1D plots
    [XX,ZZ]   = meshgrid(X,Z);
    Xfc       = (X(1:end-1)+X(2:end))./2;
    Zfc       = (Z(1:end-1)+Z(2:end))./2;
    [XXu,ZZu] = meshgrid(Xfc,Z);
    [XXw,ZZw] = meshgrid(X,Zfc);
    % prepare for plotting

    TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
    TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
    UN = {'Units','Centimeters'};
    LW = {'LineWidth',1};

    copper = colormap('copper');

    LS = {'LineStyle','-','--','-.',':'};
    % colormap(ocean(4))
    %cm = colormap(ocean(size(A,1)))
    % oceancustom = ocean(4);
%     CL = {'Color', copper(n,1:3)};
    CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};


    if Nx <= 10  % create 1D plots

        fh1 = figure(1);
        subplot(1,5,i)
        plot(mean(c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title(txt,TX{:},FS{:}); set(gca,TL{:},TS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        sgtitle('$\bar{c}/(1-f)$ [wt\%]',TX{:},FS{:});
        hold on;
        %         legend(txt,'Interpreter','latex','location','best');
        fh2 = figure(2);
        subplot(1,5,i)
        plot(mean(mu (2:end-1,2:end-1),2)*100.*(mean(mu (2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title(txt,TX{:},FS{:}); set(gca,TL{:},TS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        sgtitle(['$\mu$ [vol\%]'],TX{:},FS{:})
        hold on;
        % legend('time step 1','basaltic, time step 40','intermediate, time step 40','rhyolitic, time step 40', ...
        %     'Interpreter','latex','location','best');

        fh3 = figure(3);
        subplot(1,5,i)
        plot(mean(si(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title(txt,TX{:},FS{:}); set(gca,TL{:},TS{:});ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        sgtitle('stable isotope',TX{:},FS{:});
        hold on;
        % legend('time step 1','basaltic, time step 40','intermediate, time step 40','rhyolitic, time step 40', ...
        %     'Interpreter','latex','location','best');



    else % create 2D plots
        
        % set axis and border dimensions
        axh = 6.00; axw = axh*L/D;
        ahs = 0.40; avs = 0.2;
        axb = 1.00; axt = 0.4;
        axl = 1.20; axr = 0.4;

        % initialize figures and axes
        fh1 = figure(1); colormap(ocean);
        fh = axb + 1*axh + 0*avs + axt;
        fw = axl + 5.5*axw + 2*ahs + axr;
        set(fh1,UN{:},'Position',[3 3 fw fh]);
        set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh1,'Color','w','InvertHardcopy','off');
        set(fh1,'Resize','on');
        ax(1) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+0*axh+0*avs axw axh]);

        fh2 = figure(2); colormap(ocean);
        fh = axb + 1*axh + 0*avs + axt;
        fw = axl + 5.5*axw + 2*ahs + axr;
        set(fh2,UN{:},'Position',[3 3 fw fh]);
        set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh2,'Color','w','InvertHardcopy','off');
        set(fh2,'Resize','on');
        ax(2) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+0*axh+0*avs axw axh]);

        fh3 = figure(3);  colormap(ocean);
        fh = axb + 1*axh + 0*avs + axt;
        fw = axl + 5.5*axw + 2*ahs + axr;
        set(fh3,UN{:},'Position',[3 3 fw fh]);
        set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh3,'Color','w','InvertHardcopy','off');
        set(fh3,'Resize','off');
        ax(3) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+0*axh+0*avs axw axh]);


        fh4 = figure(4);  colormap(ocean);
        fh = axb + 1*axh + 0*avs + axt;
        fw = axl + 5*axw + 2*ahs + axr;
        set(fh4,UN{:},'Position',[3 3 fw fh]);
        set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh4,'Color','w','InvertHardcopy','off');
        set(fh4,'Resize','on');
        ax(4) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+0*axh+0*avs axw axh]);

        fh5 = figure(5);  colormap(ocean);
        fh = axb + 1*axh + 0*avs + axt;
        fw = axl + 5*axw + 2*ahs + axr;
        set(fh5,UN{:},'Position',[3 3 fw fh]);
        set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh5,'Color','w','InvertHardcopy','off');
        set(fh5,'Resize','on');
        ax(5) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+0*axh+0*avs axw axh]);


        % plot temperature in Fig. 1
        figure(1);
        
        axes(ax(1));
        imagesc(X(2:end-1),Z(2:end-1),T(2:end-1,2:end-1)-273.15); axis ij equal tight; box on;   cb = colorbar;
        set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});  hold on;
        sgtitle(['$T [^\circ$C]'],TX{:},FS{:},'Color','k');
        if i == 1
            ylabel('Depth [m]',TX{:},FS{:});
        end
        if i == 5
         
           
        end
        hold on;
       % plot  and composition in Fig. 2
        figure(2);
            sgtitle(['$\bar{c}/(1-f)$ [wt\% SiO$_2$]'],TX{:},FS{:},'Color','k');
        axes(ax(2));
        imagesc(X(2:end-1),Z(2:end-1),c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)).*100); axis ij equal tight; box on;  cb = colorbar;
        set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); hold on;
        if i == 1
            ylabel('Depth [m]',TX{:},FS{:});
        end


        figure(3);
            sgtitle(['$\chi$ [vol\%]'],TX{:},FS{:},'Color','k');
        axes(ax(3));
        imagesc(X(2:end-1),Z(2:end-1),chi(2:end-1,2:end-1).*100.*(chi(2:end-1,2:end-1)>1e-9) ); axis ij equal tight; box on;  cb = colorbar;
       set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); hold on;
       if i == 1
           ylabel('Depth [m]',TX{:},FS{:});
       end


          figure(4);
            sgtitle(['stable isotope'],TX{:},FS{:},'Color','k');
        axes(ax());
        imagesc(X(2:end-1),Z(2:end-1),si(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
       set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); hold on;
       if i == 1
           ylabel('Depth [m]',TX{:},FS{:});
       end 
%         axes(ax(2));
%         imagesc(X(2:end-1),Z(2:end-1),c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)).*100); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{c}/(1-f)$ [wt\% SiO$_2$]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
%         axes(ax(3));
%         imagesc(X(2:end-1),Z(2:end-1),v(2:end-1,2:end-1).*100.*(v(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{v}$ [wt\% H$_2$O]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

%         % plot phase fractions and reaction rates in Fig. 3
%         figure(3);
%         sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
%         axes(ax(31));
%         imagesc(X(2:end-1),Z(2:end-1),chi(2:end-1,2:end-1).*100.*(chi(2:end-1,2:end-1)>1e-9) ); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
%         axes(ax(32));
%         imagesc(X(2:end-1),Z(2:end-1),phi(2:end-1,2:end-1).*100.*(phi(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\phi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
%         axes(ax(33));
%         imagesc(X(2:end-1),Z(2:end-1),Gx(2:end-1,2:end-1)./rho(2:end-1,2:end-1)*hr*100.*(chi(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_x/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
%         axes(ax(34));
%         imagesc(X(2:end-1),Z(2:end-1),Gf(2:end-1,2:end-1)./rho(2:end-1,2:end-1)*hr*100.*(phi(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_f/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
% 
%         % plot density, rheology, and segregation speeds in Fig. 4
%         figure(4);
%         sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
%         axes(ax(41));
%         imagesc(X(2:end-1),Z(2:end-1),      rho(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\rho}$ [kg/m$^3$]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
%         axes(ax(42));
%         imagesc(X(2:end-1),Z(2:end-1),log10(eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\eta}$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
%         axes(ax(43));
%         imagesc(X(2:end-1),Z(2:end-1),-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^x$ [m/hr]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
%         axes(ax(44));
%         imagesc(X(2:end-1),Z(2:end-1),-(phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^f$ [m/hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
% 
%         % plot geochemical variables in Fig. 5
%         figure(5);
%         sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
%         axes(ax(51));
%         imagesc(X(2:end-1),Z(2:end-1),it(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['incomp. trace'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
%         axes(ax(52));
%         imagesc(X(2:end-1),Z(2:end-1),ct(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['comp. trace'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
%         axes(ax(53));
%         imagesc(X(2:end-1),Z(2:end-1),si(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['stable isotope'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
%         axes(ax(54));
%         imagesc(X(2:end-1),Z(2:end-1),rip(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['radiogen. parent'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
%         axes(ax(55));
%         imagesc(X(2:end-1),Z(2:end-1),rid(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['radiogen. daughter'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
%         axes(ax(56));
%         imagesc(X(2:end-1),Z(2:end-1),(dcy_rip(2:end-1,2:end-1)-dcy_rid(2:end-1,2:end-1))./(dcy_rip(2:end-1,2:end-1)+dcy_rid(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
%         set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['radiogen. disequilibrium'],TX{:},FS{:}); set(gca,'YTickLabel',[]);


    hold on
        %     legend('basaltic', 'intermediate', 'rhyolitic','Interpreter','latex','location','best');
        %

    end
end
% hold on
% legend('basaltic', 'intermediate', 'rhyolitic','Interpreter','latex','location','best');
% legend(txt,'Interpreter','latex','location','best')



% %save outputs
% name = 'Ta8';
% opdir = 'out/';
% if ~isfolder([opdir,'/',name])
%     mkdir([opdir,'/',name]);
% end
%      if Nx <= 10 && Nz <= 10  % print 0D plots
%         name_save = [outdir,'/',runID,'/',name,'_tch',num2str(i)];
%         print(fh1,name_save,'-dpng','-r300','-opengl');
%         name_save = [outdir,'/',runID,'/',name,'_aux',num2str(i)];
%         print(fh2,name_save,'-dpng','-r300','-opengl');
%     elseif Nx <= 10  % create 1D plots
% name_save = [opdir,'/',name,'/',name,'_compo'];
% print(fh1,name_save,'-dpng','-r600','-opengl');
% name_save = [opdir,'/',name,'/',name,'_melt'];
% print(fh2,name_save,'-dpng','-r600','-opengl');
% name_save = [opdir,'/',name,'/',name,'_SI'];
% print(fh3,name_save,'-dpng','-r600','-opengl');
%         name_save = [opdir,'/',name,'/',name,'_elements'];
%         print(fh4,name_save,'-depsc','-r600','-opengl');
%     else
%         name_save = [outdir,'/',runID,'/',name,'_vep_',num2str(i)];
%         print(fh1,name_save,'-dpng','-r600','-opengl');
%         name_save = [outdir,'/',runID,'/',name,'_tch_',num2str(i)];
%         print(fh2,name_save,'-dpng','-r600','-opengl');
%         name_save = [outdir,'/',runID,'/',name,'_phs_',num2str(i)];
%         print(fh3,name_save,'-dpng','-r600','-opengl');
%         name_save = [outdir,'/',runID,'/',name,'_sgr_',num2str(i)];
%         print(fh4,name_save,'-dpng','-r600','-opengl');
%         name_save = [outdir,'/',runID,'/',name,'_gch',num2str(i)];
%         print(fh5,name_save,'-dpng','-r600','-opengl');
% %         name_save = [runID,'_eql',num2str(i)];
% %         print(fh7,name,'-dpng','-r300','-opengl');
%     end


