%% Main plotting routine %%
% This routine enables you to plot multiple simulations together to see the
% differences.
% PAel June 2022

close all; clear all; clc;

outdir      = '../Cluster/out/';
runID       = '2D_Ta8_interm_N200';
name        = '2D_plots_rhy';
for i = 1:6
    switch i
        case 1
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_1.mat']; %change to continuum matfile
            txt = 't = 1';
        case  2
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_60.mat']; %change to continuum matfile
            n = 1;
            txt = 't = 100';
        case  3
            outdir      = '../Cluster/out/';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_120.mat']; %change to continuum matfile
            txt = 't = 200';
        case  4
            outdir      = '../Cluster/out/';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_180.mat']; %change to continuum matfile
            txt = 't = ';
        case  5
            outdir      = '../Cluster/out/';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_240.mat']; %change to continuum matfile                txt = 't = 12000';
            txt = 't = 360';
        case  6
            outdir      = '../Cluster/out/';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_300.mat']; %change to continuum matfile                txt = 't = 16000';
            txt = 't = 500';
    end


    addpath(path);


    %loading workspace
    if exist(parfile,'file'); load(parfile); end
    if exist(contfile,'file')
        load(contfile,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
    end

    addpath('src/')
    load ocean %ocean colormap

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
    CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};



    % set axis and border dimensions
    axh = 6.00; axw = axh*L/D;
    ahs = 0.40; avs = 0.2;
    axb = 1.00; axt = 0.4;
    axl = 1.20; axr = 0.4;

    % initialize figures and axes
    fh1 = figure(1); colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh1,UN{:},'Position',[3 3 fw fh]);
    set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh1,'Color','w','InvertHardcopy','off');
    set(fh1,'Resize','on');
    if i <= 3
        ax(11) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+1*axh+1*avs axw axh]);
    else
        ax(12) = axes(UN{:},'position',[axl+(i-4)*axw+(i-4)*ahs axb+0*axh+0*avs axw axh]);
    end

    fh2 = figure(2); colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh2,UN{:},'Position',[3 3 fw fh]);
    set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh2,'Color','w','InvertHardcopy','off');
    set(fh2,'Resize','on');
    if i <= 3
        ax(21) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+1*axh+1*avs axw axh]);
    else
        ax(22) = axes(UN{:},'position',[axl+(i-4)*axw+(i-4)*ahs axb+0*axh+0*avs axw axh]);
    end

    fh3 = figure(3);  colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh3,UN{:},'Position',[3 3 fw fh]);
    set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh3,'Color','w','InvertHardcopy','off');
    set(fh3,'Resize','on');
    if i <= 3
        ax(31) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+1*axh+1*avs axw axh]);
    else
        ax(32) = axes(UN{:},'position',[axl+(i-4)*axw+(i-4)*ahs axb+0*axh+0*avs axw axh]);
    end

    fh4 = figure(4);  colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh4,UN{:},'Position',[3 3 fw fh]);
    set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh4,'Color','w','InvertHardcopy','off');
    set(fh4,'Resize','on');
    if i <= 3
        ax(41) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+1*axh+1*avs axw axh]);
    else
        ax(42) = axes(UN{:},'position',[axl+(i-4)*axw+(i-4)*ahs axb+0*axh+0*avs axw axh]);
    end

    fh5 = figure(5);  colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh5,UN{:},'Position',[3 3 fw fh]);
    set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh5,'Color','w','InvertHardcopy','off');
    set(fh5,'Resize','on');
    if i <= 3
        ax(51) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+1*axh+1*avs axw axh]);
    else
        ax(52) = axes(UN{:},'position',[axl+(i-4)*axw+(i-4)*ahs axb+0*axh+0*avs axw axh]);
    end

    fh6 = figure(6);  colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh6,UN{:},'Position',[3 3 fw fh]);
    set(fh6,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh6,'Color','w','InvertHardcopy','off');
    set(fh6,'Resize','on');
    if i <= 3
        ax(61) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+1*axh+1*avs axw axh]);
    else
        ax(62) = axes(UN{:},'position',[axl+(i-4)*axw+(i-4)*ahs axb+0*axh+0*avs axw axh]);
    end

    fh7 = figure(7);  colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh7,UN{:},'Position',[3 3 fw fh]);
    set(fh7,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh7,'Color','w','InvertHardcopy','off');
    set(fh7,'Resize','on');
    if i <= 3
        ax(71) = axes(UN{:},'position',[axl+(i-1)*axw+(i-1)*ahs axb+1*axh+1*avs axw axh]);
    else
        ax(72) = axes(UN{:},'position',[axl+(i-4)*axw+(i-4)*ahs axb+0*axh+0*avs axw axh]);
    end


    % plot temperature in Fig. 1
    figure(1);
    imagesc(X(2:end-1),Z(2:end-1),T(2:end-1,2:end-1)-273.15); axis ij equal tight; box on;   cb = colorbar;
    set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); hold on;
    sgtitle(['$T [^\circ$C]'],TX{:},FS{:},'Color','k');
    if i == 1
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 4
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 5
        xlabel('Width [m]',TX{:},FS{:});
    end


    % plot  and composition in Fig. 2
    figure(2);
    sgtitle(['$\bar{c}/(1-f)$ [wt\% SiO$_2$]'],TX{:},FS{:},'Color','k');

    imagesc(X(2:end-1),Z(2:end-1),c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)).*100); axis ij equal tight; box on;  cb = colorbar;
    set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); hold on;
    if i == 1
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 4
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 5
        xlabel('Width [m]',TX{:},FS{:});
    end


    figure(3);
    sgtitle(['$\chi$ [vol\%]'],TX{:},FS{:},'Color','k');
    imagesc(X(2:end-1),Z(2:end-1),chi(2:end-1,2:end-1).*100.*(chi(2:end-1,2:end-1)>1e-9) ); axis ij equal tight; box on;  cb = colorbar;
    set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); hold on;
    if i == 1
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 4
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 5
        xlabel('Width [m]',TX{:},FS{:});
    end


    figure(4);
    sgtitle(['stable isotope $\delta^{18}O$'],TX{:},FS{:},'Color','k');
    imagesc(X(2:end-1),Z(2:end-1),si(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); hold on;
    if i == 1
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 4
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 5
        xlabel('Width [m]',TX{:},FS{:});
    end

    figure(5);
    sgtitle(['$\phi$ [vol\%]'],TX{:},FS{:},'Color','k');
    imagesc(X(2:end-1),Z(2:end-1),phi(2:end-1,2:end-1).*100.*(phi(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); hold on;
    if i == 1
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 4
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 5
        xlabel('Width [m]',TX{:},FS{:});
    end

    figure(6);
    sgtitle(['$\bar{\eta}$ [log$_{10}$ Pas]'],TX{:},FS{:},'Color','k');
    imagesc(X(2:end-1),Z(2:end-1),log10(eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); hold on;
    if i == 1
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 4
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 5
        xlabel('Width [m]',TX{:},FS{:});
    end

    figure(7);
    sgtitle(['$\bar{\rho}$ [kg/m$^3$]'],TX{:},FS{:},'Color','k');
    imagesc(X(2:end-1),Z(2:end-1),rho(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:});set(gca,TL{:},TS{:}); title(['t = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:}); hold on;
    if i == 1
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 4
        ylabel('Depth [m]',TX{:},FS{:});
    end
    if i == 5
        xlabel('Width [m]',TX{:},FS{:});
    end

    hold on
    %     legend('basaltic', 'intermediate', 'rhyolitic','Interpreter','latex','location','best');
    %

end


%save outputs

opdir = 'out/';
if ~isfolder([opdir,'/',name])
    mkdir([opdir,'/',name]);
end
if ~isfolder([opdir,'/',name])
    mkdir([name,'/',runID]);
end

name_save = [opdir,'/',name,'/',runID,'_temp'];
print(fh1,name_save,'-dpng','-r600','-opengl');
name_save = [opdir,'/',name,'/',runID,'_compo'];
print(fh2,name_save,'-dpng','-r600','-opengl');
name_save = [opdir,'/',name,'/',runID,'_xtal'];
print(fh3,name_save,'-dpng','-r600','-opengl');
name_save = [opdir,'/',name,'/',runID,'_SI'];
print(fh4,name_save,'-dpng','-r600','-opengl');
name_save = [opdir,'/',name,'/',runID,'_phi'];
print(fh5,name_save,'-dpng','-r600','-opengl');
name_save = [opdir,'/',name,'/',runID,'_eta'];
print(fh6,name_save,'-dpng','-r600','-opengl');
name_save = [opdir,'/',name,'/',runID,'_rho'];
print(fh7,name_save,'-dpng','-r600','-opengl');


