%% Main plotting routine %%
% This routine enables you to plot multiple simulations together to see the
% differences. 
% PAel June 2022

close all; clear all; clc;

outdir      = '../Cluster/out/';
name = 'Ta8'; %specify folder name where png files are saved
opdir = 'out/'; %save png's here

for i = 1:5
    switch i
        case  1
            %Files we need to plot
            outdir      = '../Cluster/out/';
            runID       = '1D_Ta4_bas';
            path        = strcat(outdir,runID);
            parfile = [path ,'/', runID, '_par.mat']; % parameter file
            contfile   = [path, '/', runID, '_0.mat']; %change to continuum matfile
            n = 1;
            txt = 't = 1';
    end
    for j = 1:3
        if j == 1
            runID = '1D_Ta8_bas';
            n = 90;
        elseif j == 2
            runID = '1D_Ta8_interm';
            n = 150;
        else

            runID = '1D_Ta8_rhy';
            n = 250;
        end

        if i == 1
            n = 1;
        end
        switch i
            case  2
                outdir      = '../Cluster/out/';
                path        = strcat(outdir,runID);
                parfile = [path ,'/', runID, '_par.mat']; % parameter file
                contfile   = [path, '/', runID, '_20.mat']; %change to continuum matfile
                txt = 't = 4000';
            case  3
                outdir      = '../Cluster/out/';
                path        = strcat(outdir,runID);
                parfile = [path ,'/', runID, '_par.mat']; % parameter file
                contfile   = [path, '/', runID, '_40.mat']; %change to continuum matfile
                txt = 't = 8000';
            case  4
                outdir      = '../Cluster/out/';
                path        = strcat(outdir,runID);
                parfile = [path ,'/', runID, '_par.mat']; % parameter file
                contfile   = [path, '/', runID, '_60.mat']; %change to continuum matfile                txt = 't = 12000';
                txt = 't = 12000';
            case  5
                outdir      = '../Cluster/out/';
                path        = strcat(outdir,runID);
                parfile = [path ,'/', runID, '_par.mat']; % parameter file
                contfile   = [path, '/', runID, '_80.mat']; %change to continuum matfile                txt = 't = 16000';
                txt = 't = 16000';

        end

        addpath(path);


        %loading workspace
        if exist(parfile,'file'); load(parfile); end
        if exist(contfile,'file')
            load(contfile,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
        end

        addpath('src\')
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
        CL = {'Color', copper(n,1:3)};



        fh1 = figure(1);
        subplot(1,5,i)
        plot(mean(c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title(txt,TX{:},FS{:}); set(gca,TL{:},TS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        sgtitle('$\bar{c}/(1-f)$ [wt\%]',TX{:},FS{:});
        hold on;
        legend(txt,'Interpreter','latex','location','best');
        
        fh2 = figure(2);
        subplot(1,5,i)
        plot(mean(mu (2:end-1,2:end-1),2)*100.*(mean(mu (2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title(txt,TX{:},FS{:}); set(gca,TL{:},TS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        sgtitle(['$\mu$ [vol\%]'],TX{:},FS{:})
        hold on;
        
        fh3 = figure(3);
        subplot(1,5,i)
        plot(mean(si(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title(txt,TX{:},FS{:}); set(gca,TL{:},TS{:});ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        sgtitle('stable isotope',TX{:},FS{:});
        hold on;


    end
end

%save outputs

if ~isfolder([opdir,'/',name])
    mkdir([opdir,'/',name]);
end

name_save = [opdir,'/',name,'/',name,'_compo'];
print(fh1,name_save,'-dpng','-r600','-opengl');
name_save = [opdir,'/',name,'/',name,'_melt'];
print(fh2,name_save,'-dpng','-r600','-opengl');
name_save = [opdir,'/',name,'/',name,'_SI'];
print(fh3,name_save,'-dpng','-r600','-opengl');



