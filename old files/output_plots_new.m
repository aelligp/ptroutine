
function [Nx,Nz,fh1,fh2,fh3,fh4] = output_plots_new(name,txt,n,runID,i,parfile, contfile)  %for 0D fh1-3 and for 1D fh1-4

%loading workspace
if exist(parfile,'file'); load(parfile); end
if exist(contfile,'file')
    load(contfile,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
end

addpath('src\')
load ocean %ocean colormap

% define Nx and Nz for output file to recognize if 0D, 1D

% % for the illustration of 2D plots in 1D use these lines
% L = 10/200;
% X         = -h/2:h:L+h/2;
% Z         = -h/2:h:D+h/2;

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
CL = {'Color', copper(n,1:3)};
%CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};



if Nx <= 10  % create 1D plots
        % Thermo chemical
        fh1 = figure(1);
        subplot(1,5,i)
        plot(mean(c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)),2)*100,Z(2:end-1).',LS{[1,2]}, 'DisplayName',txt, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title('$\bar{c}/(1-f)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,2)
%         plot(mean(cm(2:end-1,2:end-1),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%         title('$c_m$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,3)
%         plot(mean(c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)),2)*100,Z(2:end-1).',LS{[1,2]}, 'DisplayName',txt, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%         title('$\bar{c}/(1-f)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         legend(txt,'Interpreter','latex','location','best')
%         subplot(1,5,4)
% %         semilogx(geomean(min(eta(2:end-1,2:end-1),etam(2:end-1,2:end-1)),2),Z(2:end-1).', LS{[1,2]}, 'DisplayName',txt,CL{[1,2]},LW{:}); axis ij tight; box on; hold on
%         semilogx(geomean(eta(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on; 
%         title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         legend(txt,'Interpreter','latex','location','best')
       

        fh2 = figure(2);
        subplot(1,5,i)
        plot(mean(mu (2:end-1,2:end-1),2)*100.*(mean(mu (2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title(['$\mu$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,2) %change again
%          semilogx(mean(max(1e-6,v(2:end-1,2:end-1)./(1-x(2:end-1,2:end-1))*100),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on; 
%         title('$\bar{v}/(1-x)$ [wt\%]',TX{:},FS{:});
%        
%         subplot(1,5,3)
%         plot(mean(chi(2:end-1,2:end-1),2)*100.*(mean(chi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%         title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,4)
%         plot(mean(phi(2:end-1,2:end-1),2)*100.*(mean(phi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%         title(['$\phi$[vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,5)
%         plot(mean(P(2:end-1,2:end-1),2)/1e3,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%         title('$P$ [kPa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
%         fh3 = figure(3);
%         subplot(1,5,1)
%         plot(mean(rhox(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%         title('$\rho_x$ [kg/m$^3$]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,2)
%         plot(mean(rhom(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%         title('$\rho_m$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,3)
%         plot(mean(rho (2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%         title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,4)
%         plot(mean(-W(:,2:end-1),2)*hr,Zfc.',LS{[1,2]},'DisplayName',txt, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
% %         plot(mean(-(f(1:end-1,2:end-1)+f(2:end,2:end-1))/2.*wf(:,2:end-1),2)*hr,Zfc.',LS{[1,2]}, CL{[1,5]},LW{:});
%         title('$W$ [m/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% 
%         subplot(1,5,5)
%         plot(mean(-(x(1:end-1,2:end-1)+x(2:end,2:end-1))/2.*wx(:,2:end-1),2)*hr,Zfc.',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on; 
%         title('$w_\Delta^x$ [m/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         legend(txt,'Interpreter','latex','location','best')


        %         subplot(1,5,3)
        %         plot(mean(    Gx(2:end-1,2:end-1)./rho(2:end-1,2:end-1),2)*hr.*(mean(chi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        %         plot(mean(10.*Gf(2:end-1,2:end-1)./rho(2:end-1,2:end-1),2)*hr.*(mean(phi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,5]},LW{:});
        %         title('$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        %         subplot(1,5,4)
        %         plot(mean(VolSrc(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on;
        %         title('$\dot{V}$ [1/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    

        fh3 = figure(3);
        subplot(1,5,i)
        plot(mean(si(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title('stable isotope',TX{:},FS{:}); set(gca,TL{:},TS{:});ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         semilogx(mean(max(1e-3,min(1e3,itx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%         title('$it_x$',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,2)
%         semilogx(mean(max(1e-3,min(1e3,itm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
%          title('$it_m$',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,3)
       
%         subplot(1,5,4)
%         semilogx(mean(max(1e-3,min(1e3,ctx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        subplot(1,5,2)
%         semilogx(mean(max(1e-3,min(1e3,ctm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        semilogx(mean(max(1e-3,min(1e3,ct (2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title('comp. trace',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,3)
        semilogx(mean(max(1e-3,min(1e3,it (2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        title('incomp. trace',TX{:},FS{:});
       

        %         subplot(1,5,4)
        %         semilogx(mean(max(1e-3,min(1e3,ripx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        %         semilogx(mean(max(1e-3,min(1e3,ripm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:});
        %         semilogx(mean(max(1e-3,min(1e3,rip (2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
        %         title(['radiogenic parent'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        %         subplot(1,5,5)
        %         semilogx(mean(max(1e-3,min(1e3,ridx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        %         semilogx(mean(max(1e-3,min(1e3,ridm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:});
        %         semilogx(mean(max(1e-3,min(1e3,rid(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
        %         title(['radiogenic daughter'],TX{:},FS{:}); set(gca,TL{:},TS{:});

   
%         
%         fh1 = figure(1);   
% %         subplot(1,5,1)
% %         plot(mean(T(2:end-1,2:end-1),2)-273.15,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on;
% %         title('$T [^\circ$C]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,2)
%         plot(mean(cx(2:end-1,2:end-1),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
%         plot(mean(cm(2:end-1,2:end-1),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:});
%         plot(mean(c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
%         title('$\bar{c}/(1-f)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,3)
%         semilogx(mean(max(1e-6,vf(2:end-1,2:end-1)*100),2),Z(2:end-1).',LS{[1,2]}, CL{[1,5]},LW{:}); axis ij tight; box on; hold on;
%         semilogx(mean(max(1e-6,vm(2:end-1,2:end-1)*100),2),Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:});
%         semilogx(mean(max(1e-6,v(2:end-1,2:end-1)./(1-x(2:end-1,2:end-1))*100),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
%         title('$\bar{v}/(1-x)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,4) %change again
%         plot(mean(mu (2:end-1,2:end-1),2)*100.*(mean(mu (2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
%         plot(mean(phi(2:end-1,2:end-1),2)*100.*(mean(phi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,5]},LW{:});
%         plot(mean(chi(2:end-1,2:end-1),2)*100.*(mean(chi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:});
%         title(['$\mu$, $\phi$, $\chi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
% %    subplot(1,5,5)
% %         plot(mean(-W(:,2:end-1),2)*hr,Zfc.',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
% %         plot(mean(-(f(1:end-1,2:end-1)+f(2:end,2:end-1))/2.*wf(:,2:end-1),2)*hr,Zfc.',LS{[1,2]}, CL{[1,5]},LW{:});
% %         plot(mean(-(x(1:end-1,2:end-1)+x(2:end,2:end-1))/2.*wx(:,2:end-1),2)*hr,Zfc.',LS{[1,2]}, CL{[1,4]},LW{:});
% %         title('$W$, $w_\Delta^f$, $w_\Delta^x$ [m/hr]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %              
% %         fh2 = figure(3);   
% %         subplot(1,5,1)
% %         plot(mean(rhox(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
% %         plot(mean(rhom(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:});
% %         plot(mean(rho (2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
% %         title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         subplot(1,5,2)
% %         semilogx(geomean(min(eta(2:end-1,2:end-1),etam(2:end-1,2:end-1)),2),Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
% %         semilogx(geomean(eta(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
% %         title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         subplot(1,5,3)
% %         plot(mean(    Gx(2:end-1,2:end-1)./rho(2:end-1,2:end-1),2)*hr.*(mean(chi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
% %         plot(mean(10.*Gf(2:end-1,2:end-1)./rho(2:end-1,2:end-1),2)*hr.*(mean(phi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',LS{[1,2]}, CL{[1,5]},LW{:}); 
% %         title('$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         subplot(1,5,4)
% %         plot(mean(VolSrc(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on;
% %         title('$\dot{V}$ [1/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         subplot(1,5,5)
% %         plot(mean(P(2:end-1,2:end-1),2)/1e3,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on;
% %         title('$P$ [kPa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         
%         fh2 = figure(4);   
%         subplot(1,5,1)
%         semilogx(mean(max(1e-3,min(1e3,itx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
%         semilogx(mean(max(1e-3,min(1e3,itm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:});
%         semilogx(mean(max(1e-3,min(1e3,it (2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
%         title('incomp. trace',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         subplot(1,5,2)
%         semilogx(mean(max(1e-3,min(1e3,ctx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
%         semilogx(mean(max(1e-3,min(1e3,ctm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:});
%         semilogx(mean(max(1e-3,min(1e3,ct (2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
%         title('comp. trace',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         subplot(1,5,3)
% %         plot(mean(si(2:end-1,2:end-1),2),Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:}); axis ij tight; box on;
% %         title('stable isotope',TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         subplot(1,5,4)
% %         semilogx(mean(max(1e-3,min(1e3,ripx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
% %         semilogx(mean(max(1e-3,min(1e3,ripm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:});
% %         semilogx(mean(max(1e-3,min(1e3,rip (2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
% %         title(['radiogenic parent'],TX{:},FS{:}); set(gca,TL{:},TS{:});
% %         subplot(1,5,5)
% %         semilogx(mean(max(1e-3,min(1e3,ridx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
% %         semilogx(mean(max(1e-3,min(1e3,ridm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,3]},LW{:});
% %         semilogx(mean(max(1e-3,min(1e3,rid(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',LS{[1,2]}, CL{[1,2]},LW{:});
% %         title(['radiogenic daughter'],TX{:},FS{:}); set(gca,TL{:},TS{:});
% 
%    
end

% % plot phase diagram
% fh7 = figure(7);  
% TT = linspace(Tphs0+Ptop*clap,Tphs1+Ptop*clap,1e3);
% cc = [linspace(cphs1,(perCx+perCm)/2,ceil((perT-Tphs0)./(Tphs1-Tphs0)*1e3)),linspace((perCx+perCm)/2,cphs0,floor((perT-Tphs1)./(Tphs0-Tphs1)*1e3))];
% [~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,Ptop*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
% plot(CCx,TT, 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis tight; hold on; box on;
% plot(CCm,TT, 'LineStyle', LS{i},'Color', LC{1},LW{:});
% perTs  = perT;
% Tphs0s = Tphs0;
% Tphs1s = Tphs1;
% vv = 0.10*ones(size(TT));
% for i = 1:10
%     perTs  = perT -dTH2O(2)*mean(vv(abs(TT-Ptop*clap-perTs )<1)).^0.75;
%     Tphs0s = Tphs0-dTH2O(1)*mean(vv(abs(TT-Ptop*clap-Tphs0s)<1)).^0.75;
%     Tphs1s = Tphs1-dTH2O(3)*mean(vv(abs(TT-Ptop*clap-Tphs1s)<1)).^0.75;
%     TTi = Tphs0s+Ptop*clap:1:Tphs1s+Ptop*clap;
%     vv = interp1(TT,vv,TTi,'linear','extrap'); TT = TTi;
%     cc = [linspace(cphs1,(perCx+perCm)/2,round((perTs-Tphs0s)./(Tphs1s-Tphs0s)*length(TT))),linspace((perCx+perCm)/2,cphs0,round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*length(TT)))];
%     [~,CCx,CCm,~,~,vv] = equilibrium(0*TT,0*TT,TT,cc,vv,Ptop*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
% end
% plot(CCx,TT, 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis tight; hold on; box on;
% plot(CCm,TT, 'LineStyle', LS{i},'Color', LC{1},LW{:});
% 
% Tplt = T - (Pt-Ptop)*clap;
% cplt = c./(1-f);
% plot(cplt(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'k.',cx(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'b.',cm(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'r.',LW{:},'MarkerSize',15);
% 
% set(gca,'TickLabelInterpreter','latex','FontSize',15)
% title('Phase Diagram','Interpreter','latex','FontSize',22)
% xlabel('Composition','Interpreter','latex','FontSize',18)
% ylabel('Temperature','Interpreter','latex','FontSize',18)


drawnow






end


