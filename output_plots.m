
function [Nx,Nz,fh1,fh2,fh3,fh4] = output_plots(name,runID,i,parfile, contfile)  %for 0D fh1-3 and for 1D fh1-4

%loading workspace
if exist(parfile,'file'); load(parfile); end
if exist(contfile,'file'); load(contfile); end


% define Nx and Nz for output file to recognize if 0D, 1D
X         = -h/2:h:L+h/2;
Z         = -h/2:h:D+h/2;
Nx = length(X);
Nz = length(Z);

%necessary for 1D plots
[XX,ZZ]   = meshgrid(X,Z);
Xfc       = (X(1:end-1)+X(2:end))./2;
Zfc       = (Z(1:end-1)+Z(2:end))./2;
[XXu,ZZu] = meshgrid(Xfc,Z);
[XXw,ZZw] = meshgrid(X,Zfc);
% update; 
% rhom = rhom0 .* (1 - aTm.*(T-perT) - gCm.*(cm-(perCx+perCm)/2));
% rhox = rhox0 .* (1 - aTx.*(T-perT) - gCx.*(cx-(perCx+perCm)/2));
% rhof = rhof0 .* (1 - aTf.*(T-perT) + bPf.*(Pt-Ptop ));
% rhoref  = mean(mean(rho(2:end-1,2:end-1)));
% ff = max(1e-6,min(1-1e-6,permute(cat(3,chi,mu,phi),[3,1,2])));
% Csgr = ((1-ff)./[dx;1e-16;df].^2.*kv.*thtv).^-1;
% Csgr_x = squeeze(Csgr(1,:,:)) + 1e-18;
% wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-rhoref)*g0.*((Csgr_x(1:end-1,:).*Csgr_x(2:end,:)).^0.5); % crystal segregation speed

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'}; LW = {'LineWidth',1};


LS = {'-' '--' '-.' ':'};
LC = {'k' 'r' 'b' 'g'};



if Nx <= 10 && Nz <= 10  % create 0D plots

    fh1 = figure(1); 
    subplot(5,1,1)
    plot(hist.time/hr,hist.T(:,2),'LineStyle',LS{i},'Color',LC{2},LW{:}); axis xy tight; box on; hold on
    title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,2)
    plot(hist.time/hr,hist.cx(:,2)*100,'LineStyle', LS{i},'Color', LC{3},LW{:}), hold on;
    plot(hist.time/hr,hist.cm(:,2)*100,'LineStyle', LS{i},'Color', LC{2},LW{:}); hold on;
    plot(hist.time/hr,hist.c(:,2)./(1-hist.f(:,2))*100, 'LineStyle', LS{i},'Color', LC{1}, LW{:}); axis xy tight; box on;
    title('$\bar{c}/(1-f)$ [wt\% SiO$_2$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,3)
    semilogy(hist.time/hr,max(1e-6,hist.vf(:,2)*100),'LineStyle', LS{i},'Color', LC{3},LW{:}), hold on;
    semilogy(hist.time/hr,max(1e-6,hist.vm(:,2)*100),'LineStyle', LS{i},'Color', LC{2},LW{:}); hold on;
    semilogy(hist.time/hr,max(1e-6,hist.v(:,2)./(1-hist.x(:,2))*100),'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on;
    title('$\bar{v}/(1-x)$ [wt\% H$_2$O]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,4)
    plot(hist.time/hr,hist.chi(:,2)*100.*(hist.chi(:,2)>1e-9),'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on;
    title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,5)
    plot(hist.time/hr,hist.phi(:,2)*100.*(hist.phi(:,2)>1e-9),'LineStyle', LS{i},'Color', LC{2},LW{:}); axis xy tight; box on; hold on;
    title(['$\phi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [hr]',TX{:},FS{:});
    %legend(runID,'Interpreter','latex','location','best');

    fh2 = figure(2);  
    subplot(5,1,1)
    plot(hist.time/hr,hist.rho(:,2), 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on;
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,2)
    plot(hist.time/hr,log10(hist.eta(:,2)), 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on;
    title('$\bar{\eta}$ [log$_{10}$ Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,3)
    plot(hist.time/hr,hist.Gx(:,2)./hist.rho(:,2)*hr*100.*(hist.chi(:,2)>1e-9), 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on;
    title('$\Gamma_x/\bar{\rho}$ [\%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,4)
    plot(hist.time/hr,hist.Gf(:,2)./hist.rho(:,2)*hr*100.*(hist.phi(:,2)>1e-9),'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on; 
    title('$\Gamma_f/\bar{\rho}$ [\%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,5)
    plot(hist.time/hr,hist.dV(:,2)*hr, 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on;
    title('$\dot{V}$ [1/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});

    xlabel('Time [hr]',TX{:},FS{:});
    %legend(runID,'Location','bestoutside');

    fh3 = figure(3);  
    subplot(5,1,1)
    plot(hist.time/hr,hist.it(:,2), 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on;
    title('incomp. trace',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,2)
    plot(hist.time/hr,hist.ct(:,2), 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on;
    title('comp. trace',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,3)
    plot(hist.time/hr,hist.si(:,2), 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on;hold on;
    title('stable isotope',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,4)
    plot(hist.time/hr,hist.rip(:,2), 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on;
    title('radiogenic parent',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(5,1,5)
    plot(hist.time/hr,hist.rid(:,2), 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis xy tight; box on; hold on;
    title('radiogenic daughter',TX{:},FS{:}); set(gca,TL{:},TS{:});
    
    xlabel('Time [hr]',TX{:},FS{:});
    %legend(runID,'Location','bestoutside');

elseif Nx <= 10  % create 1D plots

    fh1 = figure(1); %legend(runID,'Location','bestoutside');
    subplot(1,5,1)
    plot(mean(T(2:end-1,2:end-1),2),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title('$T [^\circ$C]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,2)
    plot(mean(cx(2:end-1,2:end-1),2)*100,Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{3},LW{:}), hold on;
    plot(mean(cm(2:end-1,2:end-1),2)*100,Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{2},LW{:})
    plot(mean(c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)),2)*100,Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on;
    title('$\bar{c}/(1-f)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,3)
    semilogx(mean(vf(2:end-1,2:end-1),2)*100,Z(2:end-1).','LineStyle', LS{i},'Color', LC{3},LW{:}); hold on;
    semilogx(mean(max(1e-6,vm(2:end-1,2:end-1)*100),2),Z(2:end-1).','LineStyle', LS{i},'Color', LC{2},LW{:});
    semilogx(mean(max(1e-6,v(2:end-1,2:end-1)./(1-x(2:end-1,2:end-1))*100),2),Z(2:end-1).','LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on;
    title('$\bar{v}/(1-x)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,4)
    plot(mean(chi(2:end-1,2:end-1),2)*100.*(mean(chi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,5)
    plot(mean(phi(2:end-1,2:end-1),2)*100.*(mean(phi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title(['$\phi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
   

    fh2 = figure(2); %legend(runID,'Location','bestoutside');
    subplot(1,4,1)
    plot(mean(-W(:,2:end-1),2)*hr,Zfc.', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title('$W$ [m/hr]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,2)
%     plot(mean(-(x(1:end-1,2:end-1)+x(2:end,2:end-1))/2.*wx(:,2:end-1),2)*hr,Zfc.', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
%     title('$w_\Delta^x$ [m/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%     subplot(1,4,3)
%     plot(mean(-(f(1:end-1,2:end-1)+f(2:end,2:end-1))/2.*wf(:,2:end-1),2)*hr,Zfc.', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
%     title('$w_\Delta^f$ [m/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%     subplot(1,4,4)
    plot(mean(P(2:end-1,2:end-1),2)/1e3,Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title('$P$ [kPa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    

    fh3 = figure(3); 
    subplot(1,5,1)
    plot(mean(rho(2:end-1,2:end-1),2),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,2)
    plot(mean(log10(eta(2:end-1,2:end-1)),2),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title('$\bar{\eta}$ [log$_{10}$ Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,3)
    plot(mean(Gx(2:end-1,2:end-1)./rho(2:end-1,2:end-1),2)*100*hr.*(mean(chi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title('$\Gamma_x/\bar{\rho}$ [\%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,4)
    plot(mean(Gf(2:end-1,2:end-1)./rho(2:end-1,2:end-1),2)*100*hr.*(mean(phi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title('$\Gamma_f/\bar{\rho}$ [\%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%     subplot(1,5,5)
%     plot(mean(VolSrc(2:end-1,2:end-1),2),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
%     title('$\dot{V}$ [1/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
%     

    fh4 = figure(4);  
    subplot(1,5,1);
    plot(mean(it(2:end-1,2:end-1),2),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title('incomp. trace',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,2)
    plot(mean(ct(2:end-1,2:end-1),2),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on;hold on;
    title('comp. trace',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,3)
    plot(mean(si(2:end-1,2:end-1),2),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title('stable isotope',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,4)
    plot(mean(rip(2:end-1,2:end-1),2),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title(['radiogenic parent'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,5)
    plot(mean(rid(2:end-1,2:end-1),2),Z(2:end-1).', 'LineStyle', LS{i},'Color', LC{1},LW{:}); axis ij tight; box on; hold on;
    title(['radiogenic daughter'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    
   
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


