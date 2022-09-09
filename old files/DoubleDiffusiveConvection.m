%% Double Diffusive COnvectoin Analyses
% Use this script to analyse the different factors which influence the
% layered convection in the nakhla model 
% different density of T, c and x
% Also applicable to differences in vertical convection speed of different
% parameter sets

%define environment and load files
runID   = '1D_Ta4_bas';
% outdir  = '../Cluster/200resolution/intermediate/Ta8/out'
outdir  = '../Cluster/out/';
path    = strcat(outdir,runID);
for i = 65
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
rhom = rhom0 .* (1 - aTm.*(T-perT-273.15) - gCm.*(cm-(perCx+perCm)/2));
rhox = rhox0 .* (1 - aTx.*(T-perT-273.15) - gCx.*(cx-(perCx+perCm)/2));

%% Layer calculations
% Calculate temperature differences of layers
Ttop = mean(T(ZZ>0.5 & ZZ<4));
Tbot = mean(T(ZZ>4 & ZZ<8));
Tcp  = mean(T(ZZ>8 & ZZ<10)); % Cumulate pile

% Calculate the composition difference
% volume fraction
cmtop = mean(cm(ZZ>0.5 & ZZ<4));
cmbot = mean(cm(ZZ>4 & ZZ<8));
cxtop = mean(cx(ZZ>0.5 & ZZ<4));
cxbot = mean(cx(ZZ>4 & ZZ<8));

%Calculate melt fraction top and bot layer
mtop = mean(m(ZZ>0.5 & ZZ<4));
mbot = mean(m(ZZ>4 & ZZ<8));

%Calculate crystallinitiy difference
xtop = mean(x(ZZ>0.5 & ZZ<4));
xbot = mean(x(ZZ>4 & ZZ<8));
xcp  = mean(x(ZZ>8.5 & ZZ<10)); %CP = Cumulate pile

% xtal and melt densities of layers
rhoxbot = mean(rhox(ZZ>4 & ZZ<8)); 
rhoxtop = mean(rhox(ZZ>0.5 & ZZ<4));
rhombot = mean(rhom(ZZ>4 & ZZ<8)); 
rhomtop = mean(rhom(ZZ>0.5 & ZZ<4));

% volume fractions of layers 
mubot = mean(mu(ZZ>4 & ZZ<8)); 
mutop = mean(mu(ZZ>0.5 & ZZ<4));
chibot = mean(chi(ZZ>4 & ZZ<8)); 
chitop = mean(chi(ZZ>0.5 & ZZ<4));

% Calculated mixture density of model layers
rhotop = mean (rho(ZZ>0.5 & ZZ<4));
rhobot = mean (rho(ZZ>4 & ZZ<8));

%% Density difference analyses

drhoT = (chibot.*(-rhox0.*aTx.*(Tbot-Ttop)))+(mubot.*(-rhom0.*aTm.*(Tbot-Ttop)))


drhoC = (chibot.*(-rhox0.*gCx.*(cxbot-cxtop)))+(mubot.*(-rhom0.*gCm.*(cmbot-cmtop))) %compositional difference


drhoX = (rhoxtop -rhomtop).*(chibot-chitop)

% mixture density difference rhobar = effect on temp + compo + xtals
drhobar = chibot.*(-rhox0.*aTx.*(Tbot-Ttop))+(mubot.*(-rhom0.*aTm.*(Tbot-Ttop))) + ...
(chibot.*(-rhox0.*gCx.*(cxbot-cxtop)))+(mubot.*(-rhom0.*gCm.*(cmbot-cmtop))) + ...
(rhoxtop -rhomtop).*(chibot-chitop)


rhobar = rhobot -rhotop % model mixture density


% if drhobar <0 sprintf('unstable')
% else sprintf('stable')
% end


% % Initiate figure
% % prepare for plotting
% TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
% TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
% UN = {'Units','Centimeters'};
% LW = {'LineWidth',1};
% 
% LS = {'LineStyle','-','--','-.',':'};
% %CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
% 
% fh(1) = figure(1); 
% plot(time/hr, drhobar,'-o',LW{:},'DisplayName',txt); axis xy tight; box on; hold on
% title('stable isotope assimilation',TX{:},FS{:}); xlabel('time [hr]',TX{:},FS{:}); ylabel('Stable isotope [$\%$]',TX{:},FS{:}); set(gca,TL{:},TS{:});

 end

% xtal = (rhoxtop -rhomtop).*(chibot-chitop)
% 
% xtaleffect = rhoxtop.*chibot + rhomtop.*mubot - rhoxtop.*chitop - rhomtop.*mutop