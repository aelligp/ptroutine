%% Double Diffusive COnvectoin Analyses
% Use this script to analyse the different factors which influence the
% layered convection in the nakhla model 
% different density of T, c and x
% Also applicable to differences in vertical convection speed of different
% parameter sets

%define environment and load files
runID   = '2D_Ta8_interm_N200';
% outdir  = '../Cluster/200resolution/intermediate/Ta8/out'
outdir  = '../Cluster/out/';
path    = strcat(outdir,runID);

addpath(path);

parfile =  [path ,'/', runID, '_par.mat'];
contfile=  [path ,'/', runID, '_127.mat']; %set to desired .mat file number
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
Ttop = mean(T(ZZ>0.5 & ZZ<5));
Tbot = mean(T(ZZ>5 & ZZ<8.5));
Tcp  = mean(T(ZZ>8.5 & ZZ<10)); % Cumulate pile

% Calculate the composition difference
% volume fraction
cmtop = mean(cm(ZZ>0.5 & ZZ<5));
cmbot = mean(cm(ZZ>5 & ZZ<9));
cxtop = mean(cx(ZZ>0.5 & ZZ<5));
cxbot = mean(cx(ZZ>5 & ZZ<9));

%Calculate melt fraction top and bot layer
mtop = mean(m(ZZ>0.5 & ZZ<5));
mbot = mean(m(ZZ>5 & ZZ<9));

%Calculate crystallinitiy difference
xtop = mean(x(ZZ>0.5 & ZZ<5));
xbot = mean(x(ZZ>5 & ZZ<9));
xcp  = mean(x(ZZ>8.5 & ZZ<10)); %CP = Cumulate pile

% xtal and melt densities of layers
rhoxbot = mean(rhox(ZZ>5 & ZZ<9)); 
rhoxtop = mean(rhox(ZZ>0.5 & ZZ<5));
rhombot = mean(rhom(ZZ>5 & ZZ<9)); 
rhomtop = mean(rhom(ZZ>0.5 & ZZ<5));

% volume fractions of layers 
mubot = mean(mu(ZZ>5 & ZZ<9)); 
mutop = mean(mu(ZZ>0.5 & ZZ<5));
chibot = mean(chi(ZZ>5 & ZZ<9)); 
chitop = mean(chi(ZZ>0.5 & ZZ<5));

% Calculated mixture density of model layers
rhotop = mean (rho(ZZ>0.5 & ZZ<5));
rhobot = mean (rho(ZZ>5 & ZZ<9));

%% Density difference analyses

drhoT = (chibot.*(-rhox0.*aTx.*(Tbot-Ttop)))+(mubot.*(-rhom0.*aTm.*(Tbot-Ttop)))


drhoC = (chibot.*(-rhox0.*gCx.*(cxbot-cxtop)))+(mubot.*(-rhom0.*gCm.*(cmbot-cmtop))) %compositional difference


drhoX = (rhoxtop -rhomtop).*(chibot-chitop)

% mixture density difference rhobar = effect on temp + compo + xtals
drhobar = chibot.*(-rhox0.*aTx.*(Tbot-Ttop))+(mubot.*(-rhom0.*aTm.*(Tbot-Ttop))) + ...
(chibot.*(-rhox0.*gCx.*(cxbot-cxtop)))+(mubot.*(-rhom0.*gCm.*(cmbot-cmtop))) + ...
(rhoxtop -rhomtop).*(chibot-chitop)


rhobar = rhobot -rhotop % model mixture density


if sumrho <0 sprintf('unstable')
else sprintf('stable')
end



% xtal = (rhoxtop -rhomtop).*(chibot-chitop)
% 
% xtaleffect = rhoxtop.*chibot + rhomtop.*mubot - rhoxtop.*chitop - rhomtop.*mutop