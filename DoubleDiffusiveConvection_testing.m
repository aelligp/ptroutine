%% Double Diffusive COnvectoin Analyses
% Use this script to analyse the different factors which influence the
% layered convection in the nakhla model
% different density of T, c and x
% Also applicable to differences in vertical convection speed of different
% parameter sets

%define environment and load files
runID   = '2D_Ta8_interm_N200';
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
rhoo =  rhom0.*ones(size(T)); %starting density

%% Theramal density difference - stimmt
%Calculate the temperature difference of layers
Ttop = mean(T(ZZ>0.5 & ZZ<5));
Tbot = mean(T(ZZ>5 & ZZ<8.5));
Tcp  = mean(T(ZZ>8.5 & ZZ<10)); %CP = Cumulate pile
dT = Tbot - Ttop;
drhoT = (-rhom0.*aTm + rhox0.*aTx).*dT
%% Compositional density difference - something wrong

%Calculate the composition difference
% cmtop = mean(cm(ZZ>0.5 & ZZ<5));
% cxtop = mean(cx(ZZ>0.5 & ZZ<5));
% cmbot = mean(cm(ZZ>5 & ZZ<9.5));
% cxbot = mean(cx(ZZ>5 & ZZ<9.5));
% dcm = cmbot - cmtop;
% dcx = cxbot - cxbot;
% drhoC = -rhom.*gCm.*(dcm - (perCx+perCm)./2) + rhox.*gCx.*(dcx - (perCx+perCm)./2)
%% Compositional density difference 
ctop = mean(c(ZZ>0.5 & ZZ<5))
cbot = mean(c(ZZ>5 & ZZ<8.5));
ccp  = mean(c(ZZ>8.5 & ZZ<10)); %CP = Cumulate pile
dC = cbot - ctop;
drhoC = -dC.*gCm.*rhom0

%% Crystalline density difference
%Calculate crystallinitiy difference
xtop = mean(x(ZZ>0.5 & ZZ<5));
xbot = mean(x(ZZ>5 & ZZ<9));
xcp  = mean(x(ZZ>8.5 & ZZ<10)); %CP = Cumulate pile
dX = xbot - xtop; 

%Calculate density diff. m and x
drhoxm0 = rhox0-rhom0; % not sure if correct yet
drhotop = mean(rhox(ZZ>0.5 & ZZ<7)-rhom(ZZ>0.5 & ZZ<7));
drhobot = mean(rhox(ZZ>8.5 & ZZ<9)-rhom(ZZ>8.5 & ZZ<9));
drhocp  = mean(rhox(ZZ>8.5 & ZZ<10)-rhom(ZZ>8.5 & ZZ<10));

%calculating mixing density of layers
drhomixtop = dX*(drhotop)+mean(rhom(ZZ>0.5 & ZZ<8));
drhomixbot = dX*(drhobot)+mean(rhom(ZZ>8 & ZZ<9));

%Density difference to analyse buoyancy
drhoT = -dT.*aTm.*rhom0 

drhoX = dX.*drhoxm0
drhoX = xbot*drhobot - xtop*drhotop

drhotop = mean(rho(ZZ>0.5 & ZZ<7));
drhobot = mean(rho(ZZ>7 & ZZ<9));

drhomix = drhobot - drhotop
sumrho = drhoC + drhoX + drhoT

if sumrho <0 sprintf('unstable')
else sprintf('stable')
end
%% rho bar calc of YQW

rhoxbot = mean(rhox(ZZ>5 & ZZ<9)); 
rhoxtop = mean(rhox(ZZ>0.5 & ZZ<5));
rhombot = mean(rhom(ZZ>5 & ZZ<9)); 
rhomtop = mean(rhom(ZZ>0.5 & ZZ<5));

% rho_bar_bot - rho_bar_top
drhobar = 0




% drhoTT = [3.0, 2.8, 3.27]
% drhoCC = [-5.7, -8.18, -9.41]
% drhoXX = [2.4, 11.72, 13.42]

% figure(1)
% scatter(drhoTT,drhoCC,drhoXX);
% hold on
% scatter(X,drhox)
% scatter(X,drhoT)
% 
% legend('X','drhoC','drhox','drhoT')