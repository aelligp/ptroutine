%% Double Diffusive COnvectoin Analyses
% Use this script to analyse the different factors which influence the
% layered convection in the nakhla model
% different density of T, c and x
% Also applicable to differences in vertical convection speed of different
% parameter sets

%define environment and load files
runID   = '2D_Ta8_interm_N200';
outdir  = '../Cluster/200resolution/intermediate/Ta8/out/';
path    = strcat(outdir,runID);
addpath(path);
parfile =  [path ,'/', runID, '_par.mat'];
contfile=  [path ,'/', runID, '_51.mat']; %set to desired .mat file number
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

%Calculate the temperature difference of layers
Ttop = mean(T(ZZ>0.5 & ZZ<2));
Tbot = mean(T(ZZ>2 & ZZ<9.5));
dT = Tbot - Ttop;
%Calculate the composition difference
ctop = mean(c(ZZ>0.5 & ZZ<2));
cbot = mean(c(ZZ>2 & ZZ<9.5));
dC = cbot - ctop;
%Calculate crystallinitiy difference
xtop = mean(x(ZZ>0.5 & ZZ<2));
xbot = mean(x(ZZ>2 & ZZ<9.5));
dX = xbot - xtop; 
%Calculate density diff. m and x
drhomx0 = rhox0-rhom0; 
%Density difference to analyse buoyancy
drhoT = dT.*aTm.*rhom0 
drhoC = dC.*gCx.*rhom0
drhox = dX.*drhomx0


drhoTT = [3.0, 2.8, 3.27]
drhoCC = [-5.7, -8.18, -9.41]
drhoXX = [2.4, 11.72, 13.42]

% figure(1)
% scatter(drhoTT,drhoCC,drhoXX);
% hold on
% scatter(X,drhox)
% scatter(X,drhoT)
% 
% legend('X','drhoC','drhox','drhoT')