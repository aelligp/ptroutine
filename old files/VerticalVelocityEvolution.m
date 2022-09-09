% Velocity field analyses of different parameter set
% Use this to have an appoximation of vertical velocity in the magma
% chamber over time 
% PAel July 2022

sourcedir   = '../Cluster/';
    outdir      = 'out/';
    runID       = '2D_Ta4_rhy_N200';

    path        = strcat(sourcedir,outdir,runID);
    src         = strcat(sourcedir,'src');

    %add to path 
    addpath(path);
    addpath(src);

    %Files we need to plot
    parfile = [path ,'/', runID, '_par.mat']; % parameter file
    % ocean = load('ocean.mat'); %ocean colormap

    name = runID;
    contfile = ([runID '_' num2str(i) '.mat']);

    %loading workspace
if exist(parfile,'file'); load(parfile); end
if exist(contfile,'file')
    load(contfile,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
end

load ocean %ocean colormap

i = 155
step = i .* nop

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


maxvelo = hist.W(:,3);
VertVelo_top = mean(maxvelo(XX == 1))

% plot(hist.time/hr, hist.W(:,3).*hr); 
% legend(runID) 
hold on 