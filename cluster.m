% Plot cluster.mat files
close all, clear all

n=500;
for i= 1:n
    sourcedir   = '../Cluster/200resolution/rhyolitic/Ta4/';
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

    %for 0D fh1-2; 1D fh1-4; 2D fh1-5, fh7
    [Nx,Nz,fh1,fh2,fh3,fh4,fh5] = output_plots_cluster(name,runID,i,parfile, contfile);

    %save outputs
    if ~isfolder([outdir,'/',runID])
        mkdir([outdir,'/',runID]);
    end
    if Nx <= 10 && Nz <= 10  % print 0D plots
        name_save = [outdir,'/',runID,'/',name,'_tch'];
        print(fh1,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_aux'];
        print(fh2,name_save,'-dpng','-r300','-opengl');
    elseif Nx <= 10  % create 1D plots
        name_save = [outdir,'/',runID,'/',name,'_tch'];
        print(fh1,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_vep'];
        print(fh2,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_aux'];
        print(fh3,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_gch'];
        print(fh4,name_save,'-dpng','-r300','-opengl');
    else
        name_save = [outdir,'/',runID,'/',name,'_vep_',num2str(i)];
        print(fh1,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_tch_',num2str(i)];
        print(fh2,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_phs_',num2str(i)];
        print(fh3,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_sgr_',num2str(i)];
        print(fh4,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_gch',num2str(i)];
        print(fh5,name_save,'-dpng','-r300','-opengl');
        %name_save = [runID,'_eql',num2str(i)];
        %print(fh7,name,'-dpng','-r300','-opengl');
    end

    clear
end

% save output to file

