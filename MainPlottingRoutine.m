%% Main plotting routine %%
% This routine enables you to plot multiple simulations together to see the
% differences. It uses the output_plots.m file to plot the different 0D/1D
% plots. IF desired you can also change the line color to change with every
% run rather than the style --> change to LC{i}
% PAel June 2022

close all; clear all; clc;

%Give the plot a name to specify what parameters you plot together
name = '0D_Ta4_bas'; 

n = 1; %change to number of desired plots

for i = 1: n

    switch i
        case  1
            [outdir,runID,path,parfile,contfile] = Plot_1(); 
            name_1 = '0D_Ta4_bas';
        case  2
            [outdir,runID,path,parfile,contfile] = Plot_2(); 
            name_2 = 'TauT4 Ta8 interm wall';
        case  3
            [outdir,runID,path,parfile,contfile] = Plot_3();
            name_3 = 'TauT4 Ta8 rhyolitic wall';

    end

    %for 0D fh1-3 and for 1D fh1-4
    [Nx,Nz,fh1,fh2] = output_plots(name,runID,i,parfile, contfile); 
    
    
end
%legend(name_1,name_2,name_3, 'Interpreter','latex','location','best')
%legend

  %save outputs
    if ~isfolder([path,'/',runID])
        mkdir([outdir,'/',runID]);
    end
    if Nx <= 10 && Nz <= 10  % print 0D plots
        name = [outdir,'/',runID,'/',runID,'_tch_',num2str(i)];
        print(fh1,name,'-dpng','-r300','-opengl');
        name = [outdir,'/',runID,'/',runID,'_aux_',num2str(i)];
        print(fh2,name,'-dpng','-r300','-opengl');
        name = [opdir,'/',runID,'/',runID,'_eql',num2str(i)];
        print(fh7,name,'-dpng','-r300','-opengl');
    elseif Nx <= 10  % create 1D plots
        name = [outdir,'/',runID,'/',runID,'_sol_',num2str(i)];
        print(fh1,name,'-dpng','-r300','-opengl');
        name = [outdir,'/',runID,'/',runID,'_aux_',num2str(i)];
        print(fh2,name,'-dpng','-r300','-opengl');
        name = [outdir,'/',runID,'/',runID,'_gch_',num2str(i)];
        print(fh3,name,'-dpng','-r300','-opengl');
      
    else
        name_save = [outdir,'/',runID,'/',name,'_vep_',num2str(i)];
        print(fh1,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_tch_',num2str(i)];
        print(fh2,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_phs_',num2str(i)];
        print(fh3,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_sgr_',num2str(i)];
        print(fh4,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_gch_',num2str(i)];
        print(fh5,name_save,'-dpng','-r300','-opengl');
        %name_save = [runID,'_eql',num2str(i)];
        %print(fh7,name,'-dpng','-r300','-opengl');
    end



