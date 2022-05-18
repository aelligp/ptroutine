%% Main plotting routine %%
% This routine enables you to plot multiple simulations together to see the
% differences. It uses the output_plots.m file to plot the different 0D/1D
% plots. IF desired you can also change the line color to change with every
% run rather than the style --> change to LC{i}

close all; clear all; clc;

%Give the plot a name to specify what parameters you plot together
name = '1D_anh_TT4_Ta8'; 

n = 3; %change to number of desired plots

for i = 1: n

    switch i
        case  1
            [outdir,runID,path,parfile,contfile] = Plot_1(); 
            name_1 = 'TauT4 Ta8 basaltic wall';
        case  2
            [outdir,runID,path,parfile,contfile] = Plot_2(); 
            name_2 = 'TauT4 Ta8 interm wall';
        case  3
            [outdir,runID,path,parfile,contfile] = Plot_3();
            name_3 = 'TauT4 Ta8 rhyolitic wall';

    end

    %for 0D fh1-3 and for 1D fh1-4
    [Nx,Nz,fh1,fh2,fh3,fh4] = output_plots(name,runID,i,parfile, contfile); 
    
    
end
legend(name_1,name_2,name_3, 'Interpreter','latex','location','best')
legend

% save output to file

if Nx <= 10 && Nz <= 10  % print 0D plots
    name_save = [name,'_tch'];
    print(fh1,name_save,'-dpng','-r300','-opengl');
    name_save = [name,'_aux'];
    print(fh2,name_save,'-dpng','-r300','-opengl');
    name_save = [name,'_gch'];
    print(fh3,name_save,'-dpng','-r300','-opengl');
elseif Nx <= 10  % create 1D plots
    name_save = [name,'_tch'];
    print(fh1,name_save,'-dpng','-r300','-opengl');
    name_save = [name,'_vep'];
    print(fh2,name_save,'-dpng','-r300','-opengl');
    name_save = [name,'_aux'];
    print(fh3,name_save,'-dpng','-r300','-opengl');
    name_save = [name,'_gch'];
    print(fh4,name_save,'-dpng','-r300','-opengl');
end
