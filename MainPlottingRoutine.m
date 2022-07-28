%% Main plotting routine %%
% This routine enables you to plot multiple simulations together to see the
% differences. It uses the output_plots.m file to plot the different 0D/1D
% plots. IF desired you can also change the line color to change with every
% run rather than the style --> change to LC{i}
% PAel June 2022

 close all; clear all; clc;

%Give the plot a name to specify what parameters you plot together
name = 'Ta8_40'; 
opdir = 'out/';
% n = 3; %change to number of desired plots

for i = 1:4

    switch i
        case  1
            [outdir,runID,path,parfile,contfile] = Plot_1(); 
            txt = 'time step 1';
            n = 1;
        case  2
            [outdir,runID,path,parfile,contfile] = Plot_2(); 
            txt = 'basaltic, time step 40';
            n= 80;
        case  3
            [outdir,runID,path,parfile,contfile] = Plot_3();
            txt = 'intermediate, time step 40';
            n= 160;
        case  4
            [outdir,runID,path,parfile,contfile] = Plot_4();
            txt = 'rhyolitic, time step 40';
            n = 256;
%         case  2
%             [outdir,runID,path,parfile,contfile] = Plot_5();
%             txt = 'basaltic, time step 80';
%             n= 80;
%         case  3
%             [outdir,runID,path,parfile,contfile] = Plot_6();
%             txt = 'intermediate, time step 80';
%             n= 160;
%         case  4
%             [outdir,runID,path,parfile,contfile] = Plot_7();
%             txt = 'rhyolitic, time step 80';
%             n = 256;
%         case  8
%             [outdir,runID,path,parfile,contfile] = Plot_8();
%             txt = '$\tau_a$(8h), time step 80';
%             n = 256;
    end

    %for 0D fh1-3 and for 1D fh1-4
    [Nx,Nz,fh1,fh2,fh3,fh4] = output_plots_new(name,txt,n,runID,i,parfile, contfile); 
    hold on;
    legend('time step 1','basaltic, time step 40','intermediate, time step 40','rhyolitic, time step 40', ...
         'Interpreter','latex','location','best');
        %'basaltic, time step 80', 'intermediate, time step 80','rhyolitic, time step 80', ...
       
    
    
end
% legend(txt,'Interpreter','latex','location','best')


  %save outputs
    if ~isfolder([opdir,'/',name])
        mkdir([opdir,'/',name]);
    end
     if Nx <= 10 && Nz <= 10  % print 0D plots
        name_save = [outdir,'/',runID,'/',name,'_tch',num2str(i)];
        print(fh1,name_save,'-dpng','-r300','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_aux',num2str(i)];
        print(fh2,name_save,'-dpng','-r300','-opengl');
    elseif Nx <= 10  % create 1D plots
        name_save = [opdir,'/',name,'/',name,'_compo'];
        print(fh1,name_save,'-dpng','-r600','-opengl');
        name_save = [opdir,'/',name,'/',name,'_volume'];
        print(fh2,name_save,'-dpng','-r600','-opengl');
        name_save = [opdir,'/',name,'/',name,'_density'];
        print(fh3,name_save,'-dpng','-r600','-opengl');
%         name_save = [opdir,'/',name,'/',name,'_elements'];
%         print(fh4,name_save,'-depsc','-r600','-opengl');
    else
        name_save = [outdir,'/',runID,'/',name,'_vep_',num2str(i)];
        print(fh1,name_save,'-dpng','-r600','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_tch_',num2str(i)];
        print(fh2,name_save,'-dpng','-r600','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_phs_',num2str(i)];
        print(fh3,name_save,'-dpng','-r600','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_sgr_',num2str(i)];
        print(fh4,name_save,'-dpng','-r600','-opengl');
        name_save = [outdir,'/',runID,'/',name,'_gch',num2str(i)];
        print(fh5,name_save,'-dpng','-r600','-opengl');
%         name_save = [runID,'_eql',num2str(i)];
%         print(fh7,name,'-dpng','-r300','-opengl');
    end


