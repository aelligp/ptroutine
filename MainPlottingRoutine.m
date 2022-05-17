%% Main plotting routine %%
% After replotting add etact in output_plots again in output_plot or copy new file and
% delete clf;
close all; clear all; clc;

name = 'Name your Plots';

n = 2; %change to number of desired plots
plotLabels = cell(n,1);



for i = 1: n

    %run(['Plot' num2str(i) '()']);
    [outdir,runID,path,parfile,contfile] = Plot_1();
    output_plots(i,parfile, contfile);
    
    hold on
end

plotLabels{1} = '0D_anh_Tau_T4<<Tau_a8_bas049_wall_4H2O';
plotLabels{2} = '0D_anh_Tau_T4<<Tau_a8_interm060_wall_4H2O';
plotLabels{3} = '0D_anh_Tau_T4<<Tau_a8_rhy070_wall_4H2O';
legend(plotLabels, 'Location', 'east')
legend

% 
% for i = 1:34
%       filename = ['file',int2str(i)];
%       load([filename '.txt'],'-ascii');
%       plot( filename(:,1), filename(:,2) );
% end