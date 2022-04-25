%% Main plotting routine %%
% After replotting add etact again in output_plot or copy new file and
% delete clf;
close all; clear all; clc;

name = 'Name your Plots';

 n = 3; %
 plotLabels = cell(1,n);
for i = 1: n
 run(['Plot' num2str(i) '.m']);
 plotLabels{i} = ['Plot_' num2str(i)];
 hold on
 
end
legend(plotLabels,LS);
legend.Layout.Tile = 'east';
% % save output to file
% 
% if Nx <= 10 && Nz <= 10  % print 0D plots
%     %name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
%     print(fh1,name,'-dpng','-r300','-opengl');
%     %name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
%     print(fh2,name,'-dpng','-r300','-opengl');
%     %name = [opdir,'/',runID,'/',runID,'_gch_',num2str(floor(step/nop))];
%     print(fh3,name,'-dpng','-r300','-opengl');
%     %name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
%     print(fh7,name,'-dpng','-r300','-opengl');
% elseif Nx <= 10  % create 1D plots
%     name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
%     print(fh1,name,'-dpng','-r300','-opengl');
%     name = [opdir,'/',runID,'/',runID,'_vep_',num2str(floor(step/nop))];
%     print(fh2,name,'-dpng','-r300','-opengl');
%     name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
%     print(fh3,name,'-dpng','-r300','-opengl');
%     name = [opdir,'/',runID,'/',runID,'_gch_',num2str(floor(step/nop))];
%     print(fh4,name,'-dpng','-r300','-opengl');
%     name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
%     print(fh7,name,'-dpng','-r300','-opengl');
% else
%     name = [opdir,'/',runID,'/',runID,'_vep_',num2str(floor(step/nop))];
%     print(fh1,name,'-dpng','-r300','-opengl');
%     name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
%     print(fh2,name,'-dpng','-r300','-opengl');
%     name = [opdir,'/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
%     print(fh3,name,'-dpng','-r300','-opengl');
%     name = [opdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
%     print(fh4,name,'-dpng','-r300','-opengl');
%     name = [opdir,'/',runID,'/',runID,'_gch',num2str(floor(step/nop))];
%     print(fh5,name,'-dpng','-r300','-opengl');
%     name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
%     print(fh7,name,'-dpng','-r300','-opengl');
% end