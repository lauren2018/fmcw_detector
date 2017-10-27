%% written by Taehwan Kim (taehwan@berkeley.edu)
clear; close all; clc;
format long;
set(0,'DefaultFigureWindowStyle','normal')
addpath './frequency'
addpath './Colormaps'
addpath './distinguishable_colors/'
addpath './hex_and_rgb_v2'
addpath './fig'
addpath './legendflex'
addpath './setgetpos_V1.2'
addpath 'heatmaps'

header;

m=100;
cm_magma=magma(m);
cm_inferno=inferno(m);
cm_plasma=plasma(m);
cm_viridis=viridis(m);

% jet=colormap('jet');
% parula=fake_parula();

fig1handle=fig('units','inches','width',textwidth/2,'height',textwidth/5,'font','Times','fontsize',fs);
colormap(fig1handle,cm_viridis);


subhandle1 = subplot(1,3,1);
res = load('result_fig4a_new.mat');
std_est_periodogram_lse = res.std_est_periodogram_lse;

mincolor = 0; maxcolor = 6;


Prx = {'10^{-11}' '' '' '' '' '10^{-10}' '' '' '' '' '10^{-9}' '' '' '' '' '10^{-8}'};
fnoise_list = {'10^{0}' '' '' '' '' '10^{0.5}' '' '' '' '' '10^{1}'};

hm1 = heatmapcus(flipud(single(std_est_periodogram_lse*100)),fnoise_list,fliplr(Prx), [],...
        'ShowAllTicks', true,'TickTexInterpreter','true',...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'GridLines','-');
% get(gca)
title('Baseline','fontweight','normal');
    set(gca,'ticklength',[0 0]);
    

% set(gca,'TickLabelInterpreter', 'latex');
% xlabel('\textrm{Linewidth (MHz)}','interpreter','latex');
% ylabel('$P_\textrm{RX}$ \textrm{ (W)}','interpreter','latex');    

subhandle2 = subplot(1,3,2);
res = load('result_fig4c_new.mat');
% fnoise_list = res.fnoise_list;
% Prx = res.Prx;
std_est_periodogram_lse = res.std_est_periodogram_lse;

hm2 = heatmapcus(flipud(single(std_est_periodogram_lse*100)),fnoise_list,[], [],...
        'ShowAllTicks', true,'TickTexInterpreter','true',...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'GridLines','-');
title('150m, 1GHz','fontweight','normal');
set(gca,'ticklength',[0 0]);

subhandle3 = subplot(1,3,3);
res = load('result_fig4b_new.mat');
std_est_periodogram_lse = res.std_est_periodogram_lse;

hm3 = heatmapcus(flipud(single(std_est_periodogram_lse*100)),fnoise_list,[], [],...
        'ShowAllTicks', true,'TickTexInterpreter','true',...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'GridLines','-');
 title('75m, 0.5GHz','fontweight','normal');
   
set(gca,'ticklength',[0 0]);


h1=get(subhandle1,'Position');

plotx = h1(1);
ploty = h1(2)+0.14;

plotwidth = h1(3)*1.05;
plotheight = h1(4)*0.8;
plotpitch = plotwidth + 0.025;

set(subhandle1,'Position',[plotx ploty plotwidth plotheight]);

h1=get(subhandle1,'Position');

pos = [-0.37 h1(4)/2-0.1];
str = {'{\itP}_{RX} (W)'};
t=text(subhandle1,pos(1),pos(2), str,'fontname','times',...
    'fontsize',fs,'units','inches','rotation',90);

h2=get(subhandle2,'Position');
set(subhandle2,'Position',[plotx+plotpitch ploty plotwidth plotheight]);

h3=get(subhandle3,'Position');
set(subhandle3,'Position',[plotx+plotpitch*2 ploty plotwidth plotheight]);
    
hp4 = get(subhandle3,'Position');
cb = colorbar('Position', [plotx+plotpitch*3  ploty  0.03  hp4(4)],'fontsize',fs);
% cb.Label.String='Accuracy (m)';

pos = [1.16 0.15];
str = 'Accuracy (cm)';
t=text(subhandle3,pos(1),pos(2), str,'fontname','times',...
    'fontsize',fs,'units','inches','rotation',90);

dim1 = [plotx-0.05 0 plotwidth+0.1 0.1];
str = {'Linewidth (MHz)'};
annotation(...
    'textbox',dim1, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','middle','linestyle','none');

dim2 = [plotx+plotpitch-0.05 0.0 plotwidth+0.1 0.1];
str = {'Linewidth (MHz)'};
annotation(...
    'textbox',dim2, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','middle','linestyle','none');

dim3 = [plotx+plotpitch*2-0.05 0.0 plotwidth+0.1 0.1];
str = {'Linewidth (MHz)'};
annotation(...
    'textbox',dim3, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','middle','linestyle','none');

save2pdf('fig4.pdf',fig1handle,600);


% caxis([0 10]);

