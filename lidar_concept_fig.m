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

header;

m=100;
cm_magma=magma(m);
cm_inferno=inferno(m);
cm_plasma=plasma(m);
cm_viridis=viridis(m);

%% set constants

fig2handle=fig('units','inches','width',textwidth/2,'height',textwidth/5,'font','Times','fontsize',fs);

% target
targetw = 0.02;
targeth = 0.5;
targetx = 0.93;
targety = 0.5;
dim = [targetx targety-targeth/2 targetw targeth];
annotation('rectangle',dim,'FaceColor','k','FaceAlpha',0.5);

boxw = 0.1;
boxh = 0.25;
dim = [targetx - boxw/2 targety-targeth/2-boxh+0.05 boxw boxh];
str = {'Target'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','middle','linestyle','none');

% sense
sensew = 0.5;
senseh = 0.9;
sensex = 0.05;
sensey = 0.5;
dim = [sensex sensey-senseh/2 sensew senseh];
str = {'FMCW LIDAR'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','top','BackgroundColor','k','FaceAlpha',0.1);

% Laser
laserw = 0.15;
laserh = 0.25;
laserx = 0.3;
lasery = 0.65;
dim = [laserx lasery-laserh/2 laserw laserh];
str = {'Tunable Laser'};
annotation('textbox',dim,'String',str,'FitBoxToText','off',...
    'fontname','times','fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','middle','BackgroundColor','w');

% mod. waveform
b01w = 0.05;
b01y = lasery; %% no need to touch
b01x = 0.15;
b01margin = 0.09;

xl = [b01x - b01w b01x];
yl = [b01y-b01w/2 b01y+b01w/2];
annotation('line',xl,yl);

xl = [b01x b01x];
yl = [b01y+b01w/2 b01y-b01w/2];
annotation('line',xl,yl);

xl = [b01x b01x + b01w];
yl = [b01y-b01w/2 b01y+b01w/2];
annotation('line',xl,yl);

xl = [b01x + b01w b01x + b01w];
yl = [b01y+b01w/2 b01y-b01w/2];
annotation('line',xl,yl);

xl = [b01x+b01margin laserx-0.025];
yl = [b01y b01y];
ar=annotation('arrow',xl,yl);
ar.HeadStyle ='plain';
ar.HeadWidth =3;
ar.HeadLength =3;


% Receiver
rxw = 0.15;
rxh = 0.25;
rxx = laserx;
rxy = 0.25;
dim = [rxx rxy-rxh/2 rxw rxh];
str = {'Coherent Receiver'};
annotation('textbox',dim,'String',str,'FitBoxToText','off','fontname',...
    'times','fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','middle','BackgroundColor','w');

% sense
sensew = 0.2;
senseh = 0.25;
sensex = b01x-sensew/2;
sensey = rxy;
dim = [sensex sensey-senseh/2 sensew senseh];
str = {'Frequency Estimation'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','middle','linestyle','none');

xl = [rxx-0.025 b01x+b01margin];
yl = [sensey sensey];
ar=annotation('arrow',xl,yl);
ar.HeadStyle ='plain';
ar.HeadWidth =3;
ar.HeadLength =3;

% beam from laser
b1w = 0.07;
b1y = lasery; %% no need to touch
b1x = laserx + laserw; %% no need to touch
xl = [b1x b1x + b1w];
yl = [b1y b1y];
ar=annotation('arrow',xl,yl);
ar.HeadStyle ='plain';
ar.HeadWidth =3;
ar.HeadLength =3;
ar.Color = cg;
ar.LineWidth = 1;

% lo beam
b2y = b1y; %% no need to touch
b2x = b1x + b1w; %% no need to touch
xl = [b2x b2x];
yl = [b2y rxy;];
ar=annotation('arrow',xl,yl);
ar.HeadStyle ='plain';
ar.HeadWidth =3;
ar.HeadLength =3;
ar.Color = cr;
ar.LineWidth = 1;

boxw = 0.1;
boxh = 0.25;
dim = [b2x-boxw (b2y+rxy)/2-boxh/2 boxw boxh];
str = {'$E_\textrm{LO}$'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center','Color',cr,...
    'verticalalignment','middle','linestyle','none','Interpreter','latex');


% beam to the target
b3y = b1y; %% no need to touch
b3x = b1x + b1w; %% no need to touch
xl = [b3x targetx];
yl = [b3y targety];
ar = annotation('line',xl,yl);
ar.Color = cb;
ar.LineWidth = 1;

% beam from the target
b4y = targety; %% no need to touch
b4x = targetx; %% no need to touch
xl = [b4x b2x];
yl = [b4y rxy;];
ar=annotation('arrow',xl,yl);
ar.HeadStyle ='plain';
ar.HeadWidth =3;
ar.HeadLength =3;
ar.Color = cb;
ar.LineWidth = 1;

boxw = 0.15;
boxh = 0.2;
dim = [(b4x+b2x)/2-boxw/2 (b4y+rxy)/2-boxh boxw boxh];
str = {'$E_\textrm{RX}$'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center','Color',cb,...
    'verticalalignment','middle','linestyle','none','Interpreter','latex');

boxw = 0.4;
boxh = 0.2;
dim = [(b4x+b2x)/2-boxw/2+0.05 (b4y+lasery)/2+0.2 boxw boxh];
str = {'Round trip delay \tau=2d/c'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','middle','linestyle','none');

% distance arrow
sensew = 0.5;
sensex = 0.05;
xl = [sensex+sensew targetx];
yl = [(b4y+lasery)/2+0.1 (b4y+lasery)/2+0.1];
ar=annotation('doublearrow',xl,yl);
ar.Head1Style ='plain';
ar.Head1Width =3;
ar.Head1Length =3;
ar.Head2Style ='plain';
ar.Head2Width =3;
ar.Head2Length =3;

boxw = 0.3;
boxh = 0.2;
dim = [(b4x+b2x)/2-boxw/2 (b4y+lasery)/2+0.05 boxw boxh];
str = {'d'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',fs,'horizontalalignment','center',...
    'verticalalignment','middle','linestyle','none');

% beam to the receiver
b5y = rxy; %% no need to touch
b5x = b2x; %% no need to touch
xl = [b5x b5x-b1w];
yl = [b5y b5y];
ar=annotation('arrow',xl,yl);
ar.HeadStyle ='plain';
ar.HeadWidth =3;
ar.HeadLength =3;
ar.Color = cp;
ar.LineWidth = 1;

save2pdf('fig1a.pdf',fig2handle,600);


