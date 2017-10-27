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

%% open simulink model
% open_system('./model/symodel_simple');

%% set constants
q = 1.6e-19;
Res = 1;
c = 3e8;
timestep = 0.5e-9;
stoptime = 10e-6; %%%%%%%
Fs = 1/timestep;
chirpbw = 10e9; %%%%%%%
alpha = chirpbw/stoptime;
omega0 = 0;
phi0 = 0;
prx=1/sqrt(2);

N = round(stoptime/timestep);

fnoise = 1e6; %%%%%%%
tauc=1/(pi*fnoise);
dc = (c/(4*pi*fnoise))*log(stoptime*pi*fnoise);
% Prx = 10.^(-11:0.2:-9);
% Prx = 0.316227766016840e-9;
% Prx = [10^-3 10^-9];
Prx = 10^-9;
% Prx = 10^-9;
snoise = q/Res/Prx;
% snoise_list =0;
% snoise = 0;
amplitude=1;

%% run simulation
% distance = 20;
% distance = 20;
iter = 100;
padm = 0;

% distance=100;
distance = 10:10:200;
distance = [0.1 distance];
% distance = [11 22 33 44 55 66 77 88 99 110];
% distance = 120;

dist_res = 0.5*c/chirpbw;
% distance = [3 4];
tau_list = 2*distance/c;
reference_spower = stoptime * exp(-2*tau_list/tauc) / pi;

periods = (1./(tau_list*alpha))/timestep;
rightlength = floor(periods.*floor(0.9*N./periods));

x = zeros(length(distance), iter, (1+2*padm)*N);
x_windowed = zeros(length(distance), iter, (1+2*padm)*N);
signal_est_periodogram = zeros(length(distance), iter);

dist_est_periodogram = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_lse = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_mle = zeros(length(Prx),length(distance), iter);
dist_est_music = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_correct = zeros(length(Prx),length(distance), iter);
dist_est_rb = zeros(length(Prx),length(distance), iter);

std_est_periodogram = zeros(length(Prx),length(distance));
std_est_periodogram_lse = zeros(length(Prx),length(distance));
std_est_periodogram_mle = zeros(length(Prx),length(distance));
std_est_music = zeros(length(Prx),length(distance));
std_est_rb = zeros(length(Prx),length(distance));

max_est_periodogram = zeros(length(Prx),length(distance));

percollection = zeros(length(Prx),length(distance),N/2+1);


            
% save('data.mat','x');
    
%% plot the periodogram-estmated psd for 

tau_list = 2*distance/c;
freq = linspace(0,2*tau_list(end)*alpha,10001);

normalization = 0.5*(stoptime/(2*pi));
mainterm = 0.5*(stoptime/(2*pi))*exp(-2*tau_list/tauc);
pedestal = 0.5*(tauc/(2*pi))*(1-exp(-2*tau_list/tauc).*(1+2*tau_list/tauc));
sum = mainterm + pedestal;

tau = tau_list(2);
main = (stoptime/(4*pi))* exp(-2*tau/tauc)*(sinc(pi * stoptime * (freq-alpha*tau))).^2;
pedes = ( (tauc/(4*pi)) ./ (1 + (pi * tauc * (freq-alpha*tau) ).^2) ) .* ...
    (1 - exp(-2*tau/tauc) * ( cos(2*pi*tau*(freq-alpha*tau)) + (1./(pi*(freq-alpha*tau)*tauc)).*sin(2*pi*tau*(freq-alpha*tau))  ) );
% *(cos(2*pi*(freq-alpha*tau)*tau)+(2./(tauc*(freq-alpha*tau))).*sin(2*pi*(freq-alpha*tau)*tau)));
psd1 = main+pedes;

tau = tau_list(7);
main = (stoptime/(4*pi))* exp(-2*tau/tauc)*(sinc(pi * stoptime * (freq-alpha*tau))).^2;
pedes = ( (tauc/(4*pi)) ./ (1 + (pi * tauc * (freq-alpha*tau) ).^2) ) .* ...
    (1 - exp(-2*tau/tauc) * ( cos(2*pi*tau*(freq-alpha*tau)) + (1./(pi*(freq-alpha*tau)*tauc)).*sin(2*pi*tau*(freq-alpha*tau))  ) );
% ((0.5*(tauc/(2*pi))./(1+(pi*tauc*(freq-alpha*tau)).^2))) .* (1-exp(-2*tau/tauc)*(cos(2*pi*(freq-alpha*tau)*tau)+(2./(tauc*(freq-alpha*tau))).*sin(2*pi*(freq-alpha*tau)*tau)));
psd2 = main+pedes;

tau = tau_list(16);
main = (stoptime/(4*pi))* exp(-2*tau/tauc)*(sinc(pi * stoptime * (freq-alpha*tau))).^2;
pedes = ( (tauc/(4*pi)) ./ (1 + (pi * tauc * (freq-alpha*tau) ).^2) ) .* ...
    (1 - exp(-2*tau/tauc) * ( cos(2*pi*tau*(freq-alpha*tau)) + (1./(pi*(freq-alpha*tau)*tauc)).*sin(2*pi*tau*(freq-alpha*tau))  ) );
% ((0.5*(tauc/(2*pi))./(1+(pi*tauc*(freq-alpha*tau)).^2))) .* (1-exp(-2*tau/tauc)*(cos(2*pi*(freq-alpha*tau)*tau)+(2./(tauc*(freq-alpha*tau))).*sin(2*pi*(freq-alpha*tau)*tau)));
psd3 = main+pedes;

tau = tau_list(end);
main = (stoptime/(4*pi))* exp(-2*tau/tauc)*(sinc(pi * stoptime * (freq-alpha*tau))).^2;
pedes = ( (tauc/(4*pi)) ./ (1 + (pi * tauc * (freq-alpha*tau) ).^2) ) .* ...
    (1 - exp(-2*tau/tauc) * ( cos(2*pi*tau*(freq-alpha*tau)) + (1./(pi*(freq-alpha*tau)*tauc)).*sin(2*pi*tau*(freq-alpha*tau))  ) );
% ((0.5*(tauc/(2*pi))./(1+(pi*tauc*(freq-alpha*tau)).^2))) .* (1-exp(-2*tau/tauc)*(cos(2*pi*(freq-alpha*tau)*tau)+(2./(tauc*(freq-alpha*tau))).*sin(2*pi*(freq-alpha*tau)*tau)));
psd4 = main+pedes;

fig1handle=fig('units','inches','width',textwidth/2,'height',textwidth/4,'font','Times','fontsize',fs);
set(gcf,'DefaultAxesColorOrder',colororder);

subhandle1 = subplot(1,2,1);
hold on;
plot(distance,10*log10(mainterm/normalization),':','linewidth',4.1*lw,'color',csil,'markerfacecolor',cr,'markersize',ms);
plot(distance,10*log10(pedestal/normalization),'-.','linewidth',3*lw,'color',csil,'markerfacecolor',cr,'markersize',ms);
plot(distance,10*log10(sum/normalization),'-','linewidth',3*lw,'color',cb,'markerfacecolor',cb,'markersize',ms);
box off;
xlabel('Distance (m)');
ylabel('PSD peak (dB/Hz)');

lbl = {'Beat','Pedestal','Sum'};
legendflex(lbl,'anchor',{'n','n'},'buffer',[0.2 0],'box','off','xscale',0.5);

set(gca,'linewidth',lw);
set(gca,'yminortick','off');
set(gca,'ticklength',[0.02 0.02]);
axis([0 200 -50 25]);
xticks([0 100 200]);
yticks([-25 0 25]);

cond = main+pedes > 0;
condt = main+pedes <=0;
condt = condt*1e-9;
cond = cond+condt;


subhandle2 = subplot(1,2,2);
hold on;
plot(freq/1e9,10*log10((psd1/normalization)),'-','linewidth',lw,'markerfacecolor',cb,'markersize',ms);
plot(freq/1e9,10*log10((psd2/normalization)),'-','linewidth',lw,'markerfacecolor',cb,'markersize',ms);
plot(freq/1e9,10*log10((psd3/normalization)),'-','linewidth',lw,'markerfacecolor',cb,'markersize',ms);
% plot(freq/1e9,10*log10((psd4/normalization)),'-','linewidth',lw,'markerfacecolor',cb,'markersize',ms);

% plot(F/1e9,10*log10(squeeze(percollection(1,end,:))),'-','linewidth',lw,'markerfacecolor',cb,'markersize',ms);
box off;
xlabel('Frequency (GHz)');
ylabel('Power Density (dB/Hz)');

lbl = {'10m','70m','150m'};
legendflex(lbl,'anchor',{'n','n'},'buffer',[0.2 0],'bufferunit','inches','box','off','xscale',0.5);

set(gca,'linewidth',lw);
set(gca,'yminortick','off');
set(gca,'ticklength',[0.02 0.02]);
axis([0 1.5 -75 25]);
xticks([0 0.5 1 1.5]);
yticks([-50 -25 0 25]);

h1=get(subhandle1,'Position');
set(subhandle1,'Position',[h1(1)-0.03 h1(2)+0.02 h1(3)+0.045 h1(4)]);
h2=get(subhandle2,'Position');
set(subhandle2,'Position',[h2(1)+0.01 h2(2)+0.02 h2(3)+0.045 h2(4)]);

save2pdf('fig2.pdf',fig1handle,600);


