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
freq = 0:Fs/N:Fs/2;
fnoise = 10e6; %%%%%%%
tauc=1/(pi*fnoise);
dc = (c/(4*pi*fnoise))*log(stoptime*pi*fnoise);
% Prx = 10.^(-11:0.2:-9);
% Prx = 0.316227766016840e-9;
Prx = [10^-3 10^-9];
% Prx = 10^-9;
snoise_list = q/Res./Prx;
% snoise_list =0;
% snoise = 0;
amplitude=1;

%% run simulation


distance = [11 22 33 44 55 66 77 88 99 110];

tau_list = 2*distance/c;

T = 1.5;
Fs = 1000;
dt = 1/Fs;
t = 0:dt:T-dt;
x = 0.5*sawtooth(2*pi*(t-0.25))+0.6;
xdel = 0.5*sawtooth(2*pi*(t-0.25-100*dt))+0.6;

% save('data.mat','x');
    
%% plot the periodogram-estmated psd for 

% set(0,'defaultTextInterpreter','latex');




fig1handle=fig('units','inches','width',textwidth/2,'height',textwidth/5,'font','Times','fontsize',fs);
set(gcf,'DefaultAxesColorOrder',colororder);

% subhandle1 = subplot(1,2,1);
hold on;
plot(t,xdel,'linewidth',lw,'color',cb);
plot(t,x,'linewidth',lw,'color',cr);
plot(t,1.1*ones(1,length(x)),'k--','linewidth',lw);
plot(t,0.1*ones(1,length(x)),'k--','linewidth',lw);
% hold on;
% plot(distance,std_est_periodogram_lse(1,:),'o-','linewidth',lw,'markerfacecolor',cb,'markersize',ms);
% plot(distance,std_est_rb(1,:),'^-','linewidth',lw,'markerfacecolor',cr,'markersize',ms);
% plot(distance,std_est_music(1,:),'+-','linewidth',lw,'markerfacecolor',cg,'markersize',ms);
box off;
xlabel('Time');
ylabel('Laser frequency');

% lbl = {'LLSE','R&B','MUSIC'};
% legendflex(lbl,'anchor',{'nw','nw'},'buffer',[0.2 0],'box','off','xscale',0.5);

set(gca,'linewidth',lw);
set(gca,'yminortick','off');
set(gca,'ticklength',[0.02 0.02]);
axis([0 T 0 1.3]);
xticks([]);
yticks([0.1 1.1]);
% i=yticklabels({'$f_{min}$', 'f_{max}'});
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'YTickLabel', {'$f_\textrm{min}$', '$f_\textrm{max}$'});

txt1 = '$E_\textrm{LO}$';
text(0.15,1.2,txt1,'Interpreter','latex','color',cr)

txt2 = '$E_\textrm{RX}$';
text(0.3,1.2,txt2,'Interpreter','latex','color',cb)

X = [0.5 0.55];
Y = [0.45 0.45];
ar=annotation('arrow',X,Y);
ar.HeadStyle ='plain';
ar.HeadWidth =3;
ar.HeadLength =3;

txt3 = '$\tau$';
text(0.57, 0.5, txt3,'Interpreter','latex')

X = [0.55 0.55];
Y = [0.45 0.51];
ar=annotation('arrow',X,Y);
ar.HeadStyle ='plain';
ar.HeadWidth =3;
ar.HeadLength =3;

txt4 = '$f_\textrm{beat}$';
text(0.62, 0.7, txt4,'Interpreter','latex')
 
xl = [0.67 0.67];
yl = [0.67 0.75];
annotation('line',xl,yl)
 
xl = [0.67 0.74];
yl = [0.75 0.75];
annotation('line',xl,yl)

txt5 = '$\gamma$';
text(0.92, 1, txt5,'Interpreter','latex')

X = [0.9 0.9];
Y = [0.23 0.85];
ar=annotation('doublearrow',X,Y);
ar.Head1Style ='plain';
ar.Head1Width =3;
ar.Head1Length =3;
ar.Head2Style ='plain';
ar.Head2Width =3;
ar.Head2Length =3;

txt6 = '$f_\textrm{BW}$';
text(1.35, 0.55, txt6,'Interpreter','latex')

X = [0.305 0.85];
Y = [0.3 0.3];
ar=annotation('doublearrow',X,Y);
ar.Head1Style ='plain';
ar.Head1Width =3;
ar.Head1Length =3;
ar.Head2Style ='plain';
ar.Head2Width =3;
ar.Head2Length =3;

txt7 = '$T$';
text(T/2+0.1,0.3, txt7,'Interpreter','latex')

save2pdf('fig1b.pdf',fig1handle,600);

