%% 
%% written by Taehwan Kim (taehwan@berkeley.edu)

clear; close all; clc;

%% open simulink model
open_system('./symodel');

%% set constants
c = 3e8;
timestep = 0.001e-6;
stoptime = 10e-6;
Fs = 1/timestep;
alpha = 1e16;
omega0 = 0;
phi0 = 0;

distance = 1;
tau = distance/c;
fnoise = 2.3e6;
snoise = 0e-8;

%% run simulation
sim('symodel');

%% plot the periodogram-estmated psd for 
x=vout.Data(1:end-1);
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;
plot(freq,10*log10(psdx));

title('Photocurrent PSD');
axis([0 freq(end) -140 -40]);
xlabel('Frequency (Hz)');
grid on;
ylabel('PSD (dB/Hz)');