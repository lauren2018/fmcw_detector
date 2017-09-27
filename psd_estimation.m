%% 
%% written by Taehwan Kim (taehwan@berkeley.edu)

clear; close all; clc;

%% open simulink model
open_system('./model/symodel');

%% set constants
c = 3e8;
timestep = 0.001e-6;
stoptime = 10e-6;
Fs = 1/timestep;
alpha = 1e16;
omega0 = 0;
phi0 = 0;

distance = 6;
tau = distance/c;
fnoise = 10e6;
snoise = 0e-8;

%% run simulation
iter = 10;
padm = 4;
x = zeros(iter, (1+2*padm)*stoptime/timestep);
f_est_periodogram = zeros(1,iter);
f_est_music = zeros(1,iter);
freq = 0:Fs/length(x):Fs/2;

for i=1:iter
    fseed = 2*i;
    sseed = 2*i+1;
    sim('symodel');
    x(i,:)=padarray(vout.Data(1:end-1),padm*stoptime/timestep,0,'both');
    
    per = periodogram(x(i,:),rectwin(length(x(i,:))),length(x(i,:)),Fs);
    [Y, I] = max(per);
    f_est_periodogram(i) = freq(I);
    
    [S, F] = pmusic(x(1,:),6,freq,Fs);
    [Y, I] = max(S);
    f_est_music(i) = F(I);
    
end

save('data.mat','x');
    
%% plot the periodogram-estmated psd for 
figure(1);
% x=vout.Data(1:end-1);
N = length(x(1,:));
xdft = fft(x(1,:));
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);

plot(freq,10*log10(psdx));

title('Photocurrent PSD');
axis([0 freq(end) -140 -40]);
xlabel('Frequency (Hz)');
grid on;
ylabel('PSD (dB/Hz)');

% figure(2);
% 
% plot(freq,S);
% 
% title('MUSIC');
% xlim([0 freq(end)]);
% xlabel('Frequency (Hz)');
% grid on;